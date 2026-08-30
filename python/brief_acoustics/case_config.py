"""Public TOML case contract for the single-particle exterior solver.

Python owns this user-facing layer so invalid or contradictory inputs are
rejected before Fortran allocates dense matrices.  The validated dataclasses
are then translated into a small generated Fortran namelist.  Users should edit
the TOML file, not the generated namelist.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import re
from typing import Any, Mapping

from .mesh_io import (
    DEFAULT_ELEMENT_TYPE,
    OUTWARD_FROM_SOLID,
    exterior_free_term_sign_for,
    nodes_per_element_for,
)

try:  # Python 3.11 and newer.
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - used only on Python 3.10.
    import tomli as tomllib  # type: ignore[no-redef]


class CaseConfigurationError(ValueError):
    """Raised when a public case file is incomplete or internally inconsistent."""


@dataclass(frozen=True)
class PhysicsConfig:
    """Non-redundant acoustic inputs and their derived frequency quantities."""

    density: float
    sound_speed: float
    frequency_hz: float
    angular_frequency: float
    wavenumber: float


@dataclass(frozen=True)
class MeshConfig:
    """Resolved mesh path and readable public geometry conventions."""

    path: Path
    element_type: str = DEFAULT_ELEMENT_TYPE
    normal_orientation: str = OUTWARD_FROM_SOLID

    @property
    def nodes_per_element(self) -> int:
        """Integer connectivity width required by the Fortran mesh."""

        return nodes_per_element_for(self.element_type)

    @property
    def exterior_free_term_sign(self) -> int:
        """Internal equation sign derived from the stated mesh orientation."""

        return exterior_free_term_sign_for(self.normal_orientation)


@dataclass(frozen=True)
class SolverConfig:
    """Numerical formulation and its characteristic geometric scale."""

    formulation: str
    characteristic_length: float

    def burton_miller_coupling_length(self, wavenumber: float) -> float:
        """Return the dimensional Burton-Miller scale ``min(a, 1/k)``.

        The normal-derivative boundary integral equation contains one more
        inverse length than the ordinary equation.  Its coupling coefficient
        must therefore be a length, not an arbitrary dimensionless number.
        ``characteristic_length`` is the body scale ``a`` in mesh units, while
        ``1 / wavenumber`` supplies the high-frequency wave scale.
        """

        if not math.isfinite(wavenumber) or wavenumber <= 0.0:
            raise ValueError("Burton-Miller coupling requires a positive wavenumber.")
        return min(self.characteristic_length, 1.0 / wavenumber)


@dataclass(frozen=True)
class BoundaryConfig:
    """One constant total-field boundary condition for the closed particle."""

    kind: str
    prescribed_value: complex
    robin_a: complex
    robin_b: complex
    robin_rhs: complex


@dataclass(frozen=True)
class IncidentConfig:
    """Either no incident field or one acoustic plane wave."""

    kind: str
    potential_amplitude: complex
    direction: tuple[float, float, float]
    standing_wave: bool


@dataclass(frozen=True)
class CaseConfig:
    """Fully resolved input needed by the public Python/Fortran workflow."""

    source_path: Path
    name: str
    mesh: MeshConfig
    output_directory: Path
    physics: PhysicsConfig
    solver: SolverConfig
    boundary: BoundaryConfig
    incident: IncidentConfig


def load_case(path: Path) -> CaseConfig:
    """Parse, validate, and resolve one user TOML case file."""

    source_path = path.expanduser().resolve()
    if not source_path.is_file():
        raise FileNotFoundError(f"Case file does not exist: {source_path}")

    try:
        with source_path.open("rb") as handle:
            root = tomllib.load(handle)
    except tomllib.TOMLDecodeError as error:
        raise CaseConfigurationError(f"Invalid TOML in {source_path}: {error}") from error

    if not isinstance(root, dict):
        raise CaseConfigurationError("The case file must contain TOML tables.")
    _reject_unknown_keys(
        root,
        {"case", "mesh", "physics", "solver", "boundary", "incident", "output"},
        "top level",
    )

    case_table = _required_table(root, "case")
    _reject_unknown_keys(case_table, {"name"}, "[case]")
    name = _required_text(case_table, "name", "[case]")
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", name) is None:
        raise CaseConfigurationError(
            "[case].name must use only letters, digits, '.', '_' or '-', "
            "and must start with a letter or digit."
        )

    mesh = _parse_mesh(_required_table(root, "mesh"), source_path.parent)

    physics = _parse_physics(_required_table(root, "physics"))
    solver = _parse_solver(_required_table(root, "solver"))
    boundary = _parse_boundary(_required_table(root, "boundary"))
    incident = _parse_incident(_required_table(root, "incident"))

    output_table = _optional_table(root, "output")
    _reject_unknown_keys(output_table, {"directory"}, "[output]")
    output_value = output_table.get("directory", f"results/{name}")
    if not isinstance(output_value, str) or not output_value.strip():
        raise CaseConfigurationError("[output].directory must be a non-empty path string.")
    output_directory = _resolve_relative_path(source_path.parent, output_value)

    return CaseConfig(
        source_path=source_path,
        name=name,
        mesh=mesh,
        output_directory=output_directory,
        physics=physics,
        solver=solver,
        boundary=boundary,
        incident=incident,
    )


def write_runtime_namelist(config: CaseConfig, path: Path) -> None:
    """Write the private, generated namelist consumed by the Fortran driver."""

    boundary = config.boundary
    incident = config.incident
    lines = [
        "! Generated by brief-acoustics. Edit the source case.toml, not this file.\n",
        "&au_case_input\n",
        f"  wavenumber = {_fortran_real(config.physics.wavenumber)}\n",
        f"  angular_frequency = {_fortran_real(config.physics.angular_frequency)}\n",
        f"  exterior_density = {_fortran_real(config.physics.density)}\n",
        f"  exterior_sound_speed = {_fortran_real(config.physics.sound_speed)}\n",
        f"  exterior_free_term_sign = {config.mesh.exterior_free_term_sign}\n",
        f"  characteristic_length = {_fortran_real(config.solver.characteristic_length)}\n",
        f"  formulation_mode = '{config.solver.formulation}'\n",
        f"  boundary_mode = '{boundary.kind}'\n",
        f"  prescribed_boundary_data = {_fortran_complex(boundary.prescribed_value)}\n",
        f"  robin_a = {_fortran_complex(boundary.robin_a)}\n",
        f"  robin_b = {_fortran_complex(boundary.robin_b)}\n",
        f"  robin_rhs = {_fortran_complex(boundary.robin_rhs)}\n",
        f"  incident_mode = '{incident.kind}'\n",
        "  incident_potential_amplitude = "
        f"{_fortran_complex(incident.potential_amplitude)}\n",
        "  incident_direction = "
        + ", ".join(_fortran_real(value) for value in incident.direction)
        + "\n",
        f"  standing_wave = {'.true.' if incident.standing_wave else '.false.'}\n",
        "/\n",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines), encoding="ascii")


def _parse_mesh(table: Mapping[str, Any], base_directory: Path) -> MeshConfig:
    """Parse readable element and normal conventions next to the mesh path."""

    _reject_unknown_keys(
        table,
        {"path", "element_type", "normal_orientation"},
        "[mesh]",
    )
    mesh_path = _resolve_relative_path(
        base_directory,
        _required_text(table, "path", "[mesh]"),
    )
    if not mesh_path.is_file():
        raise FileNotFoundError(f"Mesh file does not exist: {mesh_path}")

    raw_element_type = table.get("element_type", DEFAULT_ELEMENT_TYPE)
    if not isinstance(raw_element_type, str) or not raw_element_type.strip():
        raise CaseConfigurationError(
            "[mesh].element_type must be 'linear' or 'quadratic'."
        )
    element_type = raw_element_type.strip().lower()
    try:
        nodes_per_element_for(element_type)
    except ValueError as error:
        raise CaseConfigurationError("[mesh].element_type must be 'linear' or 'quadratic'.") from error

    normal_orientation = (
        _required_text(table, "normal_orientation", "[mesh]")
        .lower()
        .replace("_", "-")
    )
    if normal_orientation != OUTWARD_FROM_SOLID:
        raise CaseConfigurationError(
            "[mesh].normal_orientation must be 'outward-from-solid' for "
            "the public exterior solver."
        )
    return MeshConfig(mesh_path, element_type, normal_orientation)


def _parse_physics(table: Mapping[str, Any]) -> PhysicsConfig:
    """Accept one frequency representation and derive the other two."""

    allowed = {
        "density",
        "sound_speed",
        "frequency_hz",
        "angular_frequency",
        "wavenumber",
    }
    _reject_unknown_keys(table, allowed, "[physics]")
    density = _number(table, "density", "[physics]", required=True)
    sound_speed = _number(table, "sound_speed", "[physics]", required=True)
    if density <= 0.0:
        raise CaseConfigurationError("[physics].density must be positive.")
    if sound_speed <= 0.0:
        raise CaseConfigurationError("[physics].sound_speed must be positive.")

    frequency_keys = [
        key
        for key in ("frequency_hz", "angular_frequency", "wavenumber")
        if key in table
    ]
    if len(frequency_keys) != 1:
        raise CaseConfigurationError(
            "[physics] must define exactly one of frequency_hz, "
            "angular_frequency, or wavenumber."
        )

    selected = frequency_keys[0]
    selected_value = _number(table, selected, "[physics]", required=True)
    if selected_value <= 0.0:
        raise CaseConfigurationError(f"[physics].{selected} must be positive.")

    if selected == "frequency_hz":
        frequency_hz = selected_value
        angular_frequency = 2.0 * math.pi * frequency_hz
        wavenumber = angular_frequency / sound_speed
    elif selected == "angular_frequency":
        angular_frequency = selected_value
        frequency_hz = angular_frequency / (2.0 * math.pi)
        wavenumber = angular_frequency / sound_speed
    else:
        wavenumber = selected_value
        angular_frequency = wavenumber * sound_speed
        frequency_hz = angular_frequency / (2.0 * math.pi)

    return PhysicsConfig(
        density=density,
        sound_speed=sound_speed,
        frequency_hz=frequency_hz,
        angular_frequency=angular_frequency,
        wavenumber=wavenumber,
    )


def _parse_solver(table: Mapping[str, Any]) -> SolverConfig:
    """Parse only the solver modes certified for the first public workflow."""

    _reject_unknown_keys(
        table,
        {"formulation", "characteristic_length"},
        "[solver]",
    )
    formulation = _required_text(table, "formulation", "[solver]").lower()
    if formulation == "bm":
        formulation = "burton-miller"
    if formulation not in {"ordinary", "burton-miller"}:
        raise CaseConfigurationError(
            "[solver].formulation must be 'ordinary' or 'burton-miller'."
        )

    characteristic_length = _number(
        table,
        "characteristic_length",
        "[solver]",
        required=True,
    )
    if characteristic_length <= 0.0:
        raise CaseConfigurationError("[solver].characteristic_length must be positive.")
    return SolverConfig(formulation, characteristic_length)


def _parse_boundary(table: Mapping[str, Any]) -> BoundaryConfig:
    """Parse one constant total-field boundary condition."""

    kind = _required_text(table, "type", "[boundary]").lower()
    if kind not in {"dirichlet", "neumann", "robin"}:
        raise CaseConfigurationError(
            "[boundary].type must be 'dirichlet', 'neumann', or 'robin'."
        )

    if kind in {"dirichlet", "neumann"}:
        _reject_unknown_keys(table, {"type", "value_real", "value_imag"}, "[boundary]")
        prescribed_value = _complex_value(table, "value", "[boundary]")
        return BoundaryConfig(kind, prescribed_value, 0.0j, 0.0j, 0.0j)

    _reject_unknown_keys(
        table,
        {
            "type",
            "a_real",
            "a_imag",
            "b_real",
            "b_imag",
            "rhs_real",
            "rhs_imag",
        },
        "[boundary]",
    )
    robin_a = _complex_value(table, "a", "[boundary]")
    robin_b = _complex_value(table, "b", "[boundary]")
    robin_rhs = _complex_value(table, "rhs", "[boundary]")
    if max(abs(robin_a), abs(robin_b)) <= 1.0e-14:
        raise CaseConfigurationError(
            "Robin coefficients a and b cannot both be zero."
        )
    return BoundaryConfig(kind, 0.0j, robin_a, robin_b, robin_rhs)


def _parse_incident(table: Mapping[str, Any]) -> IncidentConfig:
    """Parse no incident field or one travelling/standing plane wave."""

    kind = _required_text(table, "type", "[incident]").lower().replace("_", "-")
    if kind == "none":
        _reject_unknown_keys(table, {"type"}, "[incident]")
        return IncidentConfig("none", 0.0j, (1.0, 0.0, 0.0), False)
    if kind != "plane-wave":
        raise CaseConfigurationError("[incident].type must be 'none' or 'plane-wave'.")

    _reject_unknown_keys(
        table,
        {
            "type",
            "potential_amplitude_real",
            "potential_amplitude_imag",
            "direction",
            "standing_wave",
        },
        "[incident]",
    )
    potential_amplitude = _complex_value(
        table,
        "potential_amplitude",
        "[incident]",
        default_real=1.0,
    )
    raw_direction = table.get("direction", [1.0, 0.0, 0.0])
    if not isinstance(raw_direction, list) or len(raw_direction) != 3:
        raise CaseConfigurationError("[incident].direction must contain three numbers.")
    direction = tuple(
        _finite_number(value, f"[incident].direction[{index}]")
        for index, value in enumerate(raw_direction)
    )
    if math.sqrt(sum(value * value for value in direction)) <= 1.0e-14:
        raise CaseConfigurationError("[incident].direction must not be the zero vector.")
    standing_wave = table.get("standing_wave", False)
    if not isinstance(standing_wave, bool):
        raise CaseConfigurationError("[incident].standing_wave must be true or false.")
    return IncidentConfig("plane-wave", potential_amplitude, direction, standing_wave)


def _required_table(root: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    """Return a required TOML table with a useful error on type mismatch."""

    if name not in root:
        raise CaseConfigurationError(f"Missing required [{name}] table.")
    value = root[name]
    if not isinstance(value, dict):
        raise CaseConfigurationError(f"[{name}] must be a TOML table.")
    return value


def _optional_table(root: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    """Return an optional TOML table or an empty mapping."""

    if name not in root:
        return {}
    value = root[name]
    if not isinstance(value, dict):
        raise CaseConfigurationError(f"[{name}] must be a TOML table.")
    return value


def _reject_unknown_keys(
    table: Mapping[str, Any],
    allowed: set[str],
    context: str,
) -> None:
    """Reject misspelled fields rather than silently ignoring user intent."""

    unknown = sorted(set(table) - allowed)
    if unknown:
        raise CaseConfigurationError(
            f"Unknown key(s) in {context}: {', '.join(unknown)}."
        )


def _required_text(table: Mapping[str, Any], key: str, context: str) -> str:
    """Read a non-empty string value."""

    value = table.get(key)
    if not isinstance(value, str) or not value.strip():
        raise CaseConfigurationError(f"{context}.{key} must be a non-empty string.")
    return value.strip()


def _number(
    table: Mapping[str, Any],
    key: str,
    context: str,
    *,
    required: bool,
    default: float = 0.0,
) -> float:
    """Read one finite real number, explicitly rejecting booleans."""

    if key not in table:
        if required:
            raise CaseConfigurationError(f"Missing required value {context}.{key}.")
        return default
    return _finite_number(table[key], f"{context}.{key}")


def _finite_number(value: Any, label: str) -> float:
    """Convert an int/float to float and reject NaN or infinity."""

    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise CaseConfigurationError(f"{label} must be a real number.")
    result = float(value)
    if not math.isfinite(result):
        raise CaseConfigurationError(f"{label} must be finite.")
    return result


def _complex_value(
    table: Mapping[str, Any],
    prefix: str,
    context: str,
    *,
    default_real: float = 0.0,
) -> complex:
    """Build one complex number from `<prefix>_real` and `_imag` fields."""

    real_part = _number(
        table,
        f"{prefix}_real",
        context,
        required=False,
        default=default_real,
    )
    imag_part = _number(
        table,
        f"{prefix}_imag",
        context,
        required=False,
        default=0.0,
    )
    return complex(real_part, imag_part)


def _resolve_relative_path(base_directory: Path, value: str) -> Path:
    """Resolve paths relative to the case file, never the caller's cwd."""

    candidate = Path(value).expanduser()
    if not candidate.is_absolute():
        candidate = base_directory / candidate
    return candidate.resolve()


def _fortran_real(value: float) -> str:
    """Write enough decimal digits to round-trip a double-precision value."""

    return f"{value:.17e}"


def _fortran_complex(value: complex) -> str:
    """Use Fortran namelist complex syntax `(real,imaginary)`."""

    return f"({_fortran_real(value.real)},{_fortran_real(value.imag)})"
