"""End-to-end orchestration for `brief-acoustics run case.toml`."""

from __future__ import annotations

from dataclasses import dataclass
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

from .case_config import CaseConfig, load_case, write_runtime_namelist
from .mesh_io import SurfaceMesh, read_fortran_mesh, write_fortran_mesh


class CaseRunError(RuntimeError):
    """Raised when build, execution, or output validation fails."""


@dataclass(frozen=True)
class SurfaceSolution:
    """Nodal fields read from the portable CSV written by Fortran."""

    node_id: np.ndarray
    particle_id: np.ndarray
    points: np.ndarray
    normals: np.ndarray
    phi_incident: np.ndarray
    phi_scattered: np.ndarray
    phi_total: np.ndarray
    dphi_incident_dn: np.ndarray
    dphi_scattered_dn: np.ndarray
    dphi_total_dn: np.ndarray
    pressure: np.ndarray
    bm_auxiliary_dpsi_dn: np.ndarray


@dataclass(frozen=True)
class RunArtifacts:
    """Paths and metadata produced by one completed public case run."""

    output_directory: Path
    surface_csv: Path
    summary_json: Path
    surface_figure: Path | None
    mesh: SurfaceMesh


def check_case(path: Path) -> tuple[CaseConfig, SurfaceMesh]:
    """Validate the TOML contract and its complete mesh without solving."""

    config = load_case(path)
    mesh = read_fortran_mesh(config.mesh.path)
    if mesh.nodes_per_element != config.mesh.nodes_per_element:
        raise CaseRunError(
            "The declared [mesh].element_type does not match the mesh connectivity: "
            f"case declares {config.mesh.element_type!r}, but the mesh has "
            f"{mesh.nodes_per_element} nodes per element ({mesh.element_type!r})."
        )
    if mesh.particle_count != 1:
        raise CaseRunError(
            "The public case runner supports exactly one closed particle; "
            f"the mesh declares {mesh.particle_count}."
        )
    return config, mesh


def run_case(
    case_path: Path,
    *,
    output_override: Path | None = None,
    build: bool = True,
    create_plot: bool = True,
) -> RunArtifacts:
    """Validate, build, solve, and post-process one public case."""

    config, mesh = check_case(case_path)
    output_directory = (
        output_override.expanduser().resolve()
        if output_override is not None
        else config.output_directory
    )
    output_directory.mkdir(parents=True, exist_ok=True)

    # Rewriting the checked mesh gives Fortran one normalized, comment-free
    # representation regardless of how the user's source file was formatted.
    normalized_mesh = output_directory / "mesh.dat"
    runtime_case = output_directory / "runtime_case.nml"
    surface_csv = output_directory / "surface_solution.csv"
    solver_log = output_directory / "solver.log"
    case_copy = output_directory / "case.toml"
    summary_json = output_directory / "summary.json"
    surface_figure = output_directory / "surface_magnitude.png" if create_plot else None

    write_fortran_mesh(mesh, normalized_mesh)
    write_runtime_namelist(config, runtime_case)
    if config.source_path != case_copy.resolve():
        shutil.copyfile(config.source_path, case_copy)

    fortran_directory = locate_fortran_directory()
    log_parts: list[str] = []
    if build:
        build_result = _run_command(
            ["make", "-C", str(fortran_directory), "run-case"],
            purpose="Fortran build",
        )
        log_parts.append("== Build ==\n" + build_result)

    executable = fortran_directory / "bin" / "run_au_case"
    if not executable.is_file():
        raise CaseRunError(
            f"Missing Fortran executable: {executable}. "
            "Run without --no-build or build `make -C fortran run-case`."
        )

    solve_result = _run_command(
        [str(executable), str(normalized_mesh), str(runtime_case), str(surface_csv)],
        purpose="Fortran solve",
        environment=_bounded_numerical_environment(),
    )
    log_parts.append("== Solve ==\n" + solve_result)
    solver_log.write_text("\n\n".join(log_parts).rstrip() + "\n", encoding="utf-8")

    solution = read_surface_solution(surface_csv)
    if len(solution.node_id) != len(mesh.nodes):
        raise CaseRunError(
            "Fortran output node count does not match the validated input mesh."
        )
    if create_plot and surface_figure is not None:
        write_surface_figure(solution, surface_figure, config.name)

    summary = build_run_summary(config, mesh, solution, create_plot=create_plot)
    summary_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return RunArtifacts(
        output_directory=output_directory,
        surface_csv=surface_csv,
        summary_json=summary_json,
        surface_figure=surface_figure,
        mesh=mesh,
    )


def locate_fortran_directory() -> Path:
    """Locate the source-tree Fortran build, with one explicit override."""

    override = os.environ.get("AU_NSBEM_FORTRAN_DIR")
    if override:
        candidate = Path(override).expanduser().resolve()
        if (candidate / "Makefile").is_file():
            return candidate
        raise CaseRunError(
            "AU_NSBEM_FORTRAN_DIR does not point to a Fortran directory "
            f"containing Makefile: {candidate}"
        )

    source_tree_candidate = Path(__file__).resolve().parents[2] / "fortran"
    if (source_tree_candidate / "Makefile").is_file():
        return source_tree_candidate
    raise CaseRunError(
        "Unable to locate the Fortran source tree. Set AU_NSBEM_FORTRAN_DIR "
        "to the repository's `fortran` directory."
    )


def read_surface_solution(path: Path) -> SurfaceSolution:
    """Read and validate the generic Fortran surface CSV."""

    if not path.is_file():
        raise CaseRunError(f"Fortran did not create the expected output: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows or reader.fieldnames is None:
        raise CaseRunError(f"Surface output is empty: {path}")

    required = {
        "node_id",
        "particle_id",
        "x",
        "y",
        "z",
        "nx",
        "ny",
        "nz",
        "phi_inc_re",
        "phi_inc_im",
        "phi_sca_re",
        "phi_sca_im",
        "phi_total_re",
        "phi_total_im",
        "dphi_inc_re",
        "dphi_inc_im",
        "dphi_sca_re",
        "dphi_sca_im",
        "dphi_total_re",
        "dphi_total_im",
        "pressure_re",
        "pressure_im",
        "bm_aux_dpsi_dn_re",
        "bm_aux_dpsi_dn_im",
    }
    missing = sorted(required - set(reader.fieldnames))
    if missing:
        raise CaseRunError(f"Surface output is missing columns: {', '.join(missing)}")

    try:
        scalar = {
            name: np.asarray([float(row[name]) for row in rows], dtype=float)
            for name in required - {"node_id", "particle_id"}
        }
        node_id = np.asarray([int(row["node_id"]) for row in rows], dtype=int)
        particle_id = np.asarray([int(row["particle_id"]) for row in rows], dtype=int)
    except (TypeError, ValueError) as error:
        raise CaseRunError(f"Surface output contains a non-numeric value: {path}") from error
    if any(not np.all(np.isfinite(values)) for values in scalar.values()):
        raise CaseRunError("Surface output contains NaN or infinite values.")

    return SurfaceSolution(
        node_id=node_id,
        particle_id=particle_id,
        points=np.column_stack([scalar["x"], scalar["y"], scalar["z"]]),
        normals=np.column_stack([scalar["nx"], scalar["ny"], scalar["nz"]]),
        phi_incident=_complex_column(scalar, "phi_inc"),
        phi_scattered=_complex_column(scalar, "phi_sca"),
        phi_total=_complex_column(scalar, "phi_total"),
        dphi_incident_dn=_complex_column(scalar, "dphi_inc"),
        dphi_scattered_dn=_complex_column(scalar, "dphi_sca"),
        dphi_total_dn=_complex_column(scalar, "dphi_total"),
        pressure=_complex_column(scalar, "pressure"),
        bm_auxiliary_dpsi_dn=_complex_column(scalar, "bm_aux_dpsi_dn"),
    )


def write_surface_figure(solution: SurfaceSolution, path: Path, case_name: str) -> None:
    """Plot total potential and pressure magnitudes on the surface nodes."""

    # Use result-local caches so the command also works in CI, containers, and
    # read-only home directories without Matplotlib/fontconfig warnings.
    os.environ.setdefault("MPLCONFIGDIR", str(path.parent / ".matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(path.parent / ".cache"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure = plt.figure(figsize=(10, 4.6))
    fields = (
        (np.abs(solution.phi_total), "|total potential|", "viridis"),
        (np.abs(solution.pressure), "|pressure|", "magma"),
    )
    coordinate_span = np.ptp(solution.points, axis=0)
    coordinate_span[coordinate_span == 0.0] = 1.0
    marker_size = max(8.0, min(28.0, 2400.0 / len(solution.node_id)))
    for index, (field, title, color_map) in enumerate(fields, start=1):
        axis = figure.add_subplot(1, 2, index, projection="3d")
        points = axis.scatter(
            solution.points[:, 0],
            solution.points[:, 1],
            solution.points[:, 2],
            c=field,
            s=marker_size,
            cmap=color_map,
        )
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        axis.set_zlabel("z")
        axis.set_title(title)
        axis.set_box_aspect(coordinate_span)
        figure.colorbar(points, ax=axis, shrink=0.7)
    figure.suptitle(case_name)
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def build_run_summary(
    config: CaseConfig,
    mesh: SurfaceMesh,
    solution: SurfaceSolution,
    *,
    create_plot: bool,
) -> dict[str, object]:
    """Build a portable machine-readable summary without absolute host paths."""

    return {
        "status": "completed",
        "case": config.name,
        "physics": {
            "density": config.physics.density,
            "sound_speed": config.physics.sound_speed,
            "frequency_hz": config.physics.frequency_hz,
            "angular_frequency": config.physics.angular_frequency,
            "wavenumber": config.physics.wavenumber,
            "time_convention": "exp(-i*omega*t)",
        },
        "solver": {
            "formulation": config.solver.formulation,
            "characteristic_length": config.solver.characteristic_length,
        },
        "boundary": {
            "type": config.boundary.kind,
            "prescribed_value": _complex_pair(config.boundary.prescribed_value),
            "robin_a": _complex_pair(config.boundary.robin_a),
            "robin_b": _complex_pair(config.boundary.robin_b),
            "robin_rhs": _complex_pair(config.boundary.robin_rhs),
        },
        "incident": {
            "type": config.incident.kind,
            "potential_amplitude": _complex_pair(config.incident.potential_amplitude),
            "direction": list(config.incident.direction),
            "standing_wave": config.incident.standing_wave,
        },
        "mesh": {
            "element_type": mesh.element_type,
            "nodes_per_element": mesh.nodes_per_element,
            "normal_orientation": config.mesh.normal_orientation,
            "nodes": len(mesh.nodes),
            "elements": len(mesh.elements),
            "particles": mesh.particle_count,
        },
        "surface_extrema": {
            "max_abs_phi_scattered": float(np.max(np.abs(solution.phi_scattered))),
            "max_abs_phi_total": float(np.max(np.abs(solution.phi_total))),
            "max_abs_dphi_total_dn": float(np.max(np.abs(solution.dphi_total_dn))),
            "max_abs_pressure": float(np.max(np.abs(solution.pressure))),
            "minimum_normal_length": float(np.min(np.linalg.norm(solution.normals, axis=1))),
            "maximum_normal_length": float(np.max(np.linalg.norm(solution.normals, axis=1))),
        },
        "outputs": {
            "case": "case.toml",
            "normalized_mesh": "mesh.dat",
            "runtime_case": "runtime_case.nml",
            "surface_solution": "surface_solution.csv",
            "surface_figure": "surface_magnitude.png" if create_plot else None,
            "solver_log": "solver.log",
        },
    }


def _run_command(
    command: list[str],
    *,
    purpose: str,
    environment: dict[str, str] | None = None,
) -> str:
    """Run one build/solve command and preserve diagnostic output on failure."""

    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    output = (completed.stdout + completed.stderr).strip()
    if completed.returncode != 0:
        detail = f"\n{output}" if output else ""
        raise CaseRunError(f"{purpose} failed with exit code {completed.returncode}.{detail}")
    return output


def _bounded_numerical_environment() -> dict[str, str]:
    """Avoid unexpectedly using every laptop/CI core for one dense solve."""

    environment = os.environ.copy()
    environment.setdefault("OMP_NUM_THREADS", "2")
    environment.setdefault("OPENBLAS_NUM_THREADS", "2")
    environment.setdefault("VECLIB_MAXIMUM_THREADS", "2")
    return environment


def _complex_column(scalar: dict[str, np.ndarray], prefix: str) -> np.ndarray:
    """Recombine `<prefix>_re` and `_im` CSV columns."""

    return scalar[f"{prefix}_re"] + 1j * scalar[f"{prefix}_im"]


def _complex_pair(value: complex) -> list[float]:
    """Represent a complex scalar portably in JSON as `[real, imaginary]`."""

    return [float(value.real), float(value.imag)]
