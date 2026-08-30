"""Compare ordinary NSBEM and Burton-Miller near two sphere resonances.

Sun and Klaseboer, Fluids 8 (2023) 56, identify fictitious sphere frequencies
at ``ka=9.3558`` and ``ka=10.41712``.  This script performs a focused sweep
around those published values using the same rigid-sphere physical problem for
both formulations.  It compares the solved surface potential with the exact
partial-wave solution and also evaluates the scattered field at ``r=1.5a``.

Profiles
--------
``focused`` (default)
    Nine points around each resonance for figures and inspection.
``regression``
    The two exact resonance values only, suitable for the CI-style command.

This is deliberately smaller than the paper's 5762-node, ka=1..40 sweep with
step 0.001.  The summary states that limitation explicitly.
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
import importlib.util
import os
from pathlib import Path
import subprocess
import sys

THIS_DIR = Path(__file__).resolve().parent
MODERN_DIR = THIS_DIR.parents[1]
PYTHON_DIR = MODERN_DIR / "python"
ANA_DIR = MODERN_DIR / "ana"
FORTRAN_DIR = MODERN_DIR / "fortran"
SPHERE_VALIDATION_DIR = MODERN_DIR / "validation" / "sphere_scattering"
FIELD_HELPER_FILE = MODERN_DIR / "validation" / "two_sphere_gap" / "run_demo.py"
RESULTS_DIR = THIS_DIR / os.environ.get("NSBEM_RESULTS_NAME", "results")

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS_DIR / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(RESULTS_DIR / ".cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(PYTHON_DIR))
sys.path.insert(0, str(ANA_DIR))

from analytic import estimate_required_terms, rigid_sphere_scattered_field  # noqa: E402
from brief_acoustics.mesh_io import (  # noqa: E402
    OUTWARD_FROM_SOLID,
    QUADRATIC_ELEMENT,
    SurfaceMesh,
    exterior_domain_normal_sign_for,
    write_fortran_mesh,
)


def _load_module(name: str, path: Path):
    """Load one validation helper module without requiring package files."""

    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load validation helper: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


_sphere_helpers = _load_module(
    "sphere_validation_helpers",
    SPHERE_VALIDATION_DIR / "run_validation.py",
)
_field_helpers = _load_module("sphere_field_evaluation_helpers", FIELD_HELPER_FILE)

build_sphere_mesh = _sphere_helpers.build_sphere_mesh
read_nsbem_csv = _sphere_helpers.read_nsbem_csv
relative_l2 = _sphere_helpers.relative_l2
field_point_phi = _field_helpers.field_point_phi


PUBLISHED_RESONANCES = (9.3558, 10.41712)
FOCUSED_OFFSETS = (-0.020, -0.010, -0.005, -0.002, 0.0, 0.002, 0.005, 0.010, 0.020)


@dataclass(frozen=True)
class SweepConfig:
    """Numerical setup and fixed evidence thresholds."""

    radius: float = 1.0
    subdivisions: int = 2
    element_type: str = QUADRATIC_ELEMENT
    normal_orientation: str = OUTWARD_FROM_SOLID
    observation_radius_ratio: float = 1.5
    profile: str = "focused"

    # These broad thresholds test the physical claim without fitting to the
    # exact values of one machine or mesh.  Burton-Miller must remain usable,
    # ordinary NSBEM must show resonance pollution, and the improvement must be
    # material at both published frequencies.
    maximum_bm_surface_error: float = 0.15
    minimum_ordinary_surface_error: float = 0.20
    minimum_improvement_factor: float = 3.0


def main() -> int:
    """Run both formulations at each selected nondimensional frequency."""

    config = SweepConfig(profile=_env_profile())
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    surface_dir = RESULTS_DIR / "surface_solutions"
    surface_dir.mkdir(parents=True, exist_ok=True)

    mesh = build_sphere_mesh(
        radius=config.radius,
        subdivisions=config.subdivisions,
        element_type=config.element_type,
    )
    mesh_file = RESULTS_DIR / "sphere_mesh.dat"
    write_fortran_mesh(mesh, mesh_file)
    build_validation_executable()

    rows: list[dict[str, float | complex]] = []
    frequencies = selected_ka_values(config.profile)
    for index, ka in enumerate(frequencies, start=1):
        print(f"[{index:02d}/{len(frequencies):02d}] solving ka={ka:.8f}")
        rows.append(run_frequency(mesh, mesh_file, surface_dir, ka, config))

    comparison_file = RESULTS_DIR / "frequency_comparison.csv"
    summary_file = RESULTS_DIR / "summary.txt"
    write_comparison_csv(comparison_file, rows)
    write_plots(rows)
    write_summary(summary_file, rows, mesh, config)

    print(f"Wrote fictitious-frequency results to {RESULTS_DIR}")
    print(summary_file.read_text(encoding="utf-8"))
    return 0


def _env_profile() -> str:
    """Read either the focused evidence sweep or two-point regression mode."""

    profile = os.environ.get("NSBEM_SWEEP_PROFILE", SweepConfig.profile).strip().lower()
    if profile not in {"focused", "regression"}:
        raise ValueError("NSBEM_SWEEP_PROFILE must be 'focused' or 'regression'")
    return profile


def selected_ka_values(profile: str) -> tuple[float, ...]:
    """Return sorted, unique ka values for the selected profile."""

    if profile == "regression":
        return PUBLISHED_RESONANCES
    values = {
        round(resonance + offset, 8)
        for resonance in PUBLISHED_RESONANCES
        for offset in FOCUSED_OFFSETS
    }
    return tuple(sorted(values))


def build_validation_executable() -> None:
    """Build the shared rigid-sphere Fortran driver once."""

    subprocess.run(["make", "-C", str(FORTRAN_DIR), "validate-rigid-sphere"], check=True)


def run_frequency(
    mesh: SurfaceMesh,
    mesh_file: Path,
    surface_dir: Path,
    ka: float,
    config: SweepConfig,
) -> dict[str, float | complex]:
    """Solve ordinary and Burton-Miller systems for one ka value."""

    wavenumber = ka / config.radius
    n_terms = estimate_required_terms(wavenumber, config.radius, safety=35)
    analytical_surface = rigid_sphere_scattered_field(
        np.array([(node.x, node.y, node.z) for node in mesh.nodes], dtype=float),
        k=wavenumber,
        radius=config.radius,
        n_terms=n_terms,
    )
    observation_point = np.array(
        [config.observation_radius_ratio * config.radius, 0.0, 0.0],
        dtype=float,
    )
    analytical_observation = rigid_sphere_scattered_field(
        observation_point,
        k=wavenumber,
        radius=config.radius,
        n_terms=n_terms,
    )

    solutions: dict[str, dict[str, np.ndarray]] = {}
    observation_values: dict[str, complex] = {}
    surface_errors: dict[str, float] = {}
    observation_errors: dict[str, float] = {}

    for formulation in ("ordinary", "burton-miller"):
        output_file = surface_dir / f"{formulation}_{ka_tag(ka)}.csv"
        run_solver(mesh_file, output_file, wavenumber, config, formulation)
        solution = read_nsbem_csv(output_file)
        solutions[formulation] = solution
        surface_errors[formulation] = relative_l2(solution["phi_sca"], analytical_surface)

        # Eq. (4.1) is evaluated for the scattered field only.  The point is
        # half a radius away from the boundary, so near-singular integration is
        # not an issue in this comparison.
        observation_values[formulation] = field_point_phi(
            observation_point,
            mesh,
            {"phi": solution["phi_sca"], "dphi_dn": solution["dphi_sca"]},
            wavenumber,
            config.normal_orientation,
        )
        observation_errors[formulation] = abs(
            observation_values[formulation] - analytical_observation
        ) / max(abs(analytical_observation), 1.0e-14)

    bm_error = surface_errors["burton-miller"]
    improvement = surface_errors["ordinary"] / max(bm_error, 1.0e-14)
    return {
        "ka": ka,
        "ordinary_surface_error": surface_errors["ordinary"],
        "bm_surface_error": bm_error,
        "surface_improvement_factor": improvement,
        "analytic_observation": complex(analytical_observation),
        "ordinary_observation": observation_values["ordinary"],
        "bm_observation": observation_values["burton-miller"],
        "ordinary_observation_error": observation_errors["ordinary"],
        "bm_observation_error": observation_errors["burton-miller"],
    }


def run_solver(
    mesh_file: Path,
    output_file: Path,
    wavenumber: float,
    config: SweepConfig,
    formulation: str,
) -> None:
    """Run one quiet Fortran solve and include its output if it fails."""

    executable = FORTRAN_DIR / "bin" / "validate_rigid_sphere"
    environment = os.environ.copy()
    environment.setdefault("OMP_NUM_THREADS", "2")
    environment.setdefault("OPENBLAS_NUM_THREADS", "2")
    environment.setdefault("VECLIB_MAXIMUM_THREADS", "2")
    completed = subprocess.run(
        [
            str(executable),
            str(mesh_file),
            str(output_file),
            f"{wavenumber:.17g}",
            f"{config.radius:.17g}",
            str(exterior_domain_normal_sign_for(config.normal_orientation)),
            "neumann",
            "1.0",
            "0.5",
            formulation,
        ],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"{formulation} solve failed at k={wavenumber:.8g}:\n"
            f"{completed.stdout}{completed.stderr}"
        )


def ka_tag(ka: float) -> str:
    """Create a stable filename fragment without a decimal point."""

    return f"ka_{ka:.8f}".replace(".", "p")


def write_comparison_csv(path: Path, rows: list[dict[str, float | complex]]) -> None:
    """Write all scalar evidence in one machine-readable table."""

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "ka",
                "ordinary_surface_relative_l2",
                "bm_surface_relative_l2",
                "surface_improvement_factor",
                "analytic_observation_re",
                "analytic_observation_im",
                "ordinary_observation_re",
                "ordinary_observation_im",
                "bm_observation_re",
                "bm_observation_im",
                "ordinary_observation_relative_error",
                "bm_observation_relative_error",
            ]
        )
        for row in rows:
            analytic = complex(row["analytic_observation"])
            ordinary = complex(row["ordinary_observation"])
            bm = complex(row["bm_observation"])
            writer.writerow(
                [
                    row["ka"],
                    row["ordinary_surface_error"],
                    row["bm_surface_error"],
                    row["surface_improvement_factor"],
                    analytic.real,
                    analytic.imag,
                    ordinary.real,
                    ordinary.imag,
                    bm.real,
                    bm.imag,
                    row["ordinary_observation_error"],
                    row["bm_observation_error"],
                ]
            )


def write_plots(rows: list[dict[str, float | complex]]) -> None:
    """Plot surface errors and field values in two local frequency windows.

    The two resonances are about one nondimensional frequency unit apart, while
    each sweep covers only ``+/-0.02`` around its resonance.  Separate panels
    therefore show the local peaks without drawing a misleading line across
    the unsampled interval between them.
    """

    ka = np.array([float(row["ka"]) for row in rows])
    ordinary_error = np.array([float(row["ordinary_surface_error"]) for row in rows])
    bm_error = np.array([float(row["bm_surface_error"]) for row in rows])
    analytic_observation = np.array([complex(row["analytic_observation"]) for row in rows])
    ordinary_observation = np.array([complex(row["ordinary_observation"]) for row in rows])
    bm_observation = np.array([complex(row["bm_observation"]) for row in rows])

    figure, axes = local_frequency_axes(
        title="Rigid sphere near published fictitious frequencies"
    )
    for index, (axis, resonance) in enumerate(zip(axes, PUBLISHED_RESONANCES)):
        mask = local_frequency_mask(ka, resonance)
        axis.semilogy(
            ka[mask], ordinary_error[mask], "o-", label="Ordinary NSBEM", linewidth=1.3
        )
        axis.semilogy(
            ka[mask], bm_error[mask], "s-", label="Burton-Miller", linewidth=1.3
        )
        add_resonance_marker(axis, resonance)
        if index == 0:
            axis.set_ylabel("Surface relative L2 error")
            axis.legend()
    figure.tight_layout()
    figure.savefig(RESULTS_DIR / "surface_error_sweep.png", dpi=200)
    plt.close(figure)

    figure, axes = local_frequency_axes(
        title="Scattered potential at the downstream observation point"
    )
    for index, (axis, resonance) in enumerate(zip(axes, PUBLISHED_RESONANCES)):
        mask = local_frequency_mask(ka, resonance)
        axis.plot(
            ka[mask],
            np.abs(analytic_observation[mask]),
            "k-",
            label="Analytical",
            linewidth=1.5,
        )
        axis.plot(
            ka[mask],
            np.abs(ordinary_observation[mask]),
            "o-",
            label="Ordinary NSBEM",
            linewidth=1.2,
        )
        axis.plot(
            ka[mask],
            np.abs(bm_observation[mask]),
            "s-",
            label="Burton-Miller",
            linewidth=1.2,
        )
        add_resonance_marker(axis, resonance)
        if index == 0:
            axis.set_ylabel("|phi_sc| at r=1.5a")
            axis.legend()
    figure.tight_layout()
    figure.savefig(RESULTS_DIR / "observation_magnitude_sweep.png", dpi=200)
    plt.close(figure)


def local_frequency_axes(title: str):
    """Create two panels with matching scales for the resonance windows."""

    figure, axes = plt.subplots(1, 2, figsize=(10, 4.8), sharey=True)
    figure.suptitle(title)
    half_width = max(abs(offset) for offset in FOCUSED_OFFSETS) + 0.002
    for axis, resonance in zip(axes, PUBLISHED_RESONANCES):
        axis.set_xlim(resonance - half_width, resonance + half_width)
        axis.set_xlabel("ka")
        axis.set_title(f"Near ka = {resonance:.5f}")
        axis.grid(True)
    return figure, axes


def local_frequency_mask(ka: np.ndarray, resonance: float) -> np.ndarray:
    """Select only points belonging to one focused resonance window."""

    half_width = max(abs(offset) for offset in FOCUSED_OFFSETS) + 1.0e-10
    return np.abs(ka - resonance) <= half_width


def add_resonance_marker(axis, resonance: float) -> None:
    """Mark one theoretical value quoted by the 2023 paper."""

    axis.axvline(
        resonance,
        color="0.45",
        linestyle="--",
        linewidth=1.0,
        label="Published resonance",
    )


def write_summary(
    path: Path,
    rows: list[dict[str, float | complex]],
    mesh: SurfaceMesh,
    config: SweepConfig,
) -> None:
    """Apply fixed resonance checks and document the evidence scope."""

    resonance_rows = [
        min(rows, key=lambda row: abs(float(row["ka"]) - resonance))
        for resonance in PUBLISHED_RESONANCES
    ]
    checks: list[bool] = []
    detail: list[str] = []
    for resonance, row in zip(PUBLISHED_RESONANCES, resonance_rows):
        ordinary_error = float(row["ordinary_surface_error"])
        bm_error = float(row["bm_surface_error"])
        improvement = float(row["surface_improvement_factor"])
        checks.extend(
            [
                bm_error <= config.maximum_bm_surface_error,
                ordinary_error >= config.minimum_ordinary_surface_error,
                improvement >= config.minimum_improvement_factor,
            ]
        )
        detail.extend(
            [
                f"ka = {resonance:.8f}",
                f"  ordinary surface relative L2 error = {ordinary_error:.8e}",
                f"  Burton-Miller surface relative L2 error = {bm_error:.8e}",
                f"  surface error improvement factor = {improvement:.8e}",
                f"  ordinary observation relative error = {float(row['ordinary_observation_error']):.8e}",
                f"  Burton-Miller observation relative error = {float(row['bm_observation_error']):.8e}",
            ]
        )

    status = "PASS" if all(checks) else "FAIL"
    text = "\n".join(
        [
            "Rigid-sphere fictitious-frequency comparison",
            "=============================================",
            f"validation status = {status}",
            f"sweep profile = {config.profile}",
            f"frequency points = {len(rows)}",
            f"mesh element type = {mesh.element_type}",
            f"mesh nodes per element = {mesh.nodes_per_element}",
            f"mesh nodes = {len(mesh.nodes)}",
            f"mesh elements = {len(mesh.elements)}",
            f"maximum accepted Burton-Miller surface error = {config.maximum_bm_surface_error:.8e}",
            f"minimum required ordinary surface error = {config.minimum_ordinary_surface_error:.8e}",
            f"minimum required improvement factor = {config.minimum_improvement_factor:.8e}",
            "",
            *detail,
            "",
            "Evidence scope",
            "--------------",
            "Published source: Sun and Klaseboer, Fluids 8 (2023) 56.",
            "Published resonance values: ka=9.3558 and ka=10.41712.",
            "Observation point: r=1.5a along the incident-wave direction.",
            "This focused 642-node sweep is a regression, not the paper's",
            "5762-node ka=1..40 sweep with step 0.001.",
            "",
        ]
    )
    path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
