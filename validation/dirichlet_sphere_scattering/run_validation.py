"""Validate a pressure-release sphere against its partial-wave solution.

The total velocity potential is prescribed as zero on the sphere:

    phi_total = phi_incident + phi_scattered = 0.

This end-to-end case exercises the Dirichlet branch that was previously checked
only at matrix level.  Set ``NSBEM_FORMULATION`` to ``ordinary`` or
``burton-miller`` to send the same physical problem through either solver.
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
RESULTS_DIR = THIS_DIR / os.environ.get("NSBEM_RESULTS_NAME", "results")

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS_DIR / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(RESULTS_DIR / ".cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(PYTHON_DIR))
sys.path.insert(0, str(ANA_DIR))

from analytic import (  # noqa: E402
    dirichlet_sphere_scattered_field,
    dirichlet_sphere_scattered_normal_derivative,
    dirichlet_sphere_total_field,
    estimate_required_terms,
)
from brief_acoustics.mesh_io import (  # noqa: E402
    DEFAULT_ELEMENT_TYPE,
    OUTWARD_FROM_SOLID,
    SurfaceMesh,
    exterior_domain_normal_sign_for,
    write_fortran_mesh,
)

_sphere_spec = importlib.util.spec_from_file_location(
    "sphere_validation_helpers",
    SPHERE_VALIDATION_DIR / "run_validation.py",
)
if _sphere_spec is None or _sphere_spec.loader is None:
    raise ImportError("Unable to load sphere validation helpers")
_sphere_helpers = importlib.util.module_from_spec(_sphere_spec)
sys.modules[_sphere_spec.name] = _sphere_helpers
_sphere_spec.loader.exec_module(_sphere_helpers)

build_sphere_mesh = _sphere_helpers.build_sphere_mesh
plot_component = _sphere_helpers.plot_component
read_nsbem_csv = _sphere_helpers.read_nsbem_csv
relative_l2 = _sphere_helpers.relative_l2
surface_theta = _sphere_helpers.surface_theta
read_element_type_from_environment = _sphere_helpers.read_element_type_from_environment


@dataclass(frozen=True)
class ValidationConfig:
    """Inputs for one pressure-release sphere validation."""

    radius: float = 1.0
    wavenumber: float = 1.0
    subdivisions: int = 2
    element_type: str = DEFAULT_ELEMENT_TYPE
    normal_orientation: str = OUTWARD_FROM_SOLID
    formulation: str = "ordinary"


def main() -> int:
    """Generate the mesh, run Fortran, and compare every surface node."""

    config = ValidationConfig(
        subdivisions=_env_int("NSBEM_SUBDIVISIONS", ValidationConfig.subdivisions),
        element_type=read_element_type_from_environment(),
        formulation=_env_formulation(),
    )
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    mesh = build_sphere_mesh(
        radius=config.radius,
        subdivisions=config.subdivisions,
        element_type=config.element_type,
    )
    mesh_file = RESULTS_DIR / "sphere_mesh.dat"
    nsbem_file = RESULTS_DIR / "nsbem_surface.csv"
    analytic_file = RESULTS_DIR / "analytic_surface.csv"
    comparison_file = RESULTS_DIR / "comparison_surface.csv"
    summary_file = RESULTS_DIR / "summary.txt"

    write_fortran_mesh(mesh, mesh_file)
    build_validation_executable()
    run_validation_executable(mesh_file, nsbem_file, config)

    nsbem = read_nsbem_csv(nsbem_file)
    analytic = evaluate_analytical_solution(nsbem, config)
    write_analytic_csv(analytic_file, analytic)
    comparison = write_comparison_csv(comparison_file, nsbem, analytic)
    write_plots(comparison)
    write_summary(summary_file, config, mesh, comparison)

    print(f"Wrote Dirichlet validation results to {RESULTS_DIR}")
    print(summary_file.read_text(encoding="utf-8"))
    return 0


def _env_int(name: str, default: int) -> int:
    """Read an optional integer environment variable."""

    value = os.environ.get(name)
    return default if value is None else int(value)


def _env_formulation() -> str:
    """Read and validate the requested final linear-system formulation."""

    value = os.environ.get("NSBEM_FORMULATION", ValidationConfig.formulation).strip().lower()
    if value == "bm":
        value = "burton-miller"
    if value not in {"ordinary", "burton-miller"}:
        raise ValueError("NSBEM_FORMULATION must be 'ordinary' or 'burton-miller'")
    return value


def build_validation_executable() -> None:
    """Build the shared sphere validation executable."""

    subprocess.run(["make", "-C", str(FORTRAN_DIR), "validate-rigid-sphere"], check=True)


def run_validation_executable(
    mesh_file: Path,
    output_file: Path,
    config: ValidationConfig,
) -> None:
    """Run Fortran with zero prescribed total potential."""

    executable = FORTRAN_DIR / "bin" / "validate_rigid_sphere"
    subprocess.run(
        [
            str(executable),
            str(mesh_file),
            str(output_file),
            f"{config.wavenumber:.17g}",
            f"{config.radius:.17g}",
            str(exterior_domain_normal_sign_for(config.normal_orientation)),
            "dirichlet",
            "1.0",
            "0.0",
            config.formulation,
        ],
        check=True,
    )


def evaluate_analytical_solution(
    nsbem: dict[str, np.ndarray],
    config: ValidationConfig,
) -> dict[str, np.ndarray]:
    """Evaluate the pressure-release partial-wave solution at mesh nodes."""

    n_terms = estimate_required_terms(config.wavenumber, config.radius, safety=30)
    scattered_normal_derivative = dirichlet_sphere_scattered_normal_derivative(
        nsbem["points"],
        nsbem["normals"],
        k=config.wavenumber,
        radius=config.radius,
        n_terms=n_terms,
    )
    incident_normal_derivative = (
        1j
        * config.wavenumber
        * nsbem["normals"][:, 0]
        * np.exp(1j * config.wavenumber * nsbem["points"][:, 0])
    )
    return {
        "node_id": nsbem["node_id"],
        "points": nsbem["points"],
        "theta": surface_theta(nsbem["points"]),
        "phi_sca": dirichlet_sphere_scattered_field(
            nsbem["points"],
            k=config.wavenumber,
            radius=config.radius,
            n_terms=n_terms,
        ),
        "phi_total": dirichlet_sphere_total_field(
            nsbem["points"],
            k=config.wavenumber,
            radius=config.radius,
            n_terms=n_terms,
        ),
        "dphi_sca": scattered_normal_derivative,
        "dphi_total": incident_normal_derivative + scattered_normal_derivative,
    }


def write_analytic_csv(path: Path, analytic: dict[str, np.ndarray]) -> None:
    """Write the analytical nodal fields as portable CSV data."""

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "node_id",
                "x",
                "y",
                "z",
                "theta",
                "phi_sca_re",
                "phi_sca_im",
                "phi_total_re",
                "phi_total_im",
                "dphi_sca_re",
                "dphi_sca_im",
                "dphi_total_re",
                "dphi_total_im",
            ]
        )
        for node_id, point, theta, phi_sca, phi_total, dphi_sca, dphi_total in zip(
            analytic["node_id"],
            analytic["points"],
            analytic["theta"],
            analytic["phi_sca"],
            analytic["phi_total"],
            analytic["dphi_sca"],
            analytic["dphi_total"],
        ):
            writer.writerow(
                [
                    int(node_id),
                    point[0],
                    point[1],
                    point[2],
                    theta,
                    phi_sca.real,
                    phi_sca.imag,
                    phi_total.real,
                    phi_total.imag,
                    dphi_sca.real,
                    dphi_sca.imag,
                    dphi_total.real,
                    dphi_total.imag,
                ]
            )


def write_comparison_csv(
    path: Path,
    nsbem: dict[str, np.ndarray],
    analytic: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    """Write numerical errors and the directly imposed Dirichlet residual."""

    result = {
        "node_id": nsbem["node_id"],
        "points": nsbem["points"],
        "theta": analytic["theta"],
        "nsbem_phi_sca": nsbem["phi_sca"],
        "analytic_phi_sca": analytic["phi_sca"],
        "error_sca": nsbem["phi_sca"] - analytic["phi_sca"],
        "nsbem_phi_total": nsbem["phi_total"],
        "analytic_phi_total": analytic["phi_total"],
        "error_total": nsbem["phi_total"] - analytic["phi_total"],
        "nsbem_dphi_sca": nsbem["dphi_sca"],
        "analytic_dphi_sca": analytic["dphi_sca"],
        "error_dphi_sca": nsbem["dphi_sca"] - analytic["dphi_sca"],
        "nsbem_dphi_total": nsbem["dphi_total"],
        "analytic_dphi_total": analytic["dphi_total"],
        "dirichlet_residual": nsbem["phi_total"],
        "normal_radial_alignment": nsbem["normal_radial_alignment"],
    }

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "node_id",
                "x",
                "y",
                "z",
                "theta",
                "nsbem_phi_sca_re",
                "nsbem_phi_sca_im",
                "analytic_phi_sca_re",
                "analytic_phi_sca_im",
                "error_sca_abs",
                "nsbem_dphi_sca_re",
                "nsbem_dphi_sca_im",
                "analytic_dphi_sca_re",
                "analytic_dphi_sca_im",
                "error_dphi_sca_abs",
                "dirichlet_residual_abs",
            ]
        )
        for index, node_id in enumerate(result["node_id"]):
            point = result["points"][index]
            writer.writerow(
                [
                    int(node_id),
                    point[0],
                    point[1],
                    point[2],
                    result["theta"][index],
                    result["nsbem_phi_sca"][index].real,
                    result["nsbem_phi_sca"][index].imag,
                    result["analytic_phi_sca"][index].real,
                    result["analytic_phi_sca"][index].imag,
                    abs(result["error_sca"][index]),
                    result["nsbem_dphi_sca"][index].real,
                    result["nsbem_dphi_sca"][index].imag,
                    result["analytic_dphi_sca"][index].real,
                    result["analytic_dphi_sca"][index].imag,
                    abs(result["error_dphi_sca"][index]),
                    abs(result["dirichlet_residual"][index]),
                ]
            )
    return result


def write_plots(comparison: dict[str, np.ndarray]) -> None:
    """Save numerical/analytical surface comparisons."""

    order = np.argsort(comparison["theta"])
    theta = comparison["theta"][order]
    plot_component(
        theta,
        comparison["nsbem_phi_sca"][order].real,
        comparison["analytic_phi_sca"][order].real,
        "Dirichlet scattered potential real part",
        "Re(phi_sca)",
        RESULTS_DIR / "surface_scattered_real.png",
    )
    plot_component(
        theta,
        comparison["nsbem_phi_sca"][order].imag,
        comparison["analytic_phi_sca"][order].imag,
        "Dirichlet scattered potential imaginary part",
        "Im(phi_sca)",
        RESULTS_DIR / "surface_scattered_imag.png",
    )
    plt.figure(figsize=(8, 5))
    plt.semilogy(theta, np.abs(comparison["error_sca"][order]), "o", markersize=3)
    plt.xlabel("theta [rad]")
    plt.ylabel("|NSBEM - analytical|")
    plt.title("Dirichlet scattered potential absolute error")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "surface_error_abs.png", dpi=200)
    plt.close()

    plot_component(
        theta,
        comparison["nsbem_dphi_sca"][order].real,
        comparison["analytic_dphi_sca"][order].real,
        "Dirichlet scattered normal derivative real part",
        "Re(dphi_sca/dn)",
        RESULTS_DIR / "surface_dphi_scattered_real.png",
    )


def write_summary(
    path: Path,
    config: ValidationConfig,
    mesh: SurfaceMesh,
    comparison: dict[str, np.ndarray],
) -> None:
    """Write the fixed pass/fail evidence for this validation."""

    rel_l2_sca = relative_l2(comparison["nsbem_phi_sca"], comparison["analytic_phi_sca"])
    rel_l2_dphi_sca = relative_l2(
        comparison["nsbem_dphi_sca"],
        comparison["analytic_dphi_sca"],
    )
    tolerance = 5.0e-2
    # phi is prescribed by the Dirichlet condition.  The solved, independent
    # quantity is dphi/dn, so only its error can validate matrix assembly.
    status = "PASS" if rel_l2_dphi_sca <= tolerance else "FAIL"

    text = "\n".join(
        [
            "Dirichlet sphere BRIEF-Acoustics validation",
            "=====================================",
            f"validation status = {status}",
            f"validation tolerance, relative L2 dphi_sca/dn = {tolerance:.8e}",
            f"radius = {config.radius:.8g}",
            f"wavenumber = {config.wavenumber:.8g}",
            f"formulation = {config.formulation}",
            f"mesh normal orientation = {config.normal_orientation}",
            f"mesh element type = {mesh.element_type}",
            f"mesh nodes per element = {mesh.nodes_per_element}",
            f"mesh nodes = {len(mesh.nodes)}",
            f"mesh elements = {len(mesh.elements)}",
            "prescribed total potential = 0",
            f"relative L2 dphi_sca/dn error = {rel_l2_dphi_sca:.8e}",
            f"relative L2 phi_sca error = {rel_l2_sca:.8e}",
            f"max |Dirichlet residual| = {float(np.max(np.abs(comparison['dirichlet_residual']))):.8e}",
            f"normal radial alignment min = {float(np.min(comparison['normal_radial_alignment'])):.8e}",
            f"normal radial alignment max = {float(np.max(comparison['normal_radial_alignment'])):.8e}",
            "",
        ]
    )
    path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
