"""Validate the exterior solver on a non-spherical ellipsoid.

The validation uses a method of manufactured solutions.  A radiating point
source is placed at the centre of a solid triaxial ellipsoid.  Because that
source lies outside the exterior fluid domain, its field

    phi(x) = exp(i*k*r) / r

is an exact homogeneous Helmholtz solution throughout the solved domain.  We
prescribe the exact potential on the ellipsoid and compare the solved normal
derivative with the analytical gradient projected onto each mesh-node normal.

This case exercises non-spherical geometry, spatially varying Dirichlet data,
normal directions, operator assembly, and boundary-condition elimination.  It
is deliberately smooth; a sharp-edged cube needs corner-specific free-term and
normal treatment and should not be introduced as a routine regression without
that analysis.
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

from brief_acoustics.mesh_io import (  # noqa: E402
    DEFAULT_ELEMENT_TYPE,
    OUTWARD_FROM_SOLID,
    MeshNode,
    SurfaceMesh,
    exterior_free_term_sign_for,
    write_fortran_mesh,
)


def _load_sphere_helpers():
    """Load the existing sphere mesh and CSV helpers from their script."""

    spec = importlib.util.spec_from_file_location(
        "sphere_validation_helpers_for_ellipsoid",
        SPHERE_VALIDATION_DIR / "run_validation.py",
    )
    if spec is None or spec.loader is None:
        raise ImportError("Unable to load sphere validation helpers")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


_sphere_helpers = _load_sphere_helpers()
build_sphere_mesh = _sphere_helpers.build_sphere_mesh
read_nsbem_csv = _sphere_helpers.read_nsbem_csv
relative_l2 = _sphere_helpers.relative_l2
read_element_type_from_environment = _sphere_helpers.read_element_type_from_environment


@dataclass(frozen=True)
class ValidationConfig:
    """Inputs and fixed acceptance criterion for one ellipsoid solve."""

    semi_axes: tuple[float, float, float] = (1.0, 0.8, 0.6)
    source_position: tuple[float, float, float] = (0.0, 0.0, 0.0)
    wavenumber: float = 1.0
    subdivisions: int = 2
    element_type: str = DEFAULT_ELEMENT_TYPE
    normal_orientation: str = OUTWARD_FROM_SOLID
    formulation: str = "ordinary"

    # This limit is fixed after the first reproducible baselines (1.62% for the
    # ordinary linear case and 0.94% for quadratic Burton-Miller).  Five per
    # cent leaves useful cross-platform margin while still detecting material
    # sign, normal, geometry, or assembly regressions.
    maximum_relative_l2_error: float = 0.05


def main() -> int:
    """Generate the ellipsoid, run Fortran, and compare the solved derivative."""

    config = ValidationConfig(
        wavenumber=_env_float("NSBEM_WAVENUMBER", ValidationConfig.wavenumber),
        subdivisions=_env_int("NSBEM_SUBDIVISIONS", ValidationConfig.subdivisions),
        element_type=read_element_type_from_environment(),
        formulation=_env_formulation(),
    )
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    mesh = build_ellipsoid_mesh(
        semi_axes=config.semi_axes,
        subdivisions=config.subdivisions,
        element_type=config.element_type,
    )
    mesh_file = RESULTS_DIR / "ellipsoid_mesh.dat"
    numerical_file = RESULTS_DIR / "nsbem_surface.csv"
    comparison_file = RESULTS_DIR / "comparison_surface.csv"
    summary_file = RESULTS_DIR / "summary.txt"

    write_fortran_mesh(mesh, mesh_file)
    build_validation_executable()
    run_validation_executable(mesh_file, numerical_file, config)

    numerical = read_nsbem_csv(numerical_file)
    analytical = point_source_surface_fields(
        numerical["points"],
        numerical["normals"],
        wavenumber=config.wavenumber,
        source_position=np.asarray(config.source_position, dtype=float),
    )
    comparison = write_comparison_csv(
        comparison_file,
        numerical,
        analytical,
    )
    write_plots(comparison)
    passed = write_summary(summary_file, config, mesh, comparison)

    print(f"Wrote manufactured ellipsoid results to {RESULTS_DIR}")
    print(summary_file.read_text(encoding="utf-8"))
    return 0 if passed else 1


def _env_int(name: str, default: int) -> int:
    """Read an optional integer environment variable."""

    value = os.environ.get(name)
    return default if value is None else int(value)


def _env_float(name: str, default: float) -> float:
    """Read an optional floating-point environment variable."""

    value = os.environ.get(name)
    return default if value is None else float(value)


def _env_formulation() -> str:
    """Read and validate the requested final linear-system formulation."""

    value = os.environ.get("NSBEM_FORMULATION", ValidationConfig.formulation).strip().lower()
    if value == "bm":
        value = "burton-miller"
    if value not in {"ordinary", "burton-miller"}:
        raise ValueError("NSBEM_FORMULATION must be 'ordinary' or 'burton-miller'")
    return value


def build_ellipsoid_mesh(
    semi_axes: tuple[float, float, float],
    subdivisions: int,
    element_type: str,
) -> SurfaceMesh:
    """Scale a consistently oriented unit-sphere mesh into an ellipsoid.

    Positive diagonal scaling preserves element orientation.  For six-node
    elements, scaling every projected sphere midside node produces a curved
    quadratic approximation of the ellipsoid rather than a flat remesh.
    """

    axes = np.asarray(semi_axes, dtype=float)
    if axes.shape != (3,) or np.any(axes <= 0.0):
        raise ValueError("semi_axes must contain three positive values")

    sphere = build_sphere_mesh(
        radius=1.0,
        subdivisions=subdivisions,
        element_type=element_type,
    )
    nodes = tuple(
        MeshNode(
            node.node_id,
            node.x * axes[0],
            node.y * axes[1],
            node.z * axes[2],
            node.particle_id,
        )
        for node in sphere.nodes
    )
    return SurfaceMesh(
        nodes=nodes,
        elements=sphere.elements,
        nodes_per_element=sphere.nodes_per_element,
        particle_count=sphere.particle_count,
    )


def build_validation_executable() -> None:
    """Build the Fortran driver used only for the manufactured solution."""

    subprocess.run(
        ["make", "-C", str(FORTRAN_DIR), "validate-manufactured-exterior"],
        check=True,
    )


def run_validation_executable(
    mesh_file: Path,
    output_file: Path,
    config: ValidationConfig,
) -> None:
    """Run one quiet manufactured exterior solve."""

    source = config.source_position
    environment = os.environ.copy()
    environment.setdefault("OMP_NUM_THREADS", "2")
    environment.setdefault("OPENBLAS_NUM_THREADS", "2")
    environment.setdefault("VECLIB_MAXIMUM_THREADS", "2")
    subprocess.run(
        [
            str(FORTRAN_DIR / "bin" / "validate_manufactured_exterior"),
            str(mesh_file),
            str(output_file),
            f"{config.wavenumber:.17g}",
            f"{max(config.semi_axes):.17g}",
            str(exterior_free_term_sign_for(config.normal_orientation)),
            config.formulation,
            f"{source[0]:.17g}",
            f"{source[1]:.17g}",
            f"{source[2]:.17g}",
        ],
        check=True,
        env=environment,
    )


def point_source_surface_fields(
    points: np.ndarray,
    normals: np.ndarray,
    wavenumber: float,
    source_position: np.ndarray,
) -> dict[str, np.ndarray]:
    """Return exact point-source potential and mesh-normal derivative.

    For ``phi=exp(i*k*r)/r``, the Cartesian gradient is

        grad(phi) = exp(i*k*r) * (i*k*r - 1) * (x-x_s) / r**3.

    Dotting this gradient with the normals written by Fortran compares against
    the exact derivative in the same nodal-normal convention as the solver.
    """

    displacement = points - source_position[None, :]
    distance = np.linalg.norm(displacement, axis=1)
    if np.any(distance <= np.sqrt(np.finfo(float).eps)):
        raise ValueError("Point source is on or too close to a mesh node")

    phase = np.exp(1j * wavenumber * distance)
    phi = phase / distance
    gradient_factor = phase * (1j * wavenumber * distance - 1.0) / distance**3
    gradient = gradient_factor[:, None] * displacement
    dphi_dn = np.einsum("ij,ij->i", gradient, normals)
    return {"phi": phi, "dphi_dn": dphi_dn}


def write_comparison_csv(
    path: Path,
    numerical: dict[str, np.ndarray],
    analytical: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    """Write the numerical and exact nodal fields in one portable table."""

    phi_error = numerical["phi_sca"] - analytical["phi"]
    dphi_error = numerical["dphi_sca"] - analytical["dphi_dn"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "node_id",
                "x",
                "y",
                "z",
                "numerical_phi_re",
                "numerical_phi_im",
                "analytical_phi_re",
                "analytical_phi_im",
                "numerical_dphi_dn_re",
                "numerical_dphi_dn_im",
                "analytical_dphi_dn_re",
                "analytical_dphi_dn_im",
                "abs_dphi_dn_error",
            ]
        )
        for index, node_id in enumerate(numerical["node_id"]):
            point = numerical["points"][index]
            writer.writerow(
                [
                    int(node_id),
                    *point,
                    numerical["phi_sca"][index].real,
                    numerical["phi_sca"][index].imag,
                    analytical["phi"][index].real,
                    analytical["phi"][index].imag,
                    numerical["dphi_sca"][index].real,
                    numerical["dphi_sca"][index].imag,
                    analytical["dphi_dn"][index].real,
                    analytical["dphi_dn"][index].imag,
                    abs(dphi_error[index]),
                ]
            )

    return {
        "node_id": numerical["node_id"],
        "points": numerical["points"],
        "normals": numerical["normals"],
        "numerical_phi": numerical["phi_sca"],
        "analytical_phi": analytical["phi"],
        "phi_error": phi_error,
        "numerical_dphi_dn": numerical["dphi_sca"],
        "analytical_dphi_dn": analytical["dphi_dn"],
        "dphi_dn_error": dphi_error,
    }


def write_plots(comparison: dict[str, np.ndarray]) -> None:
    """Save a spatial error view and direct numerical-versus-exact comparison."""

    points = comparison["points"]
    absolute_error = np.abs(comparison["dphi_dn_error"])
    numerical = comparison["numerical_dphi_dn"]
    analytical = comparison["analytical_dphi_dn"]

    figure = plt.figure(figsize=(10, 4.6))
    spatial_axis = figure.add_subplot(1, 2, 1, projection="3d")
    scatter = spatial_axis.scatter(
        points[:, 0],
        points[:, 1],
        points[:, 2],
        c=absolute_error,
        s=12,
        cmap="viridis",
    )
    spatial_axis.set_xlabel("x")
    spatial_axis.set_ylabel("y")
    spatial_axis.set_zlabel("z")
    spatial_axis.set_title("Absolute dphi/dn error")
    figure.colorbar(scatter, ax=spatial_axis, shrink=0.7)

    comparison_axis = figure.add_subplot(1, 2, 2)
    comparison_axis.scatter(
        analytical.real,
        numerical.real,
        s=14,
        label="Real part",
        alpha=0.8,
    )
    comparison_axis.scatter(
        analytical.imag,
        numerical.imag,
        s=14,
        label="Imaginary part",
        alpha=0.8,
    )
    limits = np.array(
        [
            min(analytical.real.min(), numerical.real.min(), analytical.imag.min(), numerical.imag.min()),
            max(analytical.real.max(), numerical.real.max(), analytical.imag.max(), numerical.imag.max()),
        ]
    )
    padding = 0.05 * max(limits[1] - limits[0], 1.0)
    comparison_axis.plot(
        [limits[0] - padding, limits[1] + padding],
        [limits[0] - padding, limits[1] + padding],
        "k--",
        linewidth=1.0,
        label="Exact agreement",
    )
    comparison_axis.set_xlabel("Analytical dphi/dn")
    comparison_axis.set_ylabel("Numerical dphi/dn")
    comparison_axis.set_title("Nodal derivative comparison")
    comparison_axis.grid(True)
    comparison_axis.legend()

    figure.tight_layout()
    figure.savefig(RESULTS_DIR / "ellipsoid_derivative_comparison.png", dpi=200)
    plt.close(figure)


def write_summary(
    path: Path,
    config: ValidationConfig,
    mesh: SurfaceMesh,
    comparison: dict[str, np.ndarray],
) -> bool:
    """Apply the fixed error criterion and record the complete test context."""

    derivative_error = relative_l2(
        comparison["numerical_dphi_dn"],
        comparison["analytical_dphi_dn"],
    )
    potential_error = relative_l2(
        comparison["numerical_phi"],
        comparison["analytical_phi"],
    )
    maximum_dirichlet_residual = float(np.max(np.abs(comparison["phi_error"])))
    normal_lengths = np.linalg.norm(comparison["normals"], axis=1)
    passed = derivative_error <= config.maximum_relative_l2_error

    text = "\n".join(
        [
            "Manufactured point-source ellipsoid validation",
            "================================================",
            f"validation status = {'PASS' if passed else 'FAIL'}",
            f"validation tolerance, relative L2 dphi/dn = {config.maximum_relative_l2_error:.8e}",
            f"formulation = {config.formulation}",
            f"wavenumber = {config.wavenumber:.8g}",
            f"ellipsoid semi-axes = {config.semi_axes}",
            f"source position = {config.source_position}",
            f"mesh normal orientation = {config.normal_orientation}",
            f"mesh element type = {mesh.element_type}",
            f"mesh nodes per element = {mesh.nodes_per_element}",
            f"mesh nodes = {len(mesh.nodes)}",
            f"mesh elements = {len(mesh.elements)}",
            f"relative L2 dphi/dn error = {derivative_error:.8e}",
            f"relative L2 prescribed phi error = {potential_error:.8e}",
            f"max |Dirichlet residual| = {maximum_dirichlet_residual:.8e}",
            f"minimum nodal normal length = {normal_lengths.min():.8e}",
            f"maximum nodal normal length = {normal_lengths.max():.8e}",
            "",
        ]
    )
    path.write_text(text, encoding="utf-8")
    return passed


if __name__ == "__main__":
    raise SystemExit(main())
