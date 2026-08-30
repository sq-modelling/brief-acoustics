"""Validate BRIEF-Acoustics against rigid-sphere plane-wave scattering.

This script is deliberately written as a readable demonstration rather than a
general-purpose framework.  It performs the complete validation workflow:

1. Generate a spherical triangular surface mesh in Python.
2. Write the mesh in the plain text format read by modern Fortran.
3. Build and run the Fortran validation executable.
4. Read the numerical surface solution written by Fortran.
5. Evaluate the analytical rigid-sphere solution at the same mesh nodes.
6. Save CSV files, plots, and a compact pass/fail summary.

Environment variables for quick experiments:

    NSBEM_ELEMENT_TYPE=linear or quadratic
    NSBEM_SUBDIVISIONS=<integer>
    NSBEM_WAVENUMBER=<positive float>
    NSBEM_RADIUS=<positive float>
    NSBEM_RESULTS_NAME=<output-folder-name>
    NSBEM_FORMULATION=ordinary or burton-miller
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
import math
import os
from pathlib import Path
import subprocess
import sys

THIS_DIR = Path(__file__).resolve().parent
MODERN_DIR = THIS_DIR.parents[1]
PROJECT_DIR = MODERN_DIR.parent
PYTHON_DIR = MODERN_DIR / "python"
ANA_DIR = MODERN_DIR / "ana"
FORTRAN_DIR = MODERN_DIR / "fortran"
RESULTS_DIR = THIS_DIR / os.environ.get("NSBEM_RESULTS_NAME", "results")

os.environ.setdefault("MPLCONFIGDIR", str(RESULTS_DIR / ".matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(RESULTS_DIR / ".cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


sys.path.insert(0, str(PYTHON_DIR))
sys.path.insert(0, str(ANA_DIR))

from analytic import estimate_required_terms, rigid_sphere_scattered_field, rigid_sphere_total_field  # noqa: E402
from brief_acoustics.mesh_io import (  # noqa: E402
    DEFAULT_ELEMENT_TYPE,
    LINEAR_ELEMENT,
    OUTWARD_FROM_SOLID,
    QUADRATIC_ELEMENT,
    MeshElement,
    MeshNode,
    SurfaceMesh,
    exterior_domain_normal_sign_for,
    write_fortran_mesh,
)


@dataclass(frozen=True)
class ValidationConfig:
    """Input parameters for the rigid-sphere validation case."""

    radius: float = 1.0
    wavenumber: float = 1.0
    subdivisions: int = 2
    element_type: str = DEFAULT_ELEMENT_TYPE
    normal_orientation: str = OUTWARD_FROM_SOLID
    formulation: str = "ordinary"


def main() -> int:
    """Run the complete validation workflow."""

    # Quadratic elements are the public default. Environment variables retain
    # an explicit linear path for convergence studies and regression coverage.
    config = ValidationConfig(
        radius=_env_float("NSBEM_RADIUS", ValidationConfig.radius),
        wavenumber=_env_float("NSBEM_WAVENUMBER", ValidationConfig.wavenumber),
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

    # Python owns mesh generation for this workflow.  Fortran reads the mesh and does
    # the geometry, quadrature, boundary-integral assembly, and solve.
    write_fortran_mesh(mesh, mesh_file)
    build_validation_executable()
    run_validation_executable(mesh_file, nsbem_file, config)

    # Compare the numerical result with the analytical solution at exactly the
    # same surface nodes.  This avoids interpolation error in the validation.
    nsbem = read_nsbem_csv(nsbem_file)
    analytic = evaluate_analytical_solution(nsbem, config)
    write_analytic_csv(analytic_file, analytic)
    comparison = write_comparison_csv(comparison_file, nsbem, analytic)
    write_plots(comparison)
    write_summary(summary_file, config, mesh, comparison)

    print(f"Wrote validation results to {RESULTS_DIR}")
    print(summary_file.read_text(encoding="utf-8"))
    return 0


def _env_int(name: str, default: int) -> int:
    """Read an optional integer environment variable."""

    value = os.environ.get(name)
    if value is None:
        return default
    return int(value)


def _env_float(name: str, default: float) -> float:
    """Read an optional floating-point environment variable."""

    value = os.environ.get(name)
    if value is None:
        return default
    return float(value)


def read_element_type_from_environment() -> str:
    """Read a readable interpolation type rather than a misleading order."""

    value = os.environ.get("NSBEM_ELEMENT_TYPE", ValidationConfig.element_type).strip().lower()
    if value not in {LINEAR_ELEMENT, QUADRATIC_ELEMENT}:
        raise ValueError("NSBEM_ELEMENT_TYPE must be 'linear' or 'quadratic'")
    return value


def _env_formulation() -> str:
    """Read and validate the requested boundary-integral formulation."""

    value = os.environ.get("NSBEM_FORMULATION", ValidationConfig.formulation).strip().lower()
    if value == "bm":
        value = "burton-miller"
    if value not in {"ordinary", "burton-miller"}:
        raise ValueError("NSBEM_FORMULATION must be 'ordinary' or 'burton-miller'")
    return value


def build_validation_executable() -> None:
    """Compile the Fortran validation executable if needed."""

    subprocess.run(
        ["make", "-C", str(FORTRAN_DIR), "validate-rigid-sphere"],
        check=True,
    )


def run_validation_executable(mesh_file: Path, output_file: Path, config: ValidationConfig) -> None:
    """Run the Fortran solver for the rigid Neumann sphere case."""

    executable = FORTRAN_DIR / "bin" / "validate_rigid_sphere"
    subprocess.run(
        [
            str(executable),
            str(mesh_file),
            str(output_file),
            f"{config.wavenumber:.17g}",
            f"{config.radius:.17g}",
            str(exterior_domain_normal_sign_for(config.normal_orientation)),
            "neumann",
            "1.0",
            "0.5",
            config.formulation,
        ],
        check=True,
    )


def build_sphere_mesh(radius: float, subdivisions: int, element_type: str) -> SurfaceMesh:
    """Build either a 3-node linear sphere mesh or a 6-node quadratic one."""

    if element_type == LINEAR_ELEMENT:
        return build_icosphere(radius=radius, subdivisions=subdivisions)
    if element_type == QUADRATIC_ELEMENT:
        return build_quadratic_icosphere(radius=radius, subdivisions=subdivisions)
    raise ValueError("element_type must be 'linear' or 'quadratic'")


def build_icosphere(radius: float, subdivisions: int) -> SurfaceMesh:
    """Build a triangular sphere mesh by subdividing an icosahedron.

    An icosahedron gives a reasonably uniform starting mesh on a sphere.  Each
    subdivision splits every triangle into four smaller triangles and projects
    the new nodes back to the sphere radius.
    """

    if radius <= 0.0:
        raise ValueError("radius must be positive")
    if subdivisions < 0:
        raise ValueError("subdivisions must be non-negative")

    # The golden ratio defines the 12 vertices of a regular icosahedron.
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    raw_vertices = np.array(
        [
            (-1.0, phi, 0.0),
            (1.0, phi, 0.0),
            (-1.0, -phi, 0.0),
            (1.0, -phi, 0.0),
            (0.0, -1.0, phi),
            (0.0, 1.0, phi),
            (0.0, -1.0, -phi),
            (0.0, 1.0, -phi),
            (phi, 0.0, -1.0),
            (phi, 0.0, 1.0),
            (-phi, 0.0, -1.0),
            (-phi, 0.0, 1.0),
        ],
        dtype=float,
    )
    vertices = [_project_to_radius(vertex, radius) for vertex in raw_vertices]
    faces = [
        (0, 11, 5),
        (0, 5, 1),
        (0, 1, 7),
        (0, 7, 10),
        (0, 10, 11),
        (1, 5, 9),
        (5, 11, 4),
        (11, 10, 2),
        (10, 7, 6),
        (7, 1, 8),
        (3, 9, 4),
        (3, 4, 2),
        (3, 2, 6),
        (3, 6, 8),
        (3, 8, 9),
        (4, 9, 5),
        (2, 4, 11),
        (6, 2, 10),
        (8, 6, 7),
        (9, 8, 1),
    ]

    for _ in range(subdivisions):
        # Each old triangle becomes four new triangles.  A cache is used so that
        # two neighboring triangles share the same midpoint node on their common
        # edge instead of creating duplicate nodes.
        midpoint_cache: dict[tuple[int, int], int] = {}
        refined_faces: list[tuple[int, int, int]] = []
        for face in faces:
            a, b, c = face
            ab = _midpoint_index(vertices, midpoint_cache, a, b, radius)
            bc = _midpoint_index(vertices, midpoint_cache, b, c, radius)
            ca = _midpoint_index(vertices, midpoint_cache, c, a, radius)
            refined_faces.extend(
                [
                    (a, ab, ca),
                    (b, bc, ab),
                    (c, ca, bc),
                    (ab, bc, ca),
                ]
            )
        faces = refined_faces

    # The Fortran geometry routines assume a consistent surface orientation.
    oriented_faces = [_orient_outward(vertices, face) for face in faces]

    nodes = tuple(
        MeshNode(i + 1, float(vertex[0]), float(vertex[1]), float(vertex[2]), 1)
        for i, vertex in enumerate(vertices)
    )
    elements = tuple(
        MeshElement(i + 1, 1, tuple(index + 1 for index in face))
        for i, face in enumerate(oriented_faces)
    )
    return SurfaceMesh(nodes=nodes, elements=elements, nodes_per_element=3, particle_count=1)


def build_quadratic_icosphere(radius: float, subdivisions: int) -> SurfaceMesh:
    """Build the same sphere as build_icosphere, but with 6-node curved triangles.

    The first three nodes of each element are the triangle corners.  Nodes 4, 5,
    and 6 are the midside nodes on edges 1-2, 2-3, and 3-1.  Each midside node is
    projected back to the exact sphere radius, so the quadratic element contains
    genuine curved-surface information instead of only subdividing a flat face.
    """
    linear_mesh = build_icosphere(radius=radius, subdivisions=subdivisions)
    vertices = [
        np.array((node.x, node.y, node.z), dtype=float)
        for node in linear_mesh.nodes
    ]
    midpoint_cache: dict[tuple[int, int], int] = {}
    elements: list[MeshElement] = []

    for element in linear_mesh.elements:
        n1, n2, n3 = element.node_ids
        n4 = _quadratic_midside_node(vertices, midpoint_cache, n1, n2, radius)
        n5 = _quadratic_midside_node(vertices, midpoint_cache, n2, n3, radius)
        n6 = _quadratic_midside_node(vertices, midpoint_cache, n3, n1, radius)
        elements.append(
            MeshElement(
                element.element_id,
                element.particle_id,
                (n1, n2, n3, n4, n5, n6),
            )
        )

    nodes = tuple(
        MeshNode(i + 1, float(vertex[0]), float(vertex[1]), float(vertex[2]), 1)
        for i, vertex in enumerate(vertices)
    )
    return SurfaceMesh(
        nodes=nodes,
        elements=tuple(elements),
        nodes_per_element=6,
        particle_count=linear_mesh.particle_count,
    )


def _quadratic_midside_node(
    vertices: list[np.ndarray],
    midpoint_cache: dict[tuple[int, int], int],
    first_node_id: int,
    second_node_id: int,
    radius: float,
) -> int:
    """Create or reuse a quadratic midside node on one mesh edge.

    The input node numbers are one-based because they come from SurfaceMesh.
    Internally, Python list indices are zero-based, so we subtract 1.
    """

    first = first_node_id - 1
    second = second_node_id - 1
    key = tuple(sorted((first, second)))
    if key in midpoint_cache:
        return midpoint_cache[key] + 1

    midpoint = _project_to_radius(0.5 * (vertices[first] + vertices[second]), radius)
    vertices.append(midpoint)
    index = len(vertices) - 1
    midpoint_cache[key] = index
    return index + 1


def _project_to_radius(point: np.ndarray, radius: float) -> np.ndarray:
    """Project a point onto the sphere with the requested radius."""

    return radius * point / np.linalg.norm(point)


def _midpoint_index(
    vertices: list[np.ndarray],
    midpoint_cache: dict[tuple[int, int], int],
    first: int,
    second: int,
    radius: float,
) -> int:
    """Create or reuse a midpoint node during linear mesh refinement."""

    key = tuple(sorted((first, second)))
    if key in midpoint_cache:
        return midpoint_cache[key]

    midpoint = _project_to_radius(0.5 * (vertices[first] + vertices[second]), radius)
    vertices.append(midpoint)
    index = len(vertices) - 1
    midpoint_cache[key] = index
    return index


def _orient_outward(vertices: list[np.ndarray], face: tuple[int, int, int]) -> tuple[int, int, int]:
    """Return a triangle ordering whose normal points away from the origin."""

    a, b, c = face
    normal = np.cross(vertices[b] - vertices[a], vertices[c] - vertices[a])
    centroid = (vertices[a] + vertices[b] + vertices[c]) / 3.0
    if float(np.dot(normal, centroid)) < 0.0:
        return (a, c, b)
    return face


def read_nsbem_csv(path: Path) -> dict[str, np.ndarray]:
    """Read the surface CSV written by the Fortran validation executable."""

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows:
        raise ValueError(f"NSBEM output is empty: {path}")

    data: dict[str, np.ndarray] = {}
    for key in rows[0]:
        if key == "node_id":
            data[key] = np.array([int(row[key]) for row in rows])
        else:
            data[key] = np.array([float(row[key]) for row in rows], dtype=float)

    # Recombine scalar columns into arrays that are easier to use in NumPy.
    data["points"] = np.column_stack([data["x"], data["y"], data["z"]])
    data["phi_sca"] = data["phi_sca_re"] + 1j * data["phi_sca_im"]
    data["phi_total"] = data["phi_total_re"] + 1j * data["phi_total_im"]
    data["dphi_inc"] = data["dphi_inc_re"] + 1j * data["dphi_inc_im"]
    data["dphi_sca"] = data["dphi_sca_re"] + 1j * data["dphi_sca_im"]
    data["dphi_total"] = data["dphi_total_re"] + 1j * data["dphi_total_im"]
    # For a sphere, the outward normal should align with the radius vector.  This
    # diagnostic catches wrong element orientation or wrong normal direction.
    normals = np.column_stack([data["nx"], data["ny"], data["nz"]])
    data["normals"] = normals
    data["normal_radial_alignment"] = np.sum(data["points"] * normals, axis=1) / (
        np.linalg.norm(data["points"], axis=1) * np.linalg.norm(normals, axis=1)
    )
    return data


def evaluate_analytical_solution(nsbem: dict[str, np.ndarray], config: ValidationConfig) -> dict[str, np.ndarray]:
    """Evaluate the analytical solution at the numerical mesh nodes."""

    n_terms = estimate_required_terms(config.wavenumber, config.radius, safety=30)
    phi_sca = rigid_sphere_scattered_field(
        nsbem["points"],
        k=config.wavenumber,
        radius=config.radius,
        n_terms=n_terms,
    )
    phi_total = rigid_sphere_total_field(
        nsbem["points"],
        k=config.wavenumber,
        radius=config.radius,
        n_terms=n_terms,
    )
    return {
        "node_id": nsbem["node_id"],
        "points": nsbem["points"],
        "theta": surface_theta(nsbem["points"]),
        "phi_sca": phi_sca,
        "phi_total": phi_total,
        "n_terms": np.array([n_terms]),
    }


def write_analytic_csv(path: Path, analytic: dict[str, np.ndarray]) -> None:
    """Write analytical nodal values for inspection and plotting."""

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
            ]
        )
        for node_id, point, theta, phi_sca, phi_total in zip(
            analytic["node_id"],
            analytic["points"],
            analytic["theta"],
            analytic["phi_sca"],
            analytic["phi_total"],
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
                ]
            )


def write_comparison_csv(
    path: Path,
    nsbem: dict[str, np.ndarray],
    analytic: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    """Write node-by-node numerical/analytical differences."""

    error_sca = nsbem["phi_sca"] - analytic["phi_sca"]
    error_total = nsbem["phi_total"] - analytic["phi_total"]
    result = {
        "node_id": nsbem["node_id"],
        "points": nsbem["points"],
        "theta": analytic["theta"],
        "nsbem_phi_sca": nsbem["phi_sca"],
        "analytic_phi_sca": analytic["phi_sca"],
        "error_sca": error_sca,
        "nsbem_phi_total": nsbem["phi_total"],
        "analytic_phi_total": analytic["phi_total"],
        "error_total": error_total,
        "dphi_total": nsbem["dphi_total"],
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
                "nsbem_phi_total_re",
                "nsbem_phi_total_im",
                "analytic_phi_total_re",
                "analytic_phi_total_im",
                "error_total_abs",
                "dphi_total_abs",
            ]
        )
        for i, node_id in enumerate(result["node_id"]):
            point = result["points"][i]
            writer.writerow(
                [
                    int(node_id),
                    point[0],
                    point[1],
                    point[2],
                    result["theta"][i],
                    result["nsbem_phi_sca"][i].real,
                    result["nsbem_phi_sca"][i].imag,
                    result["analytic_phi_sca"][i].real,
                    result["analytic_phi_sca"][i].imag,
                    abs(result["error_sca"][i]),
                    result["nsbem_phi_total"][i].real,
                    result["nsbem_phi_total"][i].imag,
                    result["analytic_phi_total"][i].real,
                    result["analytic_phi_total"][i].imag,
                    abs(result["error_total"][i]),
                    abs(result["dphi_total"][i]),
                ]
            )
    return result


def surface_theta(points: np.ndarray) -> np.ndarray:
    """Return polar angle theta measured from the incident x direction."""

    radius = np.linalg.norm(points, axis=1)
    mu = np.clip(points[:, 0] / radius, -1.0, 1.0)
    return np.arccos(mu)


def write_plots(comparison: dict[str, np.ndarray]) -> None:
    """Save simple comparison plots for the sphere surface values."""

    order = np.argsort(comparison["theta"])
    theta = comparison["theta"][order]

    plot_component(
        theta,
        comparison["nsbem_phi_sca"][order].real,
        comparison["analytic_phi_sca"][order].real,
        "Scattered potential real part",
        "Re(phi_sca)",
        RESULTS_DIR / "surface_scattered_real.png",
    )
    plot_component(
        theta,
        comparison["nsbem_phi_sca"][order].imag,
        comparison["analytic_phi_sca"][order].imag,
        "Scattered potential imaginary part",
        "Im(phi_sca)",
        RESULTS_DIR / "surface_scattered_imag.png",
    )
    plt.figure(figsize=(8, 5))
    plt.semilogy(theta, np.abs(comparison["error_sca"][order]), "o", markersize=3)
    plt.xlabel("theta [rad]")
    plt.ylabel("|NSBEM - analytical|")
    plt.title("Scattered potential absolute error")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "surface_error_abs.png", dpi=200)
    plt.close()


def plot_component(
    theta: np.ndarray,
    nsbem_values: np.ndarray,
    analytic_values: np.ndarray,
    title: str,
    ylabel: str,
    path: Path,
) -> None:
    """Plot one scalar component of the numerical and analytical fields."""

    plt.figure(figsize=(8, 5))
    plt.plot(theta, analytic_values, "-", label="Analytical", linewidth=1.5)
    plt.plot(theta, nsbem_values, "o", label="NSBEM", markersize=3)
    plt.xlabel("theta [rad]")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


def write_summary(
    path: Path,
    config: ValidationConfig,
    mesh: SurfaceMesh,
    comparison: dict[str, np.ndarray],
) -> None:
    """Write a compact validation report."""

    err_sca = np.abs(comparison["error_sca"])
    err_total = np.abs(comparison["error_total"])
    ref_sca = np.abs(comparison["analytic_phi_sca"])
    ref_total = np.abs(comparison["analytic_phi_total"])
    dphi_total = np.abs(comparison["dphi_total"])
    normal_alignment = comparison["normal_radial_alignment"]

    max_err_sca = float(np.max(err_sca))
    rms_err_sca = float(np.sqrt(np.mean(err_sca**2)))
    rel_l2_sca = float(np.linalg.norm(err_sca) / max(np.linalg.norm(ref_sca), 1.0e-14))
    max_err_total = float(np.max(err_total))
    rel_l2_total = float(np.linalg.norm(err_total) / max(np.linalg.norm(ref_total), 1.0e-14))
    max_bc_residual = float(np.max(dphi_total))
    # These convention diagnostics are useful while developing BEM kernels.  A
    # wrong sign convention often shows up as agreement with -analytic or a
    # complex conjugate instead of the intended analytical field.
    rel_l2_minus = relative_l2(comparison["nsbem_phi_sca"], -comparison["analytic_phi_sca"])
    rel_l2_conjugate = relative_l2(comparison["nsbem_phi_sca"], np.conj(comparison["analytic_phi_sca"]))
    rel_l2_minus_conjugate = relative_l2(
        comparison["nsbem_phi_sca"],
        -np.conj(comparison["analytic_phi_sca"]),
    )
    alpha = np.vdot(comparison["analytic_phi_sca"], comparison["nsbem_phi_sca"]) / np.vdot(
        comparison["analytic_phi_sca"],
        comparison["analytic_phi_sca"],
    )
    best_scaled_rel_l2 = relative_l2(
        comparison["nsbem_phi_sca"],
        alpha * comparison["analytic_phi_sca"],
    )
    validation_tolerance = 5.0e-2
    validation_status = "PASS" if rel_l2_sca <= validation_tolerance else "FAIL"

    text = "\n".join(
        [
            "Rigid sphere BRIEF-Acoustics validation",
            "=================================",
            f"validation status = {validation_status}",
            f"validation tolerance, relative L2 phi_sca = {validation_tolerance:.8e}",
            f"radius = {config.radius:.8g}",
            f"wavenumber = {config.wavenumber:.8g}",
            f"formulation = {config.formulation}",
            f"mesh normal orientation = {config.normal_orientation}",
            f"mesh element type = {mesh.element_type}",
            f"mesh nodes per element = {mesh.nodes_per_element}",
            f"mesh nodes = {len(mesh.nodes)}",
            f"mesh elements = {len(mesh.elements)}",
            f"normal radial alignment min = {float(np.min(normal_alignment)):.8e}",
            f"normal radial alignment max = {float(np.max(normal_alignment)):.8e}",
            f"max |phi_sca error| = {max_err_sca:.8e}",
            f"rms |phi_sca error| = {rms_err_sca:.8e}",
            f"relative L2 phi_sca error = {rel_l2_sca:.8e}",
            f"max |phi_total error| = {max_err_total:.8e}",
            f"relative L2 phi_total error = {rel_l2_total:.8e}",
            f"max |dphi_total/dn numerical residual| = {max_bc_residual:.8e}",
            "",
            "Convention diagnostics",
            "----------------------",
            f"relative L2 against -analytic phi_sca = {rel_l2_minus:.8e}",
            f"relative L2 against conj(analytic) phi_sca = {rel_l2_conjugate:.8e}",
            f"relative L2 against -conj(analytic) phi_sca = {rel_l2_minus_conjugate:.8e}",
            f"best complex scale alpha = {alpha.real:.8e} {alpha.imag:+.8e}j",
            f"relative L2 after best scale = {best_scaled_rel_l2:.8e}",
            "",
        ]
    )
    path.write_text(text, encoding="utf-8")


def relative_l2(values: np.ndarray, reference: np.ndarray) -> float:
    """Return ||values - reference||_2 / ||reference||_2."""

    return float(np.linalg.norm(values - reference) / max(np.linalg.norm(reference), 1.0e-14))


if __name__ == "__main__":
    raise SystemExit(main())
