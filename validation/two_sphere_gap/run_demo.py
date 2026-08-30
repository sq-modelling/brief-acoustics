"""Two-sphere small-gap NSBEM demonstration.

This demo is designed to show a strength of the non-singular BEM formulation:
two boundaries can be placed very close together while the solved surface field
remains usable for post-processing.

The workflow is:

1. Generate two separate spherical surface meshes in Python.
2. Prescribe opposite Neumann values on the two spheres.
3. Run the modern Fortran NSBEM solver.
4. Evaluate the field along the narrow gap using the ordinary field-point
   boundary integral expression, Eq. (4.1) of the 2015 RSOS BRIEF acoustics
   paper.

The field evaluation here intentionally uses Eq. (4.1), as requested, because it
is enough for this first public demonstration.  Eq. (4.2) from the paper is the
more robust near-boundary field evaluation and should be added later for a more
aggressive near-surface comparison.
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
import importlib.util
import math
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
    MeshElement,
    MeshNode,
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
read_element_type_from_environment = _sphere_helpers.read_element_type_from_environment


@dataclass(frozen=True)
class DemoConfig:
    """Input parameters for the two-sphere small-gap demonstration."""

    radius: float = 1.0
    gap: float = 0.05
    wavenumber: float = math.pi / 2.0
    subdivisions: int = 2
    element_type: str = DEFAULT_ELEMENT_TYPE
    field_point_count: int = 101
    normal_orientation: str = OUTWARD_FROM_SOLID
    prescribed_normal_derivative_1: complex = 1.0 + 0.0j
    prescribed_normal_derivative_2: complex = -1.0 + 0.0j


def main() -> int:
    """Run the two-sphere demo and write data/figures."""

    config = DemoConfig(
        radius=_env_float("NSBEM_RADIUS", DemoConfig.radius),
        gap=_env_float("NSBEM_GAP", DemoConfig.gap),
        subdivisions=_env_int("NSBEM_SUBDIVISIONS", DemoConfig.subdivisions),
        element_type=read_element_type_from_environment(),
        field_point_count=_env_int("NSBEM_FIELD_POINTS", DemoConfig.field_point_count),
    )
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    mesh = build_two_sphere_mesh(config)
    mesh_file = RESULTS_DIR / "two_spheres_mesh.dat"
    surface_file = RESULTS_DIR / "surface_solution.csv"
    gap_file = RESULTS_DIR / "gap_line_field.csv"
    summary_file = RESULTS_DIR / "summary.txt"

    write_fortran_mesh(mesh, mesh_file)
    build_demo_executable()
    run_demo_executable(mesh_file, surface_file, config)

    surface = read_surface_solution(surface_file)
    gap = evaluate_gap_line(mesh, surface, config)
    write_gap_csv(gap_file, gap)
    write_plots(gap, config)
    write_summary(summary_file, config, mesh, gap)

    print(f"Wrote two-sphere gap demo results to {RESULTS_DIR}")
    print(summary_file.read_text(encoding="utf-8"))
    return 0


def _env_int(name: str, default: int) -> int:
    """Read one optional integer environment override."""

    value = os.environ.get(name)
    if value is None:
        return default
    return int(value)


def _env_float(name: str, default: float) -> float:
    """Read one optional floating-point environment override."""

    value = os.environ.get(name)
    if value is None:
        return default
    return float(value)


def build_demo_executable() -> None:
    """Compile the Fortran two-sphere demo executable if needed."""

    subprocess.run(["make", "-C", str(FORTRAN_DIR), "demo-two-spheres"], check=True)


def run_demo_executable(mesh_file: Path, output_file: Path, config: DemoConfig) -> None:
    """Run the Fortran two-sphere surface solve."""

    executable = FORTRAN_DIR / "bin" / "demo_two_spheres"
    subprocess.run(
        [
            str(executable),
            str(mesh_file),
            str(output_file),
            f"{config.wavenumber:.17g}",
            f"{config.radius:.17g}",
            str(exterior_domain_normal_sign_for(config.normal_orientation)),
            f"{config.prescribed_normal_derivative_1.real:.17g}",
            f"{config.prescribed_normal_derivative_1.imag:.17g}",
            f"{config.prescribed_normal_derivative_2.real:.17g}",
            f"{config.prescribed_normal_derivative_2.imag:.17g}",
        ],
        check=True,
    )


def build_two_sphere_mesh(config: DemoConfig) -> SurfaceMesh:
    """Create two equal spheres separated by a small gap.

    Sphere 1 is centered at x = -(radius + gap/2).  Sphere 2 is centered at
    x = +(radius + gap/2).  Therefore their closest surface points are at
    x = -gap/2 and x = +gap/2.
    """

    if config.radius <= 0.0:
        raise ValueError("radius must be positive")
    if config.gap <= 0.0:
        raise ValueError("gap must be positive")

    base = build_sphere_mesh(
        radius=config.radius,
        subdivisions=config.subdivisions,
        element_type=config.element_type,
    )
    center_offset = config.radius + 0.5 * config.gap
    left = translate_particle(base, np.array([-center_offset, 0.0, 0.0]), particle_id=1)
    right = translate_particle(
        base,
        np.array([center_offset, 0.0, 0.0]),
        particle_id=2,
        node_offset=len(left.nodes),
        element_offset=len(left.elements),
    )

    return SurfaceMesh(
        nodes=left.nodes + right.nodes,
        elements=left.elements + right.elements,
        nodes_per_element=base.nodes_per_element,
        particle_count=2,
    )


def translate_particle(
    mesh: SurfaceMesh,
    center: np.ndarray,
    particle_id: int,
    node_offset: int = 0,
    element_offset: int = 0,
) -> SurfaceMesh:
    """Translate one sphere mesh and renumber it as one particle."""

    nodes = tuple(
        MeshNode(
            node.node_id + node_offset,
            node.x + float(center[0]),
            node.y + float(center[1]),
            node.z + float(center[2]),
            particle_id,
        )
        for node in mesh.nodes
    )
    elements = tuple(
        MeshElement(
            element.element_id + element_offset,
            particle_id,
            tuple(node_id + node_offset for node_id in element.node_ids),
        )
        for element in mesh.elements
    )
    return SurfaceMesh(nodes, elements, mesh.nodes_per_element, particle_count=1)


def read_surface_solution(path: Path) -> dict[str, np.ndarray]:
    """Read nodal phi, dphi/dn, pressure, and normals from the Fortran CSV."""

    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"Surface solution is empty: {path}")

    data: dict[str, np.ndarray] = {}
    for key in rows[0]:
        if key in {"node_id", "particle_id"}:
            data[key] = np.array([int(row[key]) for row in rows])
        else:
            data[key] = np.array([float(row[key]) for row in rows], dtype=float)
    data["phi"] = data["phi_re"] + 1j * data["phi_im"]
    data["dphi_dn"] = data["dphi_dn_re"] + 1j * data["dphi_dn_im"]
    data["pressure"] = data["pressure_re"] + 1j * data["pressure_im"]
    return data


def evaluate_gap_line(
    mesh: SurfaceMesh,
    surface: dict[str, np.ndarray],
    config: DemoConfig,
) -> dict[str, np.ndarray]:
    """Evaluate potential and pressure along the line through the gap."""

    margin = min(0.1 * config.gap, 0.25 * config.gap)
    x = np.linspace(
        -0.5 * config.gap + margin,
        0.5 * config.gap - margin,
        config.field_point_count,
    )
    points = np.column_stack([x, np.zeros_like(x), np.zeros_like(x)])
    phi = np.array(
        [
            field_point_phi(point, mesh, surface, config.wavenumber, config.normal_orientation)
            for point in points
        ]
    )
    pressure = 1j * config.wavenumber * phi
    return {
        "points": points,
        "phi": phi,
        "pressure": pressure,
        "abs_pressure": np.abs(pressure),
    }


def field_point_phi(
    field_point: np.ndarray,
    mesh: SurfaceMesh,
    surface: dict[str, np.ndarray],
    wavenumber: float,
    normal_orientation: str = OUTWARD_FROM_SOLID,
) -> complex:
    """Evaluate Eq. (4.1) for one field point.

    With the Green-function convention G = exp(i*k*r)/r used in the Fortran
    kernels, Eq. (4.1) reads:

        4*pi*phi(xp) = -int phi(x) dG(xp,x)/dn dS
                       +int dphi(x)/dn G(xp,x) dS

    The element-node ordering supplies the stored mesh normal. Equation (4.1) uses
    the outward normal of the solution domain.  For an exterior problem with an
    outward-from-solid mesh these directions are opposite. The readable normal
    orientation is converted to the internal exterior-domain normal sign here.
    Both dG/dn and dphi/dn change sign together, giving one overall factor.

    The input field point must be in the exterior fluid, not on the boundary.
    """

    exterior_domain_normal_sign = exterior_domain_normal_sign_for(normal_orientation)

    total = 0.0 + 0.0j
    nodes_xyz = np.array([(node.x, node.y, node.z) for node in mesh.nodes], dtype=float)
    phi_nodes = surface["phi"]
    dphi_nodes = surface["dphi_dn"]

    for element in mesh.elements:
        node_indices = np.array(element.node_ids, dtype=int) - 1
        element_xyz = nodes_xyz[node_indices]
        element_phi = phi_nodes[node_indices]
        element_dphi = dphi_nodes[node_indices]

        for xi, eta, weight in triangle_quadrature_25():
            shape, dshape_dxi, dshape_deta = triangle_shape_functions(mesh.nodes_per_element, xi, eta)
            source_point = shape @ element_xyz
            tangent_xi = dshape_dxi @ element_xyz
            tangent_eta = dshape_deta @ element_xyz
            normal_cross = np.cross(tangent_xi, tangent_eta)
            jacobian = np.linalg.norm(normal_cross)
            if jacobian <= 0.0:
                raise ValueError("Degenerate element encountered during field evaluation")
            normal = normal_cross / jacobian

            phi_q = shape @ element_phi
            dphi_q = shape @ element_dphi
            source_minus_field = source_point - field_point
            distance = np.linalg.norm(source_minus_field)
            if distance <= 1.0e-14:
                raise ValueError("Field point lies on a boundary quadrature point")

            green = np.exp(1j * wavenumber * distance) / distance
            dgreen_dr = np.exp(1j * wavenumber * distance) * (
                1j * wavenumber * distance - 1.0
            ) / distance**2
            dgreen_dn = dgreen_dr * np.dot(source_minus_field, normal) / distance
            total += (-phi_q * dgreen_dn + dphi_q * green) * jacobian * weight

    return exterior_domain_normal_sign * total / (4.0 * math.pi)


def triangle_shape_functions(
    nodes_per_element: int,
    xi: float,
    eta: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return triangle shape functions and reference derivatives."""

    zeta = 1.0 - xi - eta
    if nodes_per_element == 3:
        shape = np.array([zeta, xi, eta], dtype=float)
        dshape_dxi = np.array([-1.0, 1.0, 0.0], dtype=float)
        dshape_deta = np.array([-1.0, 0.0, 1.0], dtype=float)
        return shape, dshape_dxi, dshape_deta
    if nodes_per_element == 6:
        shape = np.array(
            [
                zeta * (2.0 * zeta - 1.0),
                xi * (2.0 * xi - 1.0),
                eta * (2.0 * eta - 1.0),
                4.0 * zeta * xi,
                4.0 * xi * eta,
                4.0 * eta * zeta,
            ],
            dtype=float,
        )
        dshape_dxi = np.array(
            [
                1.0 - 4.0 * zeta,
                4.0 * xi - 1.0,
                0.0,
                4.0 * (zeta - xi),
                4.0 * eta,
                -4.0 * eta,
            ],
            dtype=float,
        )
        dshape_deta = np.array(
            [
                1.0 - 4.0 * zeta,
                0.0,
                4.0 * eta - 1.0,
                -4.0 * xi,
                4.0 * xi,
                4.0 * (zeta - eta),
            ],
            dtype=float,
        )
        return shape, dshape_dxi, dshape_deta
    raise ValueError("nodes_per_element must be 3 or 6")


def triangle_quadrature_25() -> tuple[tuple[float, float, float], ...]:
    """Return the same 25-point triangle rule used by the Fortran solver."""

    weight = (
        0.04087166457314298321405934249209,
        0.00667648440657478313778648919953,
        0.00667648440657478313778648919953,
        0.00667648440657478313778648919953,
        0.02297898180237236400689395481627,
        0.02297898180237236400689395481627,
        0.02297898180237236400689395481627,
        0.03195245319821202271644936688133,
        0.03195245319821202271644936688133,
        0.03195245319821202271644936688133,
        0.03195245319821202271644936688133,
        0.03195245319821202271644936688133,
        0.03195245319821202271644936688133,
        0.01709232408147971431434579202067,
        0.01709232408147971431434579202067,
        0.01709232408147971431434579202067,
        0.01709232408147971431434579202067,
        0.01709232408147971431434579202067,
        0.01709232408147971431434579202067,
        0.01264887885364419219452139534142,
        0.01264887885364419219452139534142,
        0.01264887885364419219452139534142,
        0.01264887885364419219452139534142,
        0.01264887885364419219452139534142,
        0.01264887885364419219452139534142,
    )
    xi = (
        0.33333333333333333333333333333333,
        0.03205537321694351293098458933649,
        0.93588925356611297413803082132702,
        0.03205537321694351293098458933649,
        0.14216110105656438509216210319096,
        0.71567779788687122981567579361808,
        0.14216110105656438509216210319096,
        0.32181299528883542122509756098605,
        0.53005411892734402827709567394569,
        0.14813288578382055049780676506826,
        0.53005411892734402827709567394569,
        0.14813288578382055049780676506826,
        0.32181299528883542122509756098605,
        0.02961988948872976763383626942604,
        0.60123332868345924545474289345869,
        0.36914678182781098691142083711527,
        0.60123332868345924545474289345869,
        0.36914678182781098691142083711527,
        0.02961988948872976763383626942604,
        0.02836766533993843925043575557813,
        0.80793060092287906507994990288174,
        0.16370173373718249566961434154013,
        0.80793060092287906507994990288174,
        0.16370173373718249566961434154013,
        0.02836766533993843925043575557813,
    )
    eta = (
        0.33333333333333333333333333333333,
        0.93588925356611297413803082132702,
        0.03205537321694351293098458933649,
        0.03205537321694351293098458933649,
        0.71567779788687122981567579361808,
        0.14216110105656438509216210319096,
        0.14216110105656438509216210319096,
        0.53005411892734402827709567394569,
        0.32181299528883542122509756098605,
        0.53005411892734402827709567394569,
        0.14813288578382055049780676506826,
        0.32181299528883542122509756098605,
        0.14813288578382055049780676506826,
        0.60123332868345924545474289345869,
        0.02961988948872976763383626942604,
        0.60123332868345924545474289345869,
        0.36914678182781098691142083711527,
        0.02961988948872976763383626942604,
        0.36914678182781098691142083711527,
        0.80793060092287906507994990288174,
        0.02836766533993843925043575557813,
        0.80793060092287906507994990288174,
        0.16370173373718249566961434154013,
        0.02836766533993843925043575557813,
        0.16370173373718249566961434154013,
    )
    return tuple(zip(xi, eta, weight))


def write_gap_csv(path: Path, gap: dict[str, np.ndarray]) -> None:
    """Write field values along the gap line."""

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["x", "y", "z", "phi_re", "phi_im", "pressure_re", "pressure_im", "abs_pressure"]
        )
        for point, phi, pressure, abs_pressure in zip(
            gap["points"],
            gap["phi"],
            gap["pressure"],
            gap["abs_pressure"],
        ):
            writer.writerow(
                [
                    point[0],
                    point[1],
                    point[2],
                    phi.real,
                    phi.imag,
                    pressure.real,
                    pressure.imag,
                    abs_pressure,
                ]
            )


def write_plots(gap: dict[str, np.ndarray], config: DemoConfig) -> None:
    """Save pressure and potential plots along the gap line."""

    x_over_gap = gap["points"][:, 0] / config.gap

    plt.figure(figsize=(8, 5))
    plt.plot(x_over_gap, gap["abs_pressure"], "-", linewidth=1.8)
    plt.xlabel("x / gap")
    plt.ylabel("|pressure|")
    plt.title("Two-sphere small-gap pressure magnitude")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "gap_pressure_abs.png", dpi=200)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(x_over_gap, gap["pressure"].real, "-", label="Re(pressure)", linewidth=1.8)
    plt.plot(x_over_gap, gap["pressure"].imag, "--", label="Im(pressure)", linewidth=1.8)
    plt.xlabel("x / gap")
    plt.ylabel("pressure")
    plt.title("Two-sphere small-gap complex pressure")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "gap_pressure_complex.png", dpi=200)
    plt.close()


def write_summary(
    path: Path,
    config: DemoConfig,
    mesh: SurfaceMesh,
    gap: dict[str, np.ndarray],
) -> None:
    """Write a compact summary of the demonstration run."""

    abs_pressure = gap["abs_pressure"]
    text = "\n".join(
        [
            "Two-sphere small-gap BRIEF-Acoustics demonstration",
            "============================================",
            f"radius = {config.radius:.8e}",
            f"gap = {config.gap:.8e}",
            f"gap / radius = {config.gap / config.radius:.8e}",
            f"wavenumber = {config.wavenumber:.8e}",
            f"k * radius = {config.wavenumber * config.radius:.8e}",
            f"mesh element type = {mesh.element_type}",
            f"mesh nodes per element = {mesh.nodes_per_element}",
            f"mesh normal orientation = {config.normal_orientation}",
            f"mesh nodes = {len(mesh.nodes)}",
            f"mesh elements = {len(mesh.elements)}",
            f"field points across gap = {len(gap['points'])}",
            "prescribed dphi/dn sphere 1 = "
            f"{config.prescribed_normal_derivative_1.real:.8e} "
            f"{config.prescribed_normal_derivative_1.imag:+.8e}j",
            "prescribed dphi/dn sphere 2 = "
            f"{config.prescribed_normal_derivative_2.real:.8e} "
            f"{config.prescribed_normal_derivative_2.imag:+.8e}j",
            f"min |pressure| on gap line = {float(np.min(abs_pressure)):.8e}",
            f"max |pressure| on gap line = {float(np.max(abs_pressure)):.8e}",
            "",
            "Note",
            "----",
            "Field values were evaluated with Eq. (4.1) of Sun et al. (2015).",
            "This is a workflow demonstration, not an accuracy or superiority claim.",
            "Eq. (4.2) is the preferred future upgrade for more aggressive",
            "near-boundary field evaluation.",
            "",
        ]
    )
    path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
