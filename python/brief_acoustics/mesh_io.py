"""Mesh data structures, validation, and text I/O for the Fortran solver.

Python owns the public mesh contract.  It checks closed-manifold topology,
orientation, quadratic midside-node consistency, and basic geometric quality
before the dense Fortran assembly begins.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path


LINEAR_ELEMENT = "linear"
QUADRATIC_ELEMENT = "quadratic"

# Keep the public choice in one place so case parsing, validation scripts, and
# future interfaces such as a GUI cannot drift to different defaults.
DEFAULT_ELEMENT_TYPE = QUADRATIC_ELEMENT
OUTWARD_FROM_SOLID = "outward-from-solid"

_NODES_PER_ELEMENT = {
    LINEAR_ELEMENT: 3,
    QUADRATIC_ELEMENT: 6,
}
DEFAULT_NODES_PER_ELEMENT = _NODES_PER_ELEMENT[DEFAULT_ELEMENT_TYPE]


def nodes_per_element_for(element_type: str) -> int:
    """Translate a public element label into the mesh's integer node count."""

    try:
        return _NODES_PER_ELEMENT[element_type]
    except KeyError as error:
        raise ValueError("element_type must be 'linear' or 'quadratic'.") from error


def element_type_for(nodes_per_element: int) -> str:
    """Return the public element label for a supported triangle node count."""

    for element_type, node_count in _NODES_PER_ELEMENT.items():
        if nodes_per_element == node_count:
            return element_type
    raise ValueError("nodes_per_element must be 3 or 6.")


def exterior_domain_normal_sign_for(normal_orientation: str) -> int:
    """Map the stored mesh normal to the exterior-domain normal.

    The returned sign ``sigma`` satisfies
    ``n_mesh = sigma * n_exterior_domain``.  It is an orientation convention,
    not a solid-angle value or numerical boundary coefficient.
    """

    if normal_orientation != OUTWARD_FROM_SOLID:
        raise ValueError(
            "normal_orientation must be 'outward-from-solid' for the public exterior solver."
        )
    return -1


@dataclass(frozen=True)
class MeshNode:
    """One mesh point on a boundary surface.

    The Fortran solver uses one-based node numbering, so node_id starts at 1,
    not 0.  particle_id identifies which disconnected body or material boundary
    the node belongs to.  The public examples use one particle only.
    """

    node_id: int
    x: float
    y: float
    z: float
    particle_id: int = 1


@dataclass(frozen=True)
class MeshElement:
    """One triangular surface element.

    node_ids stores the node numbers that form this element.  For a linear
    triangle it has 3 entries.  For a quadratic triangle it has 6 entries:

        (corner 1, corner 2, corner 3, midside 1-2, midside 2-3, midside 3-1)

    This ordering matches the shape functions used in the modern Fortran code.
    """

    element_id: int
    particle_id: int
    node_ids: tuple[int, ...]


@dataclass(frozen=True)
class SurfaceMesh:
    """Complete surface mesh passed from Python to Fortran.

    This is intentionally small: it stores what the current public release needs.
    Python is responsible for preparing or importing the mesh; Fortran is
    responsible for topology, geometry, integration, and solving.
    """

    nodes: tuple[MeshNode, ...]
    elements: tuple[MeshElement, ...]
    nodes_per_element: int = DEFAULT_NODES_PER_ELEMENT
    particle_count: int = 1

    @property
    def element_type(self) -> str:
        """Readable interpolation type derived from the stored node count."""

        return element_type_for(self.nodes_per_element)

    def validate(self) -> None:
        """Check the mesh before writing it to disk.

        The Fortran geometry routines assume that each particle is one closed,
        consistently oriented triangular surface.  These checks enforce that
        public contract before the numerically expensive assembly begins.  In
        addition to ids and counts, they catch common topology and ordering
        errors that would otherwise produce plausible-looking wrong normals.
        """

        if self.nodes_per_element not in {3, 6}:
            raise ValueError("SurfaceMesh.nodes_per_element must be 3 or 6.")
        if self.particle_count < 1:
            raise ValueError("SurfaceMesh.particle_count must be positive.")
        if not self.nodes:
            raise ValueError("SurfaceMesh must contain at least one node.")
        if not self.elements:
            raise ValueError("SurfaceMesh must contain at least one element.")

        node_ids = [node.node_id for node in self.nodes]
        expected_node_ids = list(range(1, len(self.nodes) + 1))
        if node_ids != expected_node_ids:
            raise ValueError("Mesh node ids must be sequential and one-based.")

        element_ids = [element.element_id for element in self.elements]
        expected_element_ids = list(range(1, len(self.elements) + 1))
        if element_ids != expected_element_ids:
            raise ValueError("Mesh element ids must be sequential and one-based.")

        node_id_set = set(node_ids)
        node_by_id = {node.node_id: node for node in self.nodes}
        for node in self.nodes:
            if not 1 <= node.particle_id <= self.particle_count:
                raise ValueError(f"Node {node.node_id} has invalid particle_id.")
            if not all(math.isfinite(value) for value in (node.x, node.y, node.z)):
                raise ValueError(f"Node {node.node_id} has a non-finite coordinate.")

        referenced_nodes: set[int] = set()
        face_keys: set[tuple[int, int, int, int]] = set()
        particle_elements: dict[int, set[int]] = {
            particle_id: set() for particle_id in range(1, self.particle_count + 1)
        }
        edge_uses: dict[
            tuple[int, int, int],
            list[tuple[int, int, int, int | None]],
        ] = {}
        for element in self.elements:
            if not 1 <= element.particle_id <= self.particle_count:
                raise ValueError(f"Element {element.element_id} has invalid particle_id.")
            if len(element.node_ids) != self.nodes_per_element:
                raise ValueError(f"Element {element.element_id} has wrong node count.")
            if any(node_id not in node_id_set for node_id in element.node_ids):
                raise ValueError(f"Element {element.element_id} references an unknown node.")
            if len(set(element.node_ids)) != len(element.node_ids):
                raise ValueError(f"Element {element.element_id} repeats a node id.")
            if any(
                node_by_id[node_id].particle_id != element.particle_id
                for node_id in element.node_ids
            ):
                raise ValueError(
                    f"Element {element.element_id} references a node owned by another particle."
                )

            referenced_nodes.update(element.node_ids)
            particle_elements[element.particle_id].add(element.element_id)
            corner_ids = element.node_ids[:3]
            face_key = (element.particle_id, *sorted(corner_ids))
            if face_key in face_keys:
                raise ValueError(f"Element {element.element_id} duplicates another triangle.")
            face_keys.add(face_key)

            _check_triangle_is_not_degenerate(element, node_by_id)

            midside_ids: tuple[int | None, int | None, int | None]
            if self.nodes_per_element == 6:
                midside_ids = (
                    element.node_ids[3],
                    element.node_ids[4],
                    element.node_ids[5],
                )
            else:
                midside_ids = (None, None, None)

            directed_edges = (
                (corner_ids[0], corner_ids[1], midside_ids[0]),
                (corner_ids[1], corner_ids[2], midside_ids[1]),
                (corner_ids[2], corner_ids[0], midside_ids[2]),
            )
            for first, second, midside in directed_edges:
                edge_key = (element.particle_id, min(first, second), max(first, second))
                edge_uses.setdefault(edge_key, []).append(
                    (element.element_id, first, second, midside)
                )

        unreferenced_nodes = node_id_set - referenced_nodes
        if unreferenced_nodes:
            first_unused = min(unreferenced_nodes)
            raise ValueError(f"Node {first_unused} is not referenced by any element.")

        element_neighbors: dict[int, set[int]] = {
            element.element_id: set() for element in self.elements
        }
        for edge_key, uses in edge_uses.items():
            particle_id, first_node, second_node = edge_key
            if len(uses) != 2:
                raise ValueError(
                    "Surface is not closed-manifold: edge "
                    f"({first_node}, {second_node}) on particle {particle_id} "
                    f"is used by {len(uses)} elements."
                )

            first_use, second_use = uses
            if first_use[1] != second_use[2] or first_use[2] != second_use[1]:
                raise ValueError(
                    "Adjacent elements have inconsistent orientation across edge "
                    f"({first_node}, {second_node}) on particle {particle_id}."
                )
            if self.nodes_per_element == 6 and first_use[3] != second_use[3]:
                raise ValueError(
                    "Quadratic elements disagree on the midside node for edge "
                    f"({first_node}, {second_node}) on particle {particle_id}."
                )

            first_element = first_use[0]
            second_element = second_use[0]
            element_neighbors[first_element].add(second_element)
            element_neighbors[second_element].add(first_element)

        for particle_id, element_set in particle_elements.items():
            if not element_set:
                raise ValueError(f"Particle {particle_id} has no elements.")
            connected = _connected_element_ids(min(element_set), element_neighbors)
            if connected != element_set:
                raise ValueError(
                    f"Particle {particle_id} contains disconnected surface components."
                )


def _check_triangle_is_not_degenerate(
    element: MeshElement,
    node_by_id: dict[int, MeshNode],
) -> None:
    """Reject a corner triangle with zero or numerically negligible area."""

    first, second, third = (
        node_by_id[node_id] for node_id in element.node_ids[:3]
    )
    edge_12 = (second.x - first.x, second.y - first.y, second.z - first.z)
    edge_13 = (third.x - first.x, third.y - first.y, third.z - first.z)
    cross = (
        edge_12[1] * edge_13[2] - edge_12[2] * edge_13[1],
        edge_12[2] * edge_13[0] - edge_12[0] * edge_13[2],
        edge_12[0] * edge_13[1] - edge_12[1] * edge_13[0],
    )
    doubled_area = math.sqrt(sum(value * value for value in cross))
    edge_23 = (third.x - second.x, third.y - second.y, third.z - second.z)
    maximum_edge_squared = max(
        sum(value * value for value in edge_12),
        sum(value * value for value in edge_13),
        sum(value * value for value in edge_23),
    )
    if maximum_edge_squared == 0.0 or doubled_area <= 1.0e-12 * maximum_edge_squared:
        raise ValueError(f"Element {element.element_id} is degenerate or nearly collinear.")


def _connected_element_ids(
    start_element: int,
    neighbors: dict[int, set[int]],
) -> set[int]:
    """Return all elements reachable through shared corner edges."""

    visited: set[int] = set()
    pending = [start_element]
    while pending:
        element_id = pending.pop()
        if element_id in visited:
            continue
        visited.add(element_id)
        pending.extend(neighbors[element_id] - visited)
    return visited


def write_fortran_mesh(mesh: SurfaceMesh, path: Path) -> None:
    """Write the simple text mesh consumed by modern `Geom_ReadMesh.f90`.

    File layout
    -----------
    First line:

        node_count element_count nodes_per_element particle_count

    Then one line per node:

        node_id x y z particle_id

    Then one line per element:

        element_id particle_id node_1 node_2 node_3 [node_4 node_5 node_6]

    The text format is deliberately plain so users can inspect it by eye and so
    that external mesh generators can produce it without a special dependency.
    """

    mesh.validate()
    lines: list[str] = [
        f"{len(mesh.nodes)} {len(mesh.elements)} {mesh.nodes_per_element} {mesh.particle_count}\n"
    ]
    for node in mesh.nodes:
        lines.append(f"{node.node_id} {node.x:.17g} {node.y:.17g} {node.z:.17g} {node.particle_id}\n")
    for element in mesh.elements:
        joined_nodes = " ".join(str(node_id) for node_id in element.node_ids)
        lines.append(f"{element.element_id} {element.particle_id} {joined_nodes}\n")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines), encoding="utf-8")


def read_fortran_mesh(path: Path) -> SurfaceMesh:
    """Read the same simple text mesh format that `write_fortran_mesh` writes."""

    # Blank lines and comment lines beginning with "#" are ignored.  This lets
    # users annotate small hand-written mesh files without confusing the parser.
    records = [
        line.split()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not records:
        raise ValueError(f"Mesh file is empty: {path}")

    node_count, element_count, nodes_per_element, particle_count = map(int, records[0])
    expected_count = 1 + node_count + element_count
    if len(records) != expected_count:
        raise ValueError(f"Mesh file has {len(records)} records, expected {expected_count}.")

    nodes: list[MeshNode] = []
    for fields in records[1 : 1 + node_count]:
        node_id = int(fields[0])
        nodes.append(MeshNode(node_id, float(fields[1]), float(fields[2]), float(fields[3]), int(fields[4])))

    elements: list[MeshElement] = []
    for fields in records[1 + node_count :]:
        element_id = int(fields[0])
        particle_id = int(fields[1])
        node_ids = tuple(int(value) for value in fields[2:])
        elements.append(MeshElement(element_id, particle_id, node_ids))

    mesh = SurfaceMesh(tuple(nodes), tuple(elements), nodes_per_element, particle_count)
    mesh.validate()
    return mesh


def tetrahedron_mesh() -> SurfaceMesh:
    """Return a tiny closed linear mesh useful for smoke tests.

    The tetrahedron is not a high-quality acoustic test mesh.  Its purpose is
    only to exercise mesh writing/reading and make sure the data path works.
    """

    nodes = (
        MeshNode(1, 0.0, 0.0, 0.0),
        MeshNode(2, 1.0, 0.0, 0.0),
        MeshNode(3, 0.0, 1.0, 0.0),
        MeshNode(4, 0.0, 0.0, 1.0),
    )
    elements = (
        MeshElement(1, 1, (1, 3, 2)),
        MeshElement(2, 1, (1, 2, 4)),
        MeshElement(3, 1, (2, 3, 4)),
        MeshElement(4, 1, (3, 1, 4)),
    )
    return SurfaceMesh(nodes, elements, nodes_per_element=3)
