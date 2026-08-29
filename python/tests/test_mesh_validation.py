from __future__ import annotations

import unittest

from brief_acoustics.mesh_io import MeshElement, MeshNode, SurfaceMesh, tetrahedron_mesh


class MeshValidationTests(unittest.TestCase):
    """Check failures that would corrupt normals or boundary topology."""

    def test_reference_tetrahedron_is_valid(self) -> None:
        tetrahedron_mesh().validate()

    def test_repeated_element_node_is_rejected(self) -> None:
        reference = tetrahedron_mesh()
        elements = list(reference.elements)
        elements[0] = MeshElement(1, 1, (1, 3, 1))
        mesh = SurfaceMesh(
            reference.nodes,
            tuple(elements),
            nodes_per_element=reference.nodes_per_element,
        )

        with self.assertRaisesRegex(ValueError, "repeats a node"):
            mesh.validate()

    def test_open_surface_edge_is_rejected(self) -> None:
        reference = tetrahedron_mesh()
        mesh = SurfaceMesh(
            reference.nodes,
            reference.elements[:-1],
            nodes_per_element=reference.nodes_per_element,
        )

        with self.assertRaisesRegex(ValueError, "not closed-manifold"):
            mesh.validate()

    def test_inconsistent_face_orientation_is_rejected(self) -> None:
        reference = tetrahedron_mesh()
        elements = list(reference.elements)
        elements[0] = MeshElement(1, 1, (1, 2, 3))
        mesh = SurfaceMesh(
            reference.nodes,
            tuple(elements),
            nodes_per_element=reference.nodes_per_element,
        )

        with self.assertRaisesRegex(ValueError, "inconsistent orientation"):
            mesh.validate()

    def test_cross_particle_node_reference_is_rejected(self) -> None:
        reference = tetrahedron_mesh()
        nodes = list(reference.nodes)
        fourth = nodes[3]
        nodes[3] = MeshNode(
            fourth.node_id,
            fourth.x,
            fourth.y,
            fourth.z,
            particle_id=2,
        )
        mesh = SurfaceMesh(
            tuple(nodes),
            reference.elements,
            nodes_per_element=3,
            particle_count=2,
        )

        with self.assertRaisesRegex(ValueError, "owned by another particle"):
            mesh.validate()

    def test_degenerate_corner_triangle_is_rejected(self) -> None:
        nodes = (
            MeshNode(1, 0.0, 0.0, 0.0),
            MeshNode(2, 1.0, 0.0, 0.0),
            MeshNode(3, 2.0, 0.0, 0.0),
            MeshNode(4, 0.0, 1.0, 0.0),
        )
        elements = (
            MeshElement(1, 1, (1, 2, 3)),
            MeshElement(2, 1, (1, 4, 2)),
            MeshElement(3, 1, (2, 4, 3)),
            MeshElement(4, 1, (3, 4, 1)),
        )
        mesh = SurfaceMesh(nodes, elements, nodes_per_element=3)

        with self.assertRaisesRegex(ValueError, "degenerate"):
            mesh.validate()


if __name__ == "__main__":
    unittest.main()
