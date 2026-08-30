from __future__ import annotations

import math
from pathlib import Path
import tempfile
import unittest

from brief_acoustics.case_config import CaseConfigurationError, SolverConfig, load_case
from brief_acoustics.mesh_io import tetrahedron_mesh, write_fortran_mesh
from brief_acoustics.runner import CaseRunError, check_case


class CaseConfigTests(unittest.TestCase):
    """Exercise the public contract without running a dense solve."""

    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary_directory.name)
        write_fortran_mesh(tetrahedron_mesh(), self.directory / "mesh.dat")

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def test_frequency_hz_derives_omega_and_wavenumber(self) -> None:
        case_path = self._write_case(frequency_line="frequency_hz = 2.0")
        config, mesh = check_case(case_path)

        self.assertEqual(config.name, "test-case")
        self.assertEqual(config.mesh.element_type, "linear")
        self.assertEqual(config.mesh.normal_orientation, "outward-from-solid")
        self.assertEqual(len(mesh.nodes), 4)
        self.assertAlmostEqual(config.physics.angular_frequency, 4.0 * math.pi)
        self.assertAlmostEqual(config.physics.wavenumber, 2.0 * math.pi)

    def test_multiple_frequency_representations_are_rejected(self) -> None:
        case_path = self._write_case(
            frequency_line="frequency_hz = 2.0\nwavenumber = 1.0"
        )
        with self.assertRaisesRegex(CaseConfigurationError, "exactly one"):
            load_case(case_path)

    def test_omitted_element_type_defaults_to_quadratic(self) -> None:
        case_path = self._write_case(
            frequency_line="wavenumber = 1.0",
            element_type=None,
        )

        config = load_case(case_path)

        self.assertEqual(config.mesh.element_type, "quadratic")
        self.assertEqual(config.mesh.nodes_per_element, 6)

    def test_linear_mesh_requires_explicit_linear_element_type(self) -> None:
        case_path = self._write_case(
            frequency_line="wavenumber = 1.0",
            element_type=None,
        )

        with self.assertRaisesRegex(CaseRunError, "does not match"):
            check_case(case_path)

    def test_unknown_key_is_rejected(self) -> None:
        case_path = self._write_case(
            frequency_line="wavenumber = 1.0",
            extra_solver_line="formulaton = 'ordinary'",
        )
        with self.assertRaisesRegex(CaseConfigurationError, "formulaton"):
            load_case(case_path)

    def test_numeric_normal_sign_key_is_rejected(self) -> None:
        case_path = self._write_case(
            frequency_line="wavenumber = 1.0",
            extra_solver_line="normal_sign = -1",
        )
        with self.assertRaisesRegex(CaseConfigurationError, "normal_sign"):
            load_case(case_path)

    def test_ambiguous_plane_wave_amplitude_key_is_rejected(self) -> None:
        case_path = self._write_case(
            frequency_line="wavenumber = 1.0",
            incident='type = "plane-wave"\namplitude_real = 1.0',
        )
        with self.assertRaisesRegex(CaseConfigurationError, "amplitude_real"):
            load_case(case_path)

    def test_unsupported_normal_orientation_is_rejected(self) -> None:
        case_path = self._write_case(
            frequency_line="wavenumber = 1.0",
            normal_orientation="inward-to-solid",
        )
        with self.assertRaisesRegex(CaseConfigurationError, "outward-from-solid"):
            load_case(case_path)

    def test_declared_element_type_must_match_mesh(self) -> None:
        case_path = self._write_case(
            frequency_line="wavenumber = 1.0",
            element_type="quadratic",
        )
        with self.assertRaisesRegex(CaseRunError, "does not match"):
            check_case(case_path)

    def test_zero_robin_coefficients_are_rejected(self) -> None:
        case_path = self._write_case(
            frequency_line="wavenumber = 1.0",
            boundary="""type = "robin"
a_real = 0.0
b_real = 0.0
rhs_real = 1.0""",
        )
        with self.assertRaisesRegex(CaseConfigurationError, "cannot both be zero"):
            load_case(case_path)

    def test_burton_miller_coupling_length_has_low_and_high_frequency_scales(self) -> None:
        solver = SolverConfig(
            formulation="burton-miller",
            characteristic_length=1.0,
        )

        self.assertAlmostEqual(solver.burton_miller_coupling_length(0.25), 1.0)
        self.assertAlmostEqual(solver.burton_miller_coupling_length(2.0), 0.5)

    def _write_case(
        self,
        *,
        frequency_line: str,
        element_type: str | None = "linear",
        normal_orientation: str = "outward-from-solid",
        extra_solver_line: str = "",
        boundary: str = 'type = "neumann"\nvalue_real = 0.0',
        incident: str = 'type = "none"',
    ) -> Path:
        path = self.directory / "case.toml"
        element_type_line = (
            "" if element_type is None else f'element_type = "{element_type}"\n'
        )
        path.write_text(
            f"""[case]
name = "test-case"

[mesh]
path = "mesh.dat"
{element_type_line}normal_orientation = "{normal_orientation}"

[physics]
density = 1.0
sound_speed = 2.0
{frequency_line}

[solver]
formulation = "ordinary"
characteristic_length = 1.0
{extra_solver_line}

[boundary]
{boundary}

[incident]
{incident}

[output]
directory = "results"
""",
            encoding="utf-8",
        )
        return path


if __name__ == "__main__":
    unittest.main()
