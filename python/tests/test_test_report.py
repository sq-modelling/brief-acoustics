"""Check the standalone feedback tool without importing or running the solver."""

from contextlib import redirect_stderr, redirect_stdout
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch


SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "make_test_report.py"
SPEC = importlib.util.spec_from_file_location("make_test_report", SCRIPT)
report_tool = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(report_tool)


class TestReportTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary.name)
        self.addCleanup(self.temporary.cleanup)

    def test_missing_tools_and_package_do_not_block_report(self) -> None:
        with patch.object(report_tool, "command_output", return_value=None), patch.object(
            report_tool.metadata, "version",
            side_effect=report_tool.metadata.PackageNotFoundError,
        ):
            report = report_tool.build_report(self.directory, None)
        self.assertIn("Not detected in this Python environment", report)
        self.assertIn("Installation succeeded: Not confirmed", report)
        self.assertIn("Minimal example succeeded: Not confirmed", report)
        self.assertIn("No results directory supplied", report)
        self.assertNotIn(str(self.directory), report)

    def test_completed_summary_is_evidence_not_automatic_success(self) -> None:
        summary = {
            "status": "completed", "case": "/private/person/confidential-case",
            "mesh": {"nodes": 162, "elements": 80, "nodes_per_element": 6,
                     "element_type": "quadratic"},
            "outputs": {"solver_log": "/private/secret.log"},
        }
        (self.directory / "summary.json").write_text(json.dumps(summary))
        (self.directory / "solver.log").write_text("SECRET log text /private/person")
        with patch.object(report_tool, "environment_rows", return_value=[]):
            report = report_tool.build_report(self.directory, self.directory)
        self.assertIn("status = completed", report)
        self.assertIn("nodes_per_element: 6", report)
        self.assertIn("surface_solution.csv: missing or empty", report)
        self.assertIn("solver.log: present and non-empty", report)
        self.assertIn("Minimal example succeeded: Not confirmed", report)
        self.assertIn("not linked to the current commit", report)
        self.assertNotIn("/private", report)
        self.assertNotIn("SECRET", report)

    def test_bad_or_missing_summaries_remain_unconfirmed(self) -> None:
        path = self.directory / "summary.json"
        for payload in (None, "{bad json", "[]", '{"status":"failed"}',
                        '{"status":"/private/token"}',
                        '{"mesh":{"element_type":{},"nodes":true}}'):
            with self.subTest(payload=payload):
                if payload is not None:
                    path.write_text(payload)
                text = "\n".join(report_tool.existing_result_lines(self.directory))
                self.assertIn("not confirmed", text)
                self.assertNotIn("status = completed", text)
                self.assertNotIn("/private/token", text)
                self.assertNotIn("nodes: True", text)

    def test_oversized_summary_is_not_loaded(self) -> None:
        (self.directory / "summary.json").write_text('{"status":"completed"}')
        with patch.object(report_tool, "MAX_SUMMARY_BYTES", 8):
            text = "\n".join(report_tool.existing_result_lines(self.directory))
        self.assertIn("too large", text)
        self.assertNotIn("status = completed", text)

    def test_unknown_compiler_output_is_not_copied(self) -> None:
        with patch.object(report_tool, "command_output", return_value="/private/person/tool"):
            rows = dict(report_tool.environment_rows(self.directory))
        self.assertEqual(rows["gfortran"], "Not detected or version unavailable")
        self.assertNotIn("/private/person", str(rows))

    def test_version_command_errors_are_handled_without_leaking_paths(self) -> None:
        for error in (FileNotFoundError("/private/person"),
                      subprocess.TimeoutExpired("gfortran", 5)):
            with self.subTest(error=type(error)), patch.object(
                report_tool.subprocess, "run", side_effect=error
            ):
                self.assertIsNone(report_tool.command_output(["gfortran", "-dumpversion"]))

    def test_git_queries_are_read_only_and_do_not_report_filenames(self) -> None:
        with patch.object(report_tool, "command_output", side_effect=[
            "13.2.0", "a" * 40, " M private-file.txt",
        ]) as query:
            rows = dict(report_tool.environment_rows(self.directory))
        self.assertEqual(rows["Current checkout commit"], "a" * 40)
        self.assertEqual(rows["Current checkout state"], "Local changes present")
        self.assertNotIn("private-file.txt", str(rows))
        for call in query.call_args_list[1:]:
            self.assertIn("--no-optional-locks", call.args[0])

    def test_cli_creates_report_and_preserves_existing_feedback(self) -> None:
        output = self.directory / "reports" / "test-report.md"
        with patch.object(report_tool, "environment_rows", return_value=[]), redirect_stdout(io.StringIO()):
            self.assertEqual(report_tool.main(["--output", str(output)]), 0)
        self.assertIn("Tester Confirmation", output.read_text())
        output.write_text("My completed feedback")
        with patch.object(report_tool, "environment_rows", return_value=[]), redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as raised:
                report_tool.main(["--output", str(output)])
        self.assertEqual(raised.exception.code, 2)
        self.assertEqual(output.read_text(), "My completed feedback")


if __name__ == "__main__":
    unittest.main()
