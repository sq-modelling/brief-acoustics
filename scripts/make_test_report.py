#!/usr/bin/env python3
"""Create a shareable installation report without running or importing the solver.

Only the Python standard library is needed. This is intentional: a tester whose
package installation failed should still be able to report their environment.
Run this script with the same Python interpreter used for the installation.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from importlib import metadata
import json
from pathlib import Path
import platform
import re
import subprocess


REPOSITORY = Path(__file__).resolve().parents[1]
UNKNOWN = "Not confirmed"
MAX_SUMMARY_BYTES = 1024 * 1024
OUTPUT_FILES = ("surface_solution.csv", "surface_magnitude.png", "solver.log")


def command_output(arguments: list[str], directory: Path | None = None) -> str | None:
    """Read version information; missing commands must not prevent a report.

    No shell is used, and each command has a short timeout. Error output is not
    copied to the report because it can contain private installation paths.
    """
    try:
        result = subprocess.run(
            arguments, cwd=directory, capture_output=True, text=True,
            check=False, timeout=5,
        )
    except (OSError, subprocess.SubprocessError, UnicodeError):
        return None
    return result.stdout.strip() if result.returncode == 0 else None


def table_text(value: str) -> str:
    """Keep a system version on one Markdown table line."""
    return " ".join(value.split()).replace("|", "/").replace("`", "'")


def environment_rows(repository: Path) -> list[tuple[str, str]]:
    """Collect versions, not usernames, hostnames, paths, or environment dumps."""
    operating_system = f"{platform.system()} {platform.release()}"
    if platform.system() == "Darwin":
        operating_system = f"macOS {platform.mac_ver()[0]}"
    elif platform.system() == "Linux":
        try:
            operating_system = platform.freedesktop_os_release().get(
                "PRETTY_NAME", operating_system
            )
        except OSError:
            pass

    compiler = command_output(["gfortran", "-dumpfullversion", "-dumpversion"])
    if compiler is None or re.fullmatch(r"\d+(?:\.\d+)*", compiler) is None:
        compiler = "Not detected or version unavailable"

    try:
        package_version = metadata.version("brief-acoustics")
    except (metadata.PackageNotFoundError, OSError, ValueError):
        package_version = "Not detected in this Python environment"

    # --no-optional-locks keeps these Git queries read-only. A source archive
    # may have no Git metadata; that is not an installation failure.
    commit = command_output(
        ["git", "--no-optional-locks", "rev-parse", "HEAD"], repository
    )
    if commit is None or re.fullmatch(r"[0-9a-f]{40,64}", commit) is None:
        commit = "Unavailable (for example, a source archive)"
        worktree = UNKNOWN
    else:
        status = command_output(
            ["git", "--no-optional-locks", "status", "--porcelain", "--untracked-files=normal"],
            repository,
        )
        worktree = UNKNOWN if status is None else (
            "Local changes present" if status else "Clean"
        )

    return [
        ("OS", operating_system),
        ("CPU architecture", platform.machine() or UNKNOWN),
        ("Python", f"{platform.python_version()} ({platform.python_implementation()})"),
        ("gfortran", compiler),
        ("Installed package metadata", package_version),
        ("Current checkout commit", commit),
        ("Current checkout state", worktree),
        ("BLAS/LAPACK implementation", "Not collected; add manually if known"),
    ]


def existing_result_lines(results: Path | None) -> list[str]:
    """Describe supplied artifacts without claiming a new test has passed.

    Existing summary files do not record the tested Git commit. They may also
    survive a later failed run. Therefore 'completed' is reported as evidence
    from a file, not as proof of a successful test of the current checkout.
    Only selected numeric and enum fields are copied; case names, paths, raw
    logs and arbitrary JSON text are deliberately excluded.
    """
    if results is None:
        return ["- No results directory supplied; example outcome is not confirmed."]

    lines = []
    try:
        with (results / "summary.json").open("rb") as handle:
            payload = handle.read(MAX_SUMMARY_BYTES + 1)
        if len(payload) > MAX_SUMMARY_BYTES:
            raise ValueError("Summary exceeds the report reader's size limit")
        summary = json.loads(payload)
        if not isinstance(summary, dict):
            raise ValueError("Summary must be an object")
    except FileNotFoundError:
        lines.append("- summary.json: missing; example outcome is not confirmed.")
    except (OSError, ValueError, RecursionError):
        lines.append("- summary.json: unreadable, invalid or too large; outcome is not confirmed.")
    else:
        if summary.get("status") == "completed":
            lines.append("- summary.json reports `status = completed` (existing artifact only).")
        else:
            lines.append("- summary.json does not report completion; outcome is not confirmed.")
        mesh = summary.get("mesh")
        if isinstance(mesh, dict):
            for key in ("nodes", "elements", "nodes_per_element"):
                value = mesh.get(key)
                if type(value) is int and 0 < value < 10**12:
                    lines.append(f"- Recorded mesh {key}: {value}.")
            if mesh.get("element_type") in ("linear", "quadratic"):
                lines.append(f"- Recorded element type: {mesh['element_type']}.")

    for name in OUTPUT_FILES:
        try:
            present = (results / name).is_file() and (results / name).stat().st_size > 0
            status = "present and non-empty" if present else "missing or empty"
        except OSError:
            status = "not readable"
        lines.append(f"- {name}: {status}; content not validated.")
    lines.append(
        "- These artifacts are not linked to the current commit. The tester must "
        "confirm which run produced them; presence alone is not a passing test."
    )
    return lines


def build_report(repository: Path, results: Path | None) -> str:
    """Assemble automatic observations followed by editable feedback fields."""
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = [
        "# BRIEF-Acoustics External Test Report", "",
        f"Report generated: {generated}", "",
        "This tool did not install the package, run the solver, or run validation.",
        "No information was uploaded. Review this report before sharing it.", "",
        "## Automatically Collected Environment", "",
        "| Item | Observation |", "| --- | --- |",
    ]
    lines.extend(
        f"| {label} | {table_text(value)} |"
        for label, value in environment_rows(repository)
    )
    lines += [
        "", "Package metadata does not prove that imports, compilation or installation worked.",
        "The environment and checkout above are those at report-generation time, not necessarily test time.",
        "", "## Existing Example Artifacts", "",
        *existing_result_lines(results),
        "", "## Tester Confirmation", "",
        "Replace 'Not confirmed' only for steps you actually checked.",
        "Use 'Not run' for skipped steps, and 'Failed' for unsuccessful steps.", "",
        "- Release tag / commit actually tested: Not confirmed",
        "- Fresh clone or source archive: Not confirmed",
        "- Clean Python environment used: Not confirmed",
        "- Installation succeeded: Not confirmed",
        "- Case check succeeded: Not confirmed",
        "- Minimal example succeeded: Not confirmed",
        "- Summary and surface output inspected: Not confirmed",
        "- Optional full validation (scripts/validate.sh): Not confirmed",
        "- Supplied artifacts belong to this test: Not confirmed",
        "- Time required (installation / example / optional validation): [fill in]",
        "- Assistance needed, and at which step: [fill in or none]",
        "", "## Confusing Or Failed Step", "",
        "[Describe the step and command, or write 'None'. Include only a short,",
        "sanitized error excerpt. Remove private paths and credentials.]",
        "", "## Potential Research Application", "",
        "[Optional: describe your intended application without unpublished or confidential details.]",
        "", "## Sharing", "",
        "Send this Markdown file as an attachment, or copy its contents into a message",
        "or the GitHub external-test issue form. A GitHub account is not required",
        "to send the report directly to the maintainer.",
        "Do not add passwords, tokens, private paths, proprietary meshes or unpublished data.", "",
    ]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results", type=Path,
        help="Existing example output directory (optional; does not run the example).",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("external-test-results/test-report.md"),
        help="New Markdown file; existing files are never overwritten.",
    )
    arguments = parser.parse_args(argv)
    report = build_report(REPOSITORY, arguments.results)
    try:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        # Exclusive creation protects any feedback already filled in by a tester.
        with arguments.output.open("x", encoding="utf-8") as handle:
            handle.write(report)
    except FileExistsError:
        parser.exit(2, "Report already exists. Use --output with a different filename.\n")
    except OSError:
        parser.exit(2, "Cannot write report. Choose a writable location with --output.\n")
    print(f"Report written: {arguments.output}")
    print("Fill in Tester Confirmation and review before sharing. Nothing was uploaded.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
