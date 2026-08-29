"""Command-line entry point for the public BRIEF-Acoustics Python workflow."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from .case_config import CaseConfigurationError
from .mesh_io import tetrahedron_mesh, write_fortran_mesh
from .runner import CaseRunError, check_case, run_case


def build_parser() -> argparse.ArgumentParser:
    """Create the small command-line interface for Python helper utilities.

    The ``run`` command is the supported public source-tree workflow.  The
    remaining commands validate cases or create a minimal test mesh.
    """

    parser = argparse.ArgumentParser(
        description="BRIEF-Acoustics preprocessing, solver, and post-processing utilities."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Validate and solve one TOML case.")
    run_parser.add_argument("case", type=Path, help="Public TOML case file.")
    run_parser.add_argument(
        "--output",
        type=Path,
        help="Override the output directory declared by the case file.",
    )
    run_parser.add_argument(
        "--no-build",
        action="store_true",
        help="Use an already-built fortran/bin/run_au_case executable.",
    )
    run_parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip the surface PNG while retaining CSV and JSON outputs.",
    )

    check_parser = subparsers.add_parser("check", help="Validate a TOML case and mesh without solving.")
    check_parser.add_argument("case", type=Path, help="Public TOML case file.")

    mesh_parser = subparsers.add_parser("mesh-smoke", help="Write a tiny BRIEF-Acoustics mesh file.")
    mesh_parser.add_argument("--mesh", required=True, type=Path, help="Mesh file to write.")

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run one command from the helper command-line interface."""

    args = build_parser().parse_args(argv)

    try:
        if args.command == "run":
            artifacts = run_case(
                args.case,
                output_override=args.output,
                build=not args.no_build,
                create_plot=not args.no_plot,
            )
            print(f"BRIEF-Acoustics case completed: {artifacts.output_directory}")
            print(f"Surface CSV: {artifacts.surface_csv}")
            print(f"Summary JSON: {artifacts.summary_json}")
            if artifacts.surface_figure is not None:
                print(f"Surface figure: {artifacts.surface_figure}")
            return 0

        if args.command == "check":
            config, mesh = check_case(args.case)
            print(f"Case '{config.name}' is valid.")
            print(
                f"Mesh: {len(mesh.nodes)} nodes, {len(mesh.elements)} elements, "
                f"{mesh.element_type} triangles ({mesh.nodes_per_element} nodes per element)"
            )
            print(
                f"Physics: k={config.physics.wavenumber:.8g}, "
                f"omega={config.physics.angular_frequency:.8g}"
            )
            print(f"Formulation: {config.solver.formulation}")
            return 0

        if args.command == "mesh-smoke":
            # Write the smallest possible closed mesh so the Fortran mesh reader can
            # be tested without needing an external mesh generator.
            write_fortran_mesh(tetrahedron_mesh(), args.mesh)
            print(f"Wrote BRIEF-Acoustics tetrahedron mesh to {args.mesh}")
            return 0

    except (CaseConfigurationError, CaseRunError, FileNotFoundError, ValueError) as error:
        parser = build_parser()
        parser.error(str(error))

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
