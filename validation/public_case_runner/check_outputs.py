"""Check that the public runner reproduces the Burton-Miller regression solve."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

THIS_DIR = Path(__file__).resolve().parent
MODERN_DIR = THIS_DIR.parents[1]


def main() -> int:
    """Validate portable outputs and compare the Burton-Miller nodal solution."""

    minimal_directory = THIS_DIR / "results_minimal"
    bm_directory = THIS_DIR / "results_burton_miller"
    for directory, expected_nodes, expected_formulation, expected_element_type in (
        (minimal_directory, 162, "ordinary", "quadratic"),
        (bm_directory, 642, "burton-miller", "quadratic"),
    ):
        summary = json.loads((directory / "summary.json").read_text(encoding="utf-8"))
        if summary["status"] != "completed":
            raise AssertionError(f"Public runner did not complete: {directory}")
        if summary["mesh"]["nodes"] != expected_nodes:
            raise AssertionError(f"Unexpected node count in {directory}")
        if summary["solver"]["formulation"] != expected_formulation:
            raise AssertionError(f"Unexpected formulation in {directory}")
        if summary["mesh"]["element_type"] != expected_element_type:
            raise AssertionError(f"Unexpected element type in {directory}")
        if summary["mesh"]["normal_orientation"] != "outward-from-solid":
            raise AssertionError(f"Unexpected normal orientation in {directory}")
        if "exterior_domain_normal_sign" in summary["solver"]:
            raise AssertionError(
                "Internal exterior-domain normal sign leaked into public result metadata"
            )
        for filename in summary["outputs"].values():
            if filename is not None and not (directory / filename).is_file():
                raise AssertionError(f"Missing public-runner output: {directory / filename}")

    generic_solution = read_scattered_potential(bm_directory / "surface_solution.csv")
    regression_solution = read_scattered_potential(
        MODERN_DIR
        / "validation"
        / "fictitious_frequency_sphere"
        / "results_regression"
        / "surface_solutions"
        / "burton-miller_ka_9p35580000.csv"
    )
    relative_difference = np.linalg.norm(generic_solution - regression_solution) / np.linalg.norm(
        regression_solution
    )
    if relative_difference > 1.0e-12:
        raise AssertionError(
            "Public runner differs from the Burton-Miller regression baseline: "
            f"relative L2 difference={relative_difference:.8e}"
        )

    print("Public case runner outputs: PASS")
    print(f"Burton-Miller workflow relative L2 difference = {relative_difference:.8e}")
    return 0


def read_scattered_potential(path: Path) -> np.ndarray:
    """Read the complex scattered potential from either compatible CSV."""

    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise AssertionError(f"Surface CSV is empty: {path}")
    return np.asarray(
        [
            complex(float(row["phi_sca_re"]), float(row["phi_sca_im"]))
            for row in rows
        ],
        dtype=complex,
    )


if __name__ == "__main__":
    raise SystemExit(main())
