# External Alpha Testing

This guide explains how to test a BRIEF-Acoustics alpha release from a fresh
copy and report installation or usability problems.

## Purpose

External testing asks whether a new user can install and run BRIEF-Acoustics from the
published instructions on an independent computer. It tests packaging,
portability, documentation, error messages, and the complete user workflow.

It does not replace the analytical and regression validations in
`validation/`, and it is not independent proof of the numerical method. The
repository's numerical claims remain tied to `bash scripts/validate.sh` and the
saved validation evidence.

## Platforms and Inputs

The supplied examples are sufficient for these tests. You do not need to
provide a private mesh or unpublished research data.

The supported alpha platforms are macOS and Ubuntu Linux. Windows is not yet a
supported native build target.

## Test Levels

### Level 1: User Workflow

Every tester should complete this level. It normally takes about 15 minutes
after the compiler and Python dependencies are available.

1. Start from a fresh clone or the source archive for the exact alpha tag.
2. Create a clean Python virtual environment.
3. Install BRIEF-Acoustics from the repository.
4. Check the default rigid-sphere case.
5. Run that case and inspect the generated summary and surface output.
6. Send a local test report or submit the external-test GitHub issue form,
   including successful reports.

On macOS or Linux:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .

brief-acoustics check examples/minimal_rigid_sphere/case.toml
brief-acoustics run examples/minimal_rigid_sphere/case.toml \
  --output external-test-results/minimal
```

The check command should identify:

```text
162 nodes
80 elements
quadratic triangles (6 nodes per element)
```

The run should create `surface_solution.csv`, `surface_magnitude.png`,
`solver.log`, and `summary.json`. The JSON status should be `completed`, and
its mesh metadata should report `element_type` as `quadratic`.

The important observation is whether the tester can reach this result from the
published instructions without changing source code. If live help is needed,
record where and why; that is useful documentation evidence rather than a
failed personal test.

### Level 2: Technical Regression

For a more detailed technical check, optionally run:

```sh
bash scripts/validate.sh
```

The final line should be:

```text
All validation checks passed.
```

This level exercises both linear and quadratic elements, ordinary and
Burton-Miller formulations, three boundary-condition families, analytical
sphere cases, the manufactured ellipsoid, and the fictitious-frequency
regression. It takes longer and is not required from every usability tester.

## What To Report

### Local Report (No GitHub Account Needed)

From the repository root, using the same Python environment as your test:

```sh
python3 scripts/make_test_report.py --results external-test-results/minimal
```

This creates `external-test-results/test-report.md`. Open it in a text editor,
fill in **Tester Confirmation**, the time required, any confusing steps and
your optional research application. Send the file to the maintainer as an
attachment, or copy its contents into a message. Nothing is uploaded automatically.
The output directory is ignored by Git.

If installation or the example failed before producing results, omit `--results`:

```sh
python3 scripts/make_test_report.py
```

The script needs only Python 3.10+ and its standard library; it does not need a
working BRIEF-Acoustics installation. It records OS, CPU architecture, Python,
`gfortran`, installed package metadata and the current Git commit when available.
For source archives, enter the tested release tag manually.

It does not run the solver or validation. Existing summary metadata and output
file presence are observations, not automatic passes: artifacts may be stale
and are not linked to the current commit. Installation and test outcomes remain
**Not confirmed** until you fill them in. Optional full-validation results and
BLAS/LAPACK details are entered manually.

An existing report is never overwritten. For another report, choose a new name:

```sh
python3 scripts/make_test_report.py --output external-test-results/test-report-second.md
```

The generated report excludes personal paths, hostnames, raw logs and mesh data.
Review anything you add manually before sharing.

### GitHub Issue Form

Alternatively, use the GitHub issue form named **External alpha test report**. Record:

- the release tag or commit tested;
- operating system and CPU architecture;
- Python and `gfortran` versions;
- BLAS/LAPACK implementation when known;
- whether installation, case checking, solving, and output inspection worked;
- whether any author assistance was needed;
- the exact command and relevant error text for a failure.

Do not include credentials, access tokens, private paths, proprietary meshes,
or unpublished input data. A short relevant log excerpt is normally enough.

## Severity

- **Release blocker:** a supported clean machine cannot install, compile, or
  complete the default example using the documentation.
- **High:** a command crashes, produces missing or inconsistent output, or the
  full validation reports a numerical failure.
- **Normal:** unclear instructions, confusing messages, portability friction,
  or a non-critical output problem.
- **Suggestion:** an improvement that is not needed for the documented alpha
  workflow.
