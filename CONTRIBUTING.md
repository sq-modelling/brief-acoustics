# Contributing To BRIEF-Acoustics

BRIEF-Acoustics is currently a research-preview acoustic solver. Contributions
should preserve the documented mathematical conventions and remain inside, or
clearly extend, the supported scope stated in `README.md`.

## Before Opening A Pull Request

1. Open an issue for changes that alter equations, sign conventions, mesh
   contracts, result formats, or public case behavior.
2. Create a virtual environment and install the source tree with
   `python -m pip install -e .`.
3. Keep numerical assembly and dense solves in Fortran; keep case validation,
   orchestration, and post-processing in Python unless there is a measured
   reason to change that boundary.
4. Add a focused test or validation for changed behavior.
5. Run `bash scripts/validate.sh` and report the result in the pull request.

Do not commit generated binaries, Fortran module files, Python caches, result
directories, publisher PDFs, private meshes, credentials, or unpublished data.
Code comments should explain equations, conventions, and non-obvious decisions;
they should not repeat the syntax of the next line.

## Scientific Changes

A change to kernels, normals, quadrature, operator assembly, boundary-condition
elimination, or solver equations needs all of the following:

- the equation or literature reference that motivates the change;
- the affected normal and time-harmonic conventions;
- a reproducible comparison against an analytical solution, manufactured
  solution, or fixed regression;
- an explanation of any changed tolerance.

New numerical results are exploratory until they pass an appropriate reference
comparison and sensitivity or convergence check. Do not broaden public claims
from a demonstration alone.

## Style

- Fortran modules use `implicit none`, explicit `intent`, descriptive names,
  and the dependency order documented in `fortran/README.md`.
- Python code supports Python 3.10 or later and uses the existing dataclass and
  validation patterns.
- Public names and messages should use physical or mathematical terminology,
  not branch history or local folder names.
