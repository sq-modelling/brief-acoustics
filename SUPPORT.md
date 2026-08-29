# BRIEF-Acoustics Support

BRIEF-Acoustics is research software maintained on a best-effort basis. The
supported public alpha scope is stated in `README.md`; unsupported geometries,
physics, platforms, and formulations may be discussed, but they are not
treated as release defects.

Use a GitHub bug report for a reproducible failure in the documented workflow.
Use the external alpha report form when testing a fresh installation. Include
the exact release or commit, operating system, Python and `gfortran` versions,
BLAS/LAPACK implementation, command, and the shortest relevant error excerpt.

Before reporting a numerical discrepancy:

1. Run `brief-acoustics check case.toml`.
2. Confirm the mesh node ordering and `normal_orientation` against
   `docs/case-file-format.md`.
3. Reproduce the problem with a supplied example or a minimal shareable mesh.
4. Run `bash scripts/validate.sh` when practical and report whether it passes.

Do not post credentials, private filesystem paths, proprietary meshes, or
unpublished research data in a public issue.
