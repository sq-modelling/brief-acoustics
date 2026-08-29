# Rigid Sphere Scattering Validation

This validation compares the BRIEF-Acoustics surface solution with the
analytical partial-wave solution in `../../ana/analytic.py`.

Case:

- acoustic exterior Helmholtz problem;
- rigid sound-hard sphere;
- incident plane wave `exp(i k x)`;
- radius `a = 1`;
- wavenumber `k = 1`;
- mesh normals point outward from the solid sphere into the acoustic fluid;
- exterior Neumann boundary condition on the total field:
  `d(phi_total)/dn = 0`.

Run from the repository root:

```sh
python3 validation/sphere_scattering/run_validation.py
```

The default mesh uses quadratic six-node curved triangles. Set
`NSBEM_ELEMENT_TYPE=linear` explicitly when reproducing a linear baseline.

To exercise the augmented Burton-Miller system, use quadratic elements:

```sh
NSBEM_FORMULATION=burton-miller NSBEM_ELEMENT_TYPE=quadratic \
  NSBEM_RESULTS_NAME=results_burton_miller \
  python3 validation/sphere_scattering/run_validation.py
```

Outputs are written to `results/`:

- `sphere_mesh.dat`
- `nsbem_surface.csv`
- `analytic_surface.csv`
- `comparison_surface.csv`
- `surface_scattered_real.png`
- `surface_scattered_imag.png`
- `surface_error_abs.png`
- `summary.txt`

Historical linear mesh-refinement diagnostics are stored in:

- `results_sub1/` for 42 nodes and 80 elements;
- `results/` for 162 nodes and 320 elements;
- `results_sub3/` for 642 nodes and 1280 elements.

The convergence table is `convergence_summary.csv`.

Current result:

```text
sub1: relative L2 phi_sca error = 1.16682810e-01
sub2: relative L2 phi_sca error = 3.10702380e-02
sub3: relative L2 phi_sca error = 7.77618603e-03
```

The quadratic Burton-Miller run has 642 nodes, 320 curved elements, and relative
L2 scattered-potential error `5.86108131e-03`.  Linear Burton-Miller elements
converge more slowly because the augmented equation contains normal derivatives;
the release regression therefore uses the quadratic mesh rather than weakening
the common 5% accuracy threshold.

This analytical workflow uses the focused validation driver. Public user cases
use the generic `brief-acoustics run case.toml` command.
