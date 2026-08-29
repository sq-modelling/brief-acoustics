# Robin Sphere Scattering Validation

This validation checks the external Robin assembly branch of BRIEF-Acoustics.

Case:

- acoustic exterior Helmholtz problem;
- sphere of radius `a = 1`;
- incident plane wave `exp(i k x)`;
- wavenumber `k = 1`;
- mesh normals point outward from the solid into the acoustic fluid;
- Robin boundary condition on the total field:
  `A * phi_total + B * d(phi_total)/dn = 0`;
- default coefficients: `A = 1`, `B = 0.5`.

Run from the repository root:

```sh
python3 validation/robin_sphere_scattering/run_validation.py
```

This command uses quadratic six-node curved triangles by default. Set
`NSBEM_ELEMENT_TYPE=linear` explicitly when reproducing the linear refinement
baseline below.

Run the same case with the augmented Burton-Miller system and quadratic elements:

```sh
NSBEM_FORMULATION=burton-miller NSBEM_ELEMENT_TYPE=quadratic \
  NSBEM_SUBDIVISIONS=2 NSBEM_RESULTS_NAME=results_burton_miller \
  python3 validation/robin_sphere_scattering/run_validation.py
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
- `results_sub2/` for 162 nodes and 320 elements;
- `results/` for 642 nodes and 1280 elements.

The convergence table is `convergence_summary.csv`.

Current result:

```text
sub1: relative L2 phi_sca error = 3.79648852e-01
sub2: relative L2 phi_sca error = 1.05238782e-01
sub3: relative L2 phi_sca error = 2.84578742e-02
```

The quadratic Burton-Miller run has 642 nodes, 320 curved elements, relative L2
scattered-potential error `1.27110176e-02`, and a maximum numerical Robin
residual below `1.0e-15`.

The validation uses the same validation-specific Fortran executable as the
rigid sphere case, but passes `robin` mode to exercise the Robin matrix assembly.
