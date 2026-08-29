# Manufactured Ellipsoid Validation

This regression validates BRIEF-Acoustics on a smooth non-spherical geometry. A
radiating point source is placed at the centre of a solid triaxial ellipsoid
with semi-axes `(1.0, 0.8, 0.6)`. Its field

```text
phi(x) = exp(i k r) / r,   r = |x - x_source|
```

is an exact outgoing Helmholtz solution in the exterior fluid because the
singularity lies inside the excluded solid. The exact potential is prescribed
as spatially varying Dirichlet data. The test then compares the normal
derivative solved by Fortran with

```text
grad(phi) = exp(i k r) (i k r - 1) (x - x_source) / r^3
```

projected onto the nodal normals written by the solver.

Run the ordinary case with the default quadratic elements from the `modern`
directory:

```sh
python3 validation/manufactured_ellipsoid/run_validation.py
```

Set `NSBEM_ELEMENT_TYPE=linear` explicitly to reproduce the linear baseline.

Run the quadratic Burton-Miller case:

```sh
NSBEM_FORMULATION=burton-miller NSBEM_ELEMENT_TYPE=quadratic \
  NSBEM_RESULTS_NAME=results_burton_miller \
  python3 validation/manufactured_ellipsoid/run_validation.py
```

Each results directory contains the mesh, raw Fortran surface solution,
node-by-node comparison, derivative-error figure, and a pass/fail summary.

The fixed baselines at `k=1` are:

```text
Ordinary linear, 162 nodes:          relative L2 dphi/dn = 1.61897423e-02
Burton-Miller quadratic, 642 nodes: relative L2 dphi/dn = 9.36468404e-03
```

Both are checked against a fixed `5%` acceptance limit.

The smooth ellipsoid is intentional. A sharp-edged cube requires explicit
corner free-term and normal treatment; adding it without that analysis would
mix a geometry-regularity question into this solver regression.
