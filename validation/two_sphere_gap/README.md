# Two-Sphere Small-Gap Demo

This demonstration solves an exterior acoustic radiation problem for two
nearby spheres using the BRIEF-Acoustics solver.

The default case uses quadratic six-node elements on two equal spheres with
opposite prescribed Neumann values:

```text
dphi/dn =  1 on sphere 1
dphi/dn = -1 on sphere 2
```

The surface solution is computed by the Fortran NSBEM core. Python then
evaluates the potential and pressure along the line through the narrow gap using
Eq. (4.1) of:

Sun, Klaseboer, Khoo & Chan, "Boundary regularized integral equation
formulation of the Helmholtz equation in acoustics", Royal Society Open Science,
2015.

The mesh normals point outward from each solid, while Eq. (4.1) is written with
the normal pointing out of the exterior fluid domain. The field evaluator
derives the required equation sign from the readable `outward-from-solid`
orientation. This is the same convention used by the analytical sphere
regressions.

## Run

From the `modern` folder:

```sh
python3 validation/two_sphere_gap/run_demo.py
```

Outputs are written to:

```text
validation/two_sphere_gap/results/
```

The main files are:

```text
surface_solution.csv       solved nodal surface potential, derivative, pressure
gap_line_field.csv         potential and pressure along the gap line
gap_pressure_abs.png       pressure magnitude along the gap line
gap_pressure_complex.png   real and imaginary pressure along the gap line
summary.txt                compact run summary
```

## Explicit Linear Run

```sh
env NSBEM_ELEMENT_TYPE=linear NSBEM_RESULTS_NAME=results_linear \
  python3 validation/two_sphere_gap/run_demo.py
```

## Adjustable Parameters

```text
NSBEM_RADIUS          sphere radius, default 1
NSBEM_GAP             surface-to-surface gap, default 0.05
NSBEM_ELEMENT_TYPE    quadratic by default; set linear explicitly to override
NSBEM_SUBDIVISIONS    sphere subdivision level, default 2
NSBEM_FIELD_POINTS    number of points across the gap, default 101
NSBEM_RESULTS_NAME    output folder name, default results
```

## Current Limitation

The field-line evaluation currently uses Eq. (4.1), the ordinary field-point
boundary integral expression. This case is explicitly a workflow
**demonstration**, not an accuracy or superiority claim. Eq. (4.2), mesh and
quadrature convergence, and an independent or conventional-BEM comparator are
required before making a near-boundary accuracy claim.
