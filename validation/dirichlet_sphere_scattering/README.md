# Dirichlet Sphere Scattering Validation

This case closes the end-to-end validation gap for the exterior Dirichlet
boundary condition.  A unit-amplitude plane wave scatters from a sphere of
radius `a=1`, with the pressure-release condition

```text
phi_total = phi_incident + phi_scattered = 0
```

on the sphere.  The Python analytical solution is the exact `A=1, B=0` special
case of the Robin partial-wave formula.

Run the ordinary formulation with the default quadratic elements from the
`modern` directory:

```sh
python3 validation/dirichlet_sphere_scattering/run_validation.py
```

Set `NSBEM_ELEMENT_TYPE=linear` explicitly to reproduce the linear baseline.

Run the quadratic Burton-Miller formulation:

```sh
NSBEM_FORMULATION=burton-miller NSBEM_ELEMENT_TYPE=quadratic \
  NSBEM_SUBDIVISIONS=2 NSBEM_RESULTS_NAME=results_burton_miller \
  python3 validation/dirichlet_sphere_scattering/run_validation.py
```

Each run writes its mesh, Fortran surface solution, analytical values,
node-by-node comparison, plots, and `summary.txt` under its results directory.

The fixed baselines at `k=1` are:

```text
Ordinary linear, 162 nodes:          relative L2 dphi_sca/dn = 2.67655932e-02
Burton-Miller quadratic, 642 nodes: relative L2 dphi_sca/dn = 8.11677778e-03
```

The error criterion uses `dphi_sca/dn`, which is the quantity actually solved
for under a Dirichlet condition. The prescribed surface potential is also
checked, but its near-machine-zero residual is not used as solver evidence.
