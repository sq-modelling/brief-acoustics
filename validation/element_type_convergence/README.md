# Rigid-Sphere Element-Type Evidence

This directory records the fixed rigid-Neumann sphere baselines for
linear three-node and quadratic six-node triangles. Every row uses `k=1`,
radius `a=1`, outward-from-solid mesh normals, and the same relative L2
surface-potential error against the partial-wave solution.

The results are stored in `baseline.csv`. They establish three points:

- ordinary linear NSBEM converges under refinement;
- linear Burton-Miller converges more slowly but reaches the common `5%`
  criterion at 2562 nodes;
- at the same 320 corner triangles, curved quadratic geometry is substantially
  more accurate for both formulations.

Reproduce the ordinary quadratic row with

```sh
NSBEM_FORMULATION=ordinary NSBEM_ELEMENT_TYPE=quadratic \
  NSBEM_SUBDIVISIONS=2 NSBEM_RESULTS_NAME=results_quadratic \
  python3 validation/sphere_scattering/run_validation.py
```

For the linear Burton-Miller refinement sequence, run

```sh
for subdivisions in 2 3 4; do
  NSBEM_FORMULATION=burton-miller NSBEM_ELEMENT_TYPE=linear \
    NSBEM_SUBDIVISIONS="$subdivisions" \
    NSBEM_RESULTS_NAME="results_bm_linear_sub${subdivisions}" \
    python3 validation/sphere_scattering/run_validation.py
done
```

The 2562-node Burton-Miller run is intentionally excluded from the minimal
CI-style command because the dense `2N x 2N` solve is much heavier than the
routine 642-node quadratic regression.
