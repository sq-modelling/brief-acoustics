#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

# Keep dense BLAS work bounded and predictable on developer laptops and CI.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-2}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-2}"

echo "==> Fortran syntax check"
make -C fortran syntax

echo "==> Build validation executable"
make -C fortran run-case
make -C fortran validate-rigid-sphere
make -C fortran validate-manufactured-exterior

echo "==> Burton-Miller matrix assembly check"
make -C fortran test-burton-miller-assembly
fortran/bin/test_burton_miller_assembly

echo "==> Python compile check"
PYTHONPATH="$ROOT_DIR/python" python3 -m compileall -q \
  "$ROOT_DIR/python/brief_acoustics" \
  "$ROOT_DIR/ana" \
  "$ROOT_DIR/validation/sphere_scattering" \
  "$ROOT_DIR/validation/robin_sphere_scattering" \
  "$ROOT_DIR/validation/dirichlet_sphere_scattering" \
  "$ROOT_DIR/validation/manufactured_ellipsoid" \
  "$ROOT_DIR/validation/fictitious_frequency_sphere" \
  "$ROOT_DIR/validation/two_sphere_gap" \
  "$ROOT_DIR/validation/public_case_runner"

echo "==> Python case and mesh contract tests"
PYTHONPATH="$ROOT_DIR/python" python3 -m unittest discover -s python/tests -v

echo "==> Rigid Neumann sphere validation (linear elements)"
NSBEM_ELEMENT_TYPE=linear \
NSBEM_SUBDIVISIONS=2 \
python3 validation/sphere_scattering/run_validation.py

echo "==> Robin sphere validation (linear elements)"
NSBEM_ELEMENT_TYPE=linear \
NSBEM_SUBDIVISIONS=3 \
python3 validation/robin_sphere_scattering/run_validation.py

echo "==> Dirichlet sphere validation (linear elements)"
NSBEM_ELEMENT_TYPE=linear \
NSBEM_SUBDIVISIONS=2 \
python3 validation/dirichlet_sphere_scattering/run_validation.py

echo "==> Ordinary rigid sphere validation (quadratic elements)"
NSBEM_ELEMENT_TYPE=quadratic \
NSBEM_SUBDIVISIONS=2 \
NSBEM_RESULTS_NAME=results_quadratic \
python3 validation/sphere_scattering/run_validation.py

echo "==> Burton-Miller rigid sphere validation (quadratic elements)"
NSBEM_FORMULATION=burton-miller \
NSBEM_ELEMENT_TYPE=quadratic \
NSBEM_SUBDIVISIONS=2 \
NSBEM_RESULTS_NAME=results_burton_miller \
python3 validation/sphere_scattering/run_validation.py

echo "==> Burton-Miller Robin sphere validation (quadratic elements)"
NSBEM_FORMULATION=burton-miller \
NSBEM_ELEMENT_TYPE=quadratic \
NSBEM_SUBDIVISIONS=2 \
NSBEM_RESULTS_NAME=results_burton_miller \
python3 validation/robin_sphere_scattering/run_validation.py

echo "==> Burton-Miller Dirichlet sphere validation (quadratic elements)"
NSBEM_FORMULATION=burton-miller \
NSBEM_ELEMENT_TYPE=quadratic \
NSBEM_SUBDIVISIONS=2 \
NSBEM_RESULTS_NAME=results_burton_miller \
python3 validation/dirichlet_sphere_scattering/run_validation.py

echo "==> Manufactured ellipsoid validation (linear elements)"
NSBEM_ELEMENT_TYPE=linear \
NSBEM_SUBDIVISIONS=2 \
python3 validation/manufactured_ellipsoid/run_validation.py

echo "==> Burton-Miller manufactured ellipsoid validation (quadratic elements)"
NSBEM_FORMULATION=burton-miller \
NSBEM_ELEMENT_TYPE=quadratic \
NSBEM_SUBDIVISIONS=2 \
NSBEM_RESULTS_NAME=results_burton_miller \
python3 validation/manufactured_ellipsoid/run_validation.py

echo "==> Fictitious-frequency regression"
NSBEM_SWEEP_PROFILE=regression \
NSBEM_RESULTS_NAME=results_regression \
python3 validation/fictitious_frequency_sphere/run_sweep.py

echo "==> Public case runner from an external working directory"
(
  cd /tmp
  PYTHONPATH="$ROOT_DIR/python" python3 -m brief_acoustics run \
    "$ROOT_DIR/examples/minimal_rigid_sphere/case.toml" \
    --output "$ROOT_DIR/validation/public_case_runner/results_minimal" \
    --no-build
  PYTHONPATH="$ROOT_DIR/python" python3 -m brief_acoustics run \
    "$ROOT_DIR/examples/burton_miller_fictitious_sphere/case.toml" \
    --output "$ROOT_DIR/validation/public_case_runner/results_burton_miller" \
    --no-build
)
PYTHONPATH="$ROOT_DIR/python" python3 validation/public_case_runner/check_outputs.py

echo "==> Check validation status"
grep -q "validation status = PASS" validation/sphere_scattering/results/summary.txt
grep -q "validation status = PASS" validation/robin_sphere_scattering/results/summary.txt
grep -q "validation status = PASS" validation/dirichlet_sphere_scattering/results/summary.txt
grep -q "validation status = PASS" validation/sphere_scattering/results_quadratic/summary.txt
grep -q "validation status = PASS" validation/sphere_scattering/results_burton_miller/summary.txt
grep -q "validation status = PASS" validation/robin_sphere_scattering/results_burton_miller/summary.txt
grep -q "validation status = PASS" validation/dirichlet_sphere_scattering/results_burton_miller/summary.txt
grep -q "validation status = PASS" validation/manufactured_ellipsoid/results/summary.txt
grep -q "validation status = PASS" validation/manufactured_ellipsoid/results_burton_miller/summary.txt
grep -q "validation status = PASS" validation/fictitious_frequency_sphere/results_regression/summary.txt

echo "All validation checks passed."
