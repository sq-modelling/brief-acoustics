# BRIEF-Acoustics

[![Validation](https://github.com/sq-modelling/brief-acoustics/actions/workflows/ci.yml/badge.svg)](https://github.com/sq-modelling/brief-acoustics/actions/workflows/ci.yml)

BRIEF-Acoustics is a research-preview non-singular boundary element solver for
frequency-domain acoustic Helmholtz problems. Its name refers to the Boundary
Regularised Integral Equation Formulation (BRIEF) developed in the cited
literature.

The public alpha combines a Fortran numerical core with a Python case runner,
mesh audit, post-processing, and validation workflows. The current public scope
is acoustic exterior Helmholtz problems for one closed particle, with
Dirichlet, Neumann, and Robin boundary conditions.

## Status

`v0.1.0-alpha.1` public prerelease and research preview.

The documented source-tree command `brief-acoustics run case.toml` validates a case and
mesh, builds the Fortran driver, solves the boundary-value problem, and writes
portable results. This alpha is not yet a prebuilt binary distribution or a
general multiple-particle acoustic solver.

## What Works

- Modern Fortran modules for geometry, quadrature, Helmholtz kernels, acoustic
  boundary conditions, surface operator assembly, and dense solves.
- Strict TOML case parsing, complete closed-surface mesh auditing, Fortran
  orchestration, CSV output, JSON metadata, and surface plots in Python.
- A one-command source-tree workflow for ordinary NSBEM and non-singular
  Burton-Miller cases.
- Rigid sound-hard sphere validation against a partial-wave analytical solution.
- Robin sphere validation against a partial-wave analytical solution.
- Pressure-release Dirichlet sphere validation against a partial-wave analytical
  solution.
- Selectable ordinary NSBEM and non-singular Burton-Miller `2N x 2N` systems
  for one exterior particle. Dirichlet, Neumann, and Robin are analytically
  validated through both solver paths.
- Matrix-level checks of all four Burton-Miller blocks and boundary-condition
  elimination.
- Focused ordinary/Burton-Miller comparison at the published rigid-sphere
  fictitious frequencies `ka=9.3558` and `ka=10.41712`.
- Quadratic six-node curved triangles as the public default, with explicit
  linear-element baselines retained for comparison and regression coverage.
- A non-spherical triaxial-ellipsoid manufactured-solution regression.
- Two-sphere small-gap radiation demonstration with field evaluation across the
  gap.
- Reproducible validation CSV files and figures.
- Default-quadratic ordinary and quadratic Burton-Miller example cases.

## Current Limitations

- Acoustic Helmholtz problems only.
- A source checkout, Python environment, Fortran compiler, and BLAS/LAPACK are
  required; no wheel or prebuilt executable is provided yet.
- The public case runner supports one closed, connected exterior particle,
  constant total-field Dirichlet/Neumann/Robin data, and either no incident
  field or one plane wave.
- Burton-Miller is currently restricted to one closed exterior particle, an
  outward-from-solid mesh, and real positive wavenumber.
- Multiple-particle, transmission, and multilayer Burton-Miller systems are not
  yet enabled; ordinary two-particle NSBEM remains available for demonstrations.
- The focused fictitious-frequency regression is not the full 5762-node,
  `ka=1..40`, step-`0.001` calculation from the 2023 paper.
- The two-sphere Eq. (4.1) field calculation is an application demonstration,
  not yet a controlled near-boundary accuracy or conventional-BEM comparison.
- Linux GitHub Actions runs the release validation on every push and pull
  request.
- Mesh generation in the validation scripts is intentionally minimal.

## Requirements

- `gfortran`
- BLAS/LAPACK
  - macOS: Accelerate framework is used by default.
  - Linux: `-llapack -lblas` is used by default.
- Python 3.10+
- Python packages: `numpy`, `scipy`, `matplotlib`

On Ubuntu, install the system toolchain with:

```sh
sudo apt-get update
sudo apt-get install -y build-essential gfortran liblapack-dev libblas-dev python3-venv
```

On macOS, install the Xcode command-line tools and GNU Fortran. With Homebrew:

```sh
xcode-select --install
brew install gcc python
```

The macOS build links against the system Accelerate framework.

## Quick Start

Clone the repository, create an isolated Python environment, install the source
tree, and run the minimal example:

```sh
git clone https://github.com/sq-modelling/brief-acoustics.git
cd brief-acoustics
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
brief-acoustics check examples/minimal_rigid_sphere/case.toml
brief-acoustics run examples/minimal_rigid_sphere/case.toml
```

The result is written to `examples/minimal_rigid_sphere/results/`. See
[`docs/case-file-format.md`](docs/case-file-format.md) for the complete input,
mesh, convention, and output contract, and [`examples/README.md`](examples/README.md)
for the quadratic Burton-Miller example. Omitting `[mesh].element_type` selects
quadratic six-node triangles; set it explicitly to `"linear"` for a linear mesh.

Run the complete release-level validation matrix with:

```sh
bash scripts/validate.sh
```

The script performs the checks used by the Linux CI job:

1. Fortran syntax check.
2. Build the public runner and validation executables.
3. Compile Python modules.
4. Run unit tests for the TOML contract and malformed-mesh rejection.
5. Check the Burton-Miller matrix blocks for Dirichlet, Neumann, and Robin data.
6. Run explicit linear Dirichlet, Neumann, and Robin sphere baselines.
7. Exercise ordinary quadratic elements and quadratic Burton-Miller versions of
   all three sphere boundary conditions.
8. Run ordinary and Burton-Miller manufactured ellipsoid validations.
9. Compare ordinary and Burton-Miller solutions at both published fictitious
   frequencies.
10. Run both public examples from outside the repository and compare the
    Burton-Miller result with its certified regression.
11. Check that every numerical validation reports `PASS`.

## Manual Commands

Validate a public case without solving:

```sh
brief-acoustics check examples/minimal_rigid_sphere/case.toml
```

Run a public case and choose an output directory:

```sh
brief-acoustics run examples/minimal_rigid_sphere/case.toml --output /tmp/brief-acoustics-run
```

Fortran syntax check:

```sh
make -C fortran syntax
```

Build the public Fortran case driver:

```sh
make -C fortran run-case
```

Build the sphere validation executable:

```sh
make -C fortran validate-rigid-sphere
```

Run the default quadratic rigid Neumann validation:

```sh
python3 validation/sphere_scattering/run_validation.py
```

Run the default quadratic Robin validation:

```sh
python3 validation/robin_sphere_scattering/run_validation.py
```

Run the default quadratic Dirichlet validation:

```sh
python3 validation/dirichlet_sphere_scattering/run_validation.py
```

Run the Burton-Miller rigid-sphere validation; quadratic elements remain the
default:

```sh
NSBEM_FORMULATION=burton-miller \
  NSBEM_RESULTS_NAME=results_burton_miller \
  python3 validation/sphere_scattering/run_validation.py
```

Run the matrix-level Burton-Miller assembly check:

```sh
make -C fortran test-burton-miller-assembly
fortran/bin/test_burton_miller_assembly
```

Run the default quadratic non-spherical manufactured ellipsoid validation:

```sh
python3 validation/manufactured_ellipsoid/run_validation.py
```

Run the 18-point comparison around both published fictitious frequencies:

```sh
python3 validation/fictitious_frequency_sphere/run_sweep.py
```

Run the two-sphere small-gap demonstration:

```sh
python3 validation/two_sphere_gap/run_demo.py
```

Run the two-sphere demonstration with explicit linear elements instead:

```sh
env NSBEM_ELEMENT_TYPE=linear NSBEM_RESULTS_NAME=results_linear \
  python3 validation/two_sphere_gap/run_demo.py
```

Run the Python mesh smoke test:

```sh
PYTHONPATH=python python3 -m brief_acoustics mesh-smoke --mesh /tmp/brief_acoustics_mesh.dat
```

## Validation Results

Rigid sound-hard sphere, explicit linear-element refinement baseline:

```text
42 nodes:   relative L2 phi_sca error = 1.16682810e-01
162 nodes:  relative L2 phi_sca error = 3.10702380e-02
642 nodes:  relative L2 phi_sca error = 7.77618603e-03
```

Robin sphere, with `phi_total + 0.5 dphi_total/dn = 0`:

```text
42 nodes:   relative L2 phi_sca error = 3.79648852e-01
162 nodes:  relative L2 phi_sca error = 1.05238782e-01
642 nodes:  relative L2 phi_sca error = 2.84578742e-02
```

Quadratic Burton-Miller sphere tests at 642 nodes and `k=1`:

```text
Rigid Neumann: relative L2 phi_sca error = 5.86108131e-03
Robin:         relative L2 phi_sca error = 1.27110176e-02
Dirichlet:     relative L2 dphi_sca/dn error = 8.11677778e-03
```

The ordinary linear Dirichlet case at 162 nodes gives
`2.67655932e-02` relative L2 error in the solved normal derivative. The full
linear/quadratic Neumann table is stored in
[`validation/element_type_convergence/baseline.csv`](validation/element_type_convergence/baseline.csv).

Manufactured point-source field on a triaxial ellipsoid:

```text
Ordinary linear, 162 nodes:          relative L2 dphi/dn error = 1.61897423e-02
Burton-Miller quadratic, 642 nodes: relative L2 dphi/dn error = 9.36468404e-03
```

At the two published fictitious frequencies on the same 642-node quadratic
mesh:

```text
ka = 9.3558:  ordinary error = 5.03276827e-01, BM error = 5.02601981e-02
ka = 10.41712: ordinary error = 7.57575174e-01, BM error = 7.22217538e-02
```

The public case file describes mesh normals directly:

```text
normal_orientation = "outward-from-solid"
```

Python derives an internal exterior-domain normal sign from this orientation.
It maps the stored mesh normal to the normal used by the exterior-domain
equations; it is not a solid-angle coefficient. The convention is validated by
the analytical sphere tests.

BRIEF cancels the local solid-angle coefficient $c(\mathbf{x}_0)$ analytically.
The implementation therefore never substitutes the smooth-surface value
$1/2$, or equivalently $2\pi$ for the unnormalised Green function. Explicit
$4\pi$ diagonal contributions remain in the unbounded-exterior equations.
They originate from the auxiliary Laplace identity on the surface at infinity,
not from a local solid-angle approximation.

### Burton--Miller coupling length

The Burton--Miller coupling is not an arbitrary dimensionless constant. The
normal-derivative boundary integral equation contains one additional inverse
length relative to the ordinary equation, so its coefficient must have units
of length. BRIEF-Acoustics distinguishes:

- `wavenumber` $k$, with units $L^{-1}$;
- `characteristic_length` $a$, a representative body scale in mesh units; and
- the coupling length $\beta=\min(a,1/k)$.

Thus $\beta=a$ for $ka\leq1$ and $\beta=1/k$ for $ka>1$. For a sphere, $a$ is
normally its radius. The complex row weight is $w=i\sigma\beta$, where
$\sigma=-1$ for the supported outward-from-solid exterior mesh, giving
$w=-i\beta$. The selected coupling length is recorded in the run log and JSON
summary.

## Repository Layout

```text
ana/          Analytical sphere solutions used by validation.
docs/         Case format, mathematical implementation notes, and testing guide.
examples/     Public ordinary and Burton-Miller TOML cases and meshes.
fortran/      Modern Fortran solver core, public driver, and validation apps.
python/       Python case runner, mesh audit, post-processing, and tests.
scripts/      Developer and CI validation commands.
validation/   Reproducible validation cases.
```

For the Burton-Miller equations and their implementation, see
[`docs/burton-miller-equation-map.md`](docs/burton-miller-equation-map.md).

## Citation

If you use this software, cite the release metadata in
[`CITATION.cff`](CITATION.cff). If you use the non-singular Burton--Miller
formulation, also cite its method paper:

> Q. Sun and E. Klaseboer, "Non-Singular Burton-Miller Boundary Element Method
> for Acoustics," *Fluids* 8(2), 56 (2023).
> <https://doi.org/10.3390/fluids8020056>

The *Fluids* article is set as `preferred-citation` so that GitHub's
**Cite this repository** menu exposes it directly. The repository metadata and
other method references are also recorded in `CITATION.cff`.

## External Alpha Testing

The alpha phase checks whether independent users can install and run the
documented default quadratic workflow on fresh macOS and Ubuntu systems. The
scope, commands, expected outputs, and reporting form are described in
[`docs/external-alpha-testing.md`](docs/external-alpha-testing.md). External
usability testing complements, but does not replace, the analytical validation
matrix.

## Contributing And Support

Before proposing a change, read [`CONTRIBUTING.md`](CONTRIBUTING.md). Use the
GitHub issue forms for reproducible bugs and external alpha reports. The
supported scope and the information needed for technical questions are listed
in [`SUPPORT.md`](SUPPORT.md); release history is recorded in
[`CHANGELOG.md`](CHANGELOG.md).

## License

BSD-3-Clause. See `LICENSE`.
