# Rigid-Sphere Fictitious-Frequency Comparison

This workflow compares ordinary NSBEM and the non-singular Burton-Miller system
near two theoretical sphere resonances reported by:

Q. Sun and E. Klaseboer, "Non-Singular Burton-Miller Boundary Element Method
for Acoustics," *Fluids* 8 (2023), 56,
<https://doi.org/10.3390/fluids8020056>.

The published values used here are `ka=9.3558` and `ka=10.41712`.  Both methods
solve the same rigid Neumann sphere on the same 642-node, 320-element quadratic
mesh.  Evidence includes:

- surface relative L2 error against the partial-wave solution;
- scattered potential at `r=1.5a` along the incident-wave direction;
- ordinary-to-Burton-Miller improvement at each resonance;
- per-frequency Fortran surface CSV files and reproducible figures.

Run the 18-point focused sweep from the `modern` directory:

```sh
python3 validation/fictitious_frequency_sphere/run_sweep.py
```

Run only the two exact resonance points for a quick regression:

```sh
NSBEM_SWEEP_PROFILE=regression \
  NSBEM_RESULTS_NAME=results_regression \
  python3 validation/fictitious_frequency_sphere/run_sweep.py
```

The default focused sweep is intentionally not presented as a full reproduction
of the paper's `ka=1..40`, step-`0.001`, 5762-node calculation.  Its purpose is
to provide a low-cost, release-level regression of the same physical failure
mechanism and Burton-Miller remedy.

The fixed baseline values at the exact published frequencies are:

```text
ka = 9.3558
  ordinary surface relative L2 error = 5.03276827e-01
  Burton-Miller surface relative L2 error = 5.02601981e-02
  improvement factor = 1.00134271e+01

ka = 10.41712
  ordinary surface relative L2 error = 7.57575174e-01
  Burton-Miller surface relative L2 error = 7.22217538e-02
  improvement factor = 1.04895704e+01
```

The focused figures use separate local frequency panels so no line is drawn
through the unsampled interval between the two resonance windows.
