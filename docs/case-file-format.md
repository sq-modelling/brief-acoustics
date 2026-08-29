# BRIEF-Acoustics Case File and Mesh Contract

This document defines the public source-tree workflow.
The supported command is:

```sh
brief-acoustics run case.toml
```

Python validates the TOML file and complete mesh, derives non-redundant physics
quantities, builds the Fortran executable when needed, runs the dense solver,
and writes portable results. The generated `runtime_case.nml` is a private
Python/Fortran interface and should not be edited by users.

## 1. Current Scope

The case runner supports:

- one closed, connected exterior particle;
- real positive acoustic wavenumber;
- ordinary NSBEM or the non-singular Burton-Miller formulation;
- linear three-node or quadratic six-node triangular elements;
- constant total-field Dirichlet, Neumann, or Robin boundary data;
- no incident field, or one travelling/standing plane wave.

It deliberately rejects multiple particles, nested layers, transmission,
complex wavenumbers, and internal-domain cases. Those data structures exist in
the research core but are not part of the validated public 1.0 contract.

## 2. Complete TOML Example

```toml
[case]
name = "rigid-sphere"

[mesh]
path = "sphere_mesh.dat"
element_type = "quadratic"                 # optional; this is the default
normal_orientation = "outward-from-solid"

[physics]
density = 1.225
sound_speed = 343.0
frequency_hz = 1000.0

[solver]
formulation = "burton-miller"
characteristic_length = 0.1

[boundary]
type = "neumann"
value_real = 0.0
value_imag = 0.0

[incident]
type = "plane-wave"
potential_amplitude_real = 1.0
potential_amplitude_imag = 0.0
direction = [1.0, 0.0, 0.0]
standing_wave = false

[output]
directory = "results"
```

Unknown or misspelled keys are errors. Relative mesh and output paths are
resolved from the directory containing `case.toml`, not from the current shell
working directory.

## 3. Case and Mesh Tables

### `[case]`

`name` is required. It may contain letters, digits, `.`, `_`, and `-`, and must
start with a letter or digit. The name appears in output metadata and figures.

### `[mesh]`

`path` and `normal_orientation` are required. `element_type` is optional:

```toml
[mesh]
path = "sphere_mesh.dat"
normal_orientation = "outward-from-solid"
# element_type = "linear"                  # explicit opt-in override
```

`path` points to the plain-text mesh described in Section 9. If
`element_type` is omitted, BRIEF-Acoustics selects quadratic six-node triangles. An
explicit value uses the mathematical interpolation family, `"linear"` or
`"quadratic"`, instead of the connectivity widths `3` and `6`. Python verifies
that the resolved element type agrees with the mesh connectivity. The current
exterior solver accepts only a consistently oriented mesh whose normals point
outward from the solid into the acoustic fluid.

Python audits the source mesh and writes a normalized copy named `mesh.dat`
into the result directory for Fortran. From `normal_orientation`, it derives
the internal exterior free-term sign. Users do not set that equation sign.

## 4. Physics Table

Required fields:

```toml
density = 1.225
sound_speed = 343.0
```

Both values must be real, finite, and positive. Define exactly one frequency
quantity:

```toml
frequency_hz = 1000.0
```

or

```toml
angular_frequency = 6283.185307179586
```

or

```toml
wavenumber = 18.31832451
```

Python derives the remaining quantities using

```text
omega = 2 pi f
k = omega / c
```

This avoids independently supplied `f`, `omega`, and `k` drifting out of
agreement.

Geometry, sound speed, and wavenumber must use one consistent unit system. For
example, metres and seconds give `c` in m/s and `k` in 1/m. Dimensionless
validation cases may use `c=1`, but should be labelled as dimensionless.

## 5. Solver Table

```toml
[solver]
formulation = "ordinary"       # or "burton-miller"
characteristic_length = 0.1
```

`characteristic_length` is a positive representative body scale in the same
length unit as the mesh. For a sphere it is normally the radius.

The corresponding internal convention is:

```text
mesh normal: points outward from the solid into the acoustic fluid
exterior-domain normal: points in the opposite direction
derived exterior free-term sign: -1
```

The case parser rejects any unsupported `normal_orientation` instead of asking
users to encode this sign convention numerically.

Use `ordinary` for controlled low-frequency work away from fictitious
frequencies. Use `burton-miller` for robust exterior frequency sweeps and when
non-uniqueness may be encountered.

## 6. Boundary Table

Boundary data applies to the **total field**, while Fortran solves and stores
the scattered field separately.

Dirichlet condition:

```toml
[boundary]
type = "dirichlet"
value_real = 0.0
value_imag = 0.0
```

This imposes

```text
phi_total = value.
```

Neumann condition:

```toml
[boundary]
type = "neumann"
value_real = 0.0
value_imag = 0.0
```

This imposes

```text
dphi_total/dn_mesh = value,
```

where `n_mesh` is the stored outward-from-solid normal. A rigid sound-hard
surface uses zero.

Robin condition:

```toml
[boundary]
type = "robin"
a_real = 1.0
a_imag = 0.0
b_real = 0.5
b_imag = 0.0
rhs_real = 0.0
rhs_imag = 0.0
```

This imposes

```text
a phi_total + b dphi_total/dn_mesh = rhs.
```

`a` and `b` cannot both be zero. Their units must make both terms compatible;
the program does not infer physical Robin scaling.

## 7. Incident Field Table

No incident field, for a radiation problem:

```toml
[incident]
type = "none"
```

Plane wave:

```toml
[incident]
type = "plane-wave"
potential_amplitude_real = 1.0
potential_amplitude_imag = 0.0
direction = [1.0, 0.0, 0.0]
standing_wave = false
```

The direction is normalized by Fortran and must not be zero. The travelling
wave is

```text
phi_inc(x) = potential_amplitude exp(i k direction.x).
```

With `standing_wave=true`, equal forward and backward waves are averaged.

The code uses the physical time convention

```text
exp(-i omega t),
```

the outgoing fundamental solution `exp(i k r)/r`, and

```text
pressure = i omega density phi.
```

`potential_amplitude` is a velocity-potential amplitude, not a pressure
amplitude. For a desired pressure amplitude, convert using the pressure
relation above.

## 8. Output Table and Files

```toml
[output]
directory = "results"
```

The table is optional. Without it, output defaults to `results/<case-name>`
relative to the case file. `--output PATH` overrides the TOML value.

A completed run contains:

```text
case.toml               copy of the public input
mesh.dat                normalized, audited mesh consumed by Fortran
runtime_case.nml        generated private Python/Fortran input
surface_solution.csv    nodal incident, scattered, total, pressure, and BM data
surface_magnitude.png   total-potential and pressure magnitudes
solver.log              build/Fortran completion output
summary.json            resolved physics, solver settings, mesh, and extrema
```

`summary.json` reports `"status": "completed"`. This means the numerical run
finished and its output is finite; it is not an analytical accuracy claim.
Accuracy claims come from the dedicated `validation/` workflows.

## 9. Mesh Format

The first record is:

```text
node_count element_count nodes_per_element particle_count
```

Each node record is:

```text
node_id x y z particle_id
```

The third header value is a connectivity width: `3` for linear elements and `6`
for quadratic elements. `nodes_per_element` describes that count directly; the
corresponding polynomial degrees are one and two. Python checks it against
`[mesh].element_type`.

Each linear element record is:

```text
element_id particle_id node_1 node_2 node_3
```

Each quadratic element record is:

```text
element_id particle_id node_1 node_2 node_3 node_4 node_5 node_6
```

Quadratic ordering is fixed:

```text
node_1, node_2, node_3: corner nodes
node_4: midside node on edge 1-2
node_5: midside node on edge 2-3
node_6: midside node on edge 3-1
```

Corner order follows the right-hand rule and must produce a normal pointing
outward from the solid. Every particle must be one connected, closed manifold.
Every corner edge must occur exactly twice and in opposite directions in its
two adjacent elements. Adjacent quadratic elements must share the same midside
node on their common edge.

The Python audit also rejects non-finite coordinates, non-sequential ids,
unknown or repeated node references, degenerate and duplicate triangles,
cross-particle references, open/non-manifold edges, inconsistent orientation,
unreferenced nodes, and disconnected components assigned one particle id.

Blank lines and lines beginning with `#` are accepted by the Python reader.
Fortran receives the normalized comment-free copy.

## 10. Commands

Install once from the repository root:

```sh
python3 -m pip install -e .
```

Check a case and mesh without compiling or solving:

```sh
brief-acoustics check case.toml
```

Run it:

```sh
brief-acoustics run case.toml
```

The command builds `fortran/bin/run_au_case` automatically. Use `--no-build`
only when that executable is already current. Use `--no-plot` for CSV/JSON-only
batch work.

Frequency and parameter sweeps belong in Python. A sweep should create or vary
single-case inputs and call the same runner; start/end/step controls are not
part of the Fortran physics core.

## 11. Common Errors

- **Open/non-manifold edge:** the input is not one watertight triangular shell.
- **Inconsistent orientation:** reverse the affected triangle corner ordering;
  do not mislabel `normal_orientation` to compensate for bad connectivity.
- **Quadratic midside disagreement:** neighboring elements do not use the same
  node for their shared corner edge.
- **More than one frequency field:** retain only one of `frequency_hz`,
  `angular_frequency`, or `wavenumber`.
- **More than one particle:** use the separate two-sphere demonstration for the
  current ordinary multi-particle research path; the public runner rejects it.
- **Missing Fortran source tree:** run from an editable source installation or
  set `AU_NSBEM_FORTRAN_DIR` to the repository's `fortran` directory.
- **Dense solve is too large:** node count controls memory and cubic direct-solve
  cost. Refine deliberately; this release does not yet include FMM/H-matrix or an
  iterative fast solver.
