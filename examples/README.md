# Examples

Run these commands from any working directory after installing the package in
editable mode from the repository root.

Minimal ordinary rigid-sphere case using the default quadratic elements:

```sh
brief-acoustics run examples/minimal_rigid_sphere/case.toml
```

Quadratic Burton-Miller case at the published sphere fictitious frequency
`ka=9.3558`:

```sh
brief-acoustics run examples/burton_miller_fictitious_sphere/case.toml
```

Each case writes `case.toml`, a normalized `mesh.dat`, the generated private
namelist, `surface_solution.csv`, `surface_magnitude.png`, `solver.log`, and
`summary.json` under its own `results/` directory.

The meshes are fixed copies of the validation meshes:

```text
meshes/sphere_linear_42.dat       42 nodes, 80 linear triangles
meshes/sphere_quadratic_162.dat  162 nodes, 80 quadratic triangles
meshes/sphere_quadratic_642.dat  642 nodes, 320 quadratic triangles
```
