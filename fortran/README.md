# Modern Fortran Core

This folder contains the Fortran solver modules used by BRIEF-Acoustics. These
files are the active numerical-core source for the public alpha. The internal
`AU_*` prefix denotes acoustics and is retained as a stable implementation
namespace.

## Module Order

Compile in dependency order:

1. `Pre_Constants.f90`
2. `Num_Quadrature.f90`
3. `Geom_Types.f90`
4. `Geom_ReadMesh.f90`
5. `Geom_MeshTopology.f90`
6. `Geom_SurfaceGeometry.f90`
7. `AU_Types.f90`
8. `AU_LayerTopology.f90`
9. `AU_HelmholtzKernels.f90`
10. `AU_BoundaryConditions.f90`
11. `AU_SurfaceOperators.f90`
12. `AU_Solver.f90`

## Applications

`app/RunAUCase.f90` is the generic single-particle exterior driver used by
`brief-acoustics run case.toml`. Python validates the public TOML and mesh,
then passes a normalized mesh and private generated namelist to this executable.

Build it with:

```sh
make run-case
```

The other programs in `app/` and `test/` support analytical regressions,
manufactured solutions, the explicitly labelled two-sphere demonstration, and
matrix-level Burton-Miller checks.
