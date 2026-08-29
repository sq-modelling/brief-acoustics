"""Public Python interface for BRIEF-Acoustics.

The public names exported here are the small Python-side API for the current
release. They cover the case contract, end-to-end source-tree runner, and mesh
input/output. The heavy boundary-integral assembly and linear solve remain in
Fortran.
"""

from .case_config import CaseConfig, CaseConfigurationError, MeshConfig, load_case
from .mesh_io import (
    DEFAULT_ELEMENT_TYPE,
    DEFAULT_NODES_PER_ELEMENT,
    MeshElement,
    MeshNode,
    SurfaceMesh,
    element_type_for,
    read_fortran_mesh,
    tetrahedron_mesh,
    write_fortran_mesh,
)
from .runner import CaseRunError, RunArtifacts, check_case, run_case

__all__ = [
    # Mesh data structures and text-file I/O.
    "DEFAULT_ELEMENT_TYPE",
    "DEFAULT_NODES_PER_ELEMENT",
    "MeshElement",
    "MeshNode",
    "SurfaceMesh",
    "element_type_for",
    "read_fortran_mesh",
    "tetrahedron_mesh",
    "write_fortran_mesh",
    # Public TOML case parsing and end-to-end execution.
    "CaseConfig",
    "CaseConfigurationError",
    "MeshConfig",
    "CaseRunError",
    "RunArtifacts",
    "check_case",
    "load_case",
    "run_case",
]
