module Geom_Types
    ! This module defines the data structures used to describe geometry.
    !
    ! In this code, a "geometry" means the surface mesh plus all geometric
    ! quantities derived from that mesh: normals, tangents, areas, element
    ! connectivity, neighbour lists, and quadrature points.  A few fields are
    ! reserved for later geometry extensions; their comments say explicitly
    ! when the current AU solver does not populate them.
    !
    ! The solver passes one geometry_type object around instead of passing many
    ! separate arrays.  This makes it much clearer which arrays belong together.
    use Pre_Constants, only: dp
    implicit none

    ! By default, names declared in this module are hidden from other modules.
    ! The public list below explicitly chooses which constants and types are
    ! part of the module interface.
    private

    public :: linear_triangle_node_count
    public :: quadratic_triangle_node_count
    public :: particle_geometry_type
    public :: surface_mesh_type
    public :: surface_differential_type
    public :: quadrature_cache_type
    public :: geometry_type

    ! Number of nodes used to describe one triangular element.
    !
    ! A linear triangle has only the three corner nodes.  Its sides are straight.
    integer, parameter :: linear_triangle_node_count = 3

    ! A quadratic triangle has three corner nodes and three midside nodes.
    ! The midside nodes allow the element to represent curved surfaces better.
    integer, parameter :: quadratic_triangle_node_count = 6

    ! New user cases default to quadratic triangles.  A mesh file still states
    ! its actual node count explicitly, and Geom_ReadMesh overwrites this
    ! placeholder after validating the file header.
    integer, parameter :: default_triangle_node_count = quadratic_triangle_node_count

    type :: particle_geometry_type
        ! Range of nodes and elements that belong to one closed surface.
        !
        ! A "particle" here means one connected boundary object, for example
        ! one sphere.  A multi-particle problem can contain several such objects.
        !
        ! The first/last values are array indices.  For example, if first_node=5
        ! and last_node=12, then nodes 5, 6, ..., 12 belong to this particle.
        integer :: first_node = 0
        integer :: last_node = 0
        integer :: first_element = 0
        integer :: last_element = 0

        ! Reserved particle summary fields.  The current AU surface preparation
        ! does not calculate either value, so both remain zero unless a caller
        ! fills them.  They must not be used as mesh-quality checks in v1.0.
        real(dp) :: surface_area = 0.0_dp
        real(dp) :: volume = 0.0_dp
    end type particle_geometry_type

    type :: surface_mesh_type
        ! nodes_per_element is the number of nodes per element.
        ! It is 3 for linear triangles and 6 for quadratic triangles.
        integer :: nodes_per_element = default_triangle_node_count

        ! Total number of nodes and triangular elements in the mesh.
        integer :: node_count = 0
        integer :: element_count = 0

        ! xyz(:, i) stores the Cartesian coordinates of node i.
        !
        ! The first index is the coordinate direction:
        !   xyz(1, i) = x coordinate of node i
        !   xyz(2, i) = y coordinate of node i
        !   xyz(3, i) = z coordinate of node i
        !
        ! The array is allocatable because its size is only known after reading
        ! the mesh file.
        real(dp), allocatable :: xyz(:, :)

        ! element_nodes(:, e) stores node ids for triangle e.
        ! Size the first dimension to nodes_per_element: 3 for linear, 6 for quadratic.
        !
        ! The required quadratic order is:
        !   1, 2, 3 = corner nodes in oriented boundary order
        !   4       = midside node on edge 1-2
        !   5       = midside node on edge 2-3
        !   6       = midside node on edge 3-1
        ! Geom_ReadMesh assumes this order; Geom_SurfaceGeometry uses it in the
        ! quadratic shape functions and normal calculation.
        integer, allocatable :: element_nodes(:, :)

        ! These arrays map each node and element back to the particle it belongs
        ! to.  For a single solid sphere, all entries are usually 1.
        integer, allocatable :: node_particle_id(:)
        integer, allocatable :: element_particle_id(:)

        ! These values store the largest list lengths needed by the connectivity
        ! arrays below.  Fortran rectangular arrays need fixed dimensions, so we
        ! store shorter lists with unused entries left as zero.
        integer :: max_node_elements = 0
        integer :: max_node_neighbors = 0

        ! node_elements(:, i) stores the elements that contain node i.
        ! node_neighbors(:, i) stores the first-ring neighboring nodes of node i.
        !
        ! Here "first-ring" has one precise meaning: two nodes occur in the same
        ! 3-node subtriangle built by Geom_MeshTopology.  For a quadratic parent
        ! triangle, the four subtriangles include its midside nodes.
        integer, allocatable :: node_elements(:, :)
        integer, allocatable :: node_neighbors(:, :)

        ! node_neighbors_2ring(:, i) stores first- and second-ring neighboring
        ! nodes of node i. This is used by local surface derivative stencils.
        !
        ! "Second-ring" means neighbours of the first-ring neighbours.  This
        ! gives a wider local patch of the surface around each node.
        integer :: max_node_neighbors_2ring = 0
        integer, allocatable :: node_neighbors_2ring(:, :)

        ! subtriangle_nodes stores a 3-node linear view of the mesh topology.
        ! Linear triangles are copied directly. Each quadratic triangle is split
        ! into four linear subtriangles using its corner and midside nodes.
        !
        ! This is convenient for calculations that only need simple flat
        ! triangles, such as estimating element areas or drawing/debugging the
        ! mesh.
        integer :: subtriangle_count = 0
        integer, allocatable :: subtriangle_nodes(:, :)
    end type surface_mesh_type

    type :: surface_differential_type
        ! Node-based local frame.
        !
        ! At each surface node we define three unit vectors:
        !   normal    = points perpendicular to the local surface
        !   tangent_1 = lies along the surface
        !   tangent_2 = lies along the surface and perpendicular to tangent_1
        !
        ! These local directions allow the solver to take derivatives along the
        ! surface rather than in the global x/y/z directions.
        real(dp), allocatable :: normal(:, :)
        real(dp), allocatable :: tangent_1(:, :)
        real(dp), allocatable :: tangent_2(:, :)

        ! Element-based area and normal.
        !
        ! These are stored per element rather than per node.  They are useful for
        ! integration over the surface and for checking element orientation.
        real(dp), allocatable :: element_area(:)
        real(dp), allocatable :: element_normal(:, :)

        ! Reserved curvature storage at each node.  Curvature would measure how
        ! strongly the surface bends, but the current AU preparation routine
        ! does not estimate it: all three arrays are allocated and set to zero.
        ! No v1.0 solver equation reads these placeholders.
        real(dp), allocatable :: curvature_1(:)
        real(dp), allocatable :: curvature_2(:)
        real(dp), allocatable :: mean_curvature(:)

        ! Surface derivative stencils at each node.
        ! Row 1 is the host node weight. Rows 2: align with
        ! mesh%node_neighbors_2ring(:, node).
        !
        ! A stencil is a small set of weights used to approximate a derivative
        ! from nearby node values.  d_dt1 and d_dt2 approximate derivatives in
        ! the two tangent directions.
        real(dp), allocatable :: d_dt1(:, :)
        real(dp), allocatable :: d_dt2(:, :)
    end type surface_differential_type

    type :: quadrature_cache_type
        ! Quadrature is numerical integration.
        !
        ! Boundary element methods integrate functions over every surface
        ! element.  A computer approximates those integrals by evaluating the
        ! function at selected quadrature points and multiplying by weights.
        !
        ! This type stores quadrature data after it has been generated, so the
        ! solver does not need to recompute point locations and weights every
        ! time it assembles a matrix entry.
        integer :: points_per_element = 0
        integer :: point_count = 0

        ! point_xyz(:, q) and point_normal(:, q) describe quadrature point q.
        !
        ! point_xyz gives the 3D position of the integration point.
        ! point_normal gives the surface normal at that point.
        real(dp), allocatable :: point_xyz(:, :)
        real(dp), allocatable :: point_normal(:, :)

        ! shape_weight(:, q) stores weighted shape functions at point q.
        !
        ! Shape functions describe how nodal values are interpolated across an
        ! element.  Multiplying them by the integration weight here saves work
        ! during matrix assembly.
        real(dp), allocatable :: shape_weight(:, :)

        ! Final scalar integration weight for each quadrature point.
        real(dp), allocatable :: integration_weight(:)
    end type quadrature_cache_type

    type :: geometry_type
        ! Top-level geometry container used by the solver.
        !
        ! This gathers all geometry-related information into one object:
        !   particles    = particle-level index ranges and reserved summaries
        !   mesh         = nodes, elements, and connectivity
        !   differential = normals, tangents, derivative stencils, and reserved
        !                  curvature arrays
        !   quadrature   = integration points and weights
        integer :: particle_count = 0
        type(particle_geometry_type), allocatable :: particles(:)
        type(surface_mesh_type) :: mesh
        type(surface_differential_type) :: differential
        type(quadrature_cache_type) :: quadrature
    contains
        ! Type-bound procedure.  This lets the code call geometry%clear().
        procedure :: clear => clear_geometry
    end type geometry_type

contains

    subroutine clear_geometry(this)
        ! Release all dynamic arrays owned by a geometry_type object and reset
        ! counters to their default empty values.
        !
        ! Fortran allocatable arrays occupy memory only after they are allocated.
        ! Deallocating them here is important when a program reads another mesh
        ! or exits a workflow cleanly.
        class(geometry_type), intent(inout) :: this

        ! Particle-level data.
        if (allocated(this%particles)) deallocate(this%particles)

        ! Mesh coordinates and connectivity.
        if (allocated(this%mesh%xyz)) deallocate(this%mesh%xyz)
        if (allocated(this%mesh%element_nodes)) deallocate(this%mesh%element_nodes)
        if (allocated(this%mesh%node_particle_id)) deallocate(this%mesh%node_particle_id)
        if (allocated(this%mesh%element_particle_id)) deallocate(this%mesh%element_particle_id)
        if (allocated(this%mesh%node_elements)) deallocate(this%mesh%node_elements)
        if (allocated(this%mesh%node_neighbors)) deallocate(this%mesh%node_neighbors)
        if (allocated(this%mesh%node_neighbors_2ring)) deallocate(this%mesh%node_neighbors_2ring)
        if (allocated(this%mesh%subtriangle_nodes)) deallocate(this%mesh%subtriangle_nodes)

        ! Surface normals, tangents, curvature, and derivative stencils.
        if (allocated(this%differential%normal)) deallocate(this%differential%normal)
        if (allocated(this%differential%tangent_1)) deallocate(this%differential%tangent_1)
        if (allocated(this%differential%tangent_2)) deallocate(this%differential%tangent_2)
        if (allocated(this%differential%element_area)) deallocate(this%differential%element_area)
        if (allocated(this%differential%element_normal)) deallocate(this%differential%element_normal)
        if (allocated(this%differential%curvature_1)) deallocate(this%differential%curvature_1)
        if (allocated(this%differential%curvature_2)) deallocate(this%differential%curvature_2)
        if (allocated(this%differential%mean_curvature)) deallocate(this%differential%mean_curvature)
        if (allocated(this%differential%d_dt1)) deallocate(this%differential%d_dt1)
        if (allocated(this%differential%d_dt2)) deallocate(this%differential%d_dt2)

        ! Cached quadrature point data.
        if (allocated(this%quadrature%point_xyz)) deallocate(this%quadrature%point_xyz)
        if (allocated(this%quadrature%point_normal)) deallocate(this%quadrature%point_normal)
        if (allocated(this%quadrature%shape_weight)) deallocate(this%quadrature%shape_weight)
        if (allocated(this%quadrature%integration_weight)) deallocate(this%quadrature%integration_weight)

        ! Reset scalar counters so the object represents an empty geometry.
        this%particle_count = 0
        this%mesh%nodes_per_element = default_triangle_node_count
        this%mesh%node_count = 0
        this%mesh%element_count = 0
        this%mesh%max_node_elements = 0
        this%mesh%max_node_neighbors = 0
        this%mesh%max_node_neighbors_2ring = 0
        this%mesh%subtriangle_count = 0
        this%quadrature%points_per_element = 0
        this%quadrature%point_count = 0
    end subroutine clear_geometry

end module Geom_Types
