module AU_Types
    ! This module defines the data structures for one acoustic problem.
    !
    ! "AU" means acoustics in this release.  The solver works in the frequency
    ! domain, so acoustic quantities are complex numbers.  The complex value
    ! stores both amplitude and phase of a time-harmonic wave.
    !
    ! The geometry module knows where the surface mesh is.  This module knows
    ! the physics to solve on that mesh:
    !   - material properties such as density and sound speed
    !   - boundary condition type for each particle
    !   - incident field values at the mesh nodes
    !   - solution arrays filled by the solver
    use Pre_Constants, only: dp, complex_zero
    use Geom_Types, only: geometry_type
    implicit none

    ! Hide helper routines and constants unless they are explicitly made public.
    private

    ! Boundary-condition identifiers.
    !
    ! The solver stores boundary-condition types as integers because they are
    ! compact and easy to use in select/case logic.  Python or a future GUI can
    ! translate user-friendly labels into these integer constants.
    public :: bc_dirichlet_external
    public :: bc_neumann_external
    public :: bc_robin_external
    public :: bc_dirichlet_internal
    public :: bc_neumann_internal
    public :: bc_robin_internal
    public :: bc_transmission
    public :: au_medium_type
    public :: au_particle_layer_type
    public :: au_boundary_condition_type
    public :: au_incident_field_type
    public :: au_surface_solution_type
    public :: au_case_type
    public :: allocate_au_case
    public :: clear_au_case
    public :: validate_au_case

    ! Dirichlet boundary condition: the potential phi is known.
    integer, parameter :: bc_dirichlet_external = 1

    ! Neumann boundary condition: the normal derivative dphi/dn is known.
    integer, parameter :: bc_neumann_external = 2

    ! Robin boundary condition: a linear combination of phi and dphi/dn is known.
    integer, parameter :: bc_robin_external = 3

    ! Internal-side identifiers are retained for the research solver path.  The
    ! public v1.0 case runner accepts external boundary conditions only.
    integer, parameter :: bc_dirichlet_internal = 4
    integer, parameter :: bc_neumann_internal = 5
    integer, parameter :: bc_robin_internal = 6

    ! Transmission means the field is coupled across an interface.  Both sides
    ! of the boundary matter, as in penetrable acoustic media.  The data value
    ! is reserved here, but AU_Solver deliberately rejects it in v1.0 because a
    ! transmission system needs two physical unknowns per interface node.
    integer, parameter :: bc_transmission = 7

    type :: au_medium_type
        ! Acoustic material data for one domain.
        !
        ! For the exterior domain, these values describe the surrounding fluid.
        ! For an interior domain, they describe the medium inside one particle.
        ! Density is stored as a real positive value.  Sound speed and
        ! wavenumber use complex storage so lossy-media extensions remain
        ! possible, but the public acoustic v1.0 path validates real positive
        ! sound speed and wavenumber.
        real(dp) :: density = 0.0_dp

        ! Sound speed c.  It is complex so lossy media can be represented later.
        complex(dp) :: sound_speed = complex_zero

        ! Acoustic wavenumber k = omega / c.
        ! Most Helmholtz kernels use k directly.
        complex(dp) :: wavenumber = complex_zero
    end type au_medium_type

    type :: au_particle_layer_type
        ! Topology information capable of describing nested acoustic domains.
        !
        ! A simple single solid sphere has one particle and parent_particle_id=0.
        ! More complex cases can have nested surfaces, where one particle is
        ! inside another.  The parent id records that nesting.
        ! A value of 0 means the outside of this particle is the infinite
        ! exterior medium.
        !
        ! This storage does not by itself imply validated multilayer support.
        ! The public v1.0 runner accepts one particle, and multilayer equations
        ! require dedicated derivation and tests before they can be released.
        integer :: parent_particle_id = 0

        ! Sign multiplying exterior-domain free terms in the boundary-integral
        ! equations. For an isolated solid whose stored mesh normals point
        ! outward into the fluid, this value is -1.  The value +1 below is only
        ! an unconfigured allocation default.  User-facing input describes the
        ! normal orientation in words, and Python derives and overwrites this
        ! internal sign before a public solve.
        integer :: exterior_free_term_sign = 1

        ! Representative particle length a, in the same unit as the mesh.
        !
        ! The normal-derivative Burton-Miller equation carries one extra
        ! inverse length compared with the ordinary boundary integral equation.
        ! Its coupling coefficient beta must therefore have units of length.
        ! The acoustic solver uses beta=min(a,1/k): the body scale a at low
        ! frequency and the wavelength-related scale 1/k at high frequency.
        ! For a sphere, a is normally its radius.
        real(dp) :: characteristic_length = 1.0_dp
    end type au_particle_layer_type

    type :: au_boundary_condition_type
        ! Boundary-condition data for one particle.
        !
        ! The public workflow stores boundary conditions at particle level first.  A helper
        ! later expands those values to every node on that particle.
        ! Python converts the public labels "dirichlet", "neumann", and
        ! "robin" into these internal integer ids.
        integer :: kind = bc_dirichlet_external

        ! Constant prescribed Dirichlet or Neumann data for the particle. If
        ! per-node data are needed, fill case%prescribed_boundary_data directly.
        complex(dp) :: prescribed_data = complex_zero

        ! Robin form: robin_a * phi + robin_b * dphi/dn = robin_rhs.
        complex(dp) :: robin_a = complex_zero
        complex(dp) :: robin_b = complex_zero
        complex(dp) :: robin_rhs = complex_zero
    end type au_boundary_condition_type

    type :: au_incident_field_type
        ! Prescribed incoming field values at every mesh node.
        !
        ! For scattering problems, the total field is:
        !   total field = incident field + scattered field
        !
        ! These arrays store the incident part.  The solver computes the
        ! scattered part.
        ! The values are intentionally precomputed outside the solver core so
        ! the solver does not mix user-input wave definitions with matrix assembly.
        complex(dp), allocatable :: external_phi(:)
        complex(dp), allocatable :: external_dphi_dn(:)
        complex(dp), allocatable :: internal_phi(:)
        complex(dp), allocatable :: internal_dphi_dn(:)
    end type au_incident_field_type

    type :: au_surface_solution_type
        ! Unknown and solved acoustic fields on the mesh nodes.
        !
        ! phi is the acoustic velocity potential.  dphi_dn means the derivative
        ! of phi in the mesh normal direction.
        ! The first group stores the solved field after incident-field terms
        ! have been removed from the boundary condition.  In an exterior
        ! scattering case this is the scattered field.  When no incident field
        ! is supplied, it is also the total solved field.
        complex(dp), allocatable :: external_phi(:)
        complex(dp), allocatable :: external_dphi_dn(:)
        complex(dp), allocatable :: internal_phi(:)
        complex(dp), allocatable :: internal_dphi_dn(:)

        ! Burton-Miller introduces an auxiliary Laplace field, usually written
        ! psi_1 in the derivation.  The augmented linear system solves for its
        ! nodal normal derivative s = d(psi_1)/dn.  This is a mathematical
        ! support variable: it is not acoustic pressure or particle velocity.
        ! We retain it so validation code can inspect both block equations.
        complex(dp), allocatable :: bm_auxiliary_dpsi_dn(:)

        ! Total field = incident field + scattered field.
        complex(dp), allocatable :: total_external_phi(:)
        complex(dp), allocatable :: total_external_dphi_dn(:)
        complex(dp), allocatable :: total_internal_phi(:)
        complex(dp), allocatable :: total_internal_dphi_dn(:)

        ! Acoustic pressure reconstructed from the velocity potential.
        complex(dp), allocatable :: external_pressure(:)
        complex(dp), allocatable :: internal_pressure(:)
    end type au_surface_solution_type

    type :: au_case_type
        ! Complete acoustic case for one prepared geometry.
        !
        ! Think of this as the main "problem setup" object.  It contains all
        ! physics input and all arrays where the solver stores output.

        ! Internal construction default.  Every public TOML case must state
        ! formulation="ordinary" or "burton-miller", and RunAUCase overwrites
        ! this value, so this initializer is not a user-facing solver default.
        logical :: use_burton_miller = .true.

        ! Angular frequency omega in radians per second.
        real(dp) :: angular_frequency = 0.0_dp

        ! Exterior medium is the infinite outside domain.
        type(au_medium_type) :: exterior_medium

        ! One interior medium, layer description, and boundary condition record
        ! per particle.
        type(au_medium_type), allocatable :: interior_medium(:)
        type(au_particle_layer_type), allocatable :: layer(:)
        type(au_boundary_condition_type), allocatable :: boundary_condition(:)

        ! Per-node known value after particle-level boundary conditions have
        ! been expanded.  A vibration/radiation caller may fill this directly
        ! or use set_normal_velocity_boundary to project a velocity vector onto
        ! the mesh normal.
        complex(dp), allocatable :: prescribed_boundary_data(:)

        type(au_incident_field_type) :: incident
        type(au_surface_solution_type) :: solution
    contains
        ! Allows object-style cleanup: call case_data%clear().
        procedure :: clear => clear_au_case_bound
    end type au_case_type

contains

    subroutine allocate_au_case(case_data, geometry, status, message)
        ! Allocate all acoustic arrays so their sizes match the geometry.
        !
        ! The number of particles decides how many particle-level records are
        ! needed.  The number of mesh nodes decides how many nodal field values
        ! are needed.
        type(au_case_type), intent(inout) :: case_data
        type(geometry_type), intent(in) :: geometry
        integer, intent(out), optional :: status
        character(len=*), intent(out), optional :: message

        integer :: error_code
        character(len=256) :: error_message
        integer :: node_count
        integer :: particle_count

        error_code = 0
        error_message = "AU case arrays allocated."

        ! Start from a clean case so old arrays cannot be mixed with the new
        ! geometry.
        call clear_au_case(case_data)

        node_count = geometry%mesh%node_count
        particle_count = geometry%particle_count

        ! A case cannot be allocated until a geometry has been read.
        if (node_count <= 0) then
            call finish(1, "Geometry has no nodes.")
            return
        end if

        if (particle_count <= 0) then
            call finish(1, "Geometry has no particles.")
            return
        end if

        ! Allocate particle-level arrays and node-level arrays in one operation.
        ! If any allocation fails, error_code becomes non-zero.
        allocate(case_data%interior_medium(particle_count), &
                 case_data%layer(particle_count), &
                 case_data%boundary_condition(particle_count), &
                 case_data%prescribed_boundary_data(node_count), &
                 case_data%incident%external_phi(node_count), &
                 case_data%incident%external_dphi_dn(node_count), &
                 case_data%incident%internal_phi(node_count), &
                 case_data%incident%internal_dphi_dn(node_count), &
                 case_data%solution%external_phi(node_count), &
                 case_data%solution%external_dphi_dn(node_count), &
                 case_data%solution%internal_phi(node_count), &
                 case_data%solution%internal_dphi_dn(node_count), &
                 case_data%solution%bm_auxiliary_dpsi_dn(node_count), &
                 case_data%solution%total_external_phi(node_count), &
                 case_data%solution%total_external_dphi_dn(node_count), &
                 case_data%solution%total_internal_phi(node_count), &
                 case_data%solution%total_internal_dphi_dn(node_count), &
                 case_data%solution%external_pressure(node_count), &
                 case_data%solution%internal_pressure(node_count), &
                 stat=error_code)

        if (error_code /= 0) then
            error_message = "Unable to allocate AU case arrays."
            call clear_au_case(case_data)
            call finish(error_code, error_message)
            return
        end if

        ! Set particle-level records to their default values.
        case_data%interior_medium = au_medium_type()
        case_data%layer = au_particle_layer_type()
        case_data%boundary_condition = au_boundary_condition_type()

        ! Initialise all nodal values to the complex zero 0 + 0i.
        case_data%prescribed_boundary_data = complex_zero

        case_data%incident%external_phi = complex_zero
        case_data%incident%external_dphi_dn = complex_zero
        case_data%incident%internal_phi = complex_zero
        case_data%incident%internal_dphi_dn = complex_zero

        case_data%solution%external_phi = complex_zero
        case_data%solution%external_dphi_dn = complex_zero
        case_data%solution%internal_phi = complex_zero
        case_data%solution%internal_dphi_dn = complex_zero
        case_data%solution%bm_auxiliary_dpsi_dn = complex_zero
        case_data%solution%total_external_phi = complex_zero
        case_data%solution%total_external_dphi_dn = complex_zero
        case_data%solution%total_internal_phi = complex_zero
        case_data%solution%total_internal_dphi_dn = complex_zero
        case_data%solution%external_pressure = complex_zero
        case_data%solution%internal_pressure = complex_zero

        call finish(error_code, error_message)

    contains

        subroutine finish(code, text)
            ! Return optional status information to the caller.
            integer, intent(in) :: code
            character(len=*), intent(in) :: text

            if (present(status)) status = code
            if (present(message)) message = text
        end subroutine finish

    end subroutine allocate_au_case

    subroutine validate_au_case(case_data, geometry, status, message)
        ! Check whether a case is complete enough for the solver to use.
        !
        ! This routine catches common setup mistakes early:
        !   - arrays not allocated
        !   - arrays allocated with the wrong size
        !   - invalid material density
        !   - unsupported boundary-condition kind
        !   - impossible layer nesting
        type(au_case_type), intent(in) :: case_data
        type(geometry_type), intent(in) :: geometry
        integer, intent(out), optional :: status
        character(len=*), intent(out), optional :: message

        integer :: code
        character(len=256) :: text
        integer :: particle_id

        code = 0
        text = "AU case is valid."

        ! Global case checks.  These are ordered so the first useful error
        ! message is returned.
        if (geometry%mesh%node_count <= 0) then
            call fail("Geometry has no nodes.")
        else if (geometry%particle_count <= 0) then
            call fail("Geometry has no particles.")
        else if (case_data%angular_frequency < 0.0_dp) then
            call fail("Angular frequency must be non-negative.")
        else if (.not. allocated(case_data%interior_medium)) then
            call fail("Interior media are not allocated.")
        else if (size(case_data%interior_medium) /= geometry%particle_count) then
            call fail("Interior media count does not match geometry particle count.")
        else if (.not. allocated(case_data%layer)) then
            call fail("Particle layers are not allocated.")
        else if (size(case_data%layer) /= geometry%particle_count) then
            call fail("Particle layer count does not match geometry particle count.")
        else if (.not. allocated(case_data%boundary_condition)) then
            call fail("Boundary conditions are not allocated.")
        else if (size(case_data%boundary_condition) /= geometry%particle_count) then
            call fail("Boundary condition count does not match geometry particle count.")
        else if (.not. allocated(case_data%prescribed_boundary_data)) then
            call fail("Prescribed boundary data are not allocated.")
        else if (size(case_data%prescribed_boundary_data) /= geometry%mesh%node_count) then
            call fail("Prescribed boundary data count does not match mesh node count.")
        else if (.not. has_node_sized_incident_arrays(case_data, geometry%mesh%node_count)) then
            call fail("Incident field arrays are not allocated with mesh node count.")
        else if (.not. has_node_sized_solution_arrays(case_data, geometry%mesh%node_count)) then
            call fail("Surface solution arrays are not allocated with mesh node count.")
        else if (case_data%exterior_medium%density <= 0.0_dp) then
            call fail("Exterior medium density must be positive.")
        end if

        if (code == 0) then
            ! Per-particle checks.
            do particle_id = 1, geometry%particle_count
                if (.not. is_valid_boundary_kind(case_data%boundary_condition(particle_id)%kind)) then
                    call fail("Boundary condition kind is not supported.")
                    exit
                end if
                if (case_data%interior_medium(particle_id)%density <= 0.0_dp) then
                    call fail("Interior medium density must be positive.")
                    exit
                end if
                if (.not. is_valid_layer(case_data, particle_id, geometry%particle_count)) then
                    call fail("Particle layer topology is not valid.")
                    exit
                end if
            end do
        end if

        if (present(status)) status = code
        if (present(message)) message = text

    contains

        subroutine fail(error_text)
            ! Store only the first validation error.  Later checks should not
            ! overwrite the message that explains the original problem.
            character(len=*), intent(in) :: error_text

            if (code == 0) then
                code = 1
                text = error_text
            end if
        end subroutine fail

    end subroutine validate_au_case

    subroutine clear_au_case_bound(this)
        ! Type-bound wrapper so an au_case_type object can call case%clear().
        class(au_case_type), intent(inout) :: this

        call clear_au_case(this)
    end subroutine clear_au_case_bound

    subroutine clear_au_case(case_data)
        ! Release all dynamic arrays in an acoustic case and reset scalar
        ! settings to their defaults.
        type(au_case_type), intent(inout) :: case_data

        ! Particle-level arrays.
        if (allocated(case_data%interior_medium)) deallocate(case_data%interior_medium)
        if (allocated(case_data%layer)) deallocate(case_data%layer)
        if (allocated(case_data%boundary_condition)) deallocate(case_data%boundary_condition)
        if (allocated(case_data%prescribed_boundary_data)) deallocate(case_data%prescribed_boundary_data)

        ! Incident-field arrays.
        if (allocated(case_data%incident%external_phi)) deallocate(case_data%incident%external_phi)
        if (allocated(case_data%incident%external_dphi_dn)) deallocate(case_data%incident%external_dphi_dn)
        if (allocated(case_data%incident%internal_phi)) deallocate(case_data%incident%internal_phi)
        if (allocated(case_data%incident%internal_dphi_dn)) deallocate(case_data%incident%internal_dphi_dn)

        ! Solution arrays.
        if (allocated(case_data%solution%external_phi)) deallocate(case_data%solution%external_phi)
        if (allocated(case_data%solution%external_dphi_dn)) deallocate(case_data%solution%external_dphi_dn)
        if (allocated(case_data%solution%internal_phi)) deallocate(case_data%solution%internal_phi)
        if (allocated(case_data%solution%internal_dphi_dn)) deallocate(case_data%solution%internal_dphi_dn)
        if (allocated(case_data%solution%bm_auxiliary_dpsi_dn)) then
            deallocate(case_data%solution%bm_auxiliary_dpsi_dn)
        end if
        if (allocated(case_data%solution%total_external_phi)) deallocate(case_data%solution%total_external_phi)
        if (allocated(case_data%solution%total_external_dphi_dn)) deallocate(case_data%solution%total_external_dphi_dn)
        if (allocated(case_data%solution%total_internal_phi)) deallocate(case_data%solution%total_internal_phi)
        if (allocated(case_data%solution%total_internal_dphi_dn)) deallocate(case_data%solution%total_internal_dphi_dn)
        if (allocated(case_data%solution%external_pressure)) deallocate(case_data%solution%external_pressure)
        if (allocated(case_data%solution%internal_pressure)) deallocate(case_data%solution%internal_pressure)

        ! Reset scalar fields to default values.
        case_data%use_burton_miller = .true.
        case_data%angular_frequency = 0.0_dp
        case_data%exterior_medium = au_medium_type()
    end subroutine clear_au_case

    pure logical function is_valid_boundary_kind(kind)
        ! Check whether an integer is one of the boundary-condition identifiers
        ! declared by this data module.  Solver-level validation separately
        ! rejects declared values that are outside the current release scope.
        integer, intent(in) :: kind

        is_valid_boundary_kind = kind >= bc_dirichlet_external .and. &
                                 kind <= bc_transmission
    end function is_valid_boundary_kind

    pure logical function is_valid_layer(case_data, particle_id, particle_count)
        ! Check that one particle's parent graph is structurally valid.
        !
        ! The parent chain must:
        !   - use valid particle ids
        !   - not point a particle to itself
        !   - not contain a cycle such as 1 inside 2 inside 1
        !
        ! This is a data-integrity check, not validation of multilayer acoustic
        ! equations.  The public v1.0 workflow does not expose nested layers.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: particle_id
        integer, intent(in) :: particle_count

        integer :: parent_id
        integer :: visited_count

        parent_id = case_data%layer(particle_id)%parent_particle_id

        ! Basic checks for this particle's immediate layer data.
        is_valid_layer = parent_id >= 0 .and. parent_id <= particle_count .and. &
                         parent_id /= particle_id .and. &
                         abs(case_data%layer(particle_id)%exterior_free_term_sign) == 1 .and. &
                         case_data%layer(particle_id)%characteristic_length > 0.0_dp
        if (.not. is_valid_layer) return

        visited_count = 0
        do while (parent_id > 0)
            ! If the chain is longer than the number of particles, it must have
            ! repeated a particle id, which means there is a cycle.
            visited_count = visited_count + 1
            if (visited_count > particle_count) then
                is_valid_layer = .false.
                return
            end if

            if (parent_id == particle_id) then
                is_valid_layer = .false.
                return
            end if

            if (case_data%layer(parent_id)%parent_particle_id < 0 .or. &
                case_data%layer(parent_id)%parent_particle_id > particle_count) then
                is_valid_layer = .false.
                return
            end if

            ! Move one level outward in the nesting tree.
            parent_id = case_data%layer(parent_id)%parent_particle_id
        end do
    end function is_valid_layer

    pure logical function has_node_sized_incident_arrays(case_data, node_count)
        ! Check that all incident-field arrays exist and have one value per node.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: node_count

        ! First check allocation.  size(...) must not be called on an
        ! unallocated array.
        has_node_sized_incident_arrays = &
            allocated(case_data%incident%external_phi) .and. &
            allocated(case_data%incident%external_dphi_dn) .and. &
            allocated(case_data%incident%internal_phi) .and. &
            allocated(case_data%incident%internal_dphi_dn)

        if (.not. has_node_sized_incident_arrays) return

        ! Then check dimensions.
        has_node_sized_incident_arrays = &
            size(case_data%incident%external_phi) == node_count .and. &
            size(case_data%incident%external_dphi_dn) == node_count .and. &
            size(case_data%incident%internal_phi) == node_count .and. &
            size(case_data%incident%internal_dphi_dn) == node_count
    end function has_node_sized_incident_arrays

    pure logical function has_node_sized_solution_arrays(case_data, node_count)
        ! Check that all solution arrays exist and have one value per node.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: node_count

        ! Allocation check first.
        has_node_sized_solution_arrays = &
            allocated(case_data%solution%external_phi) .and. &
            allocated(case_data%solution%external_dphi_dn) .and. &
            allocated(case_data%solution%internal_phi) .and. &
            allocated(case_data%solution%internal_dphi_dn) .and. &
            allocated(case_data%solution%bm_auxiliary_dpsi_dn) .and. &
            allocated(case_data%solution%total_external_phi) .and. &
            allocated(case_data%solution%total_external_dphi_dn) .and. &
            allocated(case_data%solution%total_internal_phi) .and. &
            allocated(case_data%solution%total_internal_dphi_dn) .and. &
            allocated(case_data%solution%external_pressure) .and. &
            allocated(case_data%solution%internal_pressure)

        if (.not. has_node_sized_solution_arrays) return

        ! Dimension check second.
        has_node_sized_solution_arrays = &
            size(case_data%solution%external_phi) == node_count .and. &
            size(case_data%solution%external_dphi_dn) == node_count .and. &
            size(case_data%solution%internal_phi) == node_count .and. &
            size(case_data%solution%internal_dphi_dn) == node_count .and. &
            size(case_data%solution%bm_auxiliary_dpsi_dn) == node_count .and. &
            size(case_data%solution%total_external_phi) == node_count .and. &
            size(case_data%solution%total_external_dphi_dn) == node_count .and. &
            size(case_data%solution%total_internal_phi) == node_count .and. &
            size(case_data%solution%total_internal_dphi_dn) == node_count .and. &
            size(case_data%solution%external_pressure) == node_count .and. &
            size(case_data%solution%internal_pressure) == node_count
    end function has_node_sized_solution_arrays

end module AU_Types
