module AU_BoundaryConditions
    ! This module prepares acoustic boundary-condition and incident-field data.
    !
    ! The solver itself works with arrays defined at every mesh node.  Users,
    ! however, usually describe a problem in higher-level terms:
    !   - this particle has a constant boundary value
    !   - this particle vibrates with a given velocity
    !   - the incoming wave is a plane wave
    !
    ! The routines here convert those higher-level descriptions into nodal
    ! values that the boundary-integral matrix assembly can use directly.
    ! The release uses the time convention exp(-i*omega*t), so an outgoing
    ! spatial wave is exp(i*k*r) and pressure is p=i*omega*rho*phi.
    use Pre_Constants, only: dp, complex_zero, imaginary_unit
    use Geom_Types, only: geometry_type
    use AU_Types, only: au_case_type, bc_robin_external, bc_robin_internal
    implicit none

    private

    public :: au_plane_wave_type
    public :: reset_au_fields
    public :: expand_particle_boundary_values
    public :: set_normal_velocity_boundary
    public :: set_plane_wave_incident_field
    public :: update_total_fields_and_pressure

    ! Public input tolerance used when checking that a direction vector is not
    ! effectively zero.  It matches the Python case validation threshold.
    real(dp), parameter :: tiny_direction = 1.0e-14_dp

    type :: au_plane_wave_type
        ! Compact description of an incoming acoustic plane wave.
        !
        ! A plane wave has the form:
        !   phi_inc(x) = potential_amplitude * exp(i * k * direction . x)
        !
        ! Here phi is the velocity potential, k is the wavenumber, direction is
        ! a unit vector, and "." means the dot product.
        ! Velocity-potential amplitude. Pressure amplitudes should be converted
        ! before this type is filled: phi = pressure / (i * omega * rho).
        complex(dp) :: potential_amplitude = complex_zero

        ! Propagation direction. The setter normalizes this vector.
        real(dp) :: direction(3) = [1.0_dp, 0.0_dp, 0.0_dp]

        ! If true, use 0.5 * (exp(i k d.x) + exp(-i k d.x)).
        logical :: standing_wave = .false.
    end type au_plane_wave_type

contains

    subroutine reset_au_fields(case_data, status, message)
        ! Reset all nodal acoustic arrays to zero.
        !
        ! This is useful before setting a new boundary condition or incident
        ! field on an already allocated case.  It keeps allocated memory but
        ! clears old values.
        type(au_case_type), intent(inout) :: case_data
        integer, intent(out), optional :: status
        character(len=*), intent(out), optional :: message

        integer :: error_code
        character(len=256) :: error_message

        error_code = 0
        error_message = "AU field arrays reset."

        ! The arrays must already have been allocated by allocate_au_case.
        if (.not. field_arrays_are_allocated(case_data)) then
            call finish(1, "AU field arrays are not allocated.")
            return
        end if

        ! Prescribed nodal data used by Dirichlet, Neumann, or Robin conditions.
        case_data%prescribed_boundary_data = complex_zero

        ! Incident-field arrays.
        case_data%incident%external_phi = complex_zero
        case_data%incident%external_dphi_dn = complex_zero
        case_data%incident%internal_phi = complex_zero
        case_data%incident%internal_dphi_dn = complex_zero

        ! Solver output arrays and reconstructed total fields/pressures.
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

    end subroutine reset_au_fields

    subroutine expand_particle_boundary_values(geometry, case_data, status, message)
        ! Expand particle-level boundary values to node-level values.
        !
        ! Boundary conditions are stored once per particle in case_data.  The
        ! solver needs one known value at each node, so this routine copies the
        ! relevant particle value to every node on that particle.
        type(geometry_type), intent(in) :: geometry
        type(au_case_type), intent(inout) :: case_data
        integer, intent(out), optional :: status
        character(len=*), intent(out), optional :: message

        integer :: node_id
        integer :: particle_id
        integer :: error_code
        character(len=256) :: error_message

        error_code = 0
        error_message = "Particle boundary values expanded to nodes."

        call validate_boundary_inputs(geometry, case_data, error_code, error_message)
        if (error_code /= 0) then
            call finish(error_code, error_message)
            return
        end if

        do node_id = 1, geometry%mesh%node_count
            particle_id = geometry%mesh%node_particle_id(node_id)

            ! For Robin conditions, the known quantity is the right-hand side:
            !   robin_a * phi + robin_b * dphi/dn = robin_rhs
            !
            ! For Dirichlet or Neumann conditions, the known quantity is stored
            ! in boundary_condition%prescribed_data.
            if (case_data%boundary_condition(particle_id)%kind == bc_robin_external .or. &
                case_data%boundary_condition(particle_id)%kind == bc_robin_internal) then
                case_data%prescribed_boundary_data(node_id) = &
                    case_data%boundary_condition(particle_id)%robin_rhs
            else
                case_data%prescribed_boundary_data(node_id) = &
                    case_data%boundary_condition(particle_id)%prescribed_data
            end if
        end do

        call finish(error_code, error_message)

    contains

        subroutine finish(code, text)
            ! Return optional status information to the caller.
            integer, intent(in) :: code
            character(len=*), intent(in) :: text

            if (present(status)) status = code
            if (present(message)) message = text
        end subroutine finish

    end subroutine expand_particle_boundary_values

    subroutine set_normal_velocity_boundary(geometry, case_data, particle_id, velocity, &
                                            outward_positive, status, message)
        ! Set a Neumann-style boundary value from a prescribed velocity vector.
        !
        ! A surface normal velocity is the dot product between the velocity
        ! vector and the local unit normal:
        !   normal velocity = velocity . normal
        !
        ! This routine is useful for simple vibration/radiation test cases where
        ! one particle moves with a constant complex velocity vector.
        type(geometry_type), intent(in) :: geometry
        type(au_case_type), intent(inout) :: case_data
        integer, intent(in) :: particle_id

        ! Complex velocity amplitude in x, y, z directions.
        complex(dp), intent(in) :: velocity(3)

        ! If outward_positive is false, the sign is flipped.  This lets the
        ! caller match a convention where inward normal velocity is positive.
        logical, intent(in) :: outward_positive
        integer, intent(out), optional :: status
        character(len=*), intent(out), optional :: message

        integer :: node_id
        integer :: error_code
        real(dp) :: sign_value
        character(len=256) :: error_message

        error_code = 0
        error_message = "Normal velocity boundary values set."

        call validate_boundary_inputs(geometry, case_data, error_code, error_message)
        if (error_code /= 0) then
            call finish(error_code, error_message)
            return
        end if

        if (particle_id < 1 .or. particle_id > geometry%particle_count) then
            call finish(1, "Particle id is outside the geometry.")
            return
        end if

        ! Choose the sign convention for the normal velocity.
        if (outward_positive) then
            sign_value = 1.0_dp
        else
            sign_value = -1.0_dp
        end if

        do node_id = 1, geometry%mesh%node_count
            if (geometry%mesh%node_particle_id(node_id) /= particle_id) cycle

            ! Project the constant velocity vector onto the local surface normal
            ! at this node.
            case_data%prescribed_boundary_data(node_id) = sign_value * &
                (velocity(1) * geometry%differential%normal(1, node_id) + &
                 velocity(2) * geometry%differential%normal(2, node_id) + &
                 velocity(3) * geometry%differential%normal(3, node_id))
        end do

        call finish(error_code, error_message)

    contains

        subroutine finish(code, text)
            ! Return optional status information to the caller.
            integer, intent(in) :: code
            character(len=*), intent(in) :: text

            if (present(status)) status = code
            if (present(message)) message = text
        end subroutine finish

    end subroutine set_normal_velocity_boundary

    subroutine set_plane_wave_incident_field(geometry, case_data, wave, status, message)
        ! Fill the exterior incident field arrays for a plane wave.
        !
        ! For a travelling plane wave:
        !   phi_inc = A * exp(i * k * direction . x)
        !
        ! Its normal derivative is:
        !   dphi_inc/dn = A * i * k * (direction . normal)
        !                  * exp(i * k * direction . x)
        !
        ! The internal incident arrays are not filled here because the public
        ! validation cases use an exterior incident field.
        type(geometry_type), intent(in) :: geometry
        type(au_case_type), intent(inout) :: case_data
        type(au_plane_wave_type), intent(in) :: wave
        integer, intent(out), optional :: status
        character(len=*), intent(out), optional :: message

        integer :: node_id
        integer :: error_code
        real(dp) :: direction(3)
        real(dp) :: direction_norm
        real(dp) :: phase_distance
        real(dp) :: normal_projection
        complex(dp) :: phase
        complex(dp) :: exp_forward
        complex(dp) :: exp_backward
        complex(dp) :: wavenumber
        character(len=256) :: error_message

        error_code = 0
        error_message = "Plane-wave incident field set."

        call validate_incident_inputs(geometry, case_data, error_code, error_message)
        if (error_code /= 0) then
            call finish(error_code, error_message)
            return
        end if

        ! Normalise the propagation direction.  This allows callers to provide
        ! [2, 0, 0] or [1, 0, 0] and get the same wave direction.
        direction = wave%direction
        direction_norm = sqrt(dot_product(direction, direction))
        if (direction_norm <= tiny_direction) then
            call finish(1, "Plane-wave direction has near-zero length.")
            return
        end if
        direction = direction / direction_norm

        wavenumber = case_data%exterior_medium%wavenumber

        ! Clear exterior incident arrays before filling them.
        case_data%incident%external_phi = complex_zero
        case_data%incident%external_dphi_dn = complex_zero

        do node_id = 1, geometry%mesh%node_count
            ! direction . x is the distance along the propagation direction.
            phase_distance = dot_product(direction, geometry%mesh%xyz(:, node_id))

            ! direction . normal appears when differentiating the plane wave in
            ! the surface-normal direction.
            normal_projection = dot_product(direction, geometry%differential%normal(:, node_id))
            phase = imaginary_unit * wavenumber * phase_distance
            exp_forward = exp(phase)

            if (wave%standing_wave) then
                ! Standing wave formed by adding equal forward and backward
                ! travelling waves.
                exp_backward = exp(-phase)
                case_data%incident%external_phi(node_id) = &
                    wave%potential_amplitude * 0.5_dp * (exp_forward + exp_backward)
                case_data%incident%external_dphi_dn(node_id) = &
                    wave%potential_amplitude * 0.5_dp * imaginary_unit * wavenumber * &
                    normal_projection * (exp_forward - exp_backward)
            else
                ! Travelling wave in the +direction direction.
                case_data%incident%external_phi(node_id) = wave%potential_amplitude * exp_forward
                case_data%incident%external_dphi_dn(node_id) = &
                    wave%potential_amplitude * imaginary_unit * wavenumber * normal_projection * exp_forward
            end if
        end do

        call finish(error_code, error_message)

    contains

        subroutine finish(code, text)
            ! Return optional status information to the caller.
            integer, intent(in) :: code
            character(len=*), intent(in) :: text

            if (present(status)) status = code
            if (present(message)) message = text
        end subroutine finish

    end subroutine set_plane_wave_incident_field

    subroutine update_total_fields_and_pressure(geometry, case_data, status, message)
        ! Combine incident and scattered fields, then reconstruct pressure.
        !
        ! The solver stores scattered fields directly.  For post-processing and
        ! comparison with analytical solutions, it is often more useful to have
        ! the total field:
        !   total = incident + scattered
        !
        ! This routine also computes pressure from velocity potential using the
        ! convention used in this code:
        !   time dependence = exp(-i * omega * t)
        !   pressure        = i * omega * rho * phi
        type(geometry_type), intent(in) :: geometry
        type(au_case_type), intent(inout) :: case_data
        integer, intent(out), optional :: status
        character(len=*), intent(out), optional :: message

        integer :: node_id
        integer :: particle_id
        integer :: error_code
        complex(dp) :: pressure_factor
        character(len=256) :: error_message

        error_code = 0
        error_message = "Total fields and pressure updated."

        call validate_incident_inputs(geometry, case_data, error_code, error_message)
        if (error_code /= 0) then
            call finish(error_code, error_message)
            return
        end if

        ! Add incident and scattered fields node-by-node.
        case_data%solution%total_external_phi = &
            case_data%incident%external_phi + case_data%solution%external_phi
        case_data%solution%total_external_dphi_dn = &
            case_data%incident%external_dphi_dn + case_data%solution%external_dphi_dn
        case_data%solution%total_internal_phi = &
            case_data%incident%internal_phi + case_data%solution%internal_phi
        case_data%solution%total_internal_dphi_dn = &
            case_data%incident%internal_dphi_dn + case_data%solution%internal_dphi_dn

        ! Exterior pressure uses the exterior medium density.
        pressure_factor = imaginary_unit * case_data%angular_frequency * &
                          case_data%exterior_medium%density
        case_data%solution%external_pressure = &
            pressure_factor * case_data%solution%total_external_phi

        ! Interior pressure can use a different density for each particle.
        do node_id = 1, geometry%mesh%node_count
            particle_id = geometry%mesh%node_particle_id(node_id)
            pressure_factor = imaginary_unit * case_data%angular_frequency * &
                              case_data%interior_medium(particle_id)%density
            case_data%solution%internal_pressure(node_id) = &
                pressure_factor * case_data%solution%total_internal_phi(node_id)
        end do

        call finish(error_code, error_message)

    contains

        subroutine finish(code, text)
            ! Return optional status information to the caller.
            integer, intent(in) :: code
            character(len=*), intent(in) :: text

            if (present(status)) status = code
            if (present(message)) message = text
        end subroutine finish

    end subroutine update_total_fields_and_pressure

    subroutine validate_boundary_inputs(geometry, case_data, status, message)
        ! Check inputs shared by boundary-condition routines.
        !
        ! These routines need:
        !   - a geometry with nodes and particles
        !   - node-to-particle ids
        !   - allocated particle boundary conditions
        !   - allocated nodal known-boundary array
        !   - prepared surface normals
        type(geometry_type), intent(in) :: geometry
        type(au_case_type), intent(in) :: case_data
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        status = 0
        message = "Boundary inputs are valid."

        if (geometry%mesh%node_count <= 0) then
            call fail("Geometry has no nodes.")
        else if (geometry%particle_count <= 0) then
            call fail("Geometry has no particles.")
        else if (.not. allocated(geometry%mesh%node_particle_id)) then
            call fail("Geometry node particle ids are not allocated.")
        else if (.not. allocated(case_data%boundary_condition)) then
            call fail("Boundary conditions are not allocated.")
        else if (size(case_data%boundary_condition) /= geometry%particle_count) then
            call fail("Boundary condition count does not match particle count.")
        else if (.not. allocated(case_data%prescribed_boundary_data)) then
            call fail("Prescribed boundary data are not allocated.")
        else if (size(case_data%prescribed_boundary_data) /= geometry%mesh%node_count) then
            call fail("Prescribed boundary data count does not match node count.")
        else if (.not. allocated(geometry%differential%normal)) then
            call fail("Surface normals have not been prepared.")
        end if

    contains

        subroutine fail(text)
            ! Store only the first validation error.
            character(len=*), intent(in) :: text

            if (status == 0) then
                status = 1
                message = text
            end if
        end subroutine fail

    end subroutine validate_boundary_inputs

    subroutine validate_incident_inputs(geometry, case_data, status, message)
        ! Check inputs needed by incident-field and total-field routines.
        !
        ! This builds on validate_boundary_inputs and adds checks for incident,
        ! solution, and interior-medium arrays.
        type(geometry_type), intent(in) :: geometry
        type(au_case_type), intent(in) :: case_data
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        status = 0
        message = "Incident field inputs are valid."

        call validate_boundary_inputs(geometry, case_data, status, message)
        if (status /= 0) return

        if (.not. field_arrays_are_allocated(case_data)) then
            call fail("AU field arrays are not allocated.")
        else if (.not. allocated(case_data%interior_medium)) then
            call fail("Interior media are not allocated.")
        else if (size(case_data%interior_medium) /= geometry%particle_count) then
            call fail("Interior medium count does not match particle count.")
        end if

    contains

        subroutine fail(text)
            ! Store only the first validation error.
            character(len=*), intent(in) :: text

            if (status == 0) then
                status = 1
                message = text
            end if
        end subroutine fail

    end subroutine validate_incident_inputs

    pure logical function field_arrays_are_allocated(case_data)
        ! Check that all nodal field arrays have been allocated.
        !
        ! This function checks allocation only.  Size checks are handled in
        ! AU_Types%validate_au_case.
        type(au_case_type), intent(in) :: case_data

        field_arrays_are_allocated = &
            allocated(case_data%prescribed_boundary_data) .and. &
            allocated(case_data%incident%external_phi) .and. &
            allocated(case_data%incident%external_dphi_dn) .and. &
            allocated(case_data%incident%internal_phi) .and. &
            allocated(case_data%incident%internal_dphi_dn) .and. &
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
    end function field_arrays_are_allocated

end module AU_BoundaryConditions
