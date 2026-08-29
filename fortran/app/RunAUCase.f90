program RunAUCase
    ! General single-particle exterior entry point for the public workflow.
    !
    ! Users do not write this program's namelist directly.  The Python command
    ! `brief-acoustics run case.toml` validates the public TOML file and mesh, writes a
    ! normalized mesh plus a private runtime namelist, and then invokes this
    ! executable.  Keeping TOML parsing in Python gives users better validation
    ! messages while this program stays focused on the numerical workflow.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use Pre_Constants, only: dp, complex_zero
    use Geom_Types, only: geometry_type
    use Geom_ReadMesh, only: read_mesh
    use Geom_MeshTopology, only: prepare_mesh_topology
    use Geom_SurfaceGeometry, only: prepare_surface_geometry
    use AU_Types, only: au_case_type, allocate_au_case, &
                        bc_dirichlet_external, bc_neumann_external, bc_robin_external
    use AU_BoundaryConditions, only: au_plane_wave_type, reset_au_fields, &
                                     expand_particle_boundary_values, &
                                     set_plane_wave_incident_field
    use AU_Solver, only: solve_au_surface
    implicit none

    type(geometry_type) :: geometry
    type(au_case_type) :: case_data
    type(au_plane_wave_type) :: wave

    character(len=512) :: mesh_file
    character(len=512) :: runtime_case_file
    character(len=512) :: output_file
    character(len=512) :: message
    character(len=32) :: formulation_mode
    character(len=32) :: boundary_mode
    character(len=32) :: incident_mode

    real(dp) :: wavenumber
    real(dp) :: angular_frequency
    real(dp) :: exterior_density
    real(dp) :: exterior_sound_speed
    real(dp) :: characteristic_length
    real(dp) :: incident_direction(3)
    integer :: exterior_free_term_sign
    logical :: standing_wave

    complex(dp) :: prescribed_boundary_data
    complex(dp) :: robin_a
    complex(dp) :: robin_b
    complex(dp) :: robin_rhs
    complex(dp) :: incident_potential_amplitude

    integer :: runtime_unit
    integer :: status

    namelist /au_case_input/ wavenumber, angular_frequency, exterior_density, &
        exterior_sound_speed, exterior_free_term_sign, characteristic_length, formulation_mode, &
        boundary_mode, prescribed_boundary_data, robin_a, robin_b, robin_rhs, incident_mode, &
        incident_potential_amplitude, incident_direction, standing_wave

    ! Defaults make an incomplete namelist fail predictably during validation.
    wavenumber = -1.0_dp
    angular_frequency = -1.0_dp
    exterior_density = -1.0_dp
    exterior_sound_speed = -1.0_dp
    characteristic_length = -1.0_dp
    exterior_free_term_sign = 0
    formulation_mode = ""
    boundary_mode = ""
    incident_mode = ""
    prescribed_boundary_data = complex_zero
    robin_a = complex_zero
    robin_b = complex_zero
    robin_rhs = complex_zero
    incident_potential_amplitude = complex_zero
    incident_direction = [1.0_dp, 0.0_dp, 0.0_dp]
    standing_wave = .false.

    if (command_argument_count() /= 3) then
        write(*, '(A)') "Usage: run_au_case normalized_mesh runtime_case.nml output.csv"
        error stop 2
    end if
    call get_command_argument(1, mesh_file)
    call get_command_argument(2, runtime_case_file)
    call get_command_argument(3, output_file)

    open(newunit=runtime_unit, file=trim(runtime_case_file), status="old", &
         action="read", iostat=status, iomsg=message)
    call require_ok(status, "Unable to open runtime case: " // trim(message))
    read(runtime_unit, nml=au_case_input, iostat=status, iomsg=message)
    close(runtime_unit)
    call require_ok(status, "Unable to read runtime case: " // trim(message))

    call validate_runtime_inputs(status, message)
    call require_ok(status, message)

    call read_mesh(trim(mesh_file), geometry, status, message)
    call require_ok(status, message)
    if (geometry%particle_count /= 1) then
        write(*, '(A)') "The public case runner requires exactly one particle."
        error stop 3
    end if

    call prepare_mesh_topology(geometry, status, message)
    call require_ok(status, message)
    call prepare_surface_geometry(geometry, status, message)
    call require_ok(status, message)
    call allocate_au_case(case_data, geometry, status, message)
    call require_ok(status, message)

    call configure_case(case_data, wave)
    call reset_au_fields(case_data, status, message)
    call require_ok(status, message)
    call expand_particle_boundary_values(geometry, case_data, status, message)
    call require_ok(status, message)

    if (trim(incident_mode) == "plane-wave") then
        call set_plane_wave_incident_field(geometry, case_data, wave, status, message)
        call require_ok(status, message)
    end if

    call solve_au_surface(geometry, case_data, status, message)
    call require_ok(status, message)
    call write_surface_csv(trim(output_file), geometry, case_data, status, message)
    call require_ok(status, message)

    write(*, '(A)') "BRIEF-Acoustics case solve completed."
    write(*, '(A)') "formulation = " // trim(formulation_mode)
    write(*, '(A)') "boundary = " // trim(boundary_mode)
    write(*, '(A, I0)') "nodes = ", geometry%mesh%node_count
    write(*, '(A, I0)') "elements = ", geometry%mesh%element_count

contains

    subroutine validate_runtime_inputs(status, message)
        ! Defend the numerical core even though Python already validates TOML.
        ! This protects direct executable use and catches a damaged generated
        ! namelist with a clear message before matrix allocation.
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        real(dp) :: frequency_scale
        real(dp), parameter :: relation_tolerance = 1.0e-10_dp
        real(dp), parameter :: coefficient_tolerance = 1.0e-14_dp

        status = 0
        message = "Runtime case inputs are valid."

        if (.not. ieee_is_finite(wavenumber) .or. wavenumber <= 0.0_dp) then
            call set_failure(status, message, "Wavenumber must be finite and positive.")
        else if (.not. ieee_is_finite(angular_frequency) .or. angular_frequency <= 0.0_dp) then
            call set_failure(status, message, "Angular frequency must be finite and positive.")
        else if (.not. ieee_is_finite(exterior_density) .or. exterior_density <= 0.0_dp) then
            call set_failure(status, message, "Exterior density must be finite and positive.")
        else if (.not. ieee_is_finite(exterior_sound_speed) .or. exterior_sound_speed <= 0.0_dp) then
            call set_failure(status, message, "Exterior sound speed must be finite and positive.")
        else if (.not. ieee_is_finite(characteristic_length) .or. &
                 characteristic_length <= 0.0_dp) then
            call set_failure(status, message, "Characteristic length must be finite and positive.")
        else if (exterior_free_term_sign /= -1) then
            call set_failure(status, message, "Public exterior cases require exterior_free_term_sign=-1.")
        end if

        frequency_scale = max(1.0_dp, abs(angular_frequency), &
                              abs(wavenumber * exterior_sound_speed))
        if (status == 0 .and. abs(angular_frequency - wavenumber * exterior_sound_speed) > &
            relation_tolerance * frequency_scale) then
            call set_failure(status, message, &
                             "Runtime physics is inconsistent: omega must equal k*c.")
        end if

        if (status == 0) then
            select case (trim(formulation_mode))
            case ("ordinary", "burton-miller")
                continue
            case default
                call set_failure(status, message, &
                                 "Formulation must be ordinary or burton-miller.")
            end select
        end if

        if (status == 0) then
            select case (trim(boundary_mode))
            case ("dirichlet", "neumann")
                continue
            case ("robin")
                if (max(abs(robin_a), abs(robin_b)) <= coefficient_tolerance) then
                    call set_failure(status, message, &
                                     "Robin coefficients a and b cannot both be zero.")
                end if
            case default
                call set_failure(status, message, &
                                 "Boundary mode must be dirichlet, neumann, or robin.")
            end select
        end if

        if (status == 0) then
            select case (trim(incident_mode))
            case ("none", "plane-wave")
                continue
            case default
                call set_failure(status, message, &
                                 "Incident mode must be none or plane-wave.")
            end select
        end if

    end subroutine validate_runtime_inputs

    subroutine set_failure(status, message, text)
        ! Store the first input error without overwriting its useful context.
        integer, intent(inout) :: status
        character(len=*), intent(inout) :: message
        character(len=*), intent(in) :: text

        if (status == 0) then
            status = 1
            message = text
        end if
    end subroutine set_failure

    subroutine configure_case(case_data, wave)
        ! Transfer validated runtime values into the common solver data types.
        type(au_case_type), intent(inout) :: case_data
        type(au_plane_wave_type), intent(out) :: wave

        case_data%use_burton_miller = trim(formulation_mode) == "burton-miller"
        case_data%angular_frequency = angular_frequency
        case_data%exterior_medium%density = exterior_density
        case_data%exterior_medium%sound_speed = &
            cmplx(exterior_sound_speed, 0.0_dp, kind=dp)
        case_data%exterior_medium%wavenumber = cmplx(wavenumber, 0.0_dp, kind=dp)

        ! Interior values are required by the shared allocated case type but are
        ! inactive for an impenetrable exterior boundary condition.
        case_data%interior_medium(1)%density = exterior_density
        case_data%interior_medium(1)%sound_speed = &
            cmplx(exterior_sound_speed, 0.0_dp, kind=dp)
        case_data%interior_medium(1)%wavenumber = cmplx(wavenumber, 0.0_dp, kind=dp)

        case_data%layer(1)%parent_particle_id = 0
        case_data%layer(1)%exterior_free_term_sign = exterior_free_term_sign
        case_data%layer(1)%characteristic_length = characteristic_length

        select case (trim(boundary_mode))
        case ("dirichlet")
            case_data%boundary_condition(1)%kind = bc_dirichlet_external
            case_data%boundary_condition(1)%prescribed_data = prescribed_boundary_data
        case ("neumann")
            case_data%boundary_condition(1)%kind = bc_neumann_external
            case_data%boundary_condition(1)%prescribed_data = prescribed_boundary_data
        case ("robin")
            case_data%boundary_condition(1)%kind = bc_robin_external
            case_data%boundary_condition(1)%robin_a = robin_a
            case_data%boundary_condition(1)%robin_b = robin_b
            case_data%boundary_condition(1)%robin_rhs = robin_rhs
        end select

        wave%potential_amplitude = incident_potential_amplitude
        wave%direction = incident_direction
        wave%standing_wave = standing_wave
    end subroutine configure_case

    subroutine write_surface_csv(filename, geometry, case_data, status, message)
        ! Write portable comma-separated nodal data.  Python creates plots and
        ! summaries from this file.
        character(len=*), intent(in) :: filename
        type(geometry_type), intent(in) :: geometry
        type(au_case_type), intent(in) :: case_data
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: unit
        integer :: node_id

        status = 0
        message = "Surface solution CSV written."
        open(newunit=unit, file=filename, status="replace", action="write", &
             iostat=status, iomsg=message)
        if (status /= 0) return

        write(unit, '(A)') "node_id,particle_id,x,y,z,nx,ny,nz," // &
            "phi_inc_re,phi_inc_im,phi_sca_re,phi_sca_im,phi_total_re,phi_total_im," // &
            "dphi_inc_re,dphi_inc_im,dphi_sca_re,dphi_sca_im,dphi_total_re,dphi_total_im," // &
            "pressure_re,pressure_im,bm_aux_dpsi_dn_re,bm_aux_dpsi_dn_im"

        do node_id = 1, geometry%mesh%node_count
            write(unit, '(2(I0,","),22(ES24.16E3,:,","))') &
                node_id, &
                geometry%mesh%node_particle_id(node_id), &
                geometry%mesh%xyz(1, node_id), &
                geometry%mesh%xyz(2, node_id), &
                geometry%mesh%xyz(3, node_id), &
                geometry%differential%normal(1, node_id), &
                geometry%differential%normal(2, node_id), &
                geometry%differential%normal(3, node_id), &
                real(case_data%incident%external_phi(node_id)), &
                aimag(case_data%incident%external_phi(node_id)), &
                real(case_data%solution%external_phi(node_id)), &
                aimag(case_data%solution%external_phi(node_id)), &
                real(case_data%solution%total_external_phi(node_id)), &
                aimag(case_data%solution%total_external_phi(node_id)), &
                real(case_data%incident%external_dphi_dn(node_id)), &
                aimag(case_data%incident%external_dphi_dn(node_id)), &
                real(case_data%solution%external_dphi_dn(node_id)), &
                aimag(case_data%solution%external_dphi_dn(node_id)), &
                real(case_data%solution%total_external_dphi_dn(node_id)), &
                aimag(case_data%solution%total_external_dphi_dn(node_id)), &
                real(case_data%solution%external_pressure(node_id)), &
                aimag(case_data%solution%external_pressure(node_id)), &
                real(case_data%solution%bm_auxiliary_dpsi_dn(node_id)), &
                aimag(case_data%solution%bm_auxiliary_dpsi_dn(node_id))
        end do

        close(unit)
    end subroutine write_surface_csv

    subroutine require_ok(status, message)
        integer, intent(in) :: status
        character(len=*), intent(in) :: message

        if (status /= 0) then
            write(*, '(A, I0)') "BRIEF-Acoustics case failed with status ", status
            write(*, '(A)') trim(message)
            error stop 1
        end if
    end subroutine require_ok

end program RunAUCase
