program ValidateRigidSphere
    ! Solve a sphere-scattering validation case for the Python test scripts.
    !
    ! Python creates a linear or quadratic sphere mesh and passes the physical
    ! parameters on the command line.  This program prepares the common Fortran
    ! geometry and acoustic types, solves the requested exterior boundary
    ! condition, and writes nodal fields to CSV.  Python then compares those
    ! fields with the independent spherical partial-wave solution in
    ! ana/analytic.py.
    use Pre_Constants, only: dp, complex_zero
    use Geom_Types, only: geometry_type
    use Geom_ReadMesh, only: read_mesh
    use Geom_MeshTopology, only: prepare_mesh_topology
    use Geom_SurfaceGeometry, only: prepare_surface_geometry
    use AU_Types, only: au_case_type, allocate_au_case, bc_dirichlet_external, &
                        bc_neumann_external, bc_robin_external
    use AU_BoundaryConditions, only: au_plane_wave_type, expand_particle_boundary_values, &
                                     reset_au_fields, set_plane_wave_incident_field
    use AU_Solver, only: solve_au_surface
    implicit none

    type(geometry_type) :: geometry
    type(au_case_type) :: case_data
    type(au_plane_wave_type) :: wave
    character(len=512) :: mesh_file
    character(len=512) :: output_file
    character(len=512) :: argument
    character(len=512) :: message
    character(len=32) :: boundary_mode
    character(len=32) :: formulation_mode
    real(dp) :: wavenumber
    real(dp) :: radius
    real(dp) :: robin_a
    real(dp) :: robin_b
    integer :: exterior_domain_normal_sign
    integer :: status

    wavenumber = 1.0_dp
    radius = 1.0_dp
    robin_a = 1.0_dp
    robin_b = 0.5_dp
    exterior_domain_normal_sign = -1
    boundary_mode = "neumann"
    formulation_mode = "ordinary"

    if (command_argument_count() < 2) then
        write(*, '(A)') "Usage: validate_rigid_sphere mesh_file output_csv " // &
                        "[wavenumber] [radius] [exterior_domain_normal_sign] [dirichlet|neumann|robin] " // &
                        "[robin_a] [robin_b] [ordinary|burton-miller]"
        error stop 2
    end if

    call get_command_argument(1, mesh_file)
    call get_command_argument(2, output_file)

    if (command_argument_count() >= 3) then
        call get_command_argument(3, argument)
        read(argument, *) wavenumber
    end if

    if (command_argument_count() >= 4) then
        call get_command_argument(4, argument)
        read(argument, *) radius
    end if

    if (command_argument_count() >= 5) then
        call get_command_argument(5, argument)
        read(argument, *) exterior_domain_normal_sign
    end if

    if (command_argument_count() >= 6) then
        call get_command_argument(6, boundary_mode)
    end if

    if (command_argument_count() >= 7) then
        call get_command_argument(7, argument)
        read(argument, *) robin_a
    end if

    if (command_argument_count() >= 8) then
        call get_command_argument(8, argument)
        read(argument, *) robin_b
    end if

    if (command_argument_count() >= 9) then
        call get_command_argument(9, formulation_mode)
    end if

    call read_mesh(trim(mesh_file), geometry, status, message)
    call require_ok(status, message)

    call prepare_mesh_topology(geometry, status, message)
    call require_ok(status, message)

    call prepare_surface_geometry(geometry, status, message)
    call require_ok(status, message)

    call allocate_au_case(case_data, geometry, status, message)
    call require_ok(status, message)

    call configure_sphere_case(geometry, case_data, wave, wavenumber, radius, &
                               exterior_domain_normal_sign, trim(boundary_mode), robin_a, robin_b, &
                               trim(formulation_mode))

    call reset_au_fields(case_data, status, message)
    call require_ok(status, message)

    call expand_particle_boundary_values(geometry, case_data, status, message)
    call require_ok(status, message)

    call set_plane_wave_incident_field(geometry, case_data, wave, status, message)
    call require_ok(status, message)

    call solve_au_surface(geometry, case_data, status, message)
    call require_ok(status, message)

    call write_surface_validation_csv(trim(output_file), geometry, case_data, status, message)
    call require_ok(status, message)

    write(*, '(A)') "Rigid-sphere validation solve completed with " // &
                    trim(formulation_mode) // " formulation."

contains

    subroutine configure_sphere_case(geometry, case_data, wave, wavenumber, radius, &
                                     exterior_domain_normal_sign, boundary_mode, robin_a, robin_b, formulation_mode)
        ! Fill one allocated case with nondimensional sphere-test data.
        !
        ! The validation sets density and sound speed to one, so omega=k.  The
        ! boundary values apply to the total field; AU_Solver subtracts the
        ! incident plane wave when forming its scattered-field unknowns.
        type(geometry_type), intent(in) :: geometry
        type(au_case_type), intent(inout) :: case_data
        type(au_plane_wave_type), intent(inout) :: wave
        real(dp), intent(in) :: wavenumber
        real(dp), intent(in) :: radius
        integer, intent(in) :: exterior_domain_normal_sign
        character(len=*), intent(in) :: boundary_mode
        real(dp), intent(in) :: robin_a
        real(dp), intent(in) :: robin_b
        character(len=*), intent(in) :: formulation_mode

        integer :: particle_id

        select case (formulation_mode)
        case ("ordinary")
            case_data%use_burton_miller = .false.
        case ("burton-miller", "bm")
            case_data%use_burton_miller = .true.
        case default
            write(*, '(A)') "Unsupported formulation mode: " // trim(formulation_mode)
            error stop 3
        end select
        case_data%angular_frequency = wavenumber
        case_data%exterior_medium%density = 1.0_dp
        case_data%exterior_medium%sound_speed = cmplx(1.0_dp, 0.0_dp, kind=dp)
        case_data%exterior_medium%wavenumber = cmplx(wavenumber, 0.0_dp, kind=dp)

        do particle_id = 1, geometry%particle_count
            case_data%interior_medium(particle_id)%density = 1.0_dp
            case_data%interior_medium(particle_id)%sound_speed = cmplx(1.0_dp, 0.0_dp, kind=dp)
            case_data%interior_medium(particle_id)%wavenumber = cmplx(wavenumber, 0.0_dp, kind=dp)

            case_data%layer(particle_id)%parent_particle_id = 0
            case_data%layer(particle_id)%exterior_domain_normal_sign = exterior_domain_normal_sign
            case_data%layer(particle_id)%characteristic_length = radius

            select case (boundary_mode)
            case ("dirichlet")
                ! Pressure-release boundary: the total potential is prescribed
                ! as zero, so the scattered value cancels the incident value.
                case_data%boundary_condition(particle_id)%kind = bc_dirichlet_external
                case_data%boundary_condition(particle_id)%prescribed_data = complex_zero
            case ("neumann")
                case_data%boundary_condition(particle_id)%kind = bc_neumann_external
                case_data%boundary_condition(particle_id)%prescribed_data = complex_zero
            case ("robin")
                case_data%boundary_condition(particle_id)%kind = bc_robin_external
                case_data%boundary_condition(particle_id)%robin_a = cmplx(robin_a, 0.0_dp, kind=dp)
                case_data%boundary_condition(particle_id)%robin_b = cmplx(robin_b, 0.0_dp, kind=dp)
                case_data%boundary_condition(particle_id)%robin_rhs = complex_zero
            case default
                write(*, '(A)') "Unsupported boundary mode: " // trim(boundary_mode)
                error stop 3
            end select
        end do

        wave%potential_amplitude = cmplx(1.0_dp, 0.0_dp, kind=dp)
        wave%direction = [1.0_dp, 0.0_dp, 0.0_dp]
        wave%standing_wave = .false.
    end subroutine configure_sphere_case

    subroutine write_surface_validation_csv(filename, geometry, case_data, status, message)
        ! Write the numerical and incident nodal fields used by Python's
        ! analytical error calculation.
        character(len=*), intent(in) :: filename
        type(geometry_type), intent(in) :: geometry
        type(au_case_type), intent(in) :: case_data
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: unit
        integer :: node_id

        status = 0
        message = "Rigid-sphere validation CSV written."

        open(newunit=unit, file=filename, status="replace", action="write", iostat=status, iomsg=message)
        if (status /= 0) return

        write(unit, '(A)') "node_id,x,y,z,nx,ny,nz," // &
            "phi_inc_re,phi_inc_im,phi_sca_re,phi_sca_im,phi_total_re,phi_total_im," // &
            "dphi_inc_re,dphi_inc_im,dphi_sca_re,dphi_sca_im,dphi_total_re,dphi_total_im," // &
            "bm_aux_dpsi_dn_re,bm_aux_dpsi_dn_im"

        do node_id = 1, geometry%mesh%node_count
            write(unit, '(I0,20(",",ES24.16E3))') node_id, &
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
                real(case_data%solution%bm_auxiliary_dpsi_dn(node_id)), &
                aimag(case_data%solution%bm_auxiliary_dpsi_dn(node_id))
        end do

        close(unit)
    end subroutine write_surface_validation_csv

    subroutine require_ok(status, message)
        integer, intent(in) :: status
        character(len=*), intent(in) :: message

        if (status /= 0) then
            write(*, '(A, I0)') "Validation failed with status ", status
            write(*, '(A)') trim(message)
            error stop 1
        end if
    end subroutine require_ok

end program ValidateRigidSphere
