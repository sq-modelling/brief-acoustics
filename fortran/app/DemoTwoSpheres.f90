program DemoTwoSpheres
    ! Solve a two-sphere radiation demonstration case.
    !
    ! This executable is intentionally small and specific.  Python creates a
    ! two-particle surface mesh and this program solves the exterior Neumann
    ! acoustic problem:
    !
    !   dphi/dn = prescribed_normal_derivative_1 on sphere 1
    !   dphi/dn = prescribed_normal_derivative_2 on sphere 2
    !
    ! The default Python demo prescribes +1 and -1 respectively, so the spheres
    ! radiate out of phase, similar in spirit to the small-gap example in the
    ! 2015 RSOS BRIEF acoustics paper.
    use Pre_Constants, only: dp
    use Geom_Types, only: geometry_type
    use Geom_ReadMesh, only: read_mesh
    use Geom_MeshTopology, only: prepare_mesh_topology
    use Geom_SurfaceGeometry, only: prepare_surface_geometry
    use AU_Types, only: au_case_type, allocate_au_case, bc_neumann_external
    use AU_BoundaryConditions, only: expand_particle_boundary_values, reset_au_fields
    use AU_Solver, only: solve_au_surface
    implicit none

    type(geometry_type) :: geometry
    type(au_case_type) :: case_data
    character(len=512) :: mesh_file
    character(len=512) :: output_file
    character(len=512) :: argument
    character(len=512) :: message
    real(dp) :: wavenumber
    real(dp) :: radius
    real(dp) :: normal_derivative_1_real
    real(dp) :: normal_derivative_1_imag
    real(dp) :: normal_derivative_2_real
    real(dp) :: normal_derivative_2_imag
    integer :: exterior_free_term_sign
    integer :: status

    wavenumber = 1.5707963267948966_dp
    radius = 1.0_dp
    exterior_free_term_sign = -1
    normal_derivative_1_real = 1.0_dp
    normal_derivative_1_imag = 0.0_dp
    normal_derivative_2_real = -1.0_dp
    normal_derivative_2_imag = 0.0_dp

    if (command_argument_count() < 2) then
        write(*, '(A)') "Usage: demo_two_spheres mesh_file output_csv " // &
                        "[wavenumber] [radius] [exterior_free_term_sign] " // &
                        "[dphi_dn_1_re] [dphi_dn_1_im] [dphi_dn_2_re] [dphi_dn_2_im]"
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
        read(argument, *) exterior_free_term_sign
    end if
    if (command_argument_count() >= 6) then
        call get_command_argument(6, argument)
        read(argument, *) normal_derivative_1_real
    end if
    if (command_argument_count() >= 7) then
        call get_command_argument(7, argument)
        read(argument, *) normal_derivative_1_imag
    end if
    if (command_argument_count() >= 8) then
        call get_command_argument(8, argument)
        read(argument, *) normal_derivative_2_real
    end if
    if (command_argument_count() >= 9) then
        call get_command_argument(9, argument)
        read(argument, *) normal_derivative_2_imag
    end if

    call read_mesh(trim(mesh_file), geometry, status, message)
    call require_ok(status, message)

    call prepare_mesh_topology(geometry, status, message)
    call require_ok(status, message)

    call prepare_surface_geometry(geometry, status, message)
    call require_ok(status, message)

    call allocate_au_case(case_data, geometry, status, message)
    call require_ok(status, message)

    call configure_two_sphere_case(geometry, case_data, wavenumber, radius, exterior_free_term_sign, &
                                   cmplx(normal_derivative_1_real, normal_derivative_1_imag, kind=dp), &
                                   cmplx(normal_derivative_2_real, normal_derivative_2_imag, kind=dp))

    call reset_au_fields(case_data, status, message)
    call require_ok(status, message)

    call expand_particle_boundary_values(geometry, case_data, status, message)
    call require_ok(status, message)

    call solve_au_surface(geometry, case_data, status, message)
    call require_ok(status, message)

    call write_surface_csv(trim(output_file), geometry, case_data, status, message)
    call require_ok(status, message)

    write(*, '(A)') "Two-sphere radiation demo solve completed."

contains

    subroutine configure_two_sphere_case(geometry, case_data, wavenumber, radius, &
                                         exterior_free_term_sign, prescribed_normal_derivative_1, &
                                         prescribed_normal_derivative_2)
        ! Configure the ordinary exterior two-particle demonstration.
        !
        ! Unit density and sound speed make omega=k.  There is no incident
        ! field; the prescribed Neumann values directly drive radiation from
        ! the two surfaces.  Burton-Miller remains disabled because its public
        ! implementation is currently restricted to one exterior particle.
        type(geometry_type), intent(in) :: geometry
        type(au_case_type), intent(inout) :: case_data
        real(dp), intent(in) :: wavenumber
        real(dp), intent(in) :: radius
        integer, intent(in) :: exterior_free_term_sign
        complex(dp), intent(in) :: prescribed_normal_derivative_1
        complex(dp), intent(in) :: prescribed_normal_derivative_2

        integer :: particle_id

        if (geometry%particle_count /= 2) then
            write(*, '(A)') "The two-sphere demo expects exactly two particles."
            error stop 3
        end if

        case_data%use_burton_miller = .false.
        case_data%angular_frequency = wavenumber
        case_data%exterior_medium%density = 1.0_dp
        case_data%exterior_medium%sound_speed = cmplx(1.0_dp, 0.0_dp, kind=dp)
        case_data%exterior_medium%wavenumber = cmplx(wavenumber, 0.0_dp, kind=dp)

        do particle_id = 1, geometry%particle_count
            case_data%interior_medium(particle_id)%density = 1.0_dp
            case_data%interior_medium(particle_id)%sound_speed = cmplx(1.0_dp, 0.0_dp, kind=dp)
            case_data%interior_medium(particle_id)%wavenumber = cmplx(wavenumber, 0.0_dp, kind=dp)

            case_data%layer(particle_id)%parent_particle_id = 0
            case_data%layer(particle_id)%exterior_free_term_sign = exterior_free_term_sign
            case_data%layer(particle_id)%characteristic_length = radius

            case_data%boundary_condition(particle_id)%kind = bc_neumann_external
        end do

        case_data%boundary_condition(1)%prescribed_data = prescribed_normal_derivative_1
        case_data%boundary_condition(2)%prescribed_data = prescribed_normal_derivative_2
    end subroutine configure_two_sphere_case

    subroutine write_surface_csv(filename, geometry, case_data, status, message)
        ! Write portable nodal geometry, total fields, and pressure for the
        ! Python gap-line postprocessor.
        character(len=*), intent(in) :: filename
        type(geometry_type), intent(in) :: geometry
        type(au_case_type), intent(in) :: case_data
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: unit
        integer :: node_id

        status = 0
        message = "Two-sphere surface CSV written."

        open(newunit=unit, file=filename, status="replace", action="write", iostat=status, iomsg=message)
        if (status /= 0) return

        write(unit, '(A)') "node_id,particle_id,x,y,z,nx,ny,nz," // &
            "phi_re,phi_im,dphi_dn_re,dphi_dn_im,pressure_re,pressure_im"

        do node_id = 1, geometry%mesh%node_count
            write(unit, '(2(I0,","),12(ES24.16E3,:","))') &
                node_id, &
                geometry%mesh%node_particle_id(node_id), &
                geometry%mesh%xyz(1, node_id), &
                geometry%mesh%xyz(2, node_id), &
                geometry%mesh%xyz(3, node_id), &
                geometry%differential%normal(1, node_id), &
                geometry%differential%normal(2, node_id), &
                geometry%differential%normal(3, node_id), &
                real(case_data%solution%total_external_phi(node_id)), &
                aimag(case_data%solution%total_external_phi(node_id)), &
                real(case_data%solution%total_external_dphi_dn(node_id)), &
                aimag(case_data%solution%total_external_dphi_dn(node_id)), &
                real(case_data%solution%external_pressure(node_id)), &
                aimag(case_data%solution%external_pressure(node_id))
        end do

        close(unit)
    end subroutine write_surface_csv

    subroutine require_ok(status, message)
        integer, intent(in) :: status
        character(len=*), intent(in) :: message

        if (status /= 0) then
            write(*, '(A, I0)') "Two-sphere demo failed with status ", status
            write(*, '(A)') trim(message)
            error stop 1
        end if
    end subroutine require_ok

end program DemoTwoSpheres
