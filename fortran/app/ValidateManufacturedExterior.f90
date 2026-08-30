program ValidateManufacturedExterior
    ! Validate the exterior solver on a non-spherical surface.
    !
    ! A point source placed inside the solid produces the outgoing field
    !
    !     phi(x) = exp(i*k*r) / r,
    !     r      = |x - x_source|.
    !
    ! The source singularity is outside the acoustic exterior domain, so phi is
    ! an exact homogeneous Helmholtz solution everywhere that we solve.  Its
    ! value can therefore be prescribed as Dirichlet data on any closed body.
    ! The numerical unknown dphi/dn is compared with its analytical value by
    ! the accompanying Python validation script.
    use Pre_Constants, only: dp, complex_zero, imaginary_unit
    use Geom_Types, only: geometry_type
    use Geom_ReadMesh, only: read_mesh
    use Geom_MeshTopology, only: prepare_mesh_topology
    use Geom_SurfaceGeometry, only: prepare_surface_geometry
    use AU_Types, only: au_case_type, allocate_au_case, bc_dirichlet_external
    use AU_BoundaryConditions, only: reset_au_fields, expand_particle_boundary_values
    use AU_Solver, only: solve_au_surface
    implicit none

    type(geometry_type) :: geometry
    type(au_case_type) :: case_data
    character(len=512) :: mesh_file
    character(len=512) :: output_file
    character(len=512) :: argument
    character(len=512) :: message
    character(len=32) :: formulation_mode
    real(dp) :: wavenumber
    real(dp) :: characteristic_length
    real(dp) :: source_position(3)
    integer :: exterior_domain_normal_sign
    integer :: status

    wavenumber = 1.0_dp
    characteristic_length = 1.0_dp
    exterior_domain_normal_sign = -1
    formulation_mode = "ordinary"
    source_position = 0.0_dp

    if (command_argument_count() < 2) then
        write(*, '(A)') "Usage: validate_manufactured_exterior mesh_file output_csv " // &
                        "[wavenumber] [characteristic_length] [exterior_domain_normal_sign] " // &
                        "[ordinary|burton-miller] [source_x] [source_y] [source_z]"
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
        read(argument, *) characteristic_length
    end if

    if (command_argument_count() >= 5) then
        call get_command_argument(5, argument)
        read(argument, *) exterior_domain_normal_sign
    end if

    if (command_argument_count() >= 6) then
        call get_command_argument(6, formulation_mode)
    end if

    if (command_argument_count() >= 7) then
        call get_command_argument(7, argument)
        read(argument, *) source_position(1)
    end if

    if (command_argument_count() >= 8) then
        call get_command_argument(8, argument)
        read(argument, *) source_position(2)
    end if

    if (command_argument_count() >= 9) then
        call get_command_argument(9, argument)
        read(argument, *) source_position(3)
    end if

    call read_mesh(trim(mesh_file), geometry, status, message)
    call require_ok(status, message)

    if (geometry%particle_count /= 1) then
        write(*, '(A)') "Manufactured exterior validation requires exactly one particle."
        error stop 3
    end if

    call prepare_mesh_topology(geometry, status, message)
    call require_ok(status, message)

    call prepare_surface_geometry(geometry, status, message)
    call require_ok(status, message)

    call allocate_au_case(case_data, geometry, status, message)
    call require_ok(status, message)

    call configure_case(case_data, wavenumber, characteristic_length, exterior_domain_normal_sign, &
                        trim(formulation_mode))

    call reset_au_fields(case_data, status, message)
    call require_ok(status, message)

    ! Start from the particle-level zero value, then replace it with the exact
    ! spatially varying value at every mesh node.
    call expand_particle_boundary_values(geometry, case_data, status, message)
    call require_ok(status, message)
    call set_point_source_boundary(geometry, case_data, wavenumber, source_position, &
                                   status, message)
    call require_ok(status, message)

    call solve_au_surface(geometry, case_data, status, message)
    call require_ok(status, message)

    call write_surface_csv(trim(output_file), geometry, case_data, status, message)
    call require_ok(status, message)

    write(*, '(A)') "Manufactured exterior solve completed with " // &
                    trim(formulation_mode) // " formulation."

contains

    subroutine configure_case(case_data, wavenumber, characteristic_length, exterior_domain_normal_sign, &
                              formulation_mode)
        ! Fill the allocated case with nondimensional exterior test data.
        ! Density and sound speed are one, so angular frequency equals k.
        type(au_case_type), intent(inout) :: case_data
        real(dp), intent(in) :: wavenumber
        real(dp), intent(in) :: characteristic_length
        integer, intent(in) :: exterior_domain_normal_sign
        character(len=*), intent(in) :: formulation_mode

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

        ! The current case has one impenetrable particle.  Interior properties
        ! are allocated by the common case type but are not used by this solve.
        case_data%interior_medium(1)%density = 1.0_dp
        case_data%interior_medium(1)%sound_speed = cmplx(1.0_dp, 0.0_dp, kind=dp)
        case_data%interior_medium(1)%wavenumber = cmplx(wavenumber, 0.0_dp, kind=dp)

        case_data%layer(1)%parent_particle_id = 0
        case_data%layer(1)%exterior_domain_normal_sign = exterior_domain_normal_sign
        case_data%layer(1)%characteristic_length = characteristic_length

        case_data%boundary_condition(1)%kind = bc_dirichlet_external
        case_data%boundary_condition(1)%prescribed_data = complex_zero
    end subroutine configure_case

    subroutine set_point_source_boundary(geometry, case_data, wavenumber, source_position, &
                                         status, message)
        ! Evaluate exp(i*k*r)/r at every node and prescribe it as total-field
        ! Dirichlet data.  The source must remain strictly inside the solid so
        ! the solved exterior domain contains no singularity.
        type(geometry_type), intent(in) :: geometry
        type(au_case_type), intent(inout) :: case_data
        real(dp), intent(in) :: wavenumber
        real(dp), intent(in) :: source_position(3)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        real(dp) :: displacement(3)
        real(dp) :: distance
        integer :: node_id

        status = 0
        message = "Manufactured point-source boundary values set."

        do node_id = 1, geometry%mesh%node_count
            displacement = geometry%mesh%xyz(:, node_id) - source_position
            distance = sqrt(dot_product(displacement, displacement))
            if (distance <= sqrt(epsilon(1.0_dp))) then
                status = 1
                message = "Point source is on or too close to a mesh node."
                return
            end if

            case_data%prescribed_boundary_data(node_id) = &
                exp(imaginary_unit * wavenumber * distance) / distance
        end do
    end subroutine set_point_source_boundary

    subroutine write_surface_csv(filename, geometry, case_data, status, message)
        ! Write all nodal fields needed for the independent Python derivative
        ! comparison.
        character(len=*), intent(in) :: filename
        type(geometry_type), intent(in) :: geometry
        type(au_case_type), intent(in) :: case_data
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: unit
        integer :: node_id

        status = 0
        message = "Manufactured exterior validation CSV written."

        open(newunit=unit, file=filename, status="replace", action="write", &
             iostat=status, iomsg=message)
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
    end subroutine write_surface_csv

    subroutine require_ok(status, message)
        integer, intent(in) :: status
        character(len=*), intent(in) :: message

        if (status /= 0) then
            write(*, '(A, I0)') "Validation failed with status ", status
            write(*, '(A)') trim(message)
            error stop 1
        end if
    end subroutine require_ok

end program ValidateManufacturedExterior
