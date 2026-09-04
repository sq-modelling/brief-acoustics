module AU_SurfaceOperators
    ! This module assembles the global boundary-integral operator matrices.
    !
    ! AU_HelmholtzKernels integrates one element at a time for one collocation
    ! node.  This module places those local element contributions into full
    ! node-by-node matrices:
    !
    !   external_g, external_h = operators for exterior acoustic domains
    !   internal_g, internal_h = operators for interior acoustic domains
    !
    ! A row corresponds to a collocation node where the boundary integral
    ! equation is enforced.  A column corresponds to a source node whose unknown
    ! or prescribed boundary data contributes to that equation.
    !
    ! The data structures also retain ordinary interior-domain blocks for the
    ! research path.  The public runner solves one exterior particle with
    ! either the ordinary or the Burton-Miller formulation.
    use Pre_Constants, only: dp, complex_zero
    use Geom_Types, only: geometry_type
    use AU_Types, only: au_case_type
    use AU_LayerTopology, only: bem_exterior_infinity_contribution, &
                                burton_miller_infinity_contribution, &
                                burton_miller_coupling_weight, exterior_domain_contains_particle, &
                                exterior_domain_wavenumber, interior_domain_contains_particle, &
                                interior_domain_wavenumber
    use AU_HelmholtzKernels, only: clear_helmholtz_element_integral, &
                                   helmholtz_element_integral_type, integrate_helmholtz_element
    implicit none

    private

    public :: au_surface_operator_type
    public :: assemble_au_surface_operators
    public :: clear_au_surface_operators

    type :: au_surface_operator_type
        ! Full dense operator matrices for one prepared acoustic case.
        !
        ! These matrices are dense because every boundary node generally
        ! interacts with every other boundary node in a boundary element method.

        ! External-domain single- and double-layer blocks.
        complex(dp), allocatable :: external_g(:, :)
        complex(dp), allocatable :: external_h(:, :)

        ! Internal-domain single- and double-layer blocks.
        complex(dp), allocatable :: internal_g(:, :)
        complex(dp), allocatable :: internal_h(:, :)

        ! Burton-Miller row coupling parameter stored per collocation node.
        complex(dp), allocatable :: bm_coupling_weight(:)

        ! Burton-Miller auxiliary blocks; their mathematical names are below.
        ! Regularisation and exterior-at-infinity contributions are included
        ! during assembly in this module.
        ! K = normal derivative of the Helmholtz single-layer kernel.
        complex(dp), allocatable :: bm_dg_dn0(:, :)
        ! D = mixed-normal derivative of the Helmholtz kernel.
        complex(dp), allocatable :: bm_d2g_dndn0(:, :)
        ! L = static auxiliary normal-derivative operator.
        complex(dp), allocatable :: bm_minus_dg0_dn0(:, :)
        ! G0 = static auxiliary single-layer operator.
        complex(dp), allocatable :: bm_g0(:, :)
        ! H0 = static auxiliary double-layer operator.
        complex(dp), allocatable :: bm_h0(:, :)
    contains
        ! Allows object-style cleanup: call operators%clear().
        procedure :: clear => clear_au_surface_operators_bound
    end type au_surface_operator_type

contains

    subroutine assemble_au_surface_operators(geometry, case_data, operators, status, message)
        ! Build all global surface operator matrices.
        !
        ! The assembly algorithm is:
        !   1. loop over collocation nodes, which define matrix rows
        !   2. loop over source elements
        !   3. decide whether the element belongs to the exterior or interior
        !      domain associated with the row particle
        !   4. integrate that element
        !   5. add local element-node contributions into global matrix columns
        type(geometry_type), intent(in) :: geometry
        type(au_case_type), intent(in) :: case_data
        type(au_surface_operator_type), intent(inout) :: operators
        integer, intent(out), optional :: status
        character(len=*), intent(out), optional :: message

        integer :: error_code
        character(len=256) :: error_message
        integer :: node_count
        integer :: element_count
        integer :: row_node
        integer :: element_id
        integer :: host_particle_id
        integer :: source_particle_id
        integer :: local_id
        integer :: column_node
        integer :: nodes_per_element
        logical :: same_particle
        type(helmholtz_element_integral_type) :: integral

        error_code = 0
        error_message = "AU surface operators assembled."

        ! Start from an empty operator object.
        call clear_au_surface_operators(operators)

        ! Check that geometry, quadrature, and acoustic layer data exist.
        call validate_inputs(geometry, case_data, error_code, error_message)
        if (error_code /= 0) then
            call finish(error_code, error_message)
            return
        end if

        node_count = geometry%mesh%node_count
        element_count = geometry%mesh%element_count
        nodes_per_element = geometry%mesh%nodes_per_element

        ! Allocate dense node-by-node matrices.
        call allocate_operators(operators, node_count, case_data%use_burton_miller, &
                                error_code, error_message)
        if (error_code /= 0) then
            call finish(error_code, error_message)
            return
        end if

        ! Each row is one collocation node.
        do row_node = 1, node_count
            host_particle_id = geometry%mesh%node_particle_id(row_node)

            ! Store the Burton-Miller row weight for this collocation node.
            if (case_data%use_burton_miller) then
                operators%bm_coupling_weight(row_node) = burton_miller_coupling_weight(case_data, host_particle_id)
            end if

            ! Each source element may contribute to the exterior operator, the
            ! interior operator, both, or neither, depending on layer topology.
            do element_id = 1, element_count
                source_particle_id = geometry%mesh%element_particle_id(element_id)
                same_particle = source_particle_id == host_particle_id

                if (exterior_domain_contains_particle(case_data, host_particle_id, source_particle_id)) then
                    ! Integrate this element using the wavenumber of the
                    ! exterior domain immediately outside the host particle.
                    call integrate_helmholtz_element(geometry, row_node, element_id, &
                                                     exterior_domain_wavenumber(case_data, host_particle_id), &
                                                     case_data%use_burton_miller, integral, &
                                                     error_code, error_message, same_particle)
                    if (error_code /= 0) then
                        call clear_helmholtz_element_integral(integral)
                        call clear_au_surface_operators(operators)
                        call finish(error_code, error_message)
                        return
                    end if

                    ! Scatter local element-node contributions into global
                    ! columns.  Each local_id maps to one global column node.
                    do local_id = 1, nodes_per_element
                        column_node = geometry%mesh%element_nodes(local_id, element_id)
                        operators%external_g(row_node, column_node) = &
                            operators%external_g(row_node, column_node) + integral%g(local_id)
                        operators%external_h(row_node, column_node) = &
                            operators%external_h(row_node, column_node) + integral%h(local_id)

                        if (case_data%use_burton_miller) then
                            operators%bm_dg_dn0(row_node, column_node) = &
                                operators%bm_dg_dn0(row_node, column_node) + integral%bm%dg_dn0(local_id)
                            operators%bm_d2g_dndn0(row_node, column_node) = &
                                operators%bm_d2g_dndn0(row_node, column_node) + integral%bm%d2g_dndn0(local_id)
                            operators%bm_minus_dg0_dn0(row_node, column_node) = &
                                operators%bm_minus_dg0_dn0(row_node, column_node) + integral%bm%minus_dg0_dn0(local_id)
                            operators%bm_g0(row_node, column_node) = &
                                operators%bm_g0(row_node, column_node) + integral%bm%g0(local_id)
                            operators%bm_h0(row_node, column_node) = &
                                operators%bm_h0(row_node, column_node) + integral%bm%h0(local_id)
                        end if
                    end do

                    ! Add same-particle non-singular correction terms to the
                    ! collocation diagonal.
                    operators%external_g(row_node, row_node) = &
                        operators%external_g(row_node, row_node) + integral%g_regularization
                    operators%external_h(row_node, row_node) = &
                        operators%external_h(row_node, row_node) + integral%h_regularization

                    if (case_data%use_burton_miller) then
                        operators%bm_dg_dn0(row_node, row_node) = &
                            operators%bm_dg_dn0(row_node, row_node) + integral%bm%dg0_dn_regularization
                        operators%bm_d2g_dndn0(row_node, row_node) = &
                            operators%bm_d2g_dndn0(row_node, row_node) + integral%bm%h_regularization
                        operators%bm_minus_dg0_dn0(row_node, row_node) = &
                            operators%bm_minus_dg0_dn0(row_node, row_node) + integral%bm%minus_dg0_dn_regularization
                        operators%bm_g0(row_node, row_node) = &
                            operators%bm_g0(row_node, row_node) + integral%bm%g0_regularization
                        operators%bm_h0(row_node, row_node) = &
                            operators%bm_h0(row_node, row_node) + integral%bm%h0_regularization
                    end if

                    call clear_helmholtz_element_integral(integral)
                end if

                if (interior_domain_contains_particle(case_data, host_particle_id, source_particle_id)) then
                    ! Integrate this element using the wavenumber of the
                    ! interior domain immediately inside the host particle.
                    ! Burton-Miller augmentation is only assembled for the
                    ! exterior operators in the current release.
                    call integrate_helmholtz_element(geometry, row_node, element_id, &
                                                     interior_domain_wavenumber(case_data, host_particle_id), &
                                                     .false., integral, error_code, error_message, same_particle)
                    if (error_code /= 0) then
                        call clear_helmholtz_element_integral(integral)
                        call clear_au_surface_operators(operators)
                        call finish(error_code, error_message)
                        return
                    end if

                    ! Scatter local element-node contributions into the internal
                    ! operator matrices.
                    do local_id = 1, nodes_per_element
                        column_node = geometry%mesh%element_nodes(local_id, element_id)
                        operators%internal_g(row_node, column_node) = &
                            operators%internal_g(row_node, column_node) + integral%g(local_id)
                        operators%internal_h(row_node, column_node) = &
                            operators%internal_h(row_node, column_node) + integral%h(local_id)
                    end do

                    operators%internal_g(row_node, row_node) = &
                        operators%internal_g(row_node, row_node) + integral%g_regularization
                    operators%internal_h(row_node, row_node) = &
                        operators%internal_h(row_node, row_node) + integral%h_regularization

                    call clear_helmholtz_element_integral(integral)
                end if
            end do

            ! Add the exterior-at-infinity contribution to the H diagonal.
            ! The local solid-angle coefficient has already cancelled in the
            ! BRIEF subtraction; this explicit 4*pi term comes from infinity.
            operators%external_h(row_node, row_node) = &
                operators%external_h(row_node, row_node) + &
                bem_exterior_infinity_contribution(case_data, host_particle_id)

            if (case_data%use_burton_miller) then
                ! Add the corresponding exterior-at-infinity contributions to
                ! the Burton-Miller auxiliary diagonal blocks.
                operators%bm_dg_dn0(row_node, row_node) = &
                    operators%bm_dg_dn0(row_node, row_node) - &
                    burton_miller_infinity_contribution(case_data, host_particle_id)
                operators%bm_minus_dg0_dn0(row_node, row_node) = &
                    operators%bm_minus_dg0_dn0(row_node, row_node) + &
                    burton_miller_infinity_contribution(case_data, host_particle_id)
                operators%bm_h0(row_node, row_node) = &
                    operators%bm_h0(row_node, row_node) + &
                    burton_miller_infinity_contribution(case_data, host_particle_id)
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

    end subroutine assemble_au_surface_operators

    subroutine clear_au_surface_operators_bound(this)
        ! Type-bound wrapper so an operator object can call operators%clear().
        class(au_surface_operator_type), intent(inout) :: this

        call clear_au_surface_operators(this)
    end subroutine clear_au_surface_operators_bound

    subroutine clear_au_surface_operators(operators)
        ! Release all allocated operator matrices.
        type(au_surface_operator_type), intent(inout) :: operators

        ! Ordinary exterior and interior operators.
        if (allocated(operators%external_g)) deallocate(operators%external_g)
        if (allocated(operators%external_h)) deallocate(operators%external_h)
        if (allocated(operators%internal_g)) deallocate(operators%internal_g)
        if (allocated(operators%internal_h)) deallocate(operators%internal_h)

        ! Burton-Miller auxiliary operators.
        if (allocated(operators%bm_coupling_weight)) deallocate(operators%bm_coupling_weight)
        if (allocated(operators%bm_dg_dn0)) deallocate(operators%bm_dg_dn0)
        if (allocated(operators%bm_d2g_dndn0)) deallocate(operators%bm_d2g_dndn0)
        if (allocated(operators%bm_minus_dg0_dn0)) deallocate(operators%bm_minus_dg0_dn0)
        if (allocated(operators%bm_g0)) deallocate(operators%bm_g0)
        if (allocated(operators%bm_h0)) deallocate(operators%bm_h0)
    end subroutine clear_au_surface_operators

    subroutine allocate_operators(operators, node_count, use_burton_miller, status, message)
        ! Allocate dense operator matrices for a mesh with node_count nodes.
        type(au_surface_operator_type), intent(inout) :: operators
        integer, intent(in) :: node_count
        logical, intent(in) :: use_burton_miller
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        status = 0
        message = "AU surface operator arrays allocated."

        ! Ordinary BEM blocks are always needed.
        allocate(operators%external_g(node_count, node_count), &
                 operators%external_h(node_count, node_count), &
                 operators%internal_g(node_count, node_count), &
                 operators%internal_h(node_count, node_count), &
                 stat=status)
        if (status /= 0) then
            message = "Unable to allocate AU surface operator arrays."
            call clear_au_surface_operators(operators)
            return
        end if

        operators%external_g = complex_zero
        operators%external_h = complex_zero
        operators%internal_g = complex_zero
        operators%internal_h = complex_zero

        if (use_burton_miller) then
            ! Burton-Miller blocks are allocated only when requested.
            allocate(operators%bm_coupling_weight(node_count), &
                     operators%bm_dg_dn0(node_count, node_count), &
                     operators%bm_d2g_dndn0(node_count, node_count), &
                     operators%bm_minus_dg0_dn0(node_count, node_count), &
                     operators%bm_g0(node_count, node_count), &
                     operators%bm_h0(node_count, node_count), &
                     stat=status)
            if (status /= 0) then
                message = "Unable to allocate AU Burton-Miller operator arrays."
                call clear_au_surface_operators(operators)
                return
            end if

            operators%bm_coupling_weight = complex_zero
            operators%bm_dg_dn0 = complex_zero
            operators%bm_d2g_dndn0 = complex_zero
            operators%bm_minus_dg0_dn0 = complex_zero
            operators%bm_g0 = complex_zero
            operators%bm_h0 = complex_zero
        end if
    end subroutine allocate_operators

    subroutine validate_inputs(geometry, case_data, status, message)
        ! Check that the geometry and case contain the arrays needed for
        ! operator assembly.
        type(geometry_type), intent(in) :: geometry
        type(au_case_type), intent(in) :: case_data
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        status = 0
        message = "AU surface operator inputs are valid."

        if (geometry%mesh%node_count <= 0) then
            call fail("Geometry has no nodes.")
        else if (geometry%mesh%element_count <= 0) then
            call fail("Geometry has no elements.")
        else if (.not. allocated(geometry%mesh%node_particle_id)) then
            call fail("Mesh node particle ids are not allocated.")
        else if (.not. allocated(geometry%mesh%element_particle_id)) then
            call fail("Mesh element particle ids are not allocated.")
        else if (.not. allocated(geometry%mesh%element_nodes)) then
            call fail("Mesh element connectivity is not allocated.")
        else if (.not. allocated(geometry%quadrature%point_xyz)) then
            call fail("Geometry quadrature cache is not allocated.")
        else if (.not. allocated(case_data%layer)) then
            call fail("AU layer topology is not allocated.")
        else if (.not. allocated(case_data%interior_medium)) then
            call fail("AU interior media are not allocated.")
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

    end subroutine validate_inputs

end module AU_SurfaceOperators
