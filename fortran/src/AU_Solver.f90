module AU_Solver
    ! This module builds and solves the final acoustic boundary-element linear
    ! system.
    !
    ! Earlier modules prepare geometry, boundary conditions, and dense surface
    ! operator matrices.  This module decides which boundary quantity is unknown
    ! at each node, moves known terms to the right-hand side, solves the dense
    ! complex system, and writes the result back into the case_data field arrays.
    !
    ! Validated public v1.0 scope:
    !   - one closed exterior particle
    !   - Dirichlet, Neumann, or Robin boundary conditions
    !   - ordinary NSBEM or the audited non-singular Burton-Miller system
    !
    ! An ordinary two-particle path is retained for the explicitly labelled
    ! small-gap demonstration.  Transmission, nested-layer, and multiparticle
    ! Burton-Miller systems are not public capabilities.
    use Pre_Constants, only: dp, complex_zero
    use Geom_Types, only: geometry_type
    use AU_Types, only: au_case_type, validate_au_case, &
                        bc_dirichlet_external, bc_neumann_external, bc_robin_external, &
                        bc_dirichlet_internal, bc_neumann_internal, bc_robin_internal, &
                        bc_transmission
    use AU_LayerTopology, only: exterior_domain_contains_particle, interior_domain_contains_particle
    use AU_SurfaceOperators, only: au_surface_operator_type, assemble_au_surface_operators, &
                                   clear_au_surface_operators
    use AU_BoundaryConditions, only: update_total_fields_and_pressure
    implicit none

    private

    public :: solve_au_surface
    public :: assemble_au_linear_system

    ! Which side of a particle surface the equation is written for.
    integer, parameter :: side_external = 1
    integer, parameter :: side_internal = 2

    ! Which nodal field is the unknown in the final linear system.
    integer, parameter :: unknown_phi = 1
    integer, parameter :: unknown_dphi_dn = 2

    ! Public input tolerance for deciding whether a Robin coefficient is zero.
    ! This matches the Python and runtime-case validation threshold.
    real(dp), parameter :: tiny_robin = 1.0e-14_dp

    interface
        ! LAPACK dense complex linear solver.
        !
        ! zgesv solves A * x = b for a complex dense matrix A using LU
        ! factorisation with pivoting.  On macOS this is provided through the
        ! Accelerate framework in the validation build.
        subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import dp
            integer, intent(in) :: n
            integer, intent(in) :: nrhs
            integer, intent(in) :: lda
            integer, intent(in) :: ldb
            integer, intent(out) :: ipiv(*)
            integer, intent(out) :: info
            complex(dp), intent(inout) :: a(lda, *)
            complex(dp), intent(inout) :: b(ldb, *)
        end subroutine zgesv
    end interface

contains

    subroutine solve_au_surface(geometry, case_data, status, message)
        ! Solve the acoustic surface problem for one geometry and case setup.
        !
        ! High-level workflow:
        !   1. validate geometry and case data
        !   2. choose one physical unknown per node
        !   3. assemble dense boundary-integral operators
        !   4. assemble the final linear system
        !   5. solve the dense complex system
        !   6. load nodal phi, dphi/dn, and any auxiliary values into case_data
        !   7. update total fields and pressures
        type(geometry_type), intent(in) :: geometry
        type(au_case_type), intent(inout) :: case_data
        integer, intent(out), optional :: status
        character(len=*), intent(out), optional :: message

        type(au_surface_operator_type) :: operators
        complex(dp), allocatable :: system_matrix(:, :)
        complex(dp), allocatable :: rhs(:)
        complex(dp), allocatable :: solution(:)
        integer, allocatable :: node_unknown(:)
        integer :: error_code
        character(len=256) :: error_message

        error_code = 0
        error_message = "AU surface problem solved."

        ! Check that the case is complete and within the current release scope.
        call validate_solver_inputs(geometry, case_data, error_code, error_message)
        if (error_code /= 0) then
            call finish(error_code, error_message)
            return
        end if

        ! Map each mesh node to one unknown index in the dense linear system.
        call build_unknown_map(geometry, case_data, node_unknown, error_code, error_message)
        if (error_code /= 0) then
            call finish(error_code, error_message)
            return
        end if

        ! Assemble G/H operator matrices from element integrals.
        call assemble_au_surface_operators(geometry, case_data, operators, error_code, error_message)
        if (error_code /= 0) then
            call cleanup()
            call finish(error_code, error_message)
            return
        end if

        ! Convert operator matrices plus boundary conditions into A*x=b.
        call assemble_au_linear_system(geometry, case_data, operators, node_unknown, &
                                       system_matrix, rhs, error_code, error_message)
        if (error_code /= 0) then
            call cleanup()
            call finish(error_code, error_message)
            return
        end if

        ! Solve A*x=b.
        call solve_dense_system(system_matrix, rhs, solution, error_code, error_message)
        if (error_code /= 0) then
            call cleanup()
            call finish(error_code, error_message)
            return
        end if

        ! Convert the solution vector back to physical field arrays.
        call load_solution(geometry, case_data, node_unknown, solution, error_code, error_message)
        if (error_code /= 0) then
            call cleanup()
            call finish(error_code, error_message)
            return
        end if

        ! Reconstruct total field and pressure for post-processing.
        call update_total_fields_and_pressure(geometry, case_data, error_code, error_message)
        call cleanup()
        call finish(error_code, error_message)

    contains

        subroutine cleanup()
            ! Release local work arrays and operator matrices before returning.
            if (allocated(system_matrix)) deallocate(system_matrix)
            if (allocated(rhs)) deallocate(rhs)
            if (allocated(solution)) deallocate(solution)
            if (allocated(node_unknown)) deallocate(node_unknown)
            call clear_au_surface_operators(operators)
        end subroutine cleanup

        subroutine finish(code, text)
            ! Return optional status information to the caller.
            integer, intent(in) :: code
            character(len=*), intent(in) :: text

            if (present(status)) status = code
            if (present(message)) message = text
        end subroutine finish

    end subroutine solve_au_surface

    subroutine validate_solver_inputs(geometry, case_data, status, message)
        ! Validate inputs and reject cases outside the current release scope.
        type(geometry_type), intent(in) :: geometry
        type(au_case_type), intent(in) :: case_data
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: boundary_kind
        integer :: particle_id
        real(dp) :: imaginary_tolerance
        real(dp) :: real_wavenumber
        complex(dp) :: wavenumber

        call validate_au_case(case_data, geometry, status, message)
        if (status /= 0) return

        ! Transmission problems need coupled external/internal unknowns.  The
        ! current release solves one unknown per node, so transmission is rejected.
        do particle_id = 1, geometry%particle_count
            if (case_data%boundary_condition(particle_id)%kind == bc_transmission) then
                status = 1
                message = "AU solver linear-system assembly does not yet support transmission BCs."
                return
            end if
        end do

        if (.not. case_data%use_burton_miller) return

        ! The augmented equations have been audited only for one closed surface
        ! in an unbounded exterior acoustic domain.  Rejecting broader cases is
        ! deliberate: the auxiliary Laplace field has not yet been derived and
        ! validated for disconnected particles or nested layers.
        if (geometry%particle_count /= 1) then
            status = 1
            message = "The current Burton-Miller release scope requires exactly one exterior particle."
            return
        end if

        if (case_data%layer(1)%parent_particle_id /= 0) then
            status = 1
            message = "The current Burton-Miller release scope does not support nested particle layers."
            return
        end if

        ! The first certified path assumes the mesh nodes are ordered so the
        ! geometric normals point outward from the solid.  In this convention
        ! the acoustic exterior-domain normal sign is -1.
        if (case_data%layer(1)%exterior_domain_normal_sign /= -1) then
            status = 1
            message = "The current Burton-Miller release requires an outward-from-solid mesh and exterior_domain_normal_sign=-1."
            return
        end if

        boundary_kind = case_data%boundary_condition(1)%kind
        select case (boundary_kind)
        case (bc_dirichlet_external, bc_neumann_external, bc_robin_external)
            continue
        case default
            status = 1
            message = "The current Burton-Miller release supports exterior Dirichlet, Neumann, or Robin conditions only."
            return
        end select

        ! The coupling length beta=min(a,1/k) is currently validated only for a
        ! real positive acoustic wavenumber.  Complex Helmholtz and the k=0
        ! Laplace limit need their own coupling policy and tests.
        wavenumber = case_data%exterior_medium%wavenumber
        real_wavenumber = real(wavenumber, kind=dp)
        imaginary_tolerance = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(real_wavenumber))
        if (real_wavenumber <= 100.0_dp * epsilon(1.0_dp) .or. &
            abs(aimag(wavenumber)) > imaginary_tolerance) then
            status = 1
            message = "The current Burton-Miller release requires a real, positive exterior wavenumber."
            return
        end if
    end subroutine validate_solver_inputs

    subroutine build_unknown_map(geometry, case_data, node_unknown, status, message)
        ! Assign one linear-system unknown index to each node.
        !
        ! The current release has exactly one unknown per node.  The identity map
        ! node_unknown(node_id)=node_id is therefore sufficient.
        type(geometry_type), intent(in) :: geometry
        type(au_case_type), intent(in) :: case_data
        integer, allocatable, intent(out) :: node_unknown(:)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: node_id

        status = 0
        message = "AU unknown map built."

        allocate(node_unknown(geometry%mesh%node_count), stat=status)
        if (status /= 0) then
            message = "Unable to allocate AU unknown map."
            return
        end if

        do node_id = 1, geometry%mesh%node_count
            node_unknown(node_id) = node_id
            ! Confirm that the boundary condition really leaves one valid
            ! unknown field at this node.
            if (unknown_variable(case_data, geometry%mesh%node_particle_id(node_id)) == 0) then
                status = 1
                message = "Unsupported boundary condition while building AU unknown map."
                return
            end if
        end do
    end subroutine build_unknown_map

    subroutine assemble_au_linear_system(geometry, case_data, operators, node_unknown, &
                                         system_matrix, rhs, status, message)
        ! Assemble the dense complex linear system A*x=b.
        !
        ! This lower-level routine is public so focused tests and future Python
        ! bindings can inspect the assembled matrix before factorisation.  A
        ! normal solve should still enter through solve_au_surface, which first
        ! validates the supported physical scope and assembles the operators.
        !
        ! For each collocation row, we add contributions from all source nodes
        ! that belong to the same acoustic domain. Prescribed boundary data are
        ! moved to rhs, while unknown field values are placed in A.
        type(geometry_type), intent(in) :: geometry
        type(au_case_type), intent(in) :: case_data
        type(au_surface_operator_type), intent(in) :: operators
        integer, intent(in) :: node_unknown(:)
        complex(dp), allocatable, intent(out) :: system_matrix(:, :)
        complex(dp), allocatable, intent(out) :: rhs(:)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: auxiliary_column
        integer :: auxiliary_row
        integer :: node_count
        integer :: row_node
        integer :: column_node
        integer :: row_particle
        integer :: column_particle
        integer :: row_side
        integer :: column_side
        integer :: system_size
        complex(dp) :: h_value
        complex(dp) :: g_value
        complex(dp) :: weight

        status = 0
        message = "AU linear system assembled."
        node_count = geometry%mesh%node_count

        if (case_data%use_burton_miller) then
            call validate_burton_miller_operator_blocks(operators, node_count, status, message)
            if (status /= 0) return
        end if

        system_size = node_count
        if (case_data%use_burton_miller) system_size = 2 * node_count

        allocate(system_matrix(system_size, system_size), rhs(system_size), stat=status)
        if (status /= 0) then
            message = "Unable to allocate AU linear system."
            return
        end if
        system_matrix = complex_zero
        rhs = complex_zero

        ! Loop over collocation equations.
        do row_node = 1, node_count
            row_particle = geometry%mesh%node_particle_id(row_node)

            ! Choose whether this row uses the exterior or interior equation.
            row_side = equation_side(case_data%boundary_condition(row_particle)%kind)
            if (row_side == 0) then
                status = 1
                message = "Unsupported row boundary condition in AU linear system."
                return
            end if

            ! Add contributions from source nodes in the relevant domain.
            do column_node = 1, node_count
                column_particle = geometry%mesh%node_particle_id(column_node)
                if (.not. domain_contains_source(case_data, row_side, row_particle, column_particle)) cycle

                ! Decide which side of the source particle this row sees.
                column_side = source_side_for_domain(case_data, row_side, row_particle, column_particle)
                if (column_side == 0) then
                    status = 1
                    message = "Unable to determine source side for AU domain coupling."
                    return
                end if

                ! In Burton-Miller mode the first block row is the ordinary
                ! exterior equation plus weight times its normal derivative:
                !
                !   (H + W*D)*phi - (G + W*K)*q - W*L*s = 0.
                !
                ! W multiplies rows, so the collocation node supplies weight.
                if (case_data%use_burton_miller) then
                    weight = operators%bm_coupling_weight(row_node)
                    h_value = operators%external_h(row_node, column_node) + &
                              weight * operators%bm_d2g_dndn0(row_node, column_node)
                    g_value = operators%external_g(row_node, column_node) + &
                              weight * operators%bm_dg_dn0(row_node, column_node)
                else if (row_side == side_external) then
                    h_value = operators%external_h(row_node, column_node)
                    g_value = operators%external_g(row_node, column_node)
                else
                    h_value = operators%internal_h(row_node, column_node)
                    g_value = operators%internal_g(row_node, column_node)
                end if

                ! Add this source node's contribution to one row of A and b.
                call add_boundary_contribution(case_data, column_node, column_particle, column_side, &
                                               h_value, g_value, node_unknown(column_node), &
                                               system_matrix(node_unknown(row_node), :), &
                                               rhs(node_unknown(row_node)), status, message)
                if (status /= 0) return

                if (case_data%use_burton_miller) then
                    auxiliary_row = node_count + row_node
                    auxiliary_column = node_count + column_node

                    ! Top-right block: -W*L, where s is the solved auxiliary
                    ! normal derivative d(psi_1)/dn.
                    system_matrix(node_unknown(row_node), auxiliary_column) = &
                        system_matrix(node_unknown(row_node), auxiliary_column) - &
                        weight * operators%bm_minus_dg0_dn0(row_node, column_node)

                    ! Bottom-left block and right-hand side come from H0*phi.
                    ! Passing G=0 through the same boundary-condition helper
                    ! applies Dirichlet, Neumann, and Robin elimination without
                    ! a second set of hand-written sign branches.
                    call add_boundary_contribution(case_data, column_node, column_particle, column_side, &
                                                   operators%bm_h0(row_node, column_node), complex_zero, &
                                                   node_unknown(column_node), &
                                                   system_matrix(auxiliary_row, :), rhs(auxiliary_row), &
                                                   status, message)
                    if (status /= 0) return

                    ! Bottom-right block: -G0*s.  Together with the preceding
                    ! term this is the auxiliary equation H0*phi-G0*s=0.
                    system_matrix(auxiliary_row, auxiliary_column) = &
                        system_matrix(auxiliary_row, auxiliary_column) - &
                        operators%bm_g0(row_node, column_node)
                end if
            end do
        end do
    end subroutine assemble_au_linear_system

    subroutine validate_burton_miller_operator_blocks(operators, node_count, status, message)
        ! Confirm that all six augmented-system inputs have the expected size.
        ! This protects the public low-level assembler from missing or stale
        ! operator storage.
        type(au_surface_operator_type), intent(in) :: operators
        integer, intent(in) :: node_count
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        status = 0
        message = "Burton-Miller operator blocks are allocated and node-sized."

        if (.not. allocated(operators%bm_coupling_weight) .or. &
            .not. allocated(operators%bm_dg_dn0) .or. &
            .not. allocated(operators%bm_d2g_dndn0) .or. &
            .not. allocated(operators%bm_minus_dg0_dn0) .or. &
            .not. allocated(operators%bm_g0) .or. &
            .not. allocated(operators%bm_h0)) then
            status = 1
            message = "Burton-Miller assembly requires all auxiliary operator blocks."
            return
        end if

        if (size(operators%bm_coupling_weight) /= node_count .or. &
            .not. is_node_square(operators%bm_dg_dn0, node_count) .or. &
            .not. is_node_square(operators%bm_d2g_dndn0, node_count) .or. &
            .not. is_node_square(operators%bm_minus_dg0_dn0, node_count) .or. &
            .not. is_node_square(operators%bm_g0, node_count) .or. &
            .not. is_node_square(operators%bm_h0, node_count)) then
            status = 1
            message = "Burton-Miller auxiliary operator blocks must all have node-by-node dimensions."
        end if
    end subroutine validate_burton_miller_operator_blocks

    pure logical function is_node_square(matrix, node_count)
        complex(dp), intent(in) :: matrix(:, :)
        integer, intent(in) :: node_count

        is_node_square = size(matrix, 1) == node_count .and. size(matrix, 2) == node_count
    end function is_node_square

    subroutine add_boundary_contribution(case_data, node_id, particle_id, side_id, h_value, g_value, &
                                         unknown_id, matrix_row, rhs_value, status, message)
        ! Add one source-node contribution to a linear-system row.
        !
        ! The boundary integral equation contains both phi and dphi/dn.  A
        ! boundary condition tells us which one is known and how to express the
        ! other.  This routine uses that relation to decide what belongs in the
        ! matrix and what belongs on the right-hand side.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: node_id
        integer, intent(in) :: particle_id
        integer, intent(in) :: side_id
        complex(dp), intent(in) :: h_value
        complex(dp), intent(in) :: g_value
        integer, intent(in) :: unknown_id
        complex(dp), intent(inout) :: matrix_row(:)
        complex(dp), intent(inout) :: rhs_value
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: variable_id
        complex(dp) :: known_value
        complex(dp) :: coefficient

        status = 0
        message = "Boundary contribution added."

        variable_id = unknown_variable_for_side(case_data, particle_id, side_id)
        if (variable_id == 0) then
            status = 1
            message = "Boundary condition does not provide the requested AU field side."
            return
        end if

        select case (variable_id)
        case (unknown_phi)
            ! The unknown is phi.  Express dphi/dn = known + coefficient*phi.
            call phi_unknown_relation(case_data, node_id, particle_id, side_id, known_value, coefficient)
            matrix_row(unknown_id) = matrix_row(unknown_id) + h_value - g_value * coefficient
            rhs_value = rhs_value + g_value * known_value
        case (unknown_dphi_dn)
            ! The unknown is dphi/dn.  Express phi = known + coefficient*dphi/dn.
            call dphi_dn_unknown_relation(case_data, node_id, particle_id, side_id, known_value, coefficient)
            matrix_row(unknown_id) = matrix_row(unknown_id) + h_value * coefficient - g_value
            rhs_value = rhs_value - h_value * known_value
        end select
    end subroutine add_boundary_contribution

    subroutine solve_dense_system(system_matrix, rhs, solution, status, message)
        ! Solve the dense complex system using LAPACK zgesv.
        !
        ! zgesv overwrites its matrix argument.  We preserve system_matrix and
        ! factor a private copy so the solved system can be residual-checked.
        complex(dp), intent(in) :: system_matrix(:, :)
        complex(dp), intent(in) :: rhs(:)
        complex(dp), allocatable, intent(out) :: solution(:)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        complex(dp), allocatable :: factor_matrix(:, :)
        complex(dp), allocatable :: rhs_matrix(:, :)
        integer, allocatable :: pivot(:)
        integer :: n
        integer :: info

        status = 0
        message = "Dense AU linear system solved."
        n = size(rhs)

        if (size(system_matrix, 1) /= n .or. size(system_matrix, 2) /= n) then
            status = 1
            message = "AU dense system matrix and right-hand side sizes do not match."
            return
        end if

        allocate(factor_matrix(n, n), rhs_matrix(n, 1), pivot(n), solution(n), stat=status)
        if (status /= 0) then
            message = "Unable to allocate dense AU solve work arrays."
            return
        end if

        ! LAPACK expects the right-hand side as a 2D matrix with nrhs columns.
        factor_matrix = system_matrix
        rhs_matrix(:, 1) = rhs
        call zgesv(n, 1, factor_matrix, n, pivot, rhs_matrix, n, info)
        if (info /= 0) then
            status = info
            message = "LAPACK zgesv failed while solving AU linear system."
            deallocate(factor_matrix, rhs_matrix, pivot, solution)
            return
        end if

        solution = rhs_matrix(:, 1)
        call check_dense_solution(system_matrix, rhs, solution, status, message)
        deallocate(factor_matrix, rhs_matrix, pivot)
    end subroutine solve_dense_system

    subroutine check_dense_solution(system_matrix, rhs, solution, status, message)
        ! Check the componentwise backward residual of A*x=b.
        !
        ! For row i we divide |(A*x-b)_i| by
        !
        !   |b_i| + sum_j |A_ij| |x_j|.
        !
        ! This scaling works for rows with very different physical magnitudes,
        ! as can occur between the two Burton-Miller block equations.
        complex(dp), intent(in) :: system_matrix(:, :)
        complex(dp), intent(in) :: rhs(:)
        complex(dp), intent(in) :: solution(:)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        complex(dp) :: residual
        integer :: column_id
        integer :: n
        integer :: row_id
        real(dp) :: maximum_scaled_residual
        real(dp) :: row_scale
        real(dp) :: scaled_residual
        real(dp) :: tolerance

        status = 0
        message = "Dense AU solution residual checked."
        n = size(rhs)
        maximum_scaled_residual = 0.0_dp

        do row_id = 1, n
            residual = -rhs(row_id)
            row_scale = abs(rhs(row_id))
            do column_id = 1, n
                residual = residual + system_matrix(row_id, column_id) * solution(column_id)
                row_scale = row_scale + abs(system_matrix(row_id, column_id)) * abs(solution(column_id))
            end do

            if (row_scale > tiny(1.0_dp)) then
                scaled_residual = abs(residual) / row_scale
            else
                scaled_residual = abs(residual)
            end if
            maximum_scaled_residual = max(maximum_scaled_residual, scaled_residual)
        end do

        tolerance = 10000.0_dp * epsilon(1.0_dp) * real(max(1, n), dp)
        if (.not. maximum_scaled_residual <= tolerance) then
            status = 1
            write(message, '(A,ES12.4,A,ES12.4)') &
                "Dense AU solve residual ", maximum_scaled_residual, " exceeds tolerance ", tolerance
        end if
    end subroutine check_dense_solution

    subroutine load_solution(geometry, case_data, node_unknown, solution, status, message)
        ! Copy the solved unknown vector back into physical field arrays.
        !
        ! Depending on the boundary condition, the solved unknown may be phi or
        ! dphi/dn. The other field is reconstructed from the prescribed
        ! boundary relation.
        type(geometry_type), intent(in) :: geometry
        type(au_case_type), intent(inout) :: case_data
        integer, intent(in) :: node_unknown(:)
        complex(dp), intent(in) :: solution(:)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: node_id
        integer :: node_count
        integer :: particle_id
        integer :: side_id
        integer :: variable_id
        complex(dp) :: unknown_value
        complex(dp) :: known_value
        complex(dp) :: coefficient

        status = 0
        message = "AU solution loaded into field arrays."

        node_count = geometry%mesh%node_count
        if (case_data%use_burton_miller) then
            if (size(solution) /= 2 * node_count) then
                status = 1
                message = "Burton-Miller solution vector must contain 2N values."
                return
            end if
        else if (size(solution) /= node_count) then
            status = 1
            message = "Ordinary AU solution vector must contain N values."
            return
        end if

        case_data%solution%external_phi = complex_zero
        case_data%solution%external_dphi_dn = complex_zero
        case_data%solution%internal_phi = complex_zero
        case_data%solution%internal_dphi_dn = complex_zero
        case_data%solution%bm_auxiliary_dpsi_dn = complex_zero

        do node_id = 1, geometry%mesh%node_count
            particle_id = geometry%mesh%node_particle_id(node_id)
            side_id = equation_side(case_data%boundary_condition(particle_id)%kind)
            variable_id = unknown_variable_for_side(case_data, particle_id, side_id)
            unknown_value = solution(node_unknown(node_id))

            select case (variable_id)
            case (unknown_phi)
                ! Solution vector stores phi; reconstruct dphi/dn.
                call phi_unknown_relation(case_data, node_id, particle_id, side_id, known_value, coefficient)
                call set_phi_solution(case_data, node_id, side_id, unknown_value)
                call set_dphi_dn_solution(case_data, node_id, side_id, known_value + coefficient * unknown_value)
            case (unknown_dphi_dn)
                ! Solution vector stores dphi/dn; reconstruct phi.
                call dphi_dn_unknown_relation(case_data, node_id, particle_id, side_id, known_value, coefficient)
                call set_phi_solution(case_data, node_id, side_id, known_value + coefficient * unknown_value)
                call set_dphi_dn_solution(case_data, node_id, side_id, unknown_value)
            case default
                status = 1
                message = "Unable to load solution for unsupported AU boundary condition."
                return
            end select
        end do

        if (case_data%use_burton_miller) then
            case_data%solution%bm_auxiliary_dpsi_dn = solution(node_count + 1:2 * node_count)
        end if
    end subroutine load_solution

    pure integer function equation_side(boundary_kind)
        ! Decide whether a boundary-condition kind uses an exterior or interior
        ! boundary integral equation.
        integer, intent(in) :: boundary_kind

        select case (boundary_kind)
        case (bc_dirichlet_external, bc_neumann_external, bc_robin_external)
            equation_side = side_external
        case (bc_dirichlet_internal, bc_neumann_internal, bc_robin_internal)
            equation_side = side_internal
        case default
            equation_side = 0
        end select
    end function equation_side

    pure integer function unknown_variable(case_data, particle_id)
        ! Return the unknown field for a particle's own equation side.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: particle_id

        unknown_variable = unknown_variable_for_side(case_data, particle_id, &
                                                     equation_side(case_data%boundary_condition(particle_id)%kind))
    end function unknown_variable

    pure integer function unknown_variable_for_side(case_data, particle_id, side_id)
        ! Return which field is unknown on a requested side of a particle.
        !
        ! Dirichlet: phi is known, so dphi/dn is unknown.
        ! Neumann:   dphi/dn is known, so phi is unknown.
        ! Robin:     choose an algebraically available representation from its
        !            non-zero coefficients.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: particle_id
        integer, intent(in) :: side_id

        integer :: boundary_kind

        boundary_kind = case_data%boundary_condition(particle_id)%kind
        unknown_variable_for_side = 0

        if (side_id == side_external) then
            select case (boundary_kind)
            case (bc_dirichlet_external)
                unknown_variable_for_side = unknown_dphi_dn
            case (bc_neumann_external)
                unknown_variable_for_side = unknown_phi
            case (bc_robin_external)
                unknown_variable_for_side = robin_unknown_variable(case_data, particle_id)
            end select
        else if (side_id == side_internal) then
            select case (boundary_kind)
            case (bc_dirichlet_internal)
                unknown_variable_for_side = unknown_dphi_dn
            case (bc_neumann_internal)
                unknown_variable_for_side = unknown_phi
            case (bc_robin_internal)
                unknown_variable_for_side = robin_unknown_variable(case_data, particle_id)
            end select
        end if
    end function unknown_variable_for_side

    pure integer function robin_unknown_variable(case_data, particle_id)
        ! Choose the Robin unknown.
        !
        ! Robin condition:
        !   a*phi + b*dphi/dn = rhs
        !
        ! If b is non-zero, solve for phi and express dphi/dn in terms of phi.
        ! If b is zero but a is non-zero, solve for dphi/dn and express phi in
        ! terms of dphi/dn.  If both are zero, the condition is invalid.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: particle_id

        if (abs(case_data%boundary_condition(particle_id)%robin_b) > tiny_robin) then
            robin_unknown_variable = unknown_phi
        else if (abs(case_data%boundary_condition(particle_id)%robin_a) > tiny_robin) then
            robin_unknown_variable = unknown_dphi_dn
        else
            robin_unknown_variable = 0
        end if
    end function robin_unknown_variable

    pure logical function domain_contains_source(case_data, row_side, row_particle, source_particle)
        ! Wrapper around layer-topology queries for one matrix row.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: row_side
        integer, intent(in) :: row_particle
        integer, intent(in) :: source_particle

        if (row_side == side_external) then
            domain_contains_source = exterior_domain_contains_particle(case_data, row_particle, source_particle)
        else if (row_side == side_internal) then
            domain_contains_source = interior_domain_contains_particle(case_data, row_particle, source_particle)
        else
            domain_contains_source = .false.
        end if
    end function domain_contains_source

    pure integer function source_side_for_domain(case_data, row_side, row_particle, source_particle)
        ! Determine whether a source particle contributes through its exterior
        ! or interior side for the current row domain.
        !
        ! Example:
        !   If a row is in the exterior domain of a child particle, the parent
        !   surface bounds that same region from its interior side.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: row_side
        integer, intent(in) :: row_particle
        integer, intent(in) :: source_particle

        integer :: parent_id

        source_side_for_domain = 0
        if (row_side == side_external) then
            parent_id = case_data%layer(row_particle)%parent_particle_id
            if (source_particle == parent_id) then
                source_side_for_domain = side_internal
            else
                source_side_for_domain = side_external
            end if
        else if (row_side == side_internal) then
            if (source_particle == row_particle) then
                source_side_for_domain = side_internal
            else
                source_side_for_domain = side_external
            end if
        end if
    end function source_side_for_domain

    subroutine phi_unknown_relation(case_data, node_id, particle_id, side_id, known_value, coefficient)
        ! Build the relation for cases where phi is the unknown:
        !
        !   dphi_sca/dn = known_value + coefficient * phi_sca
        !
        ! Incident-field terms are subtracted because the solver unknowns are
        ! scattered fields, while boundary conditions are applied to total
        ! fields.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: node_id
        integer, intent(in) :: particle_id
        integer, intent(in) :: side_id
        complex(dp), intent(out) :: known_value
        complex(dp), intent(out) :: coefficient

        complex(dp) :: robin_a
        complex(dp) :: robin_b
        complex(dp) :: robin_rhs

        coefficient = complex_zero
        known_value = complex_zero

        select case (case_data%boundary_condition(particle_id)%kind)
        case (bc_neumann_external, bc_neumann_internal)
            known_value = case_data%prescribed_boundary_data(node_id) - incident_dphi_dn(case_data, node_id, side_id)
        case (bc_robin_external, bc_robin_internal)
            robin_a = case_data%boundary_condition(particle_id)%robin_a
            robin_b = case_data%boundary_condition(particle_id)%robin_b
            robin_rhs = case_data%boundary_condition(particle_id)%robin_rhs
            coefficient = -robin_a / robin_b
            known_value = (robin_rhs - robin_a * incident_phi(case_data, node_id, side_id)) / robin_b - &
                          incident_dphi_dn(case_data, node_id, side_id)
        end select
    end subroutine phi_unknown_relation

    subroutine dphi_dn_unknown_relation(case_data, node_id, particle_id, side_id, known_value, coefficient)
        ! Build the relation for cases where dphi/dn is the unknown:
        !
        !   phi_sca = known_value + coefficient * dphi_sca/dn
        !
        ! Incident-field terms are subtracted because the boundary condition is
        ! written for the total field.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: node_id
        integer, intent(in) :: particle_id
        integer, intent(in) :: side_id
        complex(dp), intent(out) :: known_value
        complex(dp), intent(out) :: coefficient

        complex(dp) :: robin_a
        complex(dp) :: robin_b
        complex(dp) :: robin_rhs

        coefficient = complex_zero
        known_value = complex_zero

        select case (case_data%boundary_condition(particle_id)%kind)
        case (bc_dirichlet_external, bc_dirichlet_internal)
            known_value = case_data%prescribed_boundary_data(node_id) - incident_phi(case_data, node_id, side_id)
        case (bc_robin_external, bc_robin_internal)
            robin_a = case_data%boundary_condition(particle_id)%robin_a
            robin_b = case_data%boundary_condition(particle_id)%robin_b
            robin_rhs = case_data%boundary_condition(particle_id)%robin_rhs
            coefficient = -robin_b / robin_a
            known_value = (robin_rhs - robin_b * incident_dphi_dn(case_data, node_id, side_id)) / robin_a - &
                          incident_phi(case_data, node_id, side_id)
        end select
    end subroutine dphi_dn_unknown_relation

    pure complex(dp) function incident_phi(case_data, node_id, side_id)
        ! Return incident phi on the requested side of a node.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: node_id
        integer, intent(in) :: side_id

        if (side_id == side_external) then
            incident_phi = case_data%incident%external_phi(node_id)
        else
            incident_phi = case_data%incident%internal_phi(node_id)
        end if
    end function incident_phi

    pure complex(dp) function incident_dphi_dn(case_data, node_id, side_id)
        ! Return incident normal derivative on the requested side of a node.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: node_id
        integer, intent(in) :: side_id

        if (side_id == side_external) then
            incident_dphi_dn = case_data%incident%external_dphi_dn(node_id)
        else
            incident_dphi_dn = case_data%incident%internal_dphi_dn(node_id)
        end if
    end function incident_dphi_dn

    subroutine set_phi_solution(case_data, node_id, side_id, value)
        ! Store solved scattered phi on the requested side.
        type(au_case_type), intent(inout) :: case_data
        integer, intent(in) :: node_id
        integer, intent(in) :: side_id
        complex(dp), intent(in) :: value

        if (side_id == side_external) then
            case_data%solution%external_phi(node_id) = value
        else
            case_data%solution%internal_phi(node_id) = value
        end if
    end subroutine set_phi_solution

    subroutine set_dphi_dn_solution(case_data, node_id, side_id, value)
        ! Store solved scattered dphi/dn on the requested side.
        type(au_case_type), intent(inout) :: case_data
        integer, intent(in) :: node_id
        integer, intent(in) :: side_id
        complex(dp), intent(in) :: value

        if (side_id == side_external) then
            case_data%solution%external_dphi_dn(node_id) = value
        else
            case_data%solution%internal_dphi_dn(node_id) = value
        end if
    end subroutine set_dphi_dn_solution

end module AU_Solver
