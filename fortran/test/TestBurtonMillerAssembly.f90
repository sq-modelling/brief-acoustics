program TestBurtonMillerAssembly
    ! Check the algebra of the 2N x 2N Burton-Miller system independently of
    ! mesh integration.  We provide small artificial operator matrices, ask
    ! AU_Solver to assemble the final system, and compare every entry with
    ! the block equations written below in build_expected_system.
    use Pre_Constants, only: dp, complex_zero
    use Geom_Types, only: geometry_type
    use AU_Types, only: au_case_type, allocate_au_case, &
                        bc_dirichlet_external, bc_neumann_external, bc_robin_external
    use AU_LayerTopology, only: burton_miller_coupling_length, &
                                burton_miller_coupling_weight
    use AU_SurfaceOperators, only: au_surface_operator_type
    use AU_Solver, only: assemble_au_linear_system
    implicit none

    integer, parameter :: node_count = 2
    real(dp), parameter :: comparison_tolerance = 1.0e-12_dp

    type(geometry_type) :: geometry
    type(au_case_type) :: case_data
    type(au_surface_operator_type) :: operators
    integer :: node_unknown(node_count)
    integer :: status
    character(len=256) :: message

    call prepare_minimal_case()
    call check_coupling_weight()
    call prepare_test_operators()

    call check_dirichlet_case()
    call check_neumann_case()
    call check_robin_case()

    write(*, '(A)') "Burton-Miller matrix assembly checks passed."

contains

    subroutine prepare_minimal_case()
        ! Only node-to-particle ownership is needed by the matrix assembler.
        ! No physical mesh or quadrature is involved in this algebra test.
        geometry%particle_count = 1
        geometry%mesh%node_count = node_count
        allocate(geometry%mesh%node_particle_id(node_count))
        geometry%mesh%node_particle_id = 1

        call allocate_au_case(case_data, geometry, status, message)
        call require_ok(status, message)

        case_data%use_burton_miller = .true.
        case_data%layer(1)%parent_particle_id = 0
        case_data%layer(1)%exterior_domain_normal_sign = -1
        case_data%layer(1)%characteristic_length = 1.0_dp

        ! Non-zero incident fields make the right-hand-side checks meaningful.
        case_data%incident%external_phi(1) = cmplx(0.20_dp, -0.10_dp, kind=dp)
        case_data%incident%external_phi(2) = cmplx(-0.15_dp, 0.05_dp, kind=dp)
        case_data%incident%external_dphi_dn(1) = cmplx(0.08_dp, 0.03_dp, kind=dp)
        case_data%incident%external_dphi_dn(2) = cmplx(-0.04_dp, 0.06_dp, kind=dp)

        node_unknown = [1, 2]
    end subroutine prepare_minimal_case

    subroutine check_coupling_weight()
        ! beta=min(a,1/k) has two regimes.  With a=1, k=0.25 must use the
        ! low-frequency body scale beta=1, while k=2 must use the
        ! high-frequency wave scale beta=0.5.  Reversing the normal sign must
        ! reverse the imaginary coupling weight without changing beta.
        complex(dp) :: actual_weight
        real(dp) :: actual_length

        case_data%exterior_medium%wavenumber = cmplx(0.25_dp, 0.0_dp, kind=dp)
        actual_length = burton_miller_coupling_length(case_data, 1)
        if (abs(actual_length - 1.0_dp) > comparison_tolerance) then
            write(*, '(A)') "Burton-Miller low-frequency coupling length failed."
            error stop 1
        end if

        actual_weight = burton_miller_coupling_weight(case_data, 1)
        if (abs(actual_weight - cmplx(0.0_dp, -1.0_dp, kind=dp)) > comparison_tolerance) then
            write(*, '(A)') "Burton-Miller low-frequency coupling weight failed."
            error stop 1
        end if

        case_data%exterior_medium%wavenumber = cmplx(2.0_dp, 0.0_dp, kind=dp)
        actual_length = burton_miller_coupling_length(case_data, 1)
        if (abs(actual_length - 0.5_dp) > comparison_tolerance) then
            write(*, '(A)') "Burton-Miller high-frequency coupling length failed."
            error stop 1
        end if

        case_data%layer(1)%exterior_domain_normal_sign = -1
        actual_weight = burton_miller_coupling_weight(case_data, 1)
        if (abs(actual_weight - cmplx(0.0_dp, -0.5_dp, kind=dp)) > comparison_tolerance) then
            write(*, '(A)') "Burton-Miller coupling weight failed for exterior_domain_normal_sign=-1."
            error stop 1
        end if

        case_data%layer(1)%exterior_domain_normal_sign = 1
        actual_weight = burton_miller_coupling_weight(case_data, 1)
        if (abs(actual_weight - cmplx(0.0_dp, 0.5_dp, kind=dp)) > comparison_tolerance) then
            write(*, '(A)') "Burton-Miller coupling weight failed for exterior_domain_normal_sign=+1."
            error stop 1
        end if

        ! Restore the release convention used by the matrix cases below.
        case_data%layer(1)%exterior_domain_normal_sign = -1
        write(*, '(A)') "Burton-Miller coupling weight: PASS"
    end subroutine check_coupling_weight

    subroutine prepare_test_operators()
        ! The values are arbitrary but deliberately non-symmetric and complex.
        ! Different row weights detect the common mistake W*A -> A*W.
        integer :: column_id
        integer :: row_id

        allocate(operators%external_g(node_count, node_count), &
                 operators%external_h(node_count, node_count), &
                 operators%internal_g(node_count, node_count), &
                 operators%internal_h(node_count, node_count), &
                 operators%bm_coupling_weight(node_count), &
                 operators%bm_dg_dn0(node_count, node_count), &
                 operators%bm_d2g_dndn0(node_count, node_count), &
                 operators%bm_minus_dg0_dn0(node_count, node_count), &
                 operators%bm_g0(node_count, node_count), &
                 operators%bm_h0(node_count, node_count))

        operators%internal_g = complex_zero
        operators%internal_h = complex_zero
        operators%bm_coupling_weight(1) = cmplx(0.0_dp, -0.70_dp, kind=dp)
        operators%bm_coupling_weight(2) = cmplx(0.0_dp, -1.10_dp, kind=dp)

        do column_id = 1, node_count
            do row_id = 1, node_count
                operators%external_h(row_id, column_id) = &
                    cmplx(1.1_dp * row_id + 0.4_dp * column_id, &
                          0.07_dp * (row_id - column_id), kind=dp)
                operators%external_g(row_id, column_id) = &
                    cmplx(0.6_dp * row_id - 0.2_dp * column_id, &
                          0.05_dp * (row_id + column_id), kind=dp)
                operators%bm_d2g_dndn0(row_id, column_id) = &
                    cmplx(0.3_dp * row_id + 0.1_dp * column_id, &
                          -0.04_dp * column_id, kind=dp)
                operators%bm_dg_dn0(row_id, column_id) = &
                    cmplx(-0.2_dp * row_id + 0.15_dp * column_id, &
                          0.03_dp * row_id, kind=dp)
                operators%bm_minus_dg0_dn0(row_id, column_id) = &
                    cmplx(0.25_dp * row_id + 0.05_dp * column_id, &
                          -0.02_dp * (row_id + column_id), kind=dp)
                operators%bm_g0(row_id, column_id) = &
                    cmplx(0.9_dp * row_id + 0.2_dp * column_id, &
                          0.01_dp * (row_id - column_id), kind=dp)
                operators%bm_h0(row_id, column_id) = &
                    cmplx(-0.35_dp * row_id + 0.12_dp * column_id, &
                          0.06_dp * column_id, kind=dp)
            end do
        end do
    end subroutine prepare_test_operators

    subroutine check_dirichlet_case()
        complex(dp) :: phi_bar(node_count)
        complex(dp) :: phi_coefficient(node_count)
        complex(dp) :: q_bar(node_count)
        complex(dp) :: q_coefficient(node_count)

        case_data%boundary_condition(1)%kind = bc_dirichlet_external
        case_data%prescribed_boundary_data(1) = cmplx(0.30_dp, -0.20_dp, kind=dp)
        case_data%prescribed_boundary_data(2) = cmplx(-0.10_dp, 0.25_dp, kind=dp)

        phi_bar = case_data%prescribed_boundary_data - case_data%incident%external_phi
        phi_coefficient = complex_zero
        q_bar = complex_zero
        q_coefficient = cmplx(1.0_dp, 0.0_dp, kind=dp)

        call check_case("Dirichlet", phi_bar, phi_coefficient, q_bar, q_coefficient)
    end subroutine check_dirichlet_case

    subroutine check_neumann_case()
        complex(dp) :: phi_bar(node_count)
        complex(dp) :: phi_coefficient(node_count)
        complex(dp) :: q_bar(node_count)
        complex(dp) :: q_coefficient(node_count)

        case_data%boundary_condition(1)%kind = bc_neumann_external
        case_data%prescribed_boundary_data(1) = cmplx(0.11_dp, 0.09_dp, kind=dp)
        case_data%prescribed_boundary_data(2) = cmplx(-0.07_dp, -0.02_dp, kind=dp)

        phi_bar = complex_zero
        phi_coefficient = cmplx(1.0_dp, 0.0_dp, kind=dp)
        q_bar = case_data%prescribed_boundary_data - case_data%incident%external_dphi_dn
        q_coefficient = complex_zero

        call check_case("Neumann", phi_bar, phi_coefficient, q_bar, q_coefficient)
    end subroutine check_neumann_case

    subroutine check_robin_case()
        complex(dp) :: phi_bar(node_count)
        complex(dp) :: phi_coefficient(node_count)
        complex(dp) :: q_bar(node_count)
        complex(dp) :: q_coefficient(node_count)
        complex(dp) :: robin_a
        complex(dp) :: robin_b
        complex(dp) :: robin_rhs

        robin_a = cmplx(1.2_dp, -0.1_dp, kind=dp)
        robin_b = cmplx(0.7_dp, 0.2_dp, kind=dp)
        robin_rhs = cmplx(0.15_dp, -0.05_dp, kind=dp)

        case_data%boundary_condition(1)%kind = bc_robin_external
        case_data%boundary_condition(1)%robin_a = robin_a
        case_data%boundary_condition(1)%robin_b = robin_b
        case_data%boundary_condition(1)%robin_rhs = robin_rhs

        phi_bar = complex_zero
        phi_coefficient = cmplx(1.0_dp, 0.0_dp, kind=dp)
        q_bar = (robin_rhs - robin_a * case_data%incident%external_phi) / robin_b - &
                case_data%incident%external_dphi_dn
        q_coefficient = -robin_a / robin_b

        call check_case("Robin", phi_bar, phi_coefficient, q_bar, q_coefficient)
    end subroutine check_robin_case

    subroutine check_case(label, phi_bar, phi_coefficient, q_bar, q_coefficient)
        character(len=*), intent(in) :: label
        complex(dp), intent(in) :: phi_bar(node_count)
        complex(dp), intent(in) :: phi_coefficient(node_count)
        complex(dp), intent(in) :: q_bar(node_count)
        complex(dp), intent(in) :: q_coefficient(node_count)

        complex(dp), allocatable :: actual_matrix(:, :)
        complex(dp), allocatable :: actual_rhs(:)
        complex(dp) :: expected_matrix(2 * node_count, 2 * node_count)
        complex(dp) :: expected_rhs(2 * node_count)
        real(dp) :: matrix_error
        real(dp) :: rhs_error

        call assemble_au_linear_system(geometry, case_data, operators, node_unknown, &
                                       actual_matrix, actual_rhs, status, message)
        call require_ok(status, message)

        call build_expected_system(phi_bar, phi_coefficient, q_bar, q_coefficient, &
                                   expected_matrix, expected_rhs)

        if (size(actual_matrix, 1) /= 2 * node_count .or. &
            size(actual_matrix, 2) /= 2 * node_count .or. &
            size(actual_rhs) /= 2 * node_count) then
            write(*, '(A)') trim(label) // " check produced the wrong 2N system size."
            error stop 1
        end if

        matrix_error = maxval(abs(actual_matrix - expected_matrix))
        rhs_error = maxval(abs(actual_rhs - expected_rhs))
        if (matrix_error > comparison_tolerance .or. rhs_error > comparison_tolerance) then
            write(*, '(A,ES12.4)') trim(label) // " matrix maximum error = ", matrix_error
            write(*, '(A,ES12.4)') trim(label) // " rhs maximum error = ", rhs_error
            error stop 1
        end if

        write(*, '(A)') trim(label) // " block assembly: PASS"
    end subroutine check_case

    subroutine build_expected_system(phi_bar, phi_coefficient, q_bar, q_coefficient, &
                                     expected_matrix, expected_rhs)
        ! Construct the augmented block system directly:
        !
        !   [ Hc*P-Gc*Q   -W*L ] [u] = [-Hc*phi_bar+Gc*q_bar]
        !   [ H0*P         -G0  ] [s]   [-H0*phi_bar]
        !
        ! Boundary elimination gives phi=phi_bar+P*u and q=q_bar+Q*u.
        ! W contains the row coupling weights; Hc=H+W*D and Gc=G+W*K.
        ! s is the normal derivative of the auxiliary harmonic field.
        complex(dp), intent(in) :: phi_bar(node_count)
        complex(dp), intent(in) :: phi_coefficient(node_count)
        complex(dp), intent(in) :: q_bar(node_count)
        complex(dp), intent(in) :: q_coefficient(node_count)
        complex(dp), intent(out) :: expected_matrix(2 * node_count, 2 * node_count)
        complex(dp), intent(out) :: expected_rhs(2 * node_count)

        complex(dp) :: combined_g(node_count, node_count)
        complex(dp) :: combined_h(node_count, node_count)
        complex(dp) :: weighted_l(node_count, node_count)
        integer :: column_id
        integer :: row_id

        do column_id = 1, node_count
            do row_id = 1, node_count
                combined_h(row_id, column_id) = operators%external_h(row_id, column_id) + &
                    operators%bm_coupling_weight(row_id) * operators%bm_d2g_dndn0(row_id, column_id)
                combined_g(row_id, column_id) = operators%external_g(row_id, column_id) + &
                    operators%bm_coupling_weight(row_id) * operators%bm_dg_dn0(row_id, column_id)
                weighted_l(row_id, column_id) = operators%bm_coupling_weight(row_id) * &
                    operators%bm_minus_dg0_dn0(row_id, column_id)
            end do
        end do

        expected_matrix = complex_zero
        do column_id = 1, node_count
            expected_matrix(1:node_count, column_id) = &
                combined_h(:, column_id) * phi_coefficient(column_id) - &
                combined_g(:, column_id) * q_coefficient(column_id)
            expected_matrix(1:node_count, node_count + column_id) = -weighted_l(:, column_id)
            expected_matrix(node_count + 1:2 * node_count, column_id) = &
                operators%bm_h0(:, column_id) * phi_coefficient(column_id)
            expected_matrix(node_count + 1:2 * node_count, node_count + column_id) = &
                -operators%bm_g0(:, column_id)
        end do

        expected_rhs(1:node_count) = -matmul(combined_h, phi_bar) + matmul(combined_g, q_bar)
        expected_rhs(node_count + 1:2 * node_count) = -matmul(operators%bm_h0, phi_bar)
    end subroutine build_expected_system

    subroutine require_ok(error_code, error_message)
        integer, intent(in) :: error_code
        character(len=*), intent(in) :: error_message

        if (error_code /= 0) then
            write(*, '(A)') trim(error_message)
            error stop 1
        end if
    end subroutine require_ok

end program TestBurtonMillerAssembly
