module AU_HelmholtzKernels
    ! This module evaluates acoustic Helmholtz Green-function integrals over
    ! one triangular surface element.
    !
    ! Boundary element methods use Green functions to express the acoustic field
    ! at a collocation point in terms of integrals over the boundary surface.
    ! In this file:
    !   - the collocation point is a mesh node where the equation is enforced
    !   - the source point is a quadrature point on an integration element
    !   - n  means the source-point normal
    !   - n0 means the collocation-node normal
    !
    ! The Green function convention used here omits the common 1/(4*pi) factor.
    ! That convention is inherited by the matrix assembly routines, so it must
    ! be used consistently throughout the solver.
    use Pre_Constants, only: dp, complex_zero, imaginary_unit
    use Geom_Types, only: geometry_type, linear_triangle_node_count, quadratic_triangle_node_count
    implicit none

    ! Only expose the result types and the main integration routine.
    private

    public :: burton_miller_element_integral_type
    public :: helmholtz_element_integral_type
    public :: clear_helmholtz_element_integral
    public :: integrate_helmholtz_element

    ! Absolute underflow guard used before direct 1/r evaluation.  The
    ! non-singular formulation removes the analytical same-surface singularity,
    ! but an exactly coincident numerical point would still divide by zero.
    ! This constant is not a near-singular quadrature accuracy criterion.
    real(dp), parameter :: tiny_distance = 100.0_dp * tiny(1.0_dp)

    type :: burton_miller_element_integral_type
        ! Additional element integrals used by Burton-Miller coupling.
        !
        ! Burton-Miller combines the ordinary boundary integral equation with
        ! its normal derivative.  This helps remove non-unique exterior
        ! Helmholtz solutions at irregular frequencies.
        !
        ! Each allocatable array has one entry per local element node.  The
        ! scalar regularization values are same-particle non-singular correction
        ! terms used when the collocation node and source element belong to the
        ! same particle.

        ! Integral of N_a * dG_k/dn0 over the element (operator K).
        complex(dp), allocatable :: dg_dn0(:)

        ! Integral of N_a * d2G_k/(dn dn0) over the element (operator D).  For same-particle
        ! integrations the Laplace singular part d2G_0/(dn dn0) is subtracted.
        complex(dp), allocatable :: d2g_dndn0(:)

        ! Same-particle NSBEM correction paired with dg_dn0.
        complex(dp) :: dg0_dn_regularization = complex_zero

        ! Same-particle NSBEM correction paired with d2g_dndn0.
        complex(dp) :: h_regularization = complex_zero

        ! Integral of -N_a * dG_0/dn0 over the element (operator L).
        complex(dp), allocatable :: minus_dg0_dn0(:)

        ! Integral of N_a * G_0 over the element (operator G0).
        complex(dp), allocatable :: g0(:)

        ! Integral of N_a * dG_0/dn over the element (operator H0).
        complex(dp), allocatable :: h0(:)

        ! Same-particle NSBEM correction paired with minus_dg0_dn0.
        complex(dp) :: minus_dg0_dn_regularization = complex_zero

        ! Same-particle NSBEM correction paired with g0.
        complex(dp) :: g0_regularization = complex_zero

        ! Same-particle NSBEM correction paired with h0.
        complex(dp) :: h0_regularization = complex_zero
    end type burton_miller_element_integral_type

    type :: helmholtz_element_integral_type
        ! Result of integrating one source element for one collocation node.
        !
        ! N_a is the shape function associated with local element node a.
        ! The arrays g(:) and h(:) therefore store the element contribution to
        ! each local nodal unknown.

        ! Integral of N_a * G_k over the element. This is the BEM single-layer
        ! block contribution for each local element node a.
        complex(dp), allocatable :: g(:)

        ! Integral of N_a * dG_k/dn over the element. This is the BEM
        ! double-layer block contribution for each local element node a.
        complex(dp), allocatable :: h(:)

        ! Non-singular BEM correction that is added to the collocation diagonal
        ! when the collocation node and integration element belong to the same
        ! particle.
        complex(dp) :: g_regularization = complex_zero
        complex(dp) :: h_regularization = complex_zero

        ! Burton-Miller auxiliary integrals. These arrays are allocated only
        ! when integrate_helmholtz_element is called with use_burton_miller=.true.
        type(burton_miller_element_integral_type) :: bm
    end type helmholtz_element_integral_type

    type :: green_values_type
        ! Green-function values at one source/collocation point pair.
        !
        ! G_k is the Helmholtz kernel:
        !   G_k = exp(i*k*r) / r
        !
        ! G_0 is the Laplace/static kernel:
        !   G_0 = 1 / r
        !
        ! The derivatives are directional derivatives with respect to either the
        ! source normal n or the collocation normal n0.
        complex(dp) :: gk = complex_zero
        complex(dp) :: dgk_dn = complex_zero
        complex(dp) :: dgk_dn0 = complex_zero
        complex(dp) :: d2gk_dndn0 = complex_zero

        complex(dp) :: g0 = complex_zero
        complex(dp) :: dg0_dn = complex_zero
        complex(dp) :: dg0_dn0 = complex_zero
        complex(dp) :: d2g0_dndn0 = complex_zero

        real(dp) :: n0_dot_r = 0.0_dp
        real(dp) :: n_dot_n0 = 0.0_dp
    end type green_values_type

contains

    subroutine integrate_helmholtz_element(geometry, collocation_node_id, element_id, &
                                          wavenumber, use_burton_miller, result, &
                                          status, message, same_particle)
        ! Integrate one element contribution for one collocation node.
        !
        ! This is called many times during matrix assembly.  For each
        ! collocation node and each source element, it loops over cached
        ! quadrature points and accumulates the Green-function integrals.
        !
        ! The optional same_particle flag lets higher-level code override the
        ! default test.  Internal ordinary multi-surface assembly can therefore
        ! state explicitly whether a regularisation term belongs here.  This
        ! hook does not expand the public v1.0 single-particle scope.
        type(geometry_type), intent(in) :: geometry
        integer, intent(in) :: collocation_node_id
        integer, intent(in) :: element_id
        complex(dp), intent(in) :: wavenumber
        logical, intent(in) :: use_burton_miller
        type(helmholtz_element_integral_type), intent(inout) :: result
        integer, intent(out), optional :: status
        character(len=*), intent(out), optional :: message
        logical, intent(in), optional :: same_particle

        integer :: error_code
        character(len=256) :: error_message
        integer :: nodes_per_element
        integer :: q
        integer :: point_id
        logical :: local_same_particle
        real(dp) :: collocation_point(3)
        real(dp) :: collocation_normal(3)
        real(dp) :: source_point(3)
        real(dp) :: source_normal(3)
        real(dp) :: r_vector(3)
        real(dp) :: shape_weight(quadratic_triangle_node_count)
        real(dp) :: integration_weight
        complex(dp) :: regularized_g0
        type(green_values_type) :: values

        error_code = 0
        error_message = "Helmholtz element integral evaluated."

        ! Start from an empty result object.  This avoids accidentally reusing
        ! arrays or values from a previous element integration.
        call clear_helmholtz_element_integral(result)

        ! Check that geometry, normals, and quadrature cache are available.
        call validate_inputs(geometry, collocation_node_id, element_id, error_code, error_message)
        if (error_code /= 0) then
            call finish(error_code, error_message)
            return
        end if

        nodes_per_element = geometry%mesh%nodes_per_element

        ! Allocate result arrays with one entry per local element node.
        call allocate_result(result, nodes_per_element, use_burton_miller, error_code, error_message)
        if (error_code /= 0) then
            call finish(error_code, error_message)
            return
        end if

        ! Decide whether this source element is on the same particle as the
        ! collocation node.  Same-particle terms need NSBEM regularisation.
        if (present(same_particle)) then
            local_same_particle = same_particle
        else
            local_same_particle = &
                geometry%mesh%node_particle_id(collocation_node_id) == &
                geometry%mesh%element_particle_id(element_id)
        end if

        ! Collocation point x0 and collocation normal n0.
        collocation_point = geometry%mesh%xyz(:, collocation_node_id)
        collocation_normal = geometry%differential%normal(:, collocation_node_id)

        ! Loop over quadrature points on this source element.
        do q = 1, geometry%quadrature%points_per_element
            point_id = geometry%quadrature%points_per_element * (element_id - 1) + q

            ! Source point x and source normal n.
            source_point = geometry%quadrature%point_xyz(:, point_id)
            source_normal = geometry%quadrature%point_normal(:, point_id)

            ! shape_weight already includes the quadrature weight and surface
            ! Jacobian.  integration_weight is the scalar weight without shape
            ! functions, used by non-singular correction terms.
            shape_weight = 0.0_dp
            shape_weight(1:nodes_per_element) = geometry%quadrature%shape_weight(:, point_id)
            integration_weight = geometry%quadrature%integration_weight(point_id)

            ! r_vector points from the collocation point x0 to the source point x.
            r_vector = source_point - collocation_point

            ! Evaluate G_k, G_0, and their normal derivatives for this point pair.
            call evaluate_green_values(r_vector, source_normal, collocation_normal, &
                                       wavenumber, values, error_code, error_message)
            if (error_code /= 0) then
                call clear_helmholtz_element_integral(result)
                call finish(error_code, error_message)
                return
            end if

            ! Ordinary Helmholtz BEM element contributions.
            result%g = result%g + shape_weight(1:nodes_per_element) * values%gk
            result%h = result%h + shape_weight(1:nodes_per_element) * values%dgk_dn

            if (local_same_particle) then
                ! Non-singular correction based on the Laplace kernel G_0.
                ! This is only applied when the source element belongs to the
                ! same particle as the collocation point.
                regularized_g0 = values%n0_dot_r * values%dg0_dn - &
                                 values%n_dot_n0 * values%g0
                result%g_regularization = result%g_regularization + &
                    integration_weight * regularized_g0
                result%h_regularization = result%h_regularization - &
                    integration_weight * values%dg0_dn
            end if

            if (use_burton_miller) then
                ! Burton-Miller normal-derivative equation contributions.
                result%bm%dg_dn0 = result%bm%dg_dn0 + &
                    shape_weight(1:nodes_per_element) * values%dgk_dn0
                result%bm%d2g_dndn0 = result%bm%d2g_dndn0 + &
                    shape_weight(1:nodes_per_element) * values%d2gk_dndn0

                if (local_same_particle) then
                    ! Same-particle Burton-Miller regularisation terms.  The
                    ! Laplace singular part is subtracted from the hypersingular
                    ! kernel before assembly.
                    regularized_g0 = values%n0_dot_r * values%dg0_dn - &
                                     values%n_dot_n0 * values%g0

                    result%bm%d2g_dndn0 = result%bm%d2g_dndn0 - &
                        shape_weight(1:nodes_per_element) * values%d2g0_dndn0
                    result%bm%dg0_dn_regularization = &
                        result%bm%dg0_dn_regularization + integration_weight * values%dg0_dn
                    result%bm%h_regularization = result%bm%h_regularization + &
                        integration_weight * regularized_g0 * 0.5_dp * wavenumber * wavenumber

                    result%bm%minus_dg0_dn0 = result%bm%minus_dg0_dn0 - &
                        shape_weight(1:nodes_per_element) * values%dg0_dn0
                    result%bm%g0 = result%bm%g0 + shape_weight(1:nodes_per_element) * values%g0
                    result%bm%h0 = result%bm%h0 + shape_weight(1:nodes_per_element) * values%dg0_dn

                    result%bm%minus_dg0_dn_regularization = &
                        result%bm%minus_dg0_dn_regularization - integration_weight * values%dg0_dn
                    result%bm%g0_regularization = result%bm%g0_regularization + &
                        integration_weight * regularized_g0
                    result%bm%h0_regularization = result%bm%h0_regularization - &
                        integration_weight * values%dg0_dn
                end if
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

    end subroutine integrate_helmholtz_element

    subroutine clear_helmholtz_element_integral(result)
        ! Release all arrays inside an element-integral result and reset scalar
        ! regularisation terms to zero.
        type(helmholtz_element_integral_type), intent(inout) :: result

        ! Ordinary BEM arrays.
        if (allocated(result%g)) deallocate(result%g)
        if (allocated(result%h)) deallocate(result%h)

        ! Burton-Miller arrays.
        if (allocated(result%bm%dg_dn0)) deallocate(result%bm%dg_dn0)
        if (allocated(result%bm%d2g_dndn0)) deallocate(result%bm%d2g_dndn0)
        if (allocated(result%bm%minus_dg0_dn0)) deallocate(result%bm%minus_dg0_dn0)
        if (allocated(result%bm%g0)) deallocate(result%bm%g0)
        if (allocated(result%bm%h0)) deallocate(result%bm%h0)

        result%g_regularization = complex_zero
        result%h_regularization = complex_zero
        result%bm%dg0_dn_regularization = complex_zero
        result%bm%h_regularization = complex_zero
        result%bm%minus_dg0_dn_regularization = complex_zero
        result%bm%g0_regularization = complex_zero
        result%bm%h0_regularization = complex_zero
    end subroutine clear_helmholtz_element_integral

    subroutine validate_inputs(geometry, collocation_node_id, element_id, status, message)
        ! Check that all geometry data needed by the kernel integrator exists.
        type(geometry_type), intent(in) :: geometry
        integer, intent(in) :: collocation_node_id
        integer, intent(in) :: element_id
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        status = 0
        message = "Helmholtz kernel inputs are valid."

        if (geometry%mesh%nodes_per_element /= linear_triangle_node_count .and. &
            geometry%mesh%nodes_per_element /= quadratic_triangle_node_count) then
            call fail("Helmholtz kernels support only 3-node or 6-node triangles.")
            return
        end if

        if (collocation_node_id < 1 .or. collocation_node_id > geometry%mesh%node_count) then
            call fail("Collocation node id is outside the mesh.")
            return
        end if

        if (element_id < 1 .or. element_id > geometry%mesh%element_count) then
            call fail("Element id is outside the mesh.")
            return
        end if

        if (.not. allocated(geometry%mesh%xyz)) then
            call fail("Mesh coordinates are not allocated.")
            return
        end if

        if (.not. allocated(geometry%mesh%node_particle_id)) then
            call fail("Node particle ids are not allocated.")
            return
        end if

        if (.not. allocated(geometry%mesh%element_particle_id)) then
            call fail("Element particle ids are not allocated.")
            return
        end if

        if (.not. allocated(geometry%differential%normal)) then
            call fail("Surface normals have not been prepared.")
            return
        end if

        if (.not. allocated(geometry%quadrature%point_xyz) .or. &
            .not. allocated(geometry%quadrature%point_normal) .or. &
            .not. allocated(geometry%quadrature%shape_weight) .or. &
            .not. allocated(geometry%quadrature%integration_weight)) then
            call fail("Quadrature cache has not been prepared.")
            return
        end if

        if (geometry%quadrature%points_per_element <= 0) then
            call fail("Quadrature cache has no points per element.")
            return
        end if

    contains

        subroutine fail(text)
            ! Store one validation failure message.
            character(len=*), intent(in) :: text

            status = 1
            message = text
        end subroutine fail

    end subroutine validate_inputs

    subroutine allocate_result(result, nodes_per_element, use_burton_miller, status, message)
        ! Allocate and initialise result arrays for one element integration.
        type(helmholtz_element_integral_type), intent(inout) :: result
        integer, intent(in) :: nodes_per_element
        logical, intent(in) :: use_burton_miller
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        status = 0
        message = "Helmholtz element result allocated."

        ! Ordinary single-layer and double-layer element contributions.
        allocate(result%g(nodes_per_element), result%h(nodes_per_element), stat=status)
        if (status /= 0) then
            message = "Unable to allocate Helmholtz element result arrays."
            return
        end if

        result%g = complex_zero
        result%h = complex_zero

        if (.not. use_burton_miller) return

        ! Extra arrays needed only when Burton-Miller coupling is enabled.
        allocate(result%bm%dg_dn0(nodes_per_element), &
                 result%bm%d2g_dndn0(nodes_per_element), &
                 result%bm%minus_dg0_dn0(nodes_per_element), &
                 result%bm%g0(nodes_per_element), &
                 result%bm%h0(nodes_per_element), &
                 stat=status)
        if (status /= 0) then
            message = "Unable to allocate Burton-Miller element result arrays."
            call clear_helmholtz_element_integral(result)
            return
        end if

        result%bm%dg_dn0 = complex_zero
        result%bm%d2g_dndn0 = complex_zero
        result%bm%minus_dg0_dn0 = complex_zero
        result%bm%g0 = complex_zero
        result%bm%h0 = complex_zero
    end subroutine allocate_result

    subroutine evaluate_green_values(r_vector, source_normal, collocation_normal, &
                                     wavenumber, values, status, message)
        ! Evaluate Helmholtz and Laplace Green-function quantities for one
        ! source/collocation pair.
        !
        ! r_vector = x - x0
        !   x  = source quadrature point
        !   x0 = collocation node
        !
        ! source_normal is n, and collocation_normal is n0.
        real(dp), intent(in) :: r_vector(3)
        real(dp), intent(in) :: source_normal(3)
        real(dp), intent(in) :: collocation_normal(3)
        complex(dp), intent(in) :: wavenumber
        type(green_values_type), intent(out) :: values
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        real(dp) :: radius
        real(dp) :: inverse_radius
        real(dp) :: inverse_radius_2
        real(dp) :: inverse_radius_3
        real(dp) :: source_dot_r
        real(dp) :: source_collocation_dot
        real(dp) :: source_collocation_r_product
        complex(dp) :: ik
        complex(dp) :: ikr
        complex(dp) :: exp_ikr
        complex(dp) :: hessian_term_1
        complex(dp) :: hessian_term_2

        status = 0
        message = "Green-function values evaluated."
        values = green_values_type()

        ! Distance between source point and collocation point.
        radius = sqrt(dot_product(r_vector, r_vector))
        if (radius <= tiny_distance) then
            status = 2
            message = "Green-function source and collocation points are too close."
            return
        end if

        ! Reuse powers of 1/r instead of recomputing them many times.
        inverse_radius = 1.0_dp / radius
        inverse_radius_2 = inverse_radius * inverse_radius
        inverse_radius_3 = inverse_radius_2 * inverse_radius

        ! Geometric dot products used by the normal derivatives.
        source_dot_r = dot_product(source_normal, r_vector)
        values%n0_dot_r = dot_product(collocation_normal, r_vector)
        source_collocation_dot = dot_product(source_normal, collocation_normal)
        values%n_dot_n0 = source_collocation_dot
        source_collocation_r_product = source_dot_r * values%n0_dot_r

        ik = imaginary_unit * wavenumber
        ikr = ik * radius
        exp_ikr = exp(ikr)

        ! Helmholtz Green function and first normal derivatives.
        !
        ! G_k = exp(i*k*r) / r
        ! dgk_dn differentiates with respect to the source normal n.
        ! dgk_dn0 differentiates with respect to the collocation normal n0.
        values%gk = exp_ikr * inverse_radius
        values%dgk_dn = -source_dot_r * (inverse_radius - ik) * exp_ikr * inverse_radius_2
        values%dgk_dn0 = values%n0_dot_r * (inverse_radius - ik) * exp_ikr * inverse_radius_2

        ! Mixed second normal derivative d2G_k/(dn dn0).  This is the
        ! hypersingular quantity required by Burton-Miller coupling.
        hessian_term_1 = ikr - 1.0_dp
        hessian_term_2 = ikr * ikr - 3.0_dp * ikr + 3.0_dp
        values%d2gk_dndn0 = exp_ikr * inverse_radius_3 * &
            (-hessian_term_1 * source_collocation_dot - &
             hessian_term_2 * source_collocation_r_product * inverse_radius_2)

        ! Laplace/static Green function and derivatives.  These are used for
        ! same-particle non-singular regularisation.
        values%g0 = cmplx(inverse_radius, 0.0_dp, kind=dp)
        values%dg0_dn = cmplx(-source_dot_r * inverse_radius_3, 0.0_dp, kind=dp)
        values%dg0_dn0 = cmplx(values%n0_dot_r * inverse_radius_3, 0.0_dp, kind=dp)
        values%d2g0_dndn0 = cmplx((source_collocation_dot - &
            3.0_dp * source_collocation_r_product * inverse_radius_2) * inverse_radius_3, &
            0.0_dp, kind=dp)
    end subroutine evaluate_green_values

end module AU_HelmholtzKernels
