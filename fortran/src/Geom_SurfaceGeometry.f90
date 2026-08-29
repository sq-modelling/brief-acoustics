module Geom_SurfaceGeometry
    ! This module turns mesh connectivity into geometric information.
    !
    ! Previous geometry modules read the mesh and build topology, such as which
    ! nodes are neighbours.  This module computes quantities that depend on the
    ! actual 3D positions of the nodes:
    !   - element areas and element normals
    !   - node normals and tangent directions
    !   - quadrature point positions, normals, and integration weights
    !   - simple local derivative stencils on the surface
    !
    ! These quantities are prepared once and then reused by the acoustic
    ! boundary-integral assembly routines.
    use Pre_Constants, only: dp
    use Geom_Types, only: geometry_type, linear_triangle_node_count, quadratic_triangle_node_count
    use Num_Quadrature, only: triangle_quadrature_25
    implicit none

    private

    public :: prepare_surface_geometry

    ! Number of integration points used on each triangular element.
    ! The actual coordinates and weights are stored in Num_Quadrature.
    integer, parameter :: triangle_quadrature_points = 25

    ! Absolute underflow guard used before divisions.  This is not a
    ! scale-aware mesh-quality tolerance: the public Python mesh audit performs
    ! the practical relative degeneracy checks before Fortran is called.
    ! Direct Fortran callers must provide a non-degenerate, sensibly scaled mesh.
    real(dp), parameter :: tiny_geometry = 100.0_dp * tiny(1.0_dp)

contains

    subroutine prepare_surface_geometry(geometry, status, message)
        ! Main entry point for this module.
        !
        ! Expected input:
        !   geometry%mesh%xyz
        !   geometry%mesh%element_nodes
        !   geometry%mesh%subtriangle_nodes
        !   geometry%mesh%node_neighbors_2ring
        !
        ! Produced output:
        !   geometry%differential%element_area
        !   geometry%differential%element_normal
        !   geometry%differential%normal
        !   geometry%differential%tangent_1
        !   geometry%differential%tangent_2
        !   geometry%differential%d_dt1
        !   geometry%differential%d_dt2
        !   geometry%quadrature%...
        !
        ! The order matters: quadrature needs element geometry, and derivative
        ! stencils need node normals and tangent directions.
        type(geometry_type), intent(inout) :: geometry
        integer, intent(out), optional :: status
        character(len=*), intent(out), optional :: message

        integer :: error_code
        character(len=256) :: error_message

        error_code = 0
        error_message = "Surface geometry prepared successfully."

        ! Remove old surface geometry before rebuilding it.
        call clear_surface_geometry(geometry)

        ! Check that the mesh and topology arrays needed here exist.
        call validate_surface_inputs(geometry, error_code, error_message)
        if (error_code /= 0) then
            call finish(error_code, error_message)
            return
        end if

        ! Allocate all arrays that this module will fill.
        call allocate_surface_geometry(geometry, error_code, error_message)
        if (error_code /= 0) then
            call finish(error_code, error_message)
            return
        end if

        ! Integrate over each element to estimate its area and average normal.
        call build_element_geometry(geometry, error_code, error_message)
        if (error_code /= 0) then
            call finish(error_code, error_message)
            return
        end if

        ! Build one normal and two tangent vectors at each node.
        call build_node_frames(geometry, error_code, error_message)
        if (error_code /= 0) then
            call finish(error_code, error_message)
            return
        end if

        ! Precompute quadrature point data for later matrix assembly.
        call build_quadrature_cache(geometry, error_code, error_message)
        if (error_code /= 0) then
            call finish(error_code, error_message)
            return
        end if

        ! Build local finite-difference-like weights for surface derivatives.
        call build_derivative_stencils(geometry, error_code, error_message)
        call finish(error_code, error_message)

    contains

        subroutine finish(code, text)
            ! Return status information.  If preparation failed, remove any
            ! partly-built arrays so the geometry object is not half-valid.
            integer, intent(in) :: code
            character(len=*), intent(in) :: text

            if (code /= 0) call clear_surface_geometry(geometry)
            if (present(status)) status = code
            if (present(message)) message = text
        end subroutine finish

    end subroutine prepare_surface_geometry

    subroutine clear_surface_geometry(geometry)
        ! Deallocate geometry arrays built by this module.
        !
        ! This does not remove the raw mesh or topology arrays.  It only clears
        ! derived surface geometry and cached quadrature data.
        type(geometry_type), intent(inout) :: geometry

        if (allocated(geometry%differential%normal)) deallocate(geometry%differential%normal)
        if (allocated(geometry%differential%tangent_1)) deallocate(geometry%differential%tangent_1)
        if (allocated(geometry%differential%tangent_2)) deallocate(geometry%differential%tangent_2)
        if (allocated(geometry%differential%element_area)) deallocate(geometry%differential%element_area)
        if (allocated(geometry%differential%element_normal)) deallocate(geometry%differential%element_normal)
        if (allocated(geometry%differential%curvature_1)) deallocate(geometry%differential%curvature_1)
        if (allocated(geometry%differential%curvature_2)) deallocate(geometry%differential%curvature_2)
        if (allocated(geometry%differential%mean_curvature)) deallocate(geometry%differential%mean_curvature)
        if (allocated(geometry%differential%d_dt1)) deallocate(geometry%differential%d_dt1)
        if (allocated(geometry%differential%d_dt2)) deallocate(geometry%differential%d_dt2)

        if (allocated(geometry%quadrature%point_xyz)) deallocate(geometry%quadrature%point_xyz)
        if (allocated(geometry%quadrature%point_normal)) deallocate(geometry%quadrature%point_normal)
        if (allocated(geometry%quadrature%shape_weight)) deallocate(geometry%quadrature%shape_weight)
        if (allocated(geometry%quadrature%integration_weight)) deallocate(geometry%quadrature%integration_weight)

        geometry%quadrature%points_per_element = 0
        geometry%quadrature%point_count = 0
    end subroutine clear_surface_geometry

    subroutine validate_surface_inputs(geometry, status, message)
        ! Check that the previous pipeline stages have been completed.
        !
        ! Surface geometry needs both the raw mesh and mesh topology.  For
        ! example, node tangent frames need coordinates, while derivative
        ! stencils need two-ring neighbours.
        type(geometry_type), intent(in) :: geometry
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        status = 0
        message = "Surface geometry inputs are valid."

        if (geometry%mesh%node_count <= 0) then
            call fail("Mesh has no nodes.")
            return
        end if

        if (geometry%mesh%element_count <= 0) then
            call fail("Mesh has no elements.")
            return
        end if

        if (geometry%mesh%nodes_per_element /= linear_triangle_node_count .and. &
            geometry%mesh%nodes_per_element /= quadratic_triangle_node_count) then
            call fail("Surface geometry supports only 3-node or 6-node triangles.")
            return
        end if

        if (.not. allocated(geometry%mesh%xyz)) then
            call fail("Mesh coordinates are not allocated.")
            return
        end if

        if (.not. allocated(geometry%mesh%element_nodes)) then
            call fail("Mesh element connectivity is not allocated.")
            return
        end if

        if (.not. allocated(geometry%mesh%subtriangle_nodes)) then
            call fail("Mesh topology has not been prepared.")
            return
        end if

        if (.not. allocated(geometry%mesh%node_neighbors_2ring)) then
            call fail("Two-ring node topology has not been prepared.")
            return
        end if

    contains

        subroutine fail(text)
            ! Store one validation failure message.
            character(len=*), intent(in) :: text

            status = 1
            message = text
        end subroutine fail

    end subroutine validate_surface_inputs

    subroutine allocate_surface_geometry(geometry, status, message)
        ! Allocate arrays for surface geometry and quadrature cache.
        !
        ! Most arrays have either one value per node, one value per element, or
        ! one value per quadrature point.
        type(geometry_type), intent(inout) :: geometry
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: node_count
        integer :: element_count
        integer :: nodes_per_element
        integer :: point_count
        integer :: stencil_size

        status = 0
        message = "Surface geometry arrays allocated."

        node_count = geometry%mesh%node_count
        element_count = geometry%mesh%element_count
        nodes_per_element = geometry%mesh%nodes_per_element
        point_count = triangle_quadrature_points * element_count

        ! d_dt1 and d_dt2 store one host-node weight plus one weight for each
        ! two-ring neighbour.  Row 1 is the host node; rows 2: are neighbours.
        stencil_size = geometry%mesh%max_node_neighbors_2ring + 1

        allocate(geometry%differential%normal(3, node_count), &
                 geometry%differential%tangent_1(3, node_count), &
                 geometry%differential%tangent_2(3, node_count), &
                 geometry%differential%element_area(element_count), &
                 geometry%differential%element_normal(3, element_count), &
                 geometry%differential%curvature_1(node_count), &
                 geometry%differential%curvature_2(node_count), &
                 geometry%differential%mean_curvature(node_count), &
                 geometry%differential%d_dt1(stencil_size, node_count), &
                 geometry%differential%d_dt2(stencil_size, node_count), &
                 geometry%quadrature%point_xyz(3, point_count), &
                 geometry%quadrature%point_normal(3, point_count), &
                 geometry%quadrature%shape_weight(nodes_per_element, point_count), &
                 geometry%quadrature%integration_weight(point_count), &
                 stat=status)

        if (status /= 0) then
            message = "Unable to allocate surface geometry arrays."
            return
        end if

        geometry%differential%normal = 0.0_dp
        geometry%differential%tangent_1 = 0.0_dp
        geometry%differential%tangent_2 = 0.0_dp
        geometry%differential%element_area = 0.0_dp
        geometry%differential%element_normal = 0.0_dp

        ! The current release allocates curvature arrays for API consistency, but
        ! it does not yet estimate curvature.  These values intentionally remain
        ! zero unless a future curvature routine fills them.
        geometry%differential%curvature_1 = 0.0_dp
        geometry%differential%curvature_2 = 0.0_dp
        geometry%differential%mean_curvature = 0.0_dp
        geometry%differential%d_dt1 = 0.0_dp
        geometry%differential%d_dt2 = 0.0_dp

        geometry%quadrature%points_per_element = triangle_quadrature_points
        geometry%quadrature%point_count = point_count
        geometry%quadrature%point_xyz = 0.0_dp
        geometry%quadrature%point_normal = 0.0_dp
        geometry%quadrature%shape_weight = 0.0_dp
        geometry%quadrature%integration_weight = 0.0_dp
    end subroutine allocate_surface_geometry

    subroutine build_element_geometry(geometry, status, message)
        ! Compute area and an average normal for each element.
        !
        ! For a curved quadratic triangle, the surface Jacobian can vary across
        ! the element.  Therefore area is computed by quadrature rather than by
        ! a single flat-triangle formula.
        type(geometry_type), intent(inout) :: geometry
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: element_id
        integer :: q
        real(dp) :: xi(triangle_quadrature_points)
        real(dp) :: eta(triangle_quadrature_points)
        real(dp) :: weight(triangle_quadrature_points)
        real(dp) :: jacobian
        real(dp) :: area_sum
        real(dp) :: normal_sum(3)
        real(dp) :: point(3)
        real(dp) :: normal(3)
        real(dp) :: shape(quadratic_triangle_node_count)

        status = 0
        message = "Element geometry built."

        ! Get reference-triangle quadrature points and weights.
        call triangle_quadrature_25(xi, eta, weight)

        do element_id = 1, geometry%mesh%element_count
            area_sum = 0.0_dp
            normal_sum = 0.0_dp

            do q = 1, triangle_quadrature_points
                ! element_point_normal maps a reference point (xi, eta) to the
                ! real 3D surface element and returns the local Jacobian.
                call element_point_normal(geometry, element_id, xi(q), eta(q), &
                                          shape, point, normal, jacobian, status, message)
                if (status /= 0) return

                ! Surface area integral:
                !   area ~= sum weight(q) * jacobian(q)
                area_sum = area_sum + weight(q) * jacobian

                ! Weighted average normal.  The final normal is normalised
                ! after all quadrature points have contributed.
                normal_sum = normal_sum + weight(q) * jacobian * normal
            end do

            if (area_sum <= tiny_geometry) then
                status = 2
                message = "Element has near-zero area."
                return
            end if

            geometry%differential%element_area(element_id) = area_sum

            ! Store a unit normal for the whole element.
            call normalize_vector(normal_sum, geometry%differential%element_normal(:, element_id), &
                                  status, message)
            if (status /= 0) return
        end do
    end subroutine build_element_geometry

    subroutine build_node_frames(geometry, status, message)
        ! Compute one local coordinate frame at each node.
        !
        ! The frame consists of:
        !   normal    = unit vector perpendicular to the surface
        !   tangent_1 = first unit vector along the surface
        !   tangent_2 = second unit vector along the surface
        !
        ! Node normals are computed by angle-weighted averaging of surrounding
        ! subtriangle normals.  The method follows the geometric weighting idea
        ! in G. Thurmer and C. A. Wuthrich, "Computing Vertex Normals from
        ! Polygonal Facets," Journal of Graphics Tools 3(1), 1998,
        ! doi:10.1080/10867651.1998.10487487.
        !
        ! These nodal normals define collocation directions and tangent frames.
        ! Integration-point normals are computed separately from the
        ! isoparametric element map in element_point_normal.
        type(geometry_type), intent(inout) :: geometry
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: node_id
        integer :: subtriangle_id
        integer :: local_id
        integer :: nodes(linear_triangle_node_count)
        real(dp) :: face_normal(3)
        real(dp) :: angle_at_node
        real(dp) :: normal_sum(3)

        status = 0
        message = "Node normal and tangent frames built."

        do node_id = 1, geometry%mesh%node_count
            normal_sum = 0.0_dp

            ! Scan all subtriangles and collect contributions from the ones
            ! that contain this node.
            do subtriangle_id = 1, geometry%mesh%subtriangle_count
                nodes = geometry%mesh%subtriangle_nodes(:, subtriangle_id)

                do local_id = 1, linear_triangle_node_count
                    if (nodes(local_id) == node_id) then
                        ! Direction contribution from this subtriangle.
                        call subtriangle_normal(geometry, nodes, face_normal, status, message)
                        if (status /= 0) return

                        ! Weight by the angle at the node.  Larger local angles
                        ! should contribute more to the averaged node normal.
                        call subtriangle_angle_at_node(geometry, nodes, local_id, angle_at_node, &
                                                       status, message)
                        if (status /= 0) return

                        normal_sum = normal_sum + angle_at_node * face_normal
                    end if
                end do
            end do

            ! Convert the accumulated vector into a unit normal.
            call normalize_vector(normal_sum, geometry%differential%normal(:, node_id), &
                                  status, message)
            if (status /= 0) then
                message = "Unable to compute node normal."
                return
            end if

            ! Once the normal is known, choose two perpendicular tangent vectors
            ! that span the local tangent plane.
            call build_tangent_frame(geometry%differential%normal(:, node_id), &
                                     geometry%differential%tangent_1(:, node_id), &
                                     geometry%differential%tangent_2(:, node_id), &
                                     status, message)
            if (status /= 0) return
        end do
    end subroutine build_node_frames

    subroutine build_quadrature_cache(geometry, status, message)
        ! Precompute quadrature data for all elements.
        !
        ! Matrix assembly later loops over quadrature points many times.  It is
        ! cheaper and clearer to compute point positions, normals, and weighted
        ! shape functions once here.
        type(geometry_type), intent(inout) :: geometry
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: element_id
        integer :: q
        integer :: point_id
        integer :: nodes_per_element
        real(dp) :: xi(triangle_quadrature_points)
        real(dp) :: eta(triangle_quadrature_points)
        real(dp) :: weight(triangle_quadrature_points)
        real(dp) :: jacobian
        real(dp) :: point(3)
        real(dp) :: normal(3)
        real(dp) :: shape(quadratic_triangle_node_count)
        real(dp) :: integration_weight

        status = 0
        message = "Quadrature cache built."

        nodes_per_element = geometry%mesh%nodes_per_element
        call triangle_quadrature_25(xi, eta, weight)

        do element_id = 1, geometry%mesh%element_count
            do q = 1, triangle_quadrature_points
                ! Flatten the two indices (element_id, q) into one global
                ! quadrature-point index.
                point_id = triangle_quadrature_points * (element_id - 1) + q

                call element_point_normal(geometry, element_id, xi(q), eta(q), &
                                          shape, point, normal, jacobian, status, message)
                if (status /= 0) return

                integration_weight = weight(q) * jacobian

                ! Store point data in global arrays.  shape_weight already
                ! includes the integration weight, which saves multiplication
                ! during boundary-integral assembly.
                geometry%quadrature%point_xyz(:, point_id) = point
                geometry%quadrature%point_normal(:, point_id) = normal
                geometry%quadrature%integration_weight(point_id) = integration_weight
                geometry%quadrature%shape_weight(:, point_id) = &
                    shape(1:nodes_per_element) * integration_weight
            end do
        end do
    end subroutine build_quadrature_cache

    subroutine build_derivative_stencils(geometry, status, message)
        ! Build local weights for first derivatives along tangent_1 and tangent_2.
        !
        ! Suppose f is a scalar value stored at mesh nodes.  Around each node,
        ! this routine creates weights so that:
        !   df/dt1 ~= sum stencil_weight_1(j) * f(j)
        !   df/dt2 ~= sum stencil_weight_2(j) * f(j)
        !
        ! The neighbour coordinates are first projected onto the local tangent
        ! plane.  A small least-squares system is then used to produce derivative
        ! weights in the two tangent directions.
        type(geometry_type), intent(inout) :: geometry
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: node_id
        integer :: neighbor_row
        integer :: neighbor_id
        integer :: active_count
        real(dp) :: delta(3)
        real(dp) :: local_u
        real(dp) :: local_v
        real(dp) :: normal_offset
        real(dp) :: matrix_11
        real(dp) :: matrix_12
        real(dp) :: matrix_22
        real(dp) :: determinant
        real(dp) :: inv_11
        real(dp) :: inv_12
        real(dp) :: inv_22
        real(dp) :: weight_u
        real(dp) :: weight_v

        status = 0
        message = "Surface derivative stencils built."

        do node_id = 1, geometry%mesh%node_count
            matrix_11 = 0.0_dp
            matrix_12 = 0.0_dp
            matrix_22 = 0.0_dp
            active_count = 0

            ! First pass: build the 2x2 normal-equation matrix from the
            ! projected neighbour coordinates.
            do neighbor_row = 1, geometry%mesh%max_node_neighbors_2ring
                neighbor_id = geometry%mesh%node_neighbors_2ring(neighbor_row, node_id)
                if (neighbor_id == 0) cycle

                ! Vector from the host node to the neighbour.
                delta = geometry%mesh%xyz(:, neighbor_id) - geometry%mesh%xyz(:, node_id)

                ! Project that vector onto the tangent plane by removing its
                ! normal component.
                normal_offset = dot_product(delta, geometry%differential%normal(:, node_id))
                delta = delta - normal_offset * geometry%differential%normal(:, node_id)

                ! Coordinates of the projected neighbour in the local tangent
                ! basis.
                local_u = dot_product(delta, geometry%differential%tangent_1(:, node_id))
                local_v = dot_product(delta, geometry%differential%tangent_2(:, node_id))

                matrix_11 = matrix_11 + local_u * local_u
                matrix_12 = matrix_12 + local_u * local_v
                matrix_22 = matrix_22 + local_v * local_v
                active_count = active_count + 1
            end do

            if (active_count < 2) then
                status = 3
                message = "Node has too few neighbors for derivative stencils."
                return
            end if

            ! The determinant checks whether the local neighbours span a real
            ! two-dimensional patch.  If all neighbours are nearly collinear,
            ! derivative weights cannot be built reliably.
            determinant = matrix_11 * matrix_22 - matrix_12 * matrix_12
            if (abs(determinant) <= tiny_geometry) then
                status = 3
                message = "Node derivative stencil is singular."
                return
            end if

            inv_11 = matrix_22 / determinant
            inv_12 = -matrix_12 / determinant
            inv_22 = matrix_11 / determinant

            ! Second pass: convert each neighbour coordinate into derivative
            ! weights.  The host-node weight is the negative sum of neighbour
            ! weights, so the derivative of a constant field is exactly zero.
            do neighbor_row = 1, geometry%mesh%max_node_neighbors_2ring
                neighbor_id = geometry%mesh%node_neighbors_2ring(neighbor_row, node_id)
                if (neighbor_id == 0) cycle

                delta = geometry%mesh%xyz(:, neighbor_id) - geometry%mesh%xyz(:, node_id)
                normal_offset = dot_product(delta, geometry%differential%normal(:, node_id))
                delta = delta - normal_offset * geometry%differential%normal(:, node_id)
                local_u = dot_product(delta, geometry%differential%tangent_1(:, node_id))
                local_v = dot_product(delta, geometry%differential%tangent_2(:, node_id))

                weight_u = inv_11 * local_u + inv_12 * local_v
                weight_v = inv_12 * local_u + inv_22 * local_v

                geometry%differential%d_dt1(neighbor_row + 1, node_id) = weight_u
                geometry%differential%d_dt2(neighbor_row + 1, node_id) = weight_v
                geometry%differential%d_dt1(1, node_id) = &
                    geometry%differential%d_dt1(1, node_id) - weight_u
                geometry%differential%d_dt2(1, node_id) = &
                    geometry%differential%d_dt2(1, node_id) - weight_v
            end do
        end do
    end subroutine build_derivative_stencils

    subroutine element_point_normal(geometry, element_id, xi, eta, shape, point, normal, &
                                    jacobian, status, message)
        ! Map one reference-triangle point to a real surface element.
        !
        ! Inputs xi and eta are coordinates on the reference triangle.  The
        ! shape functions interpolate the physical node coordinates to produce:
        !   point    = 3D position on the actual element
        !   normal   = unit surface normal at that point
        !   jacobian = local area scaling from reference triangle to 3D surface
        !
        ! The normal intentionally comes from cross(dx/dxi, dx/deta), not from
        ! interpolating the angle-averaged nodal normals.  This keeps the normal
        ! and Jacobian consistent with the same isoparametric geometry map.  A
        ! linear element therefore has one flat normal, whereas a curved
        ! quadratic element can have a different normal at every quadrature
        ! point.
        type(geometry_type), intent(in) :: geometry
        integer, intent(in) :: element_id
        real(dp), intent(in) :: xi
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: shape(quadratic_triangle_node_count)
        real(dp), intent(out) :: point(3)
        real(dp), intent(out) :: normal(3)
        real(dp), intent(out) :: jacobian
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: local_id
        integer :: node_id
        integer :: nodes_per_element
        real(dp) :: dshape_dxi(quadratic_triangle_node_count)
        real(dp) :: dshape_deta(quadratic_triangle_node_count)
        real(dp) :: dx_dxi(3)
        real(dp) :: dx_deta(3)
        real(dp) :: normal_times_jacobian(3)

        status = 0
        message = "Element point and normal computed."

        nodes_per_element = geometry%mesh%nodes_per_element

        ! Shape functions tell us how much each element node contributes to the
        ! interpolated point.  Their derivatives give tangent vectors.
        call triangle_shape_functions(nodes_per_element, xi, eta, shape, dshape_dxi, dshape_deta)

        point = 0.0_dp
        dx_dxi = 0.0_dp
        dx_deta = 0.0_dp

        ! Interpolate the point and its two reference-coordinate derivatives.
        do local_id = 1, nodes_per_element
            node_id = geometry%mesh%element_nodes(local_id, element_id)
            point = point + shape(local_id) * geometry%mesh%xyz(:, node_id)
            dx_dxi = dx_dxi + dshape_dxi(local_id) * geometry%mesh%xyz(:, node_id)
            dx_deta = dx_deta + dshape_deta(local_id) * geometry%mesh%xyz(:, node_id)
        end do

        ! The cross-product order follows the local node orientation.  Reversing
        ! the corner order reverses this normal.  Its length is the surface
        ! Jacobian and is independent of that sign.
        normal_times_jacobian = cross_product(dx_dxi, dx_deta)
        jacobian = vector_norm(normal_times_jacobian)

        if (jacobian <= tiny_geometry) then
            status = 2
            message = "Element has near-zero quadrature Jacobian."
            return
        end if

        ! Divide by the Jacobian to convert normal_times_jacobian into a unit
        ! normal vector.
        normal = normal_times_jacobian / jacobian
    end subroutine element_point_normal

    subroutine triangle_shape_functions(nodes_per_element, xi, eta, shape, dshape_dxi, dshape_deta)
        ! Evaluate triangle interpolation shape functions and their derivatives.
        !
        ! Linear triangle:
        !   3 nodes, straight sides, first-order interpolation.
        !
        ! Quadratic triangle:
        !   6 nodes, curved sides possible, second-order interpolation.
        !
        ! The local node order for quadratic triangles is assumed to be:
        !   1, 2, 3 = corner nodes
        !   4       = midside node between 1 and 2
        !   5       = midside node between 2 and 3
        !   6       = midside node between 3 and 1
        integer, intent(in) :: nodes_per_element
        real(dp), intent(in) :: xi
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: shape(quadratic_triangle_node_count)
        real(dp), intent(out) :: dshape_dxi(quadratic_triangle_node_count)
        real(dp), intent(out) :: dshape_deta(quadratic_triangle_node_count)

        real(dp) :: l1
        real(dp) :: l2
        real(dp) :: l3

        shape = 0.0_dp
        dshape_dxi = 0.0_dp
        dshape_deta = 0.0_dp

        ! Barycentric coordinates on the reference triangle.
        ! They always satisfy l1 + l2 + l3 = 1.
        l1 = 1.0_dp - xi - eta
        l2 = xi
        l3 = eta

        if (nodes_per_element == linear_triangle_node_count) then
            ! Linear interpolation is exactly the barycentric coordinates.
            shape(1:3) = [l1, l2, l3]
            dshape_dxi(1:3) = [-1.0_dp, 1.0_dp, 0.0_dp]
            dshape_deta(1:3) = [-1.0_dp, 0.0_dp, 1.0_dp]
            return
        end if

        ! Quadratic Lagrange shape functions.
        shape(1) = l1 * (2.0_dp * l1 - 1.0_dp)
        shape(2) = l2 * (2.0_dp * l2 - 1.0_dp)
        shape(3) = l3 * (2.0_dp * l3 - 1.0_dp)
        shape(4) = 4.0_dp * l1 * l2
        shape(5) = 4.0_dp * l2 * l3
        shape(6) = 4.0_dp * l3 * l1

        dshape_dxi(1) = 1.0_dp - 4.0_dp * l1
        dshape_dxi(2) = 4.0_dp * l2 - 1.0_dp
        dshape_dxi(3) = 0.0_dp
        dshape_dxi(4) = 4.0_dp * (l1 - l2)
        dshape_dxi(5) = 4.0_dp * l3
        dshape_dxi(6) = -4.0_dp * l3

        dshape_deta(1) = 1.0_dp - 4.0_dp * l1
        dshape_deta(2) = 0.0_dp
        dshape_deta(3) = 4.0_dp * l3 - 1.0_dp
        dshape_deta(4) = -4.0_dp * l2
        dshape_deta(5) = 4.0_dp * l2
        dshape_deta(6) = 4.0_dp * (l1 - l3)
    end subroutine triangle_shape_functions

    subroutine subtriangle_normal(geometry, nodes, normal, status, message)
        ! Compute the unit normal of a 3-node subtriangle.
        !
        ! This helper is used when averaging subtriangle normals to build node
        ! normals.  It uses the order of nodes in the subtriangle, so mesh
        ! orientation must be consistent.
        type(geometry_type), intent(in) :: geometry
        integer, intent(in) :: nodes(linear_triangle_node_count)
        real(dp), intent(out) :: normal(3)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        real(dp) :: edge_a(3)
        real(dp) :: edge_b(3)

        status = 0
        message = "Subtriangle normal computed."

        ! Two edges from node 1 define the triangle plane.
        edge_a = geometry%mesh%xyz(:, nodes(2)) - geometry%mesh%xyz(:, nodes(1))
        edge_b = geometry%mesh%xyz(:, nodes(3)) - geometry%mesh%xyz(:, nodes(1))

        ! cross(edge_a, edge_b) points perpendicular to the triangle.
        call normalize_vector(cross_product(edge_a, edge_b), normal, status, message)
    end subroutine subtriangle_normal

    subroutine subtriangle_angle_at_node(geometry, nodes, local_id, angle, status, message)
        ! Compute the corner angle at one node of a subtriangle.
        !
        ! This angle is used as the weight in angle-weighted normal averaging.
        type(geometry_type), intent(in) :: geometry
        integer, intent(in) :: nodes(linear_triangle_node_count)
        integer, intent(in) :: local_id
        real(dp), intent(out) :: angle
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: node_center
        integer :: node_a
        integer :: node_b
        real(dp) :: edge_a(3)
        real(dp) :: edge_b(3)
        real(dp) :: norm_a
        real(dp) :: norm_b
        real(dp) :: cosine_angle

        status = 0
        message = "Subtriangle angle computed."

        ! Pick the central node and the other two nodes of the subtriangle.
        node_center = nodes(local_id)
        node_a = nodes(mod(local_id, linear_triangle_node_count) + 1)
        node_b = nodes(mod(local_id + 1, linear_triangle_node_count) + 1)

        ! Build the two edges that meet at the central node.
        edge_a = geometry%mesh%xyz(:, node_a) - geometry%mesh%xyz(:, node_center)
        edge_b = geometry%mesh%xyz(:, node_b) - geometry%mesh%xyz(:, node_center)
        norm_a = vector_norm(edge_a)
        norm_b = vector_norm(edge_b)

        if (norm_a <= tiny_geometry .or. norm_b <= tiny_geometry) then
            status = 2
            message = "Subtriangle has a near-zero edge."
            return
        end if

        ! Dot product formula:
        !   cos(angle) = (a dot b) / (|a| |b|)
        !
        ! The max/min clamp protects acos from tiny roundoff errors such as
        ! 1.0000000000000002.
        cosine_angle = dot_product(edge_a, edge_b) / (norm_a * norm_b)
        cosine_angle = max(-1.0_dp, min(1.0_dp, cosine_angle))
        angle = acos(cosine_angle)
    end subroutine subtriangle_angle_at_node

    subroutine build_tangent_frame(normal, tangent_1, tangent_2, status, message)
        ! Build two unit tangent vectors perpendicular to the supplied normal.
        !
        ! There are infinitely many valid tangent directions.  This routine
        ! chooses a stable reference axis and constructs an orthonormal frame by
        ! cross products.
        real(dp), intent(in) :: normal(3)
        real(dp), intent(out) :: tangent_1(3)
        real(dp), intent(out) :: tangent_2(3)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        real(dp) :: reference_axis(3)
        real(dp) :: tangent_2_candidate(3)

        status = 0
        message = "Tangent frame built."

        ! Choose the coordinate axis least aligned with the normal.  This avoids
        ! taking a cross product of nearly parallel vectors.
        if (abs(normal(1)) <= abs(normal(2)) .and. abs(normal(1)) <= abs(normal(3))) then
            reference_axis = [1.0_dp, 0.0_dp, 0.0_dp]
        else if (abs(normal(2)) <= abs(normal(3))) then
            reference_axis = [0.0_dp, 1.0_dp, 0.0_dp]
        else
            reference_axis = [0.0_dp, 0.0_dp, 1.0_dp]
        end if

        ! First tangent is perpendicular to both the reference axis and normal.
        call normalize_vector(cross_product(reference_axis, normal), tangent_1, status, message)
        if (status /= 0) return

        ! Second tangent completes the right-handed local frame.
        tangent_2_candidate = cross_product(normal, tangent_1)
        call normalize_vector(tangent_2_candidate, tangent_2, status, message)
    end subroutine build_tangent_frame

    subroutine normalize_vector(vector, unit_vector, status, message)
        ! Convert a vector into a unit vector.
        !
        ! If the input vector is too small, normalisation would divide by a
        ! near-zero number, so the routine returns an error instead.
        real(dp), intent(in) :: vector(3)
        real(dp), intent(out) :: unit_vector(3)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        real(dp) :: magnitude

        status = 0
        message = "Vector normalized."

        magnitude = vector_norm(vector)
        if (magnitude <= tiny_geometry) then
            unit_vector = 0.0_dp
            status = 2
            message = "Cannot normalize a near-zero vector."
            return
        end if

        unit_vector = vector / magnitude
    end subroutine normalize_vector

    pure function vector_norm(vector) result(magnitude)
        ! Euclidean length of a 3D vector:
        !   sqrt(x^2 + y^2 + z^2)
        real(dp), intent(in) :: vector(3)
        real(dp) :: magnitude

        magnitude = sqrt(dot_product(vector, vector))
    end function vector_norm

    pure function cross_product(a, b) result(c)
        ! 3D cross product.
        !
        ! The result is perpendicular to both input vectors.  In geometry code
        ! this is commonly used to compute a surface normal from two tangents.
        real(dp), intent(in) :: a(3)
        real(dp), intent(in) :: b(3)
        real(dp) :: c(3)

        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
    end function cross_product

end module Geom_SurfaceGeometry
