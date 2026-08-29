module Geom_MeshTopology
    ! This module builds mesh topology.
    !
    ! The mesh reader gives us the basic information:
    !   - node coordinates
    !   - element node lists
    !
    ! The solver also needs derived connectivity information:
    !   - which elements touch each node
    !   - which nodes are neighbours of each node
    !   - how a 6-node quadratic triangle can be viewed as smaller 3-node
    !     triangles
    !
    ! These derived lists are called "topology" because they describe how the
    ! mesh is connected, not where the mesh is located in space.
    use Geom_Types, only: geometry_type, linear_triangle_node_count, quadratic_triangle_node_count
    implicit none

    ! Only expose the main preparation routine.  The helper routines below are
    ! implementation details used to build the topology step by step.
    private

    public :: prepare_mesh_topology

contains

    subroutine prepare_mesh_topology(geometry, status, message)
        ! Build all topology arrays needed after a mesh has been read.
        !
        ! Expected input:
        !   geometry%mesh%xyz
        !   geometry%mesh%element_nodes
        !
        ! Produced output:
        !   geometry%mesh%subtriangle_nodes
        !   geometry%mesh%node_elements
        !   geometry%mesh%node_neighbors
        !   geometry%mesh%node_neighbors_2ring
        !
        ! The steps are ordered carefully.  For example, two-ring neighbours
        ! require first-ring neighbours, so first-ring neighbours must be built
        ! first.
        type(geometry_type), intent(inout) :: geometry
        integer, intent(out), optional :: status
        character(len=*), intent(out), optional :: message

        ! Internal error state passed between helper routines.
        integer :: error_code
        character(len=256) :: error_message

        error_code = 0
        error_message = "Mesh topology prepared successfully."

        ! Remove any previous topology arrays before rebuilding them.
        call clear_topology(geometry)

        ! Check that the mesh data is present and has consistent sizes.
        call validate_mesh(geometry, error_code, error_message)
        if (error_code /= 0) then
            call finish(error_code, error_message)
            return
        end if

        ! Build a 3-node-triangle view of the mesh.  This is especially useful
        ! for quadratic elements, which are split into four subtriangles.
        call build_subtriangles(geometry, error_code, error_message)
        if (error_code /= 0) then
            call finish(error_code, error_message)
            return
        end if

        ! For each node, store the elements that contain that node.
        call build_node_element_links(geometry, error_code, error_message)
        if (error_code /= 0) then
            call finish(error_code, error_message)
            return
        end if

        ! For each node, store nearby nodes that share a subtriangle with it.
        call build_first_ring_neighbors(geometry, error_code, error_message)
        if (error_code /= 0) then
            call finish(error_code, error_message)
            return
        end if

        ! For each node, store first-ring neighbours plus neighbours of those
        ! neighbours.  This wider local patch is used by derivative stencils.
        call build_two_ring_neighbors(geometry, error_code, error_message)
        call finish(error_code, error_message)

    contains

        subroutine finish(code, text)
            ! Return status information to the caller.
            !
            ! If any step failed, discard partly-built topology arrays so the
            ! geometry object is not left in a confusing half-prepared state.
            integer, intent(in) :: code
            character(len=*), intent(in) :: text

            if (code /= 0) call clear_topology(geometry)
            if (present(status)) status = code
            if (present(message)) message = text
        end subroutine finish

    end subroutine prepare_mesh_topology

    subroutine clear_topology(geometry)
        ! Remove topology arrays and reset their counters.
        !
        ! This does not remove the original mesh coordinates or element
        ! connectivity read from the mesh file.
        type(geometry_type), intent(inout) :: geometry

        if (allocated(geometry%mesh%node_elements)) deallocate(geometry%mesh%node_elements)
        if (allocated(geometry%mesh%node_neighbors)) deallocate(geometry%mesh%node_neighbors)
        if (allocated(geometry%mesh%node_neighbors_2ring)) deallocate(geometry%mesh%node_neighbors_2ring)
        if (allocated(geometry%mesh%subtriangle_nodes)) deallocate(geometry%mesh%subtriangle_nodes)

        geometry%mesh%max_node_elements = 0
        geometry%mesh%max_node_neighbors = 0
        geometry%mesh%max_node_neighbors_2ring = 0
        geometry%mesh%subtriangle_count = 0
    end subroutine clear_topology

    subroutine validate_mesh(geometry, status, message)
        ! Check that the mesh has enough valid information to build topology.
        !
        ! This routine does not judge mesh quality.  It only checks basic
        ! consistency, such as whether arrays are allocated and whether element
        ! node ids point to existing nodes.
        type(geometry_type), intent(in) :: geometry
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        ! element_id loops over elements, local_id loops over nodes inside one
        ! element, and node_id stores the global node number found there.
        integer :: element_id
        integer :: local_id
        integer :: node_id
        integer :: nodes_per_element

        status = 0
        message = "Mesh is valid."

        if (geometry%mesh%node_count <= 0) then
            call fail(1, "Mesh has no nodes.")
            return
        end if

        if (geometry%mesh%element_count <= 0) then
            call fail(1, "Mesh has no elements.")
            return
        end if

        nodes_per_element = geometry%mesh%nodes_per_element
        if (nodes_per_element /= linear_triangle_node_count .and. nodes_per_element /= quadratic_triangle_node_count) then
            call fail(1, "Mesh nodes_per_element must be 3 or 6.")
            return
        end if

        if (.not. allocated(geometry%mesh%xyz)) then
            call fail(1, "Mesh coordinates are not allocated.")
            return
        end if

        if (.not. allocated(geometry%mesh%element_nodes)) then
            call fail(1, "Mesh element connectivity is not allocated.")
            return
        end if

        ! Coordinates must be stored as xyz(1:3, 1:node_count).
        if (size(geometry%mesh%xyz, 1) /= 3 .or. &
            size(geometry%mesh%xyz, 2) /= geometry%mesh%node_count) then
            call fail(1, "Mesh coordinates have inconsistent dimensions.")
            return
        end if

        ! Connectivity must be stored as element_nodes(1:nodes_per_element,
        ! 1:element_count).
        if (size(geometry%mesh%element_nodes, 1) /= nodes_per_element .or. &
            size(geometry%mesh%element_nodes, 2) /= geometry%mesh%element_count) then
            call fail(1, "Mesh element connectivity has inconsistent dimensions.")
            return
        end if

        ! Every node id mentioned by every element must be a valid node number.
        do element_id = 1, geometry%mesh%element_count
            do local_id = 1, nodes_per_element
                node_id = geometry%mesh%element_nodes(local_id, element_id)
                if (node_id < 1 .or. node_id > geometry%mesh%node_count) then
                    call fail(1, "Mesh element connectivity references an invalid node id.")
                    return
                end if
            end do
        end do

    contains

        subroutine fail(code, text)
            ! Store a validation failure code and message.
            integer, intent(in) :: code
            character(len=*), intent(in) :: text

            status = code
            message = text
        end subroutine fail

    end subroutine validate_mesh

    subroutine build_subtriangles(geometry, status, message)
        ! Build a simple 3-node triangle representation of the mesh.
        !
        ! Why this is useful:
        !   - Many geometric operations are easiest on 3-node triangles.
        !   - First-ring neighbours can be found by asking which nodes share a
        !     subtriangle.
        !   - A quadratic triangle still has curved interpolation, but its
        !     topology can be split into four small linear pieces.
        type(geometry_type), intent(inout) :: geometry
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: element_id
        integer :: subtriangle_id
        integer :: nodes(quadratic_triangle_node_count)

        status = 0
        message = "Subtriangle topology built."

        ! Linear elements already are 3-node triangles.  Quadratic elements are
        ! split into four 3-node subtriangles.
        if (geometry%mesh%nodes_per_element == linear_triangle_node_count) then
            geometry%mesh%subtriangle_count = geometry%mesh%element_count
        else
            geometry%mesh%subtriangle_count = 4 * geometry%mesh%element_count
        end if

        allocate(geometry%mesh%subtriangle_nodes(linear_triangle_node_count, geometry%mesh%subtriangle_count), &
                 stat=status)
        if (status /= 0) then
            message = "Unable to allocate subtriangle topology."
            return
        end if

        ! For a linear mesh, the subtriangle view is exactly the original
        ! element connectivity.
        if (geometry%mesh%nodes_per_element == linear_triangle_node_count) then
            geometry%mesh%subtriangle_nodes = geometry%mesh%element_nodes
            return
        end if

        ! For a quadratic triangle, nodes 4, 5, and 6 are midside nodes.  Using
        ! those midside nodes divides the original triangle into four smaller
        ! triangles that cover the same element topology.
        do element_id = 1, geometry%mesh%element_count
            nodes = geometry%mesh%element_nodes(:, element_id)
            subtriangle_id = 4 * (element_id - 1)

            ! For quadratic triangle nodes ordered as:
            !   1, 2, 3 = corner nodes
            !   4, 5, 6 = midside nodes on edges 1-2, 2-3, 3-1
            geometry%mesh%subtriangle_nodes(:, subtriangle_id + 1) = [nodes(1), nodes(4), nodes(6)]
            geometry%mesh%subtriangle_nodes(:, subtriangle_id + 2) = [nodes(4), nodes(2), nodes(5)]
            geometry%mesh%subtriangle_nodes(:, subtriangle_id + 3) = [nodes(5), nodes(3), nodes(6)]
            geometry%mesh%subtriangle_nodes(:, subtriangle_id + 4) = [nodes(4), nodes(5), nodes(6)]
        end do
    end subroutine build_subtriangles

    subroutine build_node_element_links(geometry, status, message)
        ! Build a list of elements attached to each node.
        !
        ! Example:
        !   node_elements(:, 10) stores all elements that contain node 10.
        !
        ! This is useful when calculating node-based quantities from
        ! element-based quantities, such as averaging element normals around a
        ! node.
        type(geometry_type), intent(inout) :: geometry
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: element_id
        integer :: local_id
        integer :: node_id
        integer :: nodes_per_element
        integer, allocatable :: link_count(:)

        status = 0
        message = "Node-to-element links built."
        nodes_per_element = geometry%mesh%nodes_per_element

        ! First pass: count how many elements touch each node.
        !
        ! Fortran arrays are rectangular, so we need to know the largest count
        ! before allocating node_elements.
        allocate(link_count(geometry%mesh%node_count), stat=status)
        if (status /= 0) then
            message = "Unable to allocate node-element link counts."
            return
        end if

        link_count = 0
        do element_id = 1, geometry%mesh%element_count
            do local_id = 1, nodes_per_element
                node_id = geometry%mesh%element_nodes(local_id, element_id)
                link_count(node_id) = link_count(node_id) + 1
            end do
        end do

        ! Allocate enough rows for the most connected node.  Nodes with fewer
        ! attached elements will have trailing zero entries.
        geometry%mesh%max_node_elements = maxval(link_count)
        allocate(geometry%mesh%node_elements(geometry%mesh%max_node_elements, geometry%mesh%node_count), &
                 stat=status)
        if (status /= 0) then
            message = "Unable to allocate node-element links."
            deallocate(link_count)
            return
        end if

        geometry%mesh%node_elements = 0
        link_count = 0

        ! Second pass: fill the actual element ids into the list for each node.
        do element_id = 1, geometry%mesh%element_count
            do local_id = 1, nodes_per_element
                node_id = geometry%mesh%element_nodes(local_id, element_id)
                link_count(node_id) = link_count(node_id) + 1
                geometry%mesh%node_elements(link_count(node_id), node_id) = element_id
            end do
        end do

        deallocate(link_count)
    end subroutine build_node_element_links

    subroutine build_first_ring_neighbors(geometry, status, message)
        ! Build the first-ring neighbour list for every node.
        !
        ! A first-ring neighbour is another node that shares at least one
        ! subtriangle with the current node.  This gives a local surface patch
        ! around each node.
        type(geometry_type), intent(inout) :: geometry
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: node_id
        integer :: neighbor_id
        integer :: neighbor_count
        integer, allocatable :: neighbor_counts(:)
        logical, allocatable :: is_neighbor(:)

        status = 0
        message = "First-ring node neighbors built."

        ! is_neighbor is a temporary true/false marker array.  For one selected
        ! node, is_neighbor(j)=true means node j is currently marked as a
        ! neighbour.
        allocate(neighbor_counts(geometry%mesh%node_count), &
                 is_neighbor(geometry%mesh%node_count), &
                 stat=status)
        if (status /= 0) then
            message = "Unable to allocate first-ring neighbor work arrays."
            return
        end if

        ! First pass: count neighbours for each node so the rectangular storage
        ! array can be allocated with enough rows.
        do node_id = 1, geometry%mesh%node_count
            call mark_first_ring(geometry, node_id, is_neighbor)
            neighbor_counts(node_id) = count(is_neighbor)
        end do

        ! Allocate enough rows for the node with the most first-ring neighbours.
        geometry%mesh%max_node_neighbors = maxval(neighbor_counts)
        allocate(geometry%mesh%node_neighbors(geometry%mesh%max_node_neighbors, geometry%mesh%node_count), &
                 stat=status)
        if (status /= 0) then
            message = "Unable to allocate first-ring node neighbors."
            deallocate(neighbor_counts, is_neighbor)
            return
        end if

        geometry%mesh%node_neighbors = 0

        ! Second pass: write neighbour ids into node_neighbors(:, node_id).
        do node_id = 1, geometry%mesh%node_count
            call mark_first_ring(geometry, node_id, is_neighbor)
            neighbor_count = 0
            do neighbor_id = 1, geometry%mesh%node_count
                if (is_neighbor(neighbor_id)) then
                    neighbor_count = neighbor_count + 1
                    geometry%mesh%node_neighbors(neighbor_count, node_id) = neighbor_id
                end if
            end do
        end do

        deallocate(neighbor_counts, is_neighbor)
    end subroutine build_first_ring_neighbors

    subroutine build_two_ring_neighbors(geometry, status, message)
        ! Build the two-ring neighbour list for every node.
        !
        ! A two-ring patch contains:
        !   - first-ring neighbours of the node
        !   - first-ring neighbours of those first-ring neighbours
        !
        ! This wider patch provides enough nearby points for local derivative
        ! fitting on the surface.
        type(geometry_type), intent(inout) :: geometry
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: node_id
        integer :: neighbor_id
        integer :: neighbor_count
        integer, allocatable :: neighbor_counts(:)
        logical, allocatable :: is_neighbor(:)

        status = 0
        message = "Two-ring node neighbors built."

        allocate(neighbor_counts(geometry%mesh%node_count), &
                 is_neighbor(geometry%mesh%node_count), &
                 stat=status)
        if (status /= 0) then
            message = "Unable to allocate two-ring neighbor work arrays."
            return
        end if

        ! First pass: count two-ring neighbours for each node.
        do node_id = 1, geometry%mesh%node_count
            call mark_two_ring(geometry, node_id, is_neighbor)
            neighbor_counts(node_id) = count(is_neighbor)
        end do

        ! Allocate enough rows for the node with the largest two-ring patch.
        geometry%mesh%max_node_neighbors_2ring = maxval(neighbor_counts)
        allocate(geometry%mesh%node_neighbors_2ring( &
                 geometry%mesh%max_node_neighbors_2ring, geometry%mesh%node_count), &
                 stat=status)
        if (status /= 0) then
            message = "Unable to allocate two-ring node neighbors."
            deallocate(neighbor_counts, is_neighbor)
            return
        end if

        geometry%mesh%node_neighbors_2ring = 0

        ! Second pass: write two-ring neighbour ids into the final array.
        do node_id = 1, geometry%mesh%node_count
            call mark_two_ring(geometry, node_id, is_neighbor)
            neighbor_count = 0
            do neighbor_id = 1, geometry%mesh%node_count
                if (is_neighbor(neighbor_id)) then
                    neighbor_count = neighbor_count + 1
                    geometry%mesh%node_neighbors_2ring(neighbor_count, node_id) = neighbor_id
                end if
            end do
        end do

        deallocate(neighbor_counts, is_neighbor)
    end subroutine build_two_ring_neighbors

    subroutine mark_first_ring(geometry, node_id, is_neighbor)
        ! Mark all first-ring neighbours of node_id.
        !
        ! This routine does not allocate or store the final neighbour list.  It
        ! only fills the logical array is_neighbor.  The calling routines then
        ! either count the true entries or copy their node ids into storage.
        type(geometry_type), intent(in) :: geometry
        integer, intent(in) :: node_id
        logical, intent(out) :: is_neighbor(:)

        integer :: subtriangle_id
        integer :: local_id
        integer :: candidate_id
        logical :: node_is_on_subtriangle

        is_neighbor = .false.

        ! Scan every subtriangle.  If node_id is part of that subtriangle, the
        ! other nodes of the same subtriangle are first-ring neighbours.
        do subtriangle_id = 1, geometry%mesh%subtriangle_count
            node_is_on_subtriangle = any(geometry%mesh%subtriangle_nodes(:, subtriangle_id) == node_id)

            if (node_is_on_subtriangle) then
                do local_id = 1, linear_triangle_node_count
                    candidate_id = geometry%mesh%subtriangle_nodes(local_id, subtriangle_id)
                    if (candidate_id /= node_id) is_neighbor(candidate_id) = .true.
                end do
            end if
        end do
    end subroutine mark_first_ring

    subroutine mark_two_ring(geometry, node_id, is_neighbor)
        ! Mark all two-ring neighbours of node_id.
        !
        ! This uses the first-ring neighbour list that has already been built.
        ! For each first-ring neighbour, it also marks that neighbour's
        ! first-ring neighbours.
        type(geometry_type), intent(in) :: geometry
        integer, intent(in) :: node_id
        logical, intent(out) :: is_neighbor(:)

        integer :: first_ring_index
        integer :: first_ring_node
        integer :: second_ring_index
        integer :: second_ring_node

        is_neighbor = .false.

        ! Loop over all first-ring neighbours of node_id.
        do first_ring_index = 1, geometry%mesh%max_node_neighbors
            first_ring_node = geometry%mesh%node_neighbors(first_ring_index, node_id)

            ! A zero entry means this row is unused for this node.
            if (first_ring_node == 0) cycle

            is_neighbor(first_ring_node) = .true.

            ! Add neighbours of this first-ring neighbour.  These form the
            ! second ring around node_id.
            do second_ring_index = 1, geometry%mesh%max_node_neighbors
                second_ring_node = geometry%mesh%node_neighbors(second_ring_index, first_ring_node)
                if (second_ring_node == 0) cycle
                is_neighbor(second_ring_node) = .true.
            end do
        end do

        ! The central node should not be listed as its own neighbour.
        is_neighbor(node_id) = .false.
    end subroutine mark_two_ring

end module Geom_MeshTopology
