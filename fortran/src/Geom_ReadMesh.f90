module Geom_ReadMesh
    ! This module reads a surface mesh from a simple text file and stores it in
    ! a geometry_type object.
    !
    ! In the public workflow, Python is responsible for pre-processing tasks such
    ! as creating or converting meshes.  Fortran then reads the prepared mesh
    ! and performs the numerically expensive boundary element calculations.
    !
    ! Keeping the mesh file format simple makes the interface between Python
    ! and Fortran easier to test and easier for new users to understand.
    !
    ! This reader checks record syntax, ids, counts, and array bounds.  It does
    ! not prove that a surface is closed, manifold, or consistently oriented;
    ! the public Python mesh audit performs those topology checks before writing
    ! the normalized file consumed here.
    use Pre_Constants, only: dp
    use Geom_Types, only: geometry_type, linear_triangle_node_count, quadratic_triangle_node_count
    implicit none

    ! Hide all internal helper names.  Other modules should only call read_mesh.
    private

    public :: read_mesh

contains

    subroutine read_mesh(filename, geometry, status, message)
        ! Read one mesh file into geometry.
        !
        ! Inputs:
        !   filename = path to the mesh text file
        !
        ! Output:
        !   geometry = filled geometry object containing nodes, elements, and
        !              particle id information
        !
        ! Optional outputs:
        !   status  = 0 for success, non-zero for an error
        !   message = short explanation of what happened
        !
        ! The subroutine clears geometry first.  This prevents data from an old
        ! mesh from being mixed with data from a new mesh.
        character(len=*), intent(in) :: filename
        type(geometry_type), intent(inout) :: geometry
        integer, intent(out), optional :: status
        character(len=*), intent(out), optional :: message

        ! unit is the Fortran file handle used while the mesh file is open.
        integer :: unit

        ! io_status and io_message store information returned by open/read.
        ! A non-zero io_status means the file operation failed.
        integer :: io_status

        ! Header values read from the first line of the mesh file.
        integer :: node_count
        integer :: element_count
        integer :: nodes_per_element
        integer :: particle_count

        ! Values read from one node or element record.
        integer :: node_id
        integer :: element_id
        integer :: particle_id

        ! Loop counters for reading node and element records.
        integer :: node_index
        integer :: element_index

        ! Temporary storage for the node ids that define one element.
        ! The array is always length 6 so it can hold either a linear triangle
        ! (first 3 entries used) or a quadratic triangle (all 6 entries used).
        integer :: local_nodes(quadratic_triangle_node_count)

        ! Coordinates read from one node record.
        real(dp) :: x_coord
        real(dp) :: y_coord
        real(dp) :: z_coord

        character(len=256) :: io_message

        ! -1 means "no file is currently open".  This makes it safe for fail()
        ! to decide whether it needs to close a file.
        unit = -1
        io_message = ""

        ! The public mesh file is deliberately simple because Python writes it.
        !
        ! Header:
        !   node_count element_count nodes_per_element particle_count
        !
        ! Node records, repeated node_count times:
        !   node_id x y z particle_id
        !
        ! Element records, repeated element_count times:
        !   element_id particle_id n1 n2 n3             for linear triangles
        !   element_id particle_id n1 n2 n3 n4 n5 n6    for quadratic triangles
        !
        ! For a quadratic triangle, n1-n3 are oriented corner nodes and n4,
        ! n5, n6 lie on edges n1-n2, n2-n3, and n3-n1 respectively.  Reordering
        ! these entries changes the geometry map and can reverse or corrupt
        ! surface normals.
        !
        ! The reader expects one-based, sequential node_id and element_id values.
        ! This matches Fortran array indexing and keeps the public interface clear.

        ! Remove any old mesh data before reading a new file.
        call geometry%clear()

        ! Open the mesh file for reading.  newunit asks Fortran to choose a free
        ! file unit number automatically, avoiding clashes with other files.
        open(newunit=unit, file=filename, status="old", action="read", &
             iostat=io_status, iomsg=io_message)
        if (io_status /= 0) then
            call fail(1, "Unable to open mesh file: " // trim(io_message))
            return
        end if

        ! Read the header line.  This tells the reader how much memory to
        ! allocate before reading the node and element records.
        read(unit, *, iostat=io_status, iomsg=io_message) &
            node_count, element_count, nodes_per_element, particle_count
        if (io_status /= 0) then
            call fail(2, "Unable to read mesh header: " // trim(io_message))
            return
        end if

        ! Basic checks on the header.  These catch impossible meshes early and
        ! give a clearer error message than an array allocation or indexing
        ! failure later in the code.
        if (node_count <= 0) then
            call fail(3, "Mesh header has invalid node_count.")
            return
        end if

        if (element_count <= 0) then
            call fail(3, "Mesh header has invalid element_count.")
            return
        end if

        if (nodes_per_element /= linear_triangle_node_count .and. nodes_per_element /= quadratic_triangle_node_count) then
            call fail(3, "Mesh nodes_per_element must be 3 or 6.")
            return
        end if

        if (particle_count <= 0) then
            call fail(3, "Mesh header has invalid particle_count.")
            return
        end if

        ! Store the checked header values in the geometry object.
        geometry%particle_count = particle_count
        geometry%mesh%node_count = node_count
        geometry%mesh%element_count = element_count
        geometry%mesh%nodes_per_element = nodes_per_element

        ! Allocate arrays now that the mesh size is known.
        !
        ! xyz has 3 rows because every point has x, y, and z coordinates.
        ! element_nodes has nodes_per_element rows because each element has either
        ! 3 node ids or 6 node ids.
        allocate(geometry%particles(particle_count), &
                 geometry%mesh%xyz(3, node_count), &
                 geometry%mesh%element_nodes(nodes_per_element, element_count), &
                 geometry%mesh%node_particle_id(node_count), &
                 geometry%mesh%element_particle_id(element_count), &
                 stat=io_status)
        if (io_status /= 0) then
            call fail(4, "Unable to allocate mesh arrays.")
            return
        end if

        ! Initialise arrays to known values.  This helps avoid accidental use of
        ! uninitialised memory if a later part of the code has a bug.
        geometry%mesh%xyz = 0.0_dp
        geometry%mesh%element_nodes = 0
        geometry%mesh%node_particle_id = 0
        geometry%mesh%element_particle_id = 0

        ! Read node records.
        !
        ! Each node is one point on the surface mesh.  The particle id tells us
        ! which closed surface this node belongs to.
        do node_index = 1, node_count
            read(unit, *, iostat=io_status, iomsg=io_message) &
                node_id, x_coord, y_coord, z_coord, particle_id
            if (io_status /= 0) then
                call fail(5, "Unable to read node record: " // trim(io_message))
                return
            end if

            ! The public format requires node ids to be exactly 1, 2, 3, ...
            ! This keeps the reader simple because node_id can be used directly
            ! as the Fortran array index.
            if (node_id /= node_index) then
                call fail(5, "Node ids must be sequential and one-based.")
                return
            end if

            ! particle_id must point to one of the particles declared in the
            ! header.  For a single sphere, the valid value is usually only 1.
            if (particle_id < 1 .or. particle_id > particle_count) then
                call fail(5, "Node record has invalid particle_id.")
                return
            end if

            ! Store the node coordinates and remember which particle owns this
            ! node.  The slice xyz(:, node_id) means "all three coordinates for
            ! this one node".
            geometry%mesh%xyz(:, node_id) = [x_coord, y_coord, z_coord]
            geometry%mesh%node_particle_id(node_id) = particle_id

            ! Update the first/last node range for this particle.
            call update_range(geometry%particles(particle_id)%first_node, &
                              geometry%particles(particle_id)%last_node, node_id)
        end do

        ! Read element records.
        !
        ! Each element is a triangle.  A linear triangle stores 3 node ids.  A
        ! quadratic triangle stores 6 node ids in the exact order documented
        ! above: 3 oriented corners followed by midsides 1-2, 2-3, and 3-1.
        do element_index = 1, element_count
            local_nodes = 0

            ! Read the correct number of node ids depending on nodes_per_element.
            if (nodes_per_element == linear_triangle_node_count) then
                read(unit, *, iostat=io_status, iomsg=io_message) &
                    element_id, particle_id, local_nodes(1), local_nodes(2), local_nodes(3)
            else
                read(unit, *, iostat=io_status, iomsg=io_message) &
                    element_id, particle_id, local_nodes(1), local_nodes(2), &
                    local_nodes(3), local_nodes(4), local_nodes(5), local_nodes(6)
            end if

            if (io_status /= 0) then
                call fail(6, "Unable to read element record: " // trim(io_message))
                return
            end if

            ! As for node ids, element ids must be exactly 1, 2, 3, ...
            ! This allows element_id to be used directly as an array index.
            if (element_id /= element_index) then
                call fail(6, "Element ids must be sequential and one-based.")
                return
            end if

            ! Check that the element belongs to an existing particle.
            if (particle_id < 1 .or. particle_id > particle_count) then
                call fail(6, "Element record has invalid particle_id.")
                return
            end if

            ! Check that every node referenced by the element actually exists.
            ! This prevents later geometry calculations from reading outside the
            ! node coordinate array.
            if (any(local_nodes(1:nodes_per_element) < 1) .or. &
                any(local_nodes(1:nodes_per_element) > node_count)) then
                call fail(6, "Element record references an invalid node id.")
                return
            end if

            ! Store the element connectivity and particle ownership.
            geometry%mesh%element_nodes(:, element_id) = local_nodes(1:nodes_per_element)
            geometry%mesh%element_particle_id(element_id) = particle_id

            ! Update the first/last element range for this particle.
            call update_range(geometry%particles(particle_id)%first_element, &
                              geometry%particles(particle_id)%last_element, element_id)
        end do

        ! The file is no longer needed after all nodes and elements are read.
        close(unit)
        unit = -1

        call finish(0, "Mesh read successfully.")

    contains

        subroutine update_range(first_index, last_index, new_index)
            ! Expand an index range so that it includes new_index.
            !
            ! Example:
            !   first_index=5, last_index=8, new_index=3
            ! becomes:
            !   first_index=3, last_index=8
            !
            ! A value of zero means the range has not been set yet.
            integer, intent(inout) :: first_index
            integer, intent(inout) :: last_index
            integer, intent(in) :: new_index

            if (first_index == 0 .or. new_index < first_index) first_index = new_index
            if (last_index == 0 .or. new_index > last_index) last_index = new_index
        end subroutine update_range

        subroutine fail(error_code, text)
            ! Handle an error in one place.
            !
            ! If something goes wrong, close the file if it is open, clear any
            ! partly-read geometry data, and return the error code/message to
            ! the caller.
            integer, intent(in) :: error_code
            character(len=*), intent(in) :: text

            if (unit /= -1) close(unit)
            unit = -1
            call geometry%clear()
            call finish(error_code, text)
        end subroutine fail

        subroutine finish(error_code, text)
            ! Return status information to the caller.
            !
            ! status and message are optional arguments.  present(...) checks
            ! whether the caller supplied them before we try to write to them.
            integer, intent(in) :: error_code
            character(len=*), intent(in) :: text

            if (present(status)) status = error_code
            if (present(message)) message = text
        end subroutine finish

    end subroutine read_mesh

end module Geom_ReadMesh
