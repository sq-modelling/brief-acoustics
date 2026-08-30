module AU_LayerTopology
    ! This module describes how acoustic-domain records are connected.
    !
    ! A boundary element problem can contain several closed surfaces.  Some
    ! surfaces may sit inside other surfaces.  For example:
    !
    !   particle 1: large outer shell
    !   particle 2: smaller object inside particle 1
    !
    ! The parent_particle_id stored in AU_Types tells us this nesting:
    !   parent_particle_id = 0 means the particle is in the infinite exterior
    !   parent_particle_id = p means the particle is inside particle p
    !
    ! The routines here answer questions needed during matrix assembly, such as:
    !   - which source particles belong to a given exterior domain
    !   - which source particles belong to a given interior domain
    !   - which density and wavenumber belong to that domain
    !   - what free-term constants should be added to the diagonal
    !
    ! The topology queries are retained so the ordinary research path can
    ! represent more than one surface.  The public v1.0 runner nevertheless
    ! accepts one exterior particle, and nested or multilayer Burton-Miller
    ! equations are not certified merely because these queries exist.
    use Pre_Constants, only: dp, pi, complex_zero, imaginary_unit
    use AU_Types, only: au_case_type
    implicit none

    ! Only expose small query functions.  This module stores no mutable state.
    private

    public :: particle_parent_id
    public :: particle_exterior_free_term_sign
    public :: exterior_domain_contains_particle
    public :: interior_domain_contains_particle
    public :: exterior_domain_has_infinity
    public :: exterior_domain_wavenumber
    public :: exterior_domain_density
    public :: interior_domain_wavenumber
    public :: interior_domain_density
    public :: density_ratio_parent_to_child
    public :: bem_external_h_free_term
    public :: burton_miller_free_term
    public :: burton_miller_coupling_length
    public :: burton_miller_coupling_weight

contains

    pure integer function particle_parent_id(case_data, particle_id)
        ! Return the parent particle id for one particle.
        !
        ! A return value of 0 means the particle is directly exposed to the
        ! infinite exterior medium.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: particle_id

        particle_parent_id = case_data%layer(particle_id)%parent_particle_id
    end function particle_parent_id

    pure integer function particle_exterior_free_term_sign(case_data, particle_id)
        ! Return the exterior free-term sign for one particle.
        !
        ! The mesh has a stored normal direction.  The boundary-integral free
        ! terms depend on whether that stored normal agrees with the domain
        ! normal used in the equation.  This sign carries that convention.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: particle_id

        particle_exterior_free_term_sign = case_data%layer(particle_id)%exterior_free_term_sign
    end function particle_exterior_free_term_sign

    pure logical function exterior_domain_contains_particle(case_data, host_particle_id, source_particle_id)
        ! Decide whether source_particle_id belongs to the exterior domain of
        ! host_particle_id.
        !
        ! "Exterior domain of host" means the acoustic region immediately
        ! outside the host surface.  If the host is in open space, this is the
        ! infinite exterior fluid.  If the host is inside another particle, this
        ! is the bounded region inside its parent but outside the host.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: host_particle_id
        integer, intent(in) :: source_particle_id

        integer :: host_parent
        integer :: source_parent

        host_parent = particle_parent_id(case_data, host_particle_id)
        source_parent = particle_parent_id(case_data, source_particle_id)

        ! In plain language:
        !   - siblings share the same parent and therefore sit in the same
        !     surrounding acoustic region
        !   - the parent surface itself bounds that region from outside
        exterior_domain_contains_particle = source_parent == host_parent .or. &
                                            source_particle_id == host_parent
    end function exterior_domain_contains_particle

    pure logical function interior_domain_contains_particle(case_data, host_particle_id, source_particle_id)
        ! Decide whether source_particle_id belongs to the interior domain of
        ! host_particle_id.
        !
        ! "Interior domain of host" means the acoustic region immediately inside
        ! the host surface.  It contains the host surface itself and any
        ! immediate child surfaces inside the host.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: host_particle_id
        integer, intent(in) :: source_particle_id

        ! Include the host surface and each immediate child surface that bounds
        ! the same interior region.
        interior_domain_contains_particle = source_particle_id == host_particle_id .or. &
                                            particle_parent_id(case_data, source_particle_id) == host_particle_id
    end function interior_domain_contains_particle

    pure logical function exterior_domain_has_infinity(case_data, particle_id)
        ! Return true when the exterior side of this particle is the infinite
        ! exterior acoustic medium.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: particle_id

        exterior_domain_has_infinity = particle_parent_id(case_data, particle_id) == 0
    end function exterior_domain_has_infinity

    pure complex(dp) function exterior_domain_wavenumber(case_data, particle_id)
        ! Return the wavenumber for the acoustic domain outside one particle.
        !
        ! If the particle has no parent, the outside is the global exterior
        ! medium.  If the particle has a parent, the outside is the medium inside
        ! that parent particle.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: particle_id

        integer :: parent_id

        parent_id = particle_parent_id(case_data, particle_id)
        if (parent_id == 0) then
            exterior_domain_wavenumber = case_data%exterior_medium%wavenumber
        else
            exterior_domain_wavenumber = case_data%interior_medium(parent_id)%wavenumber
        end if
    end function exterior_domain_wavenumber

    pure real(dp) function exterior_domain_density(case_data, particle_id)
        ! Return the density for the acoustic domain outside one particle.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: particle_id

        integer :: parent_id

        parent_id = particle_parent_id(case_data, particle_id)
        if (parent_id == 0) then
            exterior_domain_density = case_data%exterior_medium%density
        else
            exterior_domain_density = case_data%interior_medium(parent_id)%density
        end if
    end function exterior_domain_density

    pure complex(dp) function interior_domain_wavenumber(case_data, particle_id)
        ! Return the wavenumber for the acoustic domain inside one particle.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: particle_id

        interior_domain_wavenumber = case_data%interior_medium(particle_id)%wavenumber
    end function interior_domain_wavenumber

    pure real(dp) function interior_domain_density(case_data, particle_id)
        ! Return the density for the acoustic domain inside one particle.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: particle_id

        interior_domain_density = case_data%interior_medium(particle_id)%density
    end function interior_domain_density

    pure real(dp) function density_ratio_parent_to_child(case_data, child_particle_id)
        ! Return density(outside child) / density(inside child).
        !
        ! Transmission conditions often contain density ratios because pressure
        ! and normal velocity use different material parameters on each side of
        ! an interface.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: child_particle_id

        density_ratio_parent_to_child = exterior_domain_density(case_data, child_particle_id) / &
                                        interior_domain_density(case_data, child_particle_id)
    end function density_ratio_parent_to_child

    pure complex(dp) function bem_external_h_free_term(case_data, particle_id)
        ! Free-term contribution for the ordinary exterior BEM equation.
        !
        ! With the Green-function convention used in this code, the full solid
        ! angle coefficient for a smooth closed exterior boundary is represented
        ! as 4*pi times the exterior free-term sign. For bounded exterior domains, this
        ! ordinary exterior free term is not added here.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: particle_id

        if (exterior_domain_has_infinity(case_data, particle_id)) then
            bem_external_h_free_term = cmplx(4.0_dp * pi * particle_exterior_free_term_sign(case_data, particle_id), &
                                            0.0_dp, kind=dp)
        else
            bem_external_h_free_term = complex_zero
        end if
    end function bem_external_h_free_term

    pure complex(dp) function burton_miller_free_term(case_data, particle_id)
        ! Free-term contribution for Burton-Miller auxiliary equations.
        !
        ! The sign and 4*pi scale are the stored-mesh-normal form of the
        ! published equation.  The full derivation and diagonal placement are
        ! recorded in docs/burton-miller-equation-map.md.  The release path
        ! calls this only after enforcing its single-particle exterior scope.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: particle_id

        burton_miller_free_term = cmplx(4.0_dp * pi * particle_exterior_free_term_sign(case_data, particle_id), &
                                       0.0_dp, kind=dp)
    end function burton_miller_free_term

    pure real(dp) function burton_miller_coupling_length(case_data, particle_id)
        ! Return the dimensional length beta used to couple the two equations.
        !
        ! The ordinary boundary integral equation and its normal derivative do
        ! not have the same physical dimensions: normal differentiation adds
        ! one factor of inverse length.  Multiplying the derivative equation by
        ! beta, which has units of length, makes their sum dimensionally
        ! homogeneous.  beta is therefore not an arbitrary tuning number.
        !
        ! Let a be a representative body length and k the positive acoustic
        ! wavenumber.  This implementation uses
        !
        !   beta = min(a, 1/k).
        !
        ! Consequently beta=a for ka<=1 and beta=1/k for ka>1.  The first
        ! branch prevents 1/k from becoming unbounded at low frequency; the
        ! second follows the wavelength-related scale used in the published
        ! high-frequency calculations.  For a sphere, a is normally its radius.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: particle_id

        real(dp) :: body_length
        real(dp) :: wavenumber_magnitude

        body_length = case_data%layer(particle_id)%characteristic_length
        wavenumber_magnitude = abs(exterior_domain_wavenumber(case_data, particle_id))

        burton_miller_coupling_length = body_length
        if (wavenumber_magnitude > tiny(1.0_dp)) then
            burton_miller_coupling_length = min(body_length, 1.0_dp / wavenumber_magnitude)
        end if
    end function burton_miller_coupling_length

    pure complex(dp) function burton_miller_coupling_weight(case_data, particle_id)
        ! Return the Burton-Miller coupling weight for one particle.
        !
        ! Burton-Miller combines the ordinary boundary integral equation with a
        ! normal-derivative equation.  The coupling weight must have dimensions
        ! of length so the two equations are scaled compatibly.
        !
        ! burton_miller_coupling_length supplies the dimensional beta.  The
        ! complex code-space weight also contains the orientation of the stored
        ! mesh normal:
        !
        !   weight = i * exterior_free_term_sign * beta.
        !
        ! The standard mesh has outward-from-solid normals and exterior_free_term_sign=-1,
        ! so this becomes -i*beta.  AU_Solver rejects zero, negative, or complex
        ! wavenumbers before this function is used by the public release path.
        type(au_case_type), intent(in) :: case_data
        integer, intent(in) :: particle_id

        real(dp) :: beta

        beta = burton_miller_coupling_length(case_data, particle_id)
        burton_miller_coupling_weight = imaginary_unit * &
            real(particle_exterior_free_term_sign(case_data, particle_id), dp) * beta
    end function burton_miller_coupling_weight

end module AU_LayerTopology
