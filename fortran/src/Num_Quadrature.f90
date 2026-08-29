module Num_Quadrature
    ! This module stores fixed numerical quadrature rules.
    !
    ! Quadrature means numerical integration: instead of integrating a function
    ! exactly by hand, the code evaluates the function at selected points and
    ! multiplies those values by selected weights.
    !
    ! Boundary element methods need many surface integrals over triangular mesh
    ! elements.  A good triangle quadrature rule is therefore important for both
    ! accuracy and speed.
    use Pre_Constants, only: dp
    implicit none

    ! Keep the quadrature table itself private unless it is explicitly listed
    ! below.  The current release exposes one 25-point triangle rule.
    private

    public :: triangle_quadrature_25

contains
    subroutine triangle_quadrature_25(xi, eta, weight)
        ! Return a fixed 25-point quadrature rule for the unit right triangle.
        !
        ! Reference triangle:
        !   0 <= xi
        !   0 <= eta
        !   xi + eta <= 1
        !
        ! The third barycentric coordinate is:
        !   zeta = 1 - xi - eta
        !
        ! The quadrature approximation is:
        !   integral_over_triangle f(xi, eta) dA
        !       ~= sum_q weight(q) * f(xi(q), eta(q))
        !
        ! The unit right triangle has area 1/2, so the weights below should sum
        ! to 0.5, not 1.0.
        !
        ! Source of the rule:
        !   F. D. Witherden and P. E. Vincent,
        !   "On the identification of symmetric quadrature rules for finite
        !   element methods",
        !   Computers & Mathematics with Applications, 69 (2015), 1232-1241.
        !   DOI: 10.1016/j.camwa.2015.03.017
        !
        ! This table is the symmetric 25-point triangle rule with degree of
        ! precision 10.  In practical terms, it integrates every polynomial
        ! xi^a * eta^b exactly, up to roundoff error, when a + b <= 10.
        !
        ! Verification performed while documenting this file:
        !   - all points are inside the reference triangle
        !   - sum(weight) = 0.5 to roundoff
        !   - monomial moments are exact through total degree 10 to roundoff
        !
        ! Important: do not renormalise these weights to sum to 1.0.  The later
        ! surface integration code multiplies these reference-triangle weights
        ! by the element Jacobian.  Changing the scale here would scale every
        ! assembled boundary-integral matrix entry.
        real(dp), intent(out) ::  xi(25), &
                                & eta(25), &
                                & weight(25)

        ! Quadrature weights.  Each weight says how much the function value at
        ! the corresponding point contributes to the integral.
        weight = (/ &
        &        0.04087166457314298321405934249209d0, &
        &        0.00667648440657478313778648919953d0, &
        &        0.00667648440657478313778648919953d0, &
        &        0.00667648440657478313778648919953d0, &
        &        0.02297898180237236400689395481627d0, &
        &        0.02297898180237236400689395481627d0, &
        &        0.02297898180237236400689395481627d0, &
        &        0.03195245319821202271644936688133d0, &
        &        0.03195245319821202271644936688133d0, &
        &        0.03195245319821202271644936688133d0, &
        &        0.03195245319821202271644936688133d0, &
        &        0.03195245319821202271644936688133d0, &
        &        0.03195245319821202271644936688133d0, &
        &        0.01709232408147971431434579202067d0, &
        &        0.01709232408147971431434579202067d0, &
        &        0.01709232408147971431434579202067d0, &
        &        0.01709232408147971431434579202067d0, &
        &        0.01709232408147971431434579202067d0, &
        &        0.01709232408147971431434579202067d0, &
        &        0.01264887885364419219452139534142d0, &
        &        0.01264887885364419219452139534142d0, &
        &        0.01264887885364419219452139534142d0, &
        &        0.01264887885364419219452139534142d0, &
        &        0.01264887885364419219452139534142d0, &
        &        0.01264887885364419219452139534142d0 /)

        ! First coordinate of each quadrature point in the reference triangle.
        xi = (/ &
        &        0.33333333333333333333333333333333d0, &
        &        0.03205537321694351293098458933649d0, &
        &        0.93588925356611297413803082132702d0, &
        &        0.03205537321694351293098458933649d0, &
        &        0.14216110105656438509216210319096d0, &
        &        0.71567779788687122981567579361808d0, &
        &        0.14216110105656438509216210319096d0, &
        &        0.32181299528883542122509756098605d0, &
        &        0.53005411892734402827709567394569d0, &
        &        0.14813288578382055049780676506826d0, &
        &        0.53005411892734402827709567394569d0, &
        &        0.14813288578382055049780676506826d0, &
        &        0.32181299528883542122509756098605d0, &
        &        0.02961988948872976763383626942604d0, &
        &        0.60123332868345924545474289345869d0, &
        &        0.36914678182781098691142083711527d0, &
        &        0.60123332868345924545474289345869d0, &
        &        0.36914678182781098691142083711527d0, &
        &        0.02961988948872976763383626942604d0, &
        &        0.02836766533993843925043575557813d0, &
        &        0.80793060092287906507994990288174d0, &
        &        0.16370173373718249566961434154013d0, &
        &        0.80793060092287906507994990288174d0, &
        &        0.16370173373718249566961434154013d0, &
        &        0.02836766533993843925043575557813d0 /)

        ! Second coordinate of each quadrature point in the reference triangle.
        ! Together, xi(q), eta(q), and 1 - xi(q) - eta(q) are the three
        ! barycentric coordinates of point q.
        eta = (/ &
        &        0.33333333333333333333333333333333d0, &
        &        0.93588925356611297413803082132702d0, &
        &        0.03205537321694351293098458933649d0, &
        &        0.03205537321694351293098458933649d0, &
        &        0.71567779788687122981567579361808d0, &
        &        0.14216110105656438509216210319096d0, &
        &        0.14216110105656438509216210319096d0, &
        &        0.53005411892734402827709567394569d0, &
        &        0.32181299528883542122509756098605d0, &
        &        0.53005411892734402827709567394569d0, &
        &        0.14813288578382055049780676506826d0, &
        &        0.32181299528883542122509756098605d0, &
        &        0.14813288578382055049780676506826d0, &
        &        0.60123332868345924545474289345869d0, &
        &        0.02961988948872976763383626942604d0, &
        &        0.60123332868345924545474289345869d0, &
        &        0.36914678182781098691142083711527d0, &
        &        0.02961988948872976763383626942604d0, &
        &        0.36914678182781098691142083711527d0, &
        &        0.80793060092287906507994990288174d0, &
        &        0.02836766533993843925043575557813d0, &
        &        0.80793060092287906507994990288174d0, &
        &        0.16370173373718249566961434154013d0, &
        &        0.02836766533993843925043575557813d0, &
        &        0.16370173373718249566961434154013d0 /)

    end subroutine triangle_quadrature_25
end module Num_Quadrature
