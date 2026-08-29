module Pre_Constants
    ! This module stores numerical constants that are used by many other
    ! Fortran files in the solver.
    !
    ! Keeping these values in one place avoids small mistakes such as using
    ! different values of pi or different real-number precision in different
    ! parts of the code.
    implicit none

    ! dp means "double precision".
    !
    ! Computers can store real numbers with different levels of accuracy.
    ! The solver uses double precision because boundary element calculations
    ! involve many additions, multiplications, and divisions of small complex
    ! numbers.  Using one shared precision kind makes all real and complex
    ! variables consistent across the code.
    integer, parameter :: dp = kind(1.0d0)

    ! Mathematical constant pi.
    !
    ! atan(1) is pi/4, so multiplying it by 4 gives pi.  Calculating pi this
    ! way lets the compiler create a value that matches the selected precision.
    real(dp), parameter :: pi = atan(1.0_dp) * 4.0_dp

    ! Complex number 0 + 0i.
    !
    ! This is used when a complex variable, vector, or matrix needs to be
    ! initialised to zero before the solver starts adding contributions to it.
    complex(dp), parameter :: complex_zero = cmplx(0.0_dp, 0.0_dp, kind=dp)

    ! Complex number 0 + 1i, often written as i or sqrt(-1).
    !
    ! The frequency-domain Helmholtz equation uses complex numbers to describe
    ! wave phase.  For example, exp(i*k*r) uses this imaginary unit.
    complex(dp), parameter :: imaginary_unit = cmplx(0.0_dp, 1.0_dp, kind=dp)
end module Pre_Constants
