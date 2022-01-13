module integrate
    !! author: JP Haupt
    !! This module contains integration schemes (e.g. Numerov).
    use precision, only: rp
    use special_functions, only: sphericalJ, sphericalN
    use potential_type, only: Potential_t
    implicit none

    private
    public :: numerov, deltaPhase
contains

    pure real(rp) function deltaPhase(angularL, Gquot, R2) result(delta)
        implicit none
        integer, intent(in) :: angularL
        real(rp), intent(in) :: Gquot, R2
            !!
        delta = -1._rp ! TODO stub
    end function

    pure subroutine numerov(pot, dr, l, E, rmax, r1, r2, ulr1, ulr2)
        implicit none
        class(Potential_t), intent(in) :: pot
        real(rp), intent(in) :: dr, E, rmax ! rmin?
        integer, intent(in) :: l
        real(rp), intent(out) :: r1, r2, ulr1, ulr2

        ! NOTE these will be different depending on potential
        ! must still implement them there
        real(rp) :: r_start, u_start, r_next, u_next
        r_start = 0
        r_next = dr
        u_start = 0
        u_next = dr**(l+1)
        r1 = -1 ! TODO stub
        r2 = -1
        ulr1 = -1
        ulr2 = -1
    end subroutine

    pure real(rp) function radialRHS(pot, r, l, energy) result(V_eff)
        implicit none
        !! returns the RHS of the radial Schroedinger equation, called F in
        !! equation 2.10 in Thijssen
        !! note this is actually 2*m/hbar^2 * F but at least for now the prefactor
        !! is unity
        class(Potential_t), intent(in) :: pot
        integer, intent(in) :: l
        real(rp), intent(in) :: r, energy

        ! TODO ?? multiply pot and energy by 2*m/hbar^2
        V_eff = pot%potential(r) + l * (l+1) / (r*r) - energy

    end function radialRHS

    pure real(rp) function getU(w, fval, h2) result(U)
        implicit none
        real(rp), intent(in) :: w, fval, h2
        !! see equation 2.13
        U = w/(1._rp - h2/12. * fval)
    end function getU

end module integrate