module integrate
    !! author: JP Haupt
    !! This module contains integration schemes (e.g. Numerov).
    use precision, only: rp
    use special_functions, only: sphericalJ, sphericalN
    use potential_type, only: Potential_t
    use constants, only: pi
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

    pure subroutine numerov(pot, dr, l, E, rmax, r1, r2, ulr1, ulr2, rsep_)
        ! use a debugger instead of print statements when debugging :)
        ! you can do this by running gdb on the executable pFUnit produces
        implicit none
        class(Potential_t), intent(in) :: pot
        real(rp), intent(in) :: dr, E, rmax ! rmin?
        integer, intent(in) :: l
        real(rp), intent(out) :: r1, r2, ulr1, ulr2
        real(rp), intent(in), optional :: rsep_
        !! r2 - r1 (defaults to wavelength/2 = pi*sqrt(hbar^2/(2*m))/sqrt(E))
        !! right now I also have hbar^2/2m = 1 so this is just pi/sqrt(E)
        real(rp) :: rsep
        real(rp) :: dr2, r, w, w_prev, sol, fval

        ! NOTE these will be different depending on potential
        ! must still implement them there
        ! or maybe make them optional arguments
        real(rp) :: r_start, w_start, r_next, w_next
        r_start = 0 ! todo do I even need this?
        r_next = dr
        w_start = 0
        w_next = dr**(l+1)

        rsep = pi/sqrt(E)/2._rp ! TODO lambda/4 rn
        if (present(rsep_)) rsep = rsep_

        dr2 = dr*dr
        w_prev = w_start
        w = w_next
        r = r_next
        fval = radialRHS(pot, r, l, E)
        ! print*, fval
        sol = (1 - dr2/12*fval)*w
        do while (r < rmax)
            w_next = 2. * w - w_prev + dr2 * sol * fval
            ! print*, r, w_next
            r = r + dr
            w_prev = w
            w = w_next
            fval = radialRHS(pot, r, l, E)
            sol = getU(w, fval, dr2)
        end do
        r1 = r
        ulr1 = sol
        do while (r < rmax+rsep)
            w_next = 2. * w - w_prev + dr2 * sol * fval
            r = r + dr
            w_prev = w
            w = w_next
            fval = radialRHS(pot, r, l, E)
            sol = getU(w, fval, dr2)
        end do
        r2 = r
        ulr2 = sol
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
        ! TODO what to do if r=0??
        V_eff = pot%potential(r) + l * (l+1) / (r*r) - energy

    end function radialRHS

    pure real(rp) function getU(w, fval, h2) result(U)
        implicit none
        real(rp), intent(in) :: w, fval, h2
        !! see equation 2.13
        U = w/(1._rp - h2/12. * fval)
    end function getU

end module integrate