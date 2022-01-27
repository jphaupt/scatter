module integrate
    !! author: JP Haupt
    !! This module contains integration schemes (e.g. Numerov).
    use precision, only: rp
    use special_functions, only: sphericalJ, sphericalN
    use potential_type, only: Potential_t
    use constants, only: pi
    implicit none

    private
    public :: numerov, deltaPhase, numerov_thijssen
contains

    pure real(rp) function deltaPhase(angularL, Gquot, R2) result(delta)
        implicit none
        integer, intent(in) :: angularL
        real(rp), intent(in) :: Gquot, R2
            !!
        delta = -1._rp ! TODO stub
    end function

    pure subroutine numerov(pot, startPoint, startVal, nextPoint, nextVal, l_ang, energy, steps_tot, rsep_)
        implicit none
        class(Potential_t), intent(in) :: pot
        real(rp), intent(in) :: startPoint, startVal, nextPoint, nextVal, energy
        integer, intent(in) :: l_ang, steps_tot
        real(rp), intent(in), optional :: rsep_
        real(rp) :: rsep

        rsep = pi/sqrt(energy)/pot%alpha
        if (present(rsep_)) rsep = rsep_        
    end subroutine numerov

    pure subroutine numerov_old(pot, step, l, E, rmax, r1, r2, ulr1, ulr2, rsep_)
        ! use a debugger instead of print statements when debugging :)
        ! you can do this by running gdb on the executable pFUnit produces
        ! TODO below is notes on Thijssen's implementation
        ! initial values: phistart, phinext
        ! integration step: delta
        ! integrate from startI to endI
        ! solution is in array Solution (size maxSol)
        ! Sing is if there is a singularity at r=0 in the potential
        ! Phi''(r_i) = FArr(r_i)*Phi(r_i)
        ! they fill FArr first:
        ! FArr(r_i) = 2*(V(R)-E)+L*(L+1)/R**2
        ! (not sure why 2* -- I guess that is alpha for LJ pot)
        implicit none
        class(Potential_t), intent(in) :: pot
        real(rp), intent(in) :: step, E, rmax ! rmin?
        integer, intent(in) :: l
        real(rp), intent(out) :: r1, r2, ulr1, ulr2
        real(rp), intent(in), optional :: rsep_
        !! r2 - r1 (defaults to wavelength/2 = pi*sqrt(hbar^2/(2*m))/sqrt(E))
        !! right now I also have hbar^2/2m = 1 so this is just pi/sqrt(E)
        real(rp) :: rsep
        real(rp) :: step_sq, r, w, w_prev, sol, fval
        real(rp) :: sol_next ! solution(startI+Istep)

        ! NOTE these will be different depending on potential
        ! must still implement them there
        ! or maybe make them optional arguments
        real(rp) :: r_start, w_start, r_next, w_next
        real(rp) :: u_start, u, u_next, u_prev
        r_start = 0 ! todo do I even need this?
        r_next = step
        u_start = 0
        u_next = step**(l+1)

        rsep = pi/sqrt(E)/pot%alpha
        if (present(rsep_)) rsep = rsep_
        
        step_sq = step*step

        if (pot%hasSingularityAtR0) then
            w_prev = u_start
        else
            fval = pot%alpha * (pot%potential(r) - E) ! assuming l=0
            w_prev = getW(u_start, fval, step_sq)
            ! (1-step_sq/12.*fval)*u_start ! just 0 anyway?!
            sol = u_start
        end if

        u = u_next
        sol_next = u_next
        fval = radialRHS(pot, r_next, l, E)
        w = getW(u, fval, step_sq)

        do while (rmax-r > 1.e-8)
            w_next = 2. * w - w_prev + step_sq * u * fval
            w_prev = w
            w = w_next
            r = r + step
            fval = radialRHS(pot, r, l, E)
            u = getU(w, fval, step_sq)
            sol_next = u
            ! TODO
        end do

        r1 = r
        ulr1 = sol_next

        ! TODO !! ulr2
        r2 = -1
        ulr2 = -1

        ! w_prev = w_start
        ! w = w_next
        ! r = r_next
        ! fval = radialRHS(pot, r, l, E)
        ! ! print*, fval
        ! sol = (1 - step_sq/12*fval)*w
        ! do while (r < rmax)
        !     w_next = 2. * w - w_prev + step_sq * sol * fval
        !     ! print*, r, w_next
        !     r = r + step
        !     w_prev = w
        !     w = w_next
        !     fval = radialRHS(pot, r, l, E)
        !     sol = getU(w, fval, step_sq)
        ! end do
        ! r1 = r
        ! ulr1 = sol
        ! do while (r < rmax+rsep)
        !     w_next = 2. * w - w_prev + step_sq * sol * fval
        !     r = r + step
        !     w_prev = w
        !     w = w_next
        !     fval = radialRHS(pot, r, l, E)
        !     sol = getU(w, fval, step_sq)
        ! end do
        ! r2 = r
        ! ulr2 = sol
    end subroutine

    subroutine numerov_thijssen(Delta, StartI, EndI, MaxSol, FArr, &
    &                     Sing, PhiStart, PhiNext, Solution)
        ! copying thijssen's numerov subroutine because mine doesn't seem to work...
    IMPLICIT NONE

    INTEGER I, StartI, EndI, MaxSol, IStep

    real(rp) Phi, PhiStart, PhiNext, W, WNext, WPrev, &
  &        Delta, Fac, FArr(MaxSol), DeltaSq, Solution(MaxSol)

    LOGICAL Sing

    IF (Delta.LT.0) THEN
      IStep = -1
    ELSE
      IStep = 1
    ENDIF

    DeltaSq = Delta*Delta
    Fac = DeltaSq/12.D0

    IF (Sing) THEN
      WPrev = PhiStart
    ELSE
      WPrev = (1-Fac*FArr(StartI))*PhiStart
      Solution(StartI) = PhiStart
    ENDIF

    Phi = PhiNext
    Solution(StartI+IStep) = PhiNext
    W   = (1-Fac*FArr(StartI+IStep))*Phi

    DO I = StartI+IStep, EndI-IStep, IStep
       WNext = W*2.D0 - WPrev + DeltaSq*Phi*FArr(I)
       WPrev = W
       W     = WNext
       Phi   = W/(1-Fac*FArr(I+IStep))
       Solution(I+IStep) = Phi
    ENDDO   
    end subroutine numerov_thijssen

    pure elemental real(rp) function radialRHS(pot, r, l, energy) result(V_eff)
        implicit none
        !! returns the RHS of the radial Schroedinger equation, called F in
        !! equation 2.10 in Thijssen
        !! note this is actually 2*m/hbar^2 * F but at least for now the prefactor
        !! is unity
        class(Potential_t), intent(in) :: pot
        integer, intent(in) :: l
        real(rp), intent(in) :: r, energy

        ! TODO what to do if r=0??
        V_eff = pot%alpha * (pot%potential(r) - energy) + l * (l+1) / (r*r)

    end function radialRHS

    pure real(rp) function getU(w, fval, h2) result(U)
        implicit none
        real(rp), intent(in) :: w, fval, h2
        !! see equation 2.13
        U = w/(1._rp - h2/12. * fval)
    end function getU

    pure real(rp) function getW(u, fval, h2) result(W)
    implicit none
    real(rp), intent(in) :: u, fval, h2
    !! see equation 2.13
    W = u*(1._rp - h2/12. * fval)
    end function getW


end module integrate