module integrate
   !! author: JP Haupt
   !! This module contains integration schemes (e.g. Numerov).
   use precision, only: rp
   use special_functions, only: sphericalJ, sphericalN
   use potential_type, only: Potential_t
   use constants, only: pi
   implicit none

   private
   public :: numerov, numerov_thijssen
contains

    pure subroutine numerov(pot, startPoint, startVal, nextPoint, nextVal, l_ang, energy, steps_tot, rsep_, sols)
      implicit none
      class(Potential_t), intent(in) :: pot
      real(rp), intent(inout) :: startPoint, startVal, nextPoint, nextVal
      !! starting value of r, u(r) and the immediate successor. These are
      !! replaced by the solution values r1, u(r1), r2, u(r2)
      !! steps_tot also gets updated to be up to where the r2 is
      real(rp), intent(in) :: energy
      integer, intent(in) :: l_ang
      integer, intent(inout) :: steps_tot
      real(rp), intent(in), optional :: rsep_
      real(rp) :: rsep
      real(rp) :: hstep, hstep_sq !hstep2
      real(rp) :: wprev, wcurr, wnext, ucurr, unext
      real(rp) :: fval, position, new_steps_tot
    !   real(rp), allocatable :: sols
      real(rp), optional, allocatable, intent(out) :: sols(:)

      integer :: i

      if (present(sols)) allocate(sols(steps_tot+1))
      rsep = pi/sqrt(energy)/pot%alpha
      if (present(rsep_)) rsep = rsep_


      hstep = nextPoint - startPoint
      hstep_sq = hstep*hstep
      position = startPoint
      if (pot%hasSingularityAtR0) then
        wprev = startVal
      else
        if (abs(position)<=1.e-10_rp) then
            ! just ignore the l term... Not a great solution but hopefully it works
            fval = pot%alpha * (pot%potential(position) - energy) ! should this be 2*energy???
        else
            fval = radialRHS(pot, position, l_ang, energy)
        endif
        ! fval = radialRHS(pot, position, l_ang, energy)
        wprev = getW(startVal, fval, hstep_sq)
        if (allocated(sols)) sols(1) = startVal
      endif

      ucurr = nextVal
      if (allocated(sols)) sols(2) = nextVal
      position = position + hstep
      fval = radialRHS(pot, position, l_ang, energy)
      wcurr = getW(ucurr, fval, hstep_sq)

      do i=1, steps_tot-1
        ! fval = radialRHS(pot, position, l_ang, energy)
        wnext = wcurr*2._rp - wprev + hstep_sq*ucurr*fval
        wprev = wcurr
        wcurr = wnext
        position = position + hstep
        ! print*, position
        fval = radialRHS(pot, position, l_ang, energy)
        ucurr = getU(wcurr, fval, hstep_sq)
        if (allocated(sols)) sols(i+1) = ucurr
      enddo
      startVal = ucurr
      startPoint = position

      steps_tot = int(rsep/hstep)
      do i=1, steps_tot
        ! fval = radialRHS(pot, position, l_ang, energy)
        wnext = wcurr*2._rp - wprev + hstep_sq*ucurr*fval
        wprev = wcurr
        wcurr = wnext
        position = position + hstep
        ! print*, position
        fval = radialRHS(pot, position, l_ang, energy)
        ucurr = getU(wcurr, fval, hstep_sq)
        ! TODO sols
        ! if (allocated(sols)) sols(i+1) = ucurr
      enddo
      nextPoint = position
      nextVal = ucurr

    end subroutine numerov

   subroutine numerov_thijssen(Delta, StartI, EndI, MaxSol, FArr, &
   &                     Sing, PhiStart, PhiNext, Solution)
      ! used just for testing
      ! written by Jos Thijssen (author of _Computational Physics_), not my work
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

      V_eff = pot%alpha * (pot%potential(r) - energy) + 1_rp * l * (l+1._rp) / (r*r)

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
