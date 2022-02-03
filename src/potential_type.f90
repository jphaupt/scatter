module potential_type
    use precision, only: rp
    use constants, only: pi
    use iso_fortran_env, only: stdout => output_unit
    use special_functions, only: sphJ => sphericalJ, sphN => sphericalN
    implicit none

    private
    public :: Potential_t, Harmonic_Potential_t

    ! type, abstract :: Potential_t ! TODO abstract type and extend to LJ and others
    type, abstract :: Potential_t
    !! object containing information about the scattering potential, such as
    !! the reduced mass, and potential parameters (epsilon, r0)

    ! defaults to H-Kr collision
    real(rp) :: maxE=3.5, minE=0.1, &
            &   start_r=0.75, end_r=5.0, rstep=0.02
    integer :: lmax=6, numE=200
    ! real(rp) :: rydConst ! = 0.48*mu*mu, not initialised in constructor

    ! real(rp) :: some_val
    real(rp), allocatable :: data(:,:)
    real(rp) :: alpha=1 ! 2m/hbar^2, called alpha in Thijssen
    logical :: hasSingularityAtR0 = .false.

    contains

        procedure(pot_int), pass(this), deferred :: potential
        procedure, pass(this) :: init_potential ! helper for extended constructors
        procedure, pass(this) :: partial_cross_section

    end type Potential_t

    ! interface Potential_t
    !     module procedure :: potential_constructor
    ! end interface Potential_t

    abstract interface
        ! ??? I guess there is no way to get an abstract constructor?
        ! module procedure :: potential_constructor

        pure real(rp) elemental function pot_int(this, r) result(retval)
            import rp
            import Potential_t
            class(Potential_t), intent(in) :: this
            real(rp), intent(in) :: r
        end function pot_int

        pure real(rp) function delta_int(this, l, energy, r1, r2, u1, u2) result(delta)
            import rp, Potential_t
            class(Potential_t), intent(in) :: this
            integer, intent(in) :: l
            real(rp), intent(in) :: energy, r1, r2, u1, u2
        end function delta_int

    end interface

    type, extends(Potential_t) :: Harmonic_Potential_t

    contains
    procedure, pass(this) :: potential => potential_harmonic
    !! TODO introduce deferred stuff

    end type
contains

pure subroutine init_potential(this, maxE, minE, numE, lmax, start_r, end_r, rstep)
    !! generic constructor used inside extended constructors
    !! (this is simply to avoid repetition)
    implicit none
    class(Potential_t), intent(inout) :: this
    real(rp), intent(in), optional :: maxE, minE, start_r, end_r, rstep
    integer, intent(in), optional :: lmax, numE

    ! ??? is there some fast way to do this? Maybe with a preprocessor?
    ! if (present(mu)) this%mu=mu
    if (present(maxE)) this%maxE=maxE
    if (present(minE)) this%minE=minE
    if (present(numE)) this%numE=numE
    if (present(lmax)) this%lmax=lmax
    if (present(start_r)) this%start_r=start_r
    if (present(end_r)) this%end_r=end_r
    if (present(rstep)) this%rstep=rstep
    ! this%rydConst = 0.48*this%mu * this%mu

end subroutine init_potential

pure elemental real(rp) function potential_harmonic(this, r) result(retval)
    implicit none
    class(Harmonic_Potential_t), intent(in) :: this
    real(rp), intent(in) :: r
    retval = this%alpha * r*r
end function potential_harmonic

pure real(rp) function calc_delta(this, l, energy, r1, r2, u1, u2) result(delta)
    implicit none
    !! Thijssen eq 2.9
    class(Potential_t), intent(in) :: this
    integer, intent(in) :: l
    real(rp), intent(in) :: energy, r1, r2, u1, u2
    real(rp) :: K, wavevec

    wavevec = sqrt(this%alpha * energy)

    K = (r1*u2)/(r2*u1)
    delta = (K*sphJ(l, r1*wavevec) - sphJ(l, wavevec*r2))
    delta = delta/(K*sphN(l, wavevec*r1) - sphN(l, wavevec*r2))
    delta = atan(delta)

end function calc_delta

pure real(rp) function partial_cross_section(this, l, energy, r1, r2, u1, u2) result(sigma_l)
!     !! @note this function may be deprecated once we do loss rates,
!     !! it assumes integration over the full domain, which is not true
!     !! for the application of interest
!     !! @endnote
    implicit none
    class(Potential_t), intent(in) :: this
    integer, intent(in) :: l
    real(rp), intent(in) :: energy, r1, r2, u1, u2
    real(rp) :: delta, sdl, k2
    delta = calc_delta(this, l, energy, r1, r2, u1, u2)
    sdl = sin(delta)
    k2 = this%alpha * energy ! 2mE/hbar^2
    sigma_l = 4*pi*sdl*sdl/k2

end function partial_cross_section

end module potential_type

! ! this subroutine cannot be pure because it may also produce a file
! subroutine solve(this, verbose_, tofile_, storeData_)
!     class(Potential_t), intent(inout) :: this ! inout necessary only if storeE
!     logical, intent(in), optional :: verbose_, tofile_, storeData_

!     logical :: verbose=.true., tofile=.true., storeData=.false.

!     integer :: fileunit
!     real(rp) :: dE, E_curr, denom, sig=0., delta
!     integer :: i, l

!     if (present(verbose_)) verbose=verbose_
!     if (present(tofile_)) tofile=tofile_
!     if (present(storeData_)) storeData=storeData_

!     if (storeData) allocate(this%data(this%numE, 2))

!     dE = (this%maxE - this%minE)/this%numE
!     E_curr = this%minE

!     if (tofile) then
!         open(newunit=fileunit, file='sigma.dat')
!         write(fileunit, *) '# energy sigma_tot'
!     end if

!     if (verbose) write(stdout, *) 'iteration energy sigma tan(delta)'

!     do i=1, this%numE ! +1 ! should I include a +1 here?
!         denom = sqrt(this%rydConst * E_curr)! TODO
!         do l=0, this%lmax
!             delta = this%calc_delta(l, E_curr, denom)
!             sig = sig + 4*pi/(denom*denom)*(2*l+1)*sin(delta)**2
!         end do
!         E_curr = E_curr + dE
!         if (tofile) write(fileunit, *) E_curr, sig
!         if (verbose) write(stdout, *) i, E_curr, sig, tan(delta)
!     end do


!     print*,this%end_r ! TODO STUB

!     if (tofile) then
!         if (verbose) write(stdout, *) 'cross section written to sigma.dat'
!         close(fileunit)
!     end if

! end subroutine solve
