module potential_type
    !! @todo separate different type to different files
    use precision, only: rp
    use constants, only: pi, hbar_sq
    use iso_fortran_env, only: stdout => output_unit
    use special_functions, only: sphJ => sphericalJ, sphN => sphericalN
    implicit none

    private
    public :: Potential_t, Potential_Harmonic_t, Potential_LennardJones_t

    ! type, abstract :: Potential_t ! TODO abstract type and extend to LJ and others
    type, abstract :: Potential_t
    !! object containing information about the scattering potential, such as
    !! the reduced mass, and potential parameters (epsilon, r0)

        ! @todo I think these should be moved out of the potential object, maybe to
        !       an Option_t object or something similar, since they don't strictly
        !       have anything to do with the potential
        real(rp) :: maxE = 3.5, minE = 0.1, &
                &   start_r = 0.75, end_r = 5.0, rstep = 0.02
        integer :: lmax = 6, numE = 200
        ! real(rp) :: rydConst ! = 0.48*mu*mu, not initialised in constructor

        ! real(rp) :: some_val
        real(rp), allocatable :: data(:, :)
        real(rp) :: alpha = 1 ! 2m/hbar^2, called alpha in Thijssen
        logical :: hasSingularityAtR0 = .false.

    contains

        procedure(pot_int), pass(this), deferred :: potential
        procedure(small_r_solution_int), pass(this), deferred :: small_r_solution
        procedure(small_r_derivative_int), pass(this), deferred :: small_r_derivative
        procedure, pass(this) :: init_potential ! helper for extended constructors
        procedure, pass(this) :: partial_cross_section
        procedure, pass(this) :: calc_delta
        procedure, pass(this) :: get_wavelength

    end type Potential_t

    ! interface Potential_t
    !     module procedure :: potential_constructor
    ! end interface Potential_t

    abstract interface
        ! ??? I guess there is no way to get an abstract constructor?
        ! module procedure :: potential_constructor

        pure elemental real(rp) function small_r_solution_int(this, r) result(retval)
            import rp
            import Potential_t
    !! solution to the spherical Schrodinger equation with a Lennard-Jones
    !! potential for small r (see eq 2.17 of Thijssen)
            class(Potential_t), intent(in) :: this
            real(rp), intent(in) :: r
            real(rp) :: const
        end function small_r_solution_int

        pure elemental real(rp) function small_r_derivative_int(this, r) result(retval)
            import rp
            import Potential_t
            class(Potential_t), intent(in) :: this
            real(rp), intent(in) :: r
            real(rp) :: const
        end function small_r_derivative_int

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

    type, extends(Potential_t) :: Potential_Harmonic_t
    !! only really used for testing

    contains
        procedure, pass(this) :: potential => potential_harmonic
        procedure, pass(this) :: small_r_solution => small_r_solution_Harmonic
        procedure, pass(this) :: small_r_derivative => small_r_derivative_Harmonic
        ! procedure, pass(this) :: constructor_harmonic

    end type Potential_Harmonic_t

    interface Potential_Harmonic_t
        module procedure :: constructor_harmonic
    end interface Potential_Harmonic_t

    ! TODO Lennard Jones
    type, extends(Potential_t) :: Potential_LennardJones_t
        ! TODO also include stuff like mu, etc.
        ! maybe make the contructor simply use a json file??
        ! TODO move to own file since it's important enough
        character(len=:), allocatable :: filename
        ! default to H-Kr collision
        real(rp) :: rho = 3.57
        !! rho in angstrom
        real(rp) :: epsilon = 5.9
        !! epsilon in meV
        real(rp) :: m1 = 1.008
        !! mass of first atom in Daltons
        real(rp) :: m2 = 83.798
        !! mass of second atom in Daltons
        ! real(rp) :: alpha=6.1
        ! @note you have to call the constructor to initialise alpha!
        ! (2*(1 Dalton)*(1.008*83.798)/(83.798+1.008)/hbar^2) *meV * (3.57 angstrom)^2 = 6.07

        ! ((1 Dalton)/hbar^2) *meV * angstrom^2 ~ 0.239225334
        ! I think what Thijssen does is take (2*(1 Dalton)/hbar^2) *meV ~ 0.48 A^-2
        ! then multiplies by rho^2. This is a good approximation but I would prefer
        ! to be more exact; i.e. include reduced mass and don't truncate
        ! I will also do the calculation in units of epsilon and rho
        ! I'm not completely decided if that's a good idea

    contains
        procedure, pass(this) :: potential => potential_LennardJones
        procedure, pass(this) :: small_r_solution => small_r_solution_LennardJones
        procedure, pass(this) :: small_r_derivative => small_r_derivative_LennardJones
        ! procedure, pass(this) :: init_potential => init_potential_LennardJones
        ! TODO constructor, potential function, define alpha, ...
        ! TODO there are some functions and values I must still implement...

    end type Potential_LennardJones_t
    interface Potential_LennardJones_t
        module procedure :: constructor_LennardJones
    end interface Potential_LennardJones_t
contains

    pure elemental real(rp) function small_r_solution_LennardJones(this, r) result(retval)
    !! solution to the spherical Schrodinger equation with a Lennard-Jones
    !! potential for small r (see eq 2.17 of Thijssen)
        class(Potential_LennardJones_t), intent(in) :: this
        real(rp), intent(in) :: r
        real(rp) :: const
        const = sqrt(this%alpha*this%epsilon/25)
        retVal = exp(-const*r**(-5))
    end function small_r_solution_LennardJones

    pure elemental real(rp) function small_r_derivative_LennardJones(this, r) result(retval)
    !! derivative at small r
        class(Potential_LennardJones_t), intent(in) :: this
        real(rp), intent(in) :: r
        real(rp) :: const
        const = sqrt(this%alpha*this%epsilon/25)
        retVal = 5*const*r**(-6)*this%small_r_solution(r)

    end function small_r_derivative_LennardJones

    pure elemental real(rp) function potential_LennardJones(this, r) result(retval)
    !! Potential function for the Lennard-Jones potential
    !! $$
    !! V(r) = \epsilon(\frac{\rho^{12}}{r^{12}} - 2\frac{\rho^6}{r^6})
    !! $$
        implicit none
        class(Potential_LennardJones_t), intent(in) :: this
        real(rp), intent(in) :: r
        real(rp) :: rho_r
        rho_r = this%rho/r
        retVal = this%epsilon*(rho_r**12 - 2*rho_r**6)

    end function potential_LennardJones

    pure type(Potential_LennardJones_t) function constructor_LennardJones(m1, m2, rho, epsilon) result(this)
    !! constructor for the harmonic potential type
    !! all it does is call the function init_potential
        implicit none
        real(rp), intent(in), optional :: m1, m2, rho, epsilon
        if (present(m1)) this%m1 = m1
        if (present(m2)) this%m2 = m2
        if (present(rho)) this%rho = rho
        if (present(epsilon)) this%epsilon = epsilon
        this%alpha = 2*(this%m1*this%m2)/(this%m1 + this%m2)/hbar_sq
        this%alpha = this%alpha*this%rho*this%rho ! units of rho^-2 meV
    end function constructor_LennardJones

    pure elemental real(rp) function small_r_solution_Harmonic(this, r) result(retval)
    !! solution to the spherical Schrodinger equation with a harmonic
    !! potential for small r (not implemented)
        class(Potential_Harmonic_t), intent(in) :: this
        real(rp), intent(in) :: r
        real(rp) :: const
        retVal = -1 ! TODO
    end function small_r_solution_Harmonic

    pure elemental real(rp) function small_r_derivative_Harmonic(this, r) result(retval)
    !! derivative at small r
        class(Potential_Harmonic_t), intent(in) :: this
        real(rp), intent(in) :: r
        real(rp) :: const
        ! const = sqrt(this%alpha * this%epsilon / 25)
        retVal = -1 ! TODO

    end function small_r_derivative_Harmonic

    pure type(Potential_Harmonic_t) function constructor_harmonic(maxE, minE, numE, lmax, start_r, end_r, rstep) result(this)
    !! constructor for the harmonic potential type
    !! all it does is call the function init_potential
        implicit none
        real(rp), intent(in), optional :: maxE, minE, start_r, end_r, rstep
        integer, intent(in), optional :: lmax, numE
        call this%init_potential(maxE, minE, numE, lmax, start_r, end_r, rstep)
    end function constructor_harmonic

    pure subroutine init_potential(this, maxE, minE, numE, lmax, start_r, end_r, rstep)
    !! generic constructor used inside extended constructors
    !! (this is simply to avoid repetition)
        implicit none
        class(Potential_t), intent(inout) :: this
        real(rp), intent(in), optional :: maxE, minE, start_r, end_r, rstep
        integer, intent(in), optional :: lmax, numE

        ! ??? is there some fast way to do this? Maybe with a preprocessor?
        ! if (present(mu)) this%mu=mu
        if (present(maxE)) this%maxE = maxE
        if (present(minE)) this%minE = minE
        if (present(numE)) this%numE = numE
        if (present(lmax)) this%lmax = lmax
        if (present(start_r)) this%start_r = start_r
        if (present(end_r)) this%end_r = end_r
        if (present(rstep)) this%rstep = rstep

    end subroutine init_potential

    pure elemental real(rp) function potential_harmonic(this, r) result(retval)
        implicit none
        class(Potential_Harmonic_t), intent(in) :: this
        real(rp), intent(in) :: r
        retval = this%alpha*r*r
    end function potential_harmonic

    pure real(rp) function calc_delta(this, l, energy, r1, r2, u1, u2) result(delta)
        implicit none
    !! Thijssen eq 2.9
    !! @note I am using a matching procedure, but this is not a great
    !! approximation. Later it would be good to use the log-derivative method
    !! @endnote
        class(Potential_t), intent(in) :: this
        integer, intent(in) :: l
        real(rp), intent(in) :: energy, r1, r2, u1, u2
        real(rp) :: K, wavevec

        wavevec = sqrt(this%alpha*energy)

        K = (r1*u2)/(r2*u1)
        delta = (K*sphJ(l, r1*wavevec) - sphJ(l, wavevec*r2))
        delta = delta/(K*sphN(l, wavevec*r1) - sphN(l, wavevec*r2))
        delta = atan(delta)

    end function calc_delta

    pure real(rp) function partial_cross_section(this, l, energy, r1, r2, u1, u2) result(sigma_l)
    !! @note this function may be deprecated once we do loss rates,
    !! it assumes integration over the full domain, which is not true
    !! for the application of interest
    !! @endnote
        implicit none
        class(Potential_t), intent(in) :: this
        integer, intent(in) :: l
        real(rp), intent(in) :: energy, r1, r2, u1, u2
        real(rp) :: delta, sdl, k2
        delta = calc_delta(this, l, energy, r1, r2, u1, u2)
        sdl = sin(delta)
        k2 = this%alpha*energy ! 2mE/hbar^2
        sigma_l = 4*pi*sdl*sdl/k2

    end function partial_cross_section

! this subroutine cannot be pure because it may also produce a file
    subroutine solve(this, verbose_, tofile_, storeData_)
        class(Potential_t), intent(inout) :: this ! inout necessary only if storeE
        !! the potential for the problem
        logical, intent(in), optional :: verbose_
        !! whether or not to print to stdout (default false)
        logical, intent(in), optional :: tofile_
        !! whether or not to send data to a file (default true)
        logical, intent(in), optional :: storeData_
        !! whether ot not to keep data in memory for further calculations (default false)

        logical :: verbose = .true., tofile = .true., storeData = .false.

        integer :: fileunit
        real(rp) :: dE, E_curr, denom, sig = 0., delta
        integer :: i, l

        verbose = .false.
        tofile = .true.
        storeData = .false.
        if (present(verbose_)) verbose = verbose_
        if (present(tofile_)) tofile = tofile_
        if (present(storeData_)) storeData = storeData_

        ! TODO implement
        ! if (storeData) allocate(this%data(this%numE, 2))

        ! dE = (this%maxE - this%minE)/this%numE
        ! E_curr = this%minE

        ! if (tofile) then
        !     open(newunit=fileunit, file='sigma.dat')
        !     write(fileunit, *) '# energy sigma_tot'
        ! end if

        ! if (verbose) write(stdout, *) 'iteration energy sigma tan(delta)'

        ! do i=1, this%numE ! +1 ! should I include a +1 here?
        !     denom = sqrt(this%alpha * E_curr)! TODO
        !     do l=0, this%lmax
        !         delta = this%calc_delta(l, E_curr, denom)
        !         sig = sig + 4*pi/(denom*denom)*(2*l+1)*sin(delta)**2
        !     end do
        !     E_curr = E_curr + dE
        !     if (tofile) write(fileunit, *) E_curr, sig
        !     if (verbose) write(stdout, *) i, E_curr, sig, tan(delta)
        ! end do

        ! print*,this%end_r ! TODO STUB

        ! if (tofile) then
        !     if (verbose) write(stdout, *) 'cross section written to sigma.dat'
        !     close(fileunit)
        ! end if

    end subroutine solve

    pure real(rp) function get_wavelength(this, energy) result(wvl)
        class(Potential_t), intent(in) :: this
        real(rp), intent(in) :: energy
        wvl = 2*pi/sqrt(this%alpha*energy)
    end function get_wavelength

end module potential_type

