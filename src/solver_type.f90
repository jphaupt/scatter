module solver_type
    !! the main driver for the calculation
    !! @todo MPI drivers I think will be here
    use precision, only: rp
    use special_functions, only: sphJ => sphericalJ, sphN => sphericalN
    use potential_type, only: Potential_t, Potential_LennardJones_t
    use json_module
    use integrate, only: numerov, radialRHS
    ! TODO figure out how to make an external member function
    implicit none

    public :: Solver_t, json_get_helper

    type :: Solver_t
        !! type containing all the system data, including the potential, min and
        !! max energy, maximum l, etc.
        class(Potential_t), allocatable :: pot
        !! not sure if this should be a pointer or allocatable variable
        character(len=:), allocatable :: filename
        real(rp) :: maxE=3.5, minE=0.1, &
                &   start_r=0.75, end_r=5.0, rstep=0.02 !, &
                ! &   startPt, startVal, nextPt, nextVal
        integer :: lmax=6, numE=200
            !! @todo these two are not set properly without an input file

    contains
        procedure, pass(this) :: read_LennardJones
        procedure, pass(this) :: numerov => numerov_solver
        procedure, pass(this) :: calc_full_cross_section
        ! no idea if you can do this (creates circular dependency as written)

    end type Solver_t

    interface Solver_t
        module procedure :: constructor_solver
    end interface Solver_t

    ! generic procedure
    interface json_get_helper
        module procedure :: json_get_helper_real, json_get_helper_int
    end interface json_get_helper

    ! procedure :: json_get_helper =>
contains

! this function cannot be pure because it reads the JSON file
    type(Solver_t) function constructor_solver(filename) result(this)
        implicit none
        character(*), optional, intent(in) :: filename
        type(json_file) :: config
        logical :: found
        character(len=:), allocatable :: pot_name
        if (present(filename)) then
            this%filename = filename
        else
            this%filename = "input.json"
        end if
        call config%initialize(comment_char='//')
        ! load file
        call config%load(this%filename)
        if (config%failed()) then
            print *, "Please include an input file!"
            stop
        end if
        ! get the type of the potential
        call config%get('potential_name', pot_name, found)
        if (.not. found) then
            print *, "Please include a potential type!"
            stop
        end if
        select case (trim(pot_name))
        case ("LennardJones")
            print *, "Initialising LennardJones Potential!"
            call this%read_LennardJones(config)
        case default
            print *, "Potential type "//trim(pot_name)//" not yet implemented"
            print *, "(Remember, it is case sensitive!)"
            stop
        end select
        call read_general(this, config) ! to get maxE, lmax, etc.

    end function constructor_solver

    subroutine read_general(this, config)
        !! initalises paramers in the solver that are not potential-specific
        !! such as maximum l angular momentum, maximum energy, minimum r, etc
        implicit none
        class(Solver_t), intent(inout) :: this
        type(json_file), intent(inout) :: config
        ! real(rp) :: maxE, minE, numE, lmax, start_r, end_r, rstep
        ! TODO test this one to make sure values actually change
        call json_get_helper(config, 'maxE', this%maxE)
        call json_get_helper(config, 'minE', this%minE)
        call json_get_helper(config, 'numE', this%numE)
        call json_get_helper(config, 'lmax', this%lmax)
        call json_get_helper(config, 'start_r', this%start_r)
        call json_get_helper(config, 'end_r', this%end_r)
        call json_get_helper(config, 'maxE', this%maxE)

    end subroutine read_general

    subroutine read_LennardJones(this, config)
        !! initialises the Lennard-Jones potential from a JSON file
        implicit none
        class(Solver_t), intent(inout) :: this
        type(json_file), intent(inout) :: config
        real(rp) :: epsilon, rho, m1, m2
            !! hacky fix: no skipped variables, still not working somehow :/
        ! ensure Lennard-Jones potential
        allocate (Potential_LennardJones_t::this%pot)
        ! in order for %get to work here, json_file must be inout for some reason
        ! remember it allows default values
        call json_get_helper(config, 'm1', m1)
        call json_get_helper(config, 'm2', m2)
        call json_get_helper(config, 'epsilon', epsilon)
        call json_get_helper(config, 'rho', rho)
        ! constructor for Lennard-Jones potential
        !! TODO !! the constructor below is taking the uninitialized variables!
        this%pot = Potential_LennardJones_t(m1, m2, rho, epsilon)
    end subroutine read_LennardJones

    subroutine json_get_helper_real(config, varname, var)
    !! just a simple helper for when I need to run config%get with defaults
        implicit none
        type(json_file), intent(inout) :: config
        character(*), intent(in) :: varname
        logical :: found
        real(rp), intent(out) :: var
        ! TODO for some reason it doesn't like reals...
        ! make sure to build json-fortran with the same type as this program!
        ! (I default to real64, they default to real32)
        call config%get(varname, var, found);
        if (.not. found) then
            print *, "default values not yet enabled "//varname
            stop
        end if
    end subroutine json_get_helper_real

    subroutine json_get_helper_int(config, varname, var)
    !! there should be a better way to do this with less repetition
        implicit none
        type(json_file), intent(inout) :: config
        character(*), intent(in) :: varname
        logical :: found
        integer, intent(out) :: var
        call config%get(varname, var, found)
        if (.not. found) then
            print *, "default values not yet enabled "//varname
            stop
        end if
    end subroutine json_get_helper_int


    ! pure 
    subroutine numerov_solver(this, l, energy, endPoint1, endVal1, endPoint2, endVal2, rsep_)
        !! wrapper function for the Numerov function below for use in the Solver
        !! type (OOP approach)
        !! this does not store the solution
        implicit none
        class(Solver_t), intent(in) :: this
        real(rp), intent(out) :: endPoint1, endVal1, endPoint2, endVal2
        real(rp), intent(in), optional :: rsep_
        real(rp), intent(in) :: energy
        integer, intent(in) :: l
        real(rp) :: startVal, nextVal
        real(rp) :: h_sq, h2f, h2fnext, h2fprev, deriv
        integer :: maxI
        maxI = nint((this%end_r-this%start_r)/this%rstep)
        startVal = this%pot%small_r_solution(this%start_r)
        ! nextVal complicated to get...
        ! see equations A.53 and A.54
        h_sq = this%rstep * this%rstep
        deriv = this%pot%small_r_derivative(this%start_r) ! at r=0
        h2f = h_sq * radialRHS(this%pot, this%start_r, l, energy) ! at r=0
        ! at r=h
        h2fnext = h_sq * radialRHS(this%pot, this%start_r+this%rstep, l, energy)
        ! at r=-h
        h2fprev = h_sq * radialRHS(this%pot, this%start_r-this%rstep, l, energy)
        ! I think I actually did wnext here?!
        nextVal = ((2_rp + 5*h2f/6)*(1_rp - h2fprev/12)*startVal &
        &         + 2*this%rstep*deriv*(1_rp - h2fprev)/6) &
        &       / ((1_rp - h2fnext/12)*(1_rp - h2fprev/6) &
        &        + (1_rp - h2fprev/12)*(1_rp - h2fnext/6))
        ! nextVal = h2f*(1_rp-h2fprev/6) & ! this is completely wrong...
        ! &         + 2*this%rstep*deriv*(1_rp-h2fprev/12) &
        ! &       / ((1_rp - h2fnext/12)*(1_rp - h2fprev/6) &
        ! &        + (1_rp - h2fprev/12)*(1_rp - h2fnext/6))
        ! @todo ? maybe something is wrong with the nextVal calculation

        ! print*,"phideriv", deriv
        ! print*, "phi1, phi2"
        ! print*, startVal, nextVal

        call numerov(this%pot, this%start_r, startVal, this%start_r+this%rstep,&
        &            nextVal, l, energy, maxI, endPoint1, endVal1, endPoint2, &
        &            endVal2, rsep_=this%pot%get_wavelength(energy)/4_rp)

    end subroutine numerov_solver

    subroutine calc_full_cross_section(this)
        !! calculates the full cross section given all the parameters in Solver
        !! @todo output to file!
        implicit none
        class(Solver_t) :: this

        real(rp) :: E_curr, dE, sigma_tot
        integer :: l
        real(rp) :: r1, r2, u1, u2
        !! @todo
        !! for some reason, numE is zero here, looks like the default doesn't work
        !! properly...
        !! same thing with lmax
        !! @endtodo
        dE = (this%maxE - this%minE)/this%numE
        ! print*, this%maxE
        ! print*, this%minE
        ! print*, this%numE
        ! print*, dE
        E_curr = this%minE
        ! TODO !!!
        E_loop: do while (E_curr <= this%maxE)
            ! print*, "E", E_curr
            ! print*, "r1,u1,r2,u2"
            sigma_tot = 0
            l_loop: do l=1,this%lmax
                call this%numerov(l, E_curr, r1, u1, r2, u2)
                ! print*, r1,u1,r2,u2
                ! TODO numerov values blow the fuck up
                ! TODO numerov needs even more testing :(
                sigma_tot = sigma_tot &
                &   + this%pot%partial_cross_section(l, E_curr, r1, r2, u1, u2)
                ! print*, l
            end do l_loop
            E_curr = E_curr + dE
            print*, E_curr, sigma_tot
        end do E_loop
        ! print*, "stub"
    end subroutine


end module solver_type
