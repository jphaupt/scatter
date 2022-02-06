module solver_type
    !! the main driver for the calculation
    use precision, only: rp
    use special_functions, only: sphJ => sphericalJ, sphN => sphericalN
    use potential_type, only: Potential_t, Potential_LennardJones_t
    use json_module
    implicit none

    public :: Solver_t, json_get_helper

    type :: Solver_t
        !! type containing all the system data, including the potential, min and
        !! max energy, maximum l, etc.
        class(Potential_t), allocatable :: pot
        !! not sure if this should be a pointer or allocatable variable
        character(len=:), allocatable :: filename
        real(rp) :: maxE=3.5, minE=0.1, &
                &   start_r=0.75, end_r=5.0, rstep=0.02
        integer :: lmax=6, numE=200

    contains
        procedure, pass(this) :: read_LennardJones

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
        real(rp) :: x
        character(len=:), allocatable :: pot_name
        if (present(filename)) then
            this%filename = filename
        else
            this%filename = "input.json"
        end if
        call config%initialize(comment_char='//')
        ! load file
        call config%load_file(this%filename)
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
        end select
        ! TODO read_general(this, config) ! to get maxE, lmax, etc.
    !! TODO
        allocate (Potential_LennardJones_t::this%pot)

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
        ! TODO

    end subroutine read_general

    subroutine read_LennardJones(this, config)
        !! initialises the Lennard-Jones potential from a JSON file
        implicit none
        class(Solver_t), intent(inout) :: this
        type(json_file), intent(inout) :: config
        real(rp) :: epsilon, rho, m1, m2
        ! ensure Lennard-Jones potential
        allocate (Potential_LennardJones_t::this%pot)
        ! in order for %get to work here, json_file must be inout for some reason
        ! remember it allows default values
        call json_get_helper(config, 'm1', m1)
        call json_get_helper(config, 'm2', m2)
        call json_get_helper(config, 'epsilon', epsilon)
        call json_get_helper(config, 'rho', rho)
        ! constructor for Lennard-Jones potential
        this%pot = Potential_LennardJones_t(m1, m2, rho, epsilon)
    end subroutine read_LennardJones

    subroutine json_get_helper_real(config, varname, var)
    !! just a simple helper for when I need to run config%get with defaults
    !! right now only works for real var
        implicit none
        type(json_file), intent(inout) :: config
        character(*), intent(in) :: varname
        logical :: found
        real(rp), intent(inout) :: var
        call config%get(varname, var, found); if (.not. found) &
 &       print *, "using default "//varname
    end subroutine json_get_helper_real

    subroutine json_get_helper_int(config, varname, var)
    !! there should be a better way to do this with less repetition
        implicit none
        type(json_file), intent(inout) :: config
        character(*), intent(in) :: varname
        logical :: found
        integer, intent(inout) :: var
        call config%get(varname, var, found); if (.not. found) &
 &       print *, "using default "//varname
    end subroutine json_get_helper_int

end module solver_type
