module solver_type
    !! the main driver for the calculation
    use precision, only: rp
    use special_functions, only: sphJ => sphericalJ, sphN => sphericalN
    use potential_type, only: Potential_t, Potential_LennardJones_t
    use json_module
    implicit none

    type :: Solver_t
        !! type containing all the system data, including the potential, min and
        !! max energy, maximum l, etc.
    class(Potential_t), allocatable :: pot
        !! not sure if this should be a pointer or allocatable variable
    character(len=:), allocatable :: filename

    contains
    procedure, pass(this) :: read_LennardJones

    end type Solver_t

    interface Solver_t
        module procedure :: constructor_solver
    end interface Solver_t
    ! print*, pot%potential(5._rp)

contains

! this function cannot be pure because it reads the JSON file
type(Solver_t) function constructor_solver(filename) result(this)
    implicit none
    character(*), optional, intent(in) :: filename
    type(json_file) :: config
    logical :: found
    character(len=:), allocatable :: pot_name
    if(present(filename)) then
        this%filename=filename
    else
        this%filename="input.json"
    endif
    call config%initialize(comment_char = '//')
    ! load file
    call config%load_file(this%filename)
    if(config%failed()) then
        print*, "Please include an input file!"
        stop
    endif
    ! get the type of the potential
    call config%get('potential_name',pot_name,found)
    if (.not. found) then
        print*, "Please include a potential type!"
        stop
    endif
    select case (trim(pot_name))
        case ("LennardJones")
            call this%read_LennardJones(config)
            print*, "LennardJones Potential!"
        case default
            print*, "Potential type "//trim(pot_name)//" not yet implemented"
            print*, "(Remember, it is case sensitive!)"
    end select
    ! TODO read_general(this, config) ! to get maxE, lmax, etc.
    !! TODO
    allocate(Potential_LennardJones_t::this%pot)

end function constructor_solver

subroutine read_LennardJones(this, config)
    !! initialises the Lennard-Jones potential from a JSON file
    implicit none
    class(Solver_t), intent(inout) :: this
    type(json_file), intent(in) :: config
    real(rp) :: epsilon, rho, m1, m2
    ! ensure Lennard-Jones potential
    allocate(Potential_LennardJones_t::this%pot)
end subroutine read_LennardJones

end module solver_type