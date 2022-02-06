program loss_rate
    !! as of right now, my program is basically a refactoring of Thijssen's
    !! chapter two exercise solution
    !! however, it looks like his approach is different from what is outlined
    !! in his book. The book makes more sense to me, so I will refactor further
    !! to match that
    !! test against /home/jph/thijssen/ch2/sigmadat
    ! (generated by solution to textbook exercise, provided by the author)
    use, intrinsic :: iso_fortran_env, only: rp => real32
    use json_module
    ! use precision, only: rp
    use mpi_f08
    use potential_type, only: Potential_Harmonic_t
    ! use solver_type, only: Solver_t
    implicit none
    integer :: t0
    ! for some reason, config%get doesn't like reals...
    real(rp) :: dt, tf, mu
    real(rp),dimension(:),allocatable :: x0
    type(json_file) :: config
    logical :: found
    ! type(Potential_t) :: pot
    ! type(Solver_t) :: model

    call config%initialize(comment_char = '//')

    !load the file:
    call config%load('config.json'); if (config%failed()) stop

    !read in the data:
    call config%get('t0',t0,found); if (.not. found) stop
    ! call config%json_file_get_real('dt',dt,found); if (.not. found) stop

    ! todo read input parameters from json file
    ! for now, defaulting to some values inside Potential_t
    ! pot = Potential_t()
    ! model = Solver_t()
    print*, "test"
    ! print*, model%filename
    ! print*, pot%potential(5._rp)

    contains

    ! subroutine scatter_init()

end program loss_rate


! real(rp) :: t0, dt, tf, mu
! real(rp),dimension(:),allocatable :: x0
! type(json_file) :: config
! logical :: found
!     ! call MPI_Init_f08()
!     ! call MPI_Init_f08
!     ! call MPI_Init()


!     ! call MPI_Comm_size_f08()

!     !initialize the class:
! call config%initialize(comment_char = '//')

! !load the file:
! call config%load_file('config.json'); if (config%failed()) stop

! !read in the data:
! call config%get('t0',t0,found); if (.not. found) stop
! call config%get('dt',dt,found); if (.not. found) stop
! call config%get('tf',tf,found); if (.not. found) stop
! call config%get('mu',mu,found); if (.not. found) stop
! call config%get('x0',x0,found); if (.not. found) stop

! ! call config%
! ! call config

! print*,mu

! !propagate:
! ! ...

! !cleanup:
! call config%destroy()
! ! call MPI_Finalize()
! ! calculate the loss rate
! ! print *, "todo stub"