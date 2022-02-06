module precision
    !! gives working precisions from iso fortran
    !!
    !! inspired by (and basically copied from)
    !! https://github.com/modern-fortran/neural-fortran/blob/master/src/mod_kinds.f90
    use, intrinsic :: iso_fortran_env, only: real32, real64, real128
    implicit none

    private
    public :: rp

#ifdef REAL32
    integer, parameter :: rp = real32 ! single precision
#elif REAL128
    integer, parameter :: rp = real128 ! quad precision
#else
    integer, parameter :: rp = real64 ! double precision
#endif

contains

end module precision
