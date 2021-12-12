module special_functions
    !! provides special functions as needed
    !! I suspect the only ones needed will be the spherical Bessel functions,
    !! but I've been known to be wrong.
    use precision, only: rp
    implicit none
    private
    public :: sphericalJ
contains
    elemental recursive real(rp) function sphericalJ(l, x) result(out)
        !! spherical Bessel function j_l(x)
        !! @todo surely there is a better way to do this? Some stdlib?
        ! TODO test this function?
        implicit none
        integer, intent(in) :: l
        real(rp), intent(in) :: x
        ! real(rp), intent(out) :: out

        if ( l == 0 ) then
            out = j0(x)
        else if (l == 1 ) then
            out = j1(x)
        else
            out = (2*l-1)*sphericalJ(l-1,x)/x - sphericalJ(l-2,x)
        end if
    end function sphericalJ

    pure real(rp) function j0(x)
        ! spherical Bessel with l=0
        ! TODO ? put in "contains" keyword inside sphericalJ?
        implicit none
        real(rp), intent(in) :: x
        j0 = sin(x)/x
    end function j0

    pure real(rp) function j1(x)
        implicit none
        real(rp), intent(in) :: x
        j1 = sin(x)/(x*x) - cos(x)/x
    end function j1

    elemental recursive real(rp) function sphericalN(l, x) result(out)
        !! returns spherical Bessel function n_l(x)
        ! defined recursively (less efficient I think, but pretty)
        implicit none
        integer, intent(in) :: l
        real(rp), intent(in) :: x
        if ( l == 0 ) then
            out = n0(x)
        else if ( l == 1 ) then
            out = n1(x)
        else
            out = (2*l-1)*sphericalN(l-1,x)/x - sphericalN(l-2,x)
        end if

    end function sphericalN

    pure real(rp) function n0(x)
        implicit none
        real(rp), intent(in) :: x
        n0 = -cos(x)/x
    end function n0

    pure real(rp) function n1(x)
        implicit none
        real(rp), intent(in) :: x
        n1 = -cos(x)/(x*x) - sin(x)/x
    end function n1
end module special_functions