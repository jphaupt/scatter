module constants
    use precision, only: rp
    implicit none
    real(rp), parameter :: pi = 4_rp*atan(1._rp)
    real(rp), parameter :: hbar_sq = 4.18015929_rp
        !! \( hbar^2 \) in units of Dalton*Angstrom^2*meV
contains
    
end module constants