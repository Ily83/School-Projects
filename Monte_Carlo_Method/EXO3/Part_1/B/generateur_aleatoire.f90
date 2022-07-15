module generateur_aleatoire
    use, intrinsic :: iso_fortran_env
    implicit none
    integer, parameter :: dp    = real64
! Genérateur aléatoire
    contains

real(dp) function random() 
    integer(int64) :: ix=51477, na=16807, nmax=2_int64**31_int64-1

    ix=abs(ix*na)
    ix=mod(ix,nmax)
    random=(real(ix)/real(nmax))
end function random

real(dp) function rangerandom(xmin,xmax)
    real(kind=8)  :: xmin,xmax,xdiff

    xdiff = xmax - xmin
rangerandom = (xdiff * random() )+xmin
end function rangerandom


end module generateur_aleatoire