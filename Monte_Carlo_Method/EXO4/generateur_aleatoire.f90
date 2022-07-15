module generateur_aleatoire
    use, intrinsic :: iso_fortran_env
    implicit none
    integer, parameter :: dp  = real64
! Genérateur aléatoire
    contains

real(dp) function random() 
    integer(int64) :: ix=51477, na=16807, nmax=2_int64**31_int64-1

    ix=abs(ix*na)
    ix=mod(ix,nmax)
    random=(real(ix)/real(nmax))
end function random

real(dp) function rangerandom(xmin,xmax)
    real(kind=dp)  :: xmin,xmax,xdiff

    xdiff = xmax - xmin
rangerandom = (xdiff * random() )+xmin
end function rangerandom 



function Gaussien(U1,U2)
    real(dp), parameter   :: PI = 4_dp*atan(1.0_dp)
    real(dp) :: U1,U2, Gaussien(2)
    Gaussien(1) = (-2._dp*log(U1))**(0.5_dp)*cos(2._dp*PI*U2);
    Gaussien(2) = (-2._dp*log(U1))**(0.5_dp)*sin(2._dp*PI*U2);
end function Gaussien


end module generateur_aleatoire