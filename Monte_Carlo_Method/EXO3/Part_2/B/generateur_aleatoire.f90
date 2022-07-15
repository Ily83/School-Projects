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


real(dp) function G1()
    real(dp) :: U1,U2
    real(dp), parameter   :: PI = 4_dp*atan(1.0_dp)
    U1 = random(); U2 = random();
    G1 = sqrt(-2_dp*log(U1))*cos(2_dp*PI*U2)
end function G1 

real(dp) function G2()
    real(dp) :: U1,U2
    real(dp), parameter   :: PI = 4_dp*atan(1.0_dp)
    U1 = random(); U2 = random();
    G2 = sqrt(-2_dp*log(U1))*sin(2_dp*PI*U2)
end function G2
end module generateur_aleatoire