program MN_3
    use, intrinsic :: iso_fortran_env    
    use ogpf ! module Traçage 
    use generateur_aleatoire
    implicit none
    type(gpf)  :: gp
    real(dp) :: t1,t2
! ----------------------
    integer :: N,Nit, i,plot
    real(wp), allocatable :: maximum(:,:), minimum(:,:)
    plot = 0 ! 0 désactivé, 1 activé
    call cpu_time(t1)
    N = 50
    Nit = 100000
    allocate(maximum(2,N), minimum(2,N))
    maximum = 0._wp; minimum = 0._wp;
    call algomax(N,Nit,plot,maximum)
    print *,"maximum: (",maximum(1,1),",",maximum(2,1),")"

    call algomin(N,Nit,plot,minimum)
    print *,"minimum: (", minimum(1,1),",",minimum(2,1),")"

    call cpu_time(t2)
    print *, "Temps d'execution: ",t2-t1
    deallocate(maximum,minimum)
contains

subroutine algomax(N,Nit,plot,vect)
implicit none 
integer  :: N, Nit,xi,yi,plot
real(dp) :: vect(2,N),G(2)
real(dp), parameter :: a = 1._dp, b = 1._dp, c = 1.0_dp, d = 1._dp
real(dp), parameter :: p = 0.999

vect(1,:) = (/( 10.d0*(random()-0.5d0), i = 1,N )/)
vect(2,:) = (/( 10.d0*(random()-0.5d0), i = 1,N )/)

do i = 1,Nit
    xi = int(N*random()); yi = int(N*random());
    G = Gaussien(random(),random());
    if ( f1(vect(1,xi),vect(2,xi),a,b,c,d) < f1(vect(1,yi),vect(2,yi),a,b,c,d) ) then
        vect(1,xi) = vect(1,yi) + p**i*G(1)
        vect(2,xi) = vect(2,yi) + p**i*G(2)
    else
        vect(1,yi) = vect(1,xi) + p**i*G(1)
        vect(2,yi) = vect(2,xi) + p**i*G(2)
    end if
end do 
if (plot == 1) then
    call gp%axis([-3._wp, 3._wp, -3._wp, 3._wp])
    call gp%options('set key off')
    call gp%options('set autoscale fix')
    call gp%plot(vect(1,:),vect(2,:), 'title "Points Kouhonen" with points pt 7 ps 1.5 lc rgb "#011582"')
end if 
end subroutine algomax

subroutine algomin(N,Nit,plot,vect)
    implicit none 
    integer  :: N, Nit,xi,yi,plot
    real(dp) :: vect(2,N),G(2)
    real(dp), parameter :: a = 1._dp, b = 1._dp, c = 1._dp, d = 1._dp
    real(dp), parameter :: p = 0.999
    
    vect(1,:) = (/( 10.d0*(random()-0.5d0), i = 1,N )/)
    vect(2,:) = (/( 10.d0*(random()-0.5d0), i = 1,N )/)
    
    do i = 1,Nit
        xi = int(N*random()); yi = int(N*random());
        G = Gaussien(random(),random());
        if ( f2(vect(1,xi),vect(2,xi),a,b,c,d) > f2(vect(1,yi),vect(2,yi),a,b,c,d) ) then
            vect(1,xi) = vect(1,yi) + p**i*G(1)
            vect(2,xi) = vect(2,yi) + p**i*G(2)
        else
            vect(1,yi) = vect(1,xi) + p**i*G(1)
            vect(2,yi) = vect(2,xi) + p**i*G(2)
        end if
    end do 
    if (plot == 1) then
        call gp%axis([-3._wp, 3._wp, -3._wp, 3._wp])
        call gp%options('set key off')
        call gp%options('set autoscale fix')
        call gp%plot(vect(1,:),vect(2,:), 'title "Points Kouhonen" with points pt 7 ps 1.5 lc rgb "#011582"')
    end if
    end subroutine algomin

! ------ Fonctions --------

real(dp) function f1(x,y,a,b,c,d)
    real(dp) x,y,a,b,c,d
    f1 = a*exp(-b*((x-1.d0)**2.d0 + (y-2.d0)**2.d0))+c*exp(-d*((x+1.d0)**2.d0+(y+3.d0)**2.d0))
end function

! f(x, y) = 20 + x^2 + y^2 - 10(cos(2pix) + cos(2piy))
real(dp) function f2(x,y,a,b,c,d)
    real(dp) x,y,a,b,c,d
    real(dp), parameter   :: PI = 4_dp*atan(1.0_dp)
    f2 = 20._dp + x**2._dp + y**2._dp - 10._dp*(cos(2._dp*PI*x)+cos(2._dp*PI*y))
end function

end program MN_3