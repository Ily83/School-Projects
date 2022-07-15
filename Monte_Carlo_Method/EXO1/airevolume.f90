module airevolume 
use ogpf ! module Traçage 
use, intrinsic :: iso_fortran_env    
implicit none
type(gpf)  :: gp
integer, parameter :: dp    = real64
real(dp) :: aire,t1,t2,a,b,R,x0,y0,z0


contains    

! ---  Aire Cercle

subroutine cercle(N,x0,y0,R,aire,plot)
    integer     :: N, plot,cpt
    real(dp)                  :: R,x0,y0
    real(dp),allocatable      :: xtab(:),ytab(:) 
    real(dp)                  :: var,ecart, moyenne, x,y,aire, nd
    allocate(xtab(N+1),ytab(N+1))
    xtab = 0._dp; ytab = 0._dp;

    ! Initialisations du plot
    ! call gp%title('Monte Carlo Cercle')
    call gp%axis([x0-R, x0+R, y0-R, y0+R])
    call gp%options('set key top left')
    call gp%options('set autoscale fix')
    call gp%options('set style line 1 lc rgb "blue" lt 1 lw 2 pt 6 ps 1.5')

    nd = 0 ! Nombre de points qui sont dedans
    cpt = 0 ! Compteurs d'essais
    do while(cpt<N)
        ! x = random() ; y = random() ;
        x = rangerandom(x0-R,x0+R); y = rangerandom(y0-R,y0+R)
        if ( (x-x0)**2._dp +(y-y0)**2._dp <= R**2._dp) then
            nd = nd + 1
            xtab(cpt+1) = x
            ytab(cpt+1) = y
        endif
        cpt=cpt+1
    enddo
    moyenne = nd/real(N)
    ! var = sqrt(moyenne*(1-moyenne)/(N-1))
    var = sqrt((real(n)/(real(n)-1))*(moyenne - moyenne*moyenne))
    ecart = 1.96_dp*sqrt(var/real(N))
    ! aire = (2*R)**2*moyenne
    aire = (2_dp*R)**2_dp*moyenne
    if (plot == 1) call gp%plot(xtab, ytab, 'title "x y" with points pt 7 ps 1.2 lc rgb "#ad2060"')
    ! print *,var,ecart
    print *, " -------------------------------------- "
    print *, "Intervalle de confiance: [",aire-ecart,";",aire+ecart,"]"
    print *, " -------------------------------------- "
    deallocate(xtab,ytab)
end subroutine Cercle


! --- Ellipse

subroutine ellipse(N,x0,y0,a,b,aire,plot)
    integer     :: N, plot, nd,cpt
    real(dp),intent(in)       :: a,b,x0,y0
    real(dp),allocatable      :: xtab(:),ytab(:) 
    real(dp)                  :: var,ecart, moyenne, x,y,aire
    allocate(xtab(N+1),ytab(N+1))
    xtab = 0.0; ytab = 0.0;

    ! Initialisations du plot
    ! call gp%title('Monte Carlo Cercle')
    call gp%axis([x0-a, x0+a, y0-b, y0+b])
    call gp%options('set key top left')
    call gp%options('set autoscale fix')
    call gp%options('set style line 1 lc rgb "blue" lt 1 lw 2 pt 6 ps 1.5')

    nd = 0
    cpt = 0
	
    do while(cpt<N)
        ! x = random() ; y = random() ;
        x = rangerandom(x0-a,x0+a); y = rangerandom(y0-b,y0+b)
        if ( sqrt(((x-x0)/a)**2._dp +((y-y0)/b)**2._dp) <= 1.0 ) then
            nd = nd +1
            xtab(cpt+1) = x
            ytab(cpt+1) = y
        endif
        cpt=cpt+1
    enddo
	
    moyenne = nd/real(N)
    var = sqrt((real(n)/(real(n)-1))*(moyenne - moyenne*moyenne))
    ecart = 1.96_dp*sqrt(var/real(N))
    aire = 4*(a*b)*moyenne
	
    print *, " -------------------------------------- "
    print *, "Intervalle de confiance: [",aire-ecart,";",aire+ecart,"]"
    print *, " -------------------------------------- "
	
    if (plot == 1) call gp%plot(xtab, ytab, 'title "x y" with points pt 7 ps 1.2 lc rgb "#ad2060"')
	
    ! print *,var,ecart
    deallocate(xtab,ytab)
end subroutine ellipse

! --- Sphere

subroutine sphere(N,x0,y0,z0,R,volume,plot)
    integer     :: N, plot, nd,cpt
    real(dp),intent(in)       :: R,x0,y0,z0
    real(dp),allocatable      :: xtab(:),ytab(:),ztab(:)
    real(dp)                  :: var,ecart, moyenne, x,y,z,volume
    allocate(xtab(N+1),ytab(N+1),ztab(N+1))
    xtab = 0.0; ytab = 0.0;

    ! Initialisations du plot
    ! call gp%title('Monte Carlo Cercle')
    call gp%axis([x0-R, x0+R, y0-R, y0+R,z0-R,z0+R])
    call gp%options('set key off')
    call gp%options('set autoscale fix')
    call gp%options('set style line 1 lc rgb "blue" lt 1 lw 2 pt 6 ps 1.5')

    nd = 0
    cpt = 0
    do while(cpt<N)
        ! x = random() ; y = random() ;
        x = rangerandom(x0-R,x0+R); y = rangerandom(y0-R,y0+R); z = rangerandom(z0-R,z0+R);
        if ( sqrt((x-x0)**2._dp +(y-y0)**2._dp+(z-z0)**2._dp) <= R ) then
            nd = nd +1
            xtab(cpt+1) = x
            ytab(cpt+1) = y
            ztab(cpt+1) = z
        endif
        cpt=cpt+1
    enddo
    moyenne = nd/real(N)
    var = sqrt((real(n)/(real(n)-1))*(moyenne - moyenne*moyenne))
    ecart = 1.96_dp*sqrt(var/real(N))
    volume = (4./3.)*6*R**3*moyenne
    print *, " -------------------------------------- "
    print *, "Intervalle de confiance: [",volume-ecart,";",volume+ecart,"]"
    print *, " -------------------------------------- "
    if (plot == 1) call gp%scatter3d(xtab, ytab,ztab, lspec='using 1:2:3 with points palette pointsize 0.5 pointtype 7')
    ! print *,var,ecart
    deallocate(xtab,ytab)
end subroutine sphere

! Genérateur aléatoire
real(dp) function random() 
    integer(int64) :: ix=51477, na=16807, nmax=2_int64**31_int64-1

    ix=abs(ix*na)
    ix=mod(ix,nmax)
    random=(real(ix)/real(nmax))

end function random

real(dp) function rangerandom(xmin,xmax)
    real(kind=8)        :: xmin,xmax,xdiff

    xdiff = xmax - xmin
    rangerandom = (xdiff * random() )+xmin

end function

real(dp) function aire_cercle(R)
    real(dp)            :: R
    real(dp), parameter :: PI = 4*ATAN(1.d0)
    aire_cercle = PI*R**2
end function

real(dp) function aire_ellipse(a,b)
    real(dp)            :: a,b
    real(dp), parameter :: PI = 4*ATAN(1.d0)
    aire_ellipse = a*b*pi
end function

real(dp) function aire_sphere(R)
    real(dp)            :: R
    real(dp), parameter :: PI = 4*ATAN(1.d0)
    aire_sphere = (4./3.)*pi*R**3
end function

end module airevolume