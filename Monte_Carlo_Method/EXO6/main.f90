program MN_3
    use, intrinsic :: iso_fortran_env    
    use ogpf ! module Traçage 
    use generateur_aleatoire
    implicit none
    type(gpf)  :: gp
    real(dp)   :: t1,t2
! ----------------------
    integer    :: N, i,plot,choix
	real(dp)   :: Temp=0._dp,x0,y0,z0,dt,eps
    ! call input(x0,y0,z0)
    x0 = 2.5_dp; y0 = 2.5_dp; z0 = 1.25_dp;
	N = 10000;
	dt = 1.d-4
    call cpu_time(t1)
    call PbStudio(N,x0,y0,z0,dt,Temp)
    call cpu_time(t2)
    print *, "Temps d'éxecution (mvt brownien) : ",t2-t1
    print *, "Temperature (mvt brownien)       :", Temp
    eps = 1d-4
    call cpu_time(t1)
    call MarcheSphere(N,x0,y0,z0,eps,Temp)
    print *, "Temperature (Marche Sphere) :", Temp
    call cpu_time(t2)
    print *, "Temps d'éxecution: (Marche sphere) ",t2-t1
    contains

subroutine PbStudio(N,x0,y0,z0,dt,Temp)
    implicit none
    Integer   :: N,i,Ext
    real(dp)  :: G(3)
    real(dp)  :: x0,y0,z0,dt,Temp,xn,yn,zn, TempCoin
    real(dp), parameter   :: PI = 4._dp*DATAN(1.d0)
    Temp = 0._dp
    Ext = 0 ! On commence à l'interieur
    do i = 1,N
        xn = x0; yn = y0; zn=z0; 
        do while (Ext == 0) 
            ! Tirages Gaussien
            G = Gaussien3D(random(),random(),random(),random());
            !MAJ des coordoonees
            xn = xn + sqrt(dt)*G(1);
            yn = yn + sqrt(dt)*G(2);
            zn = zn + sqrt(dt)*G(3);
            ! Check si on est sorti d'un coin 
            call Exits(xn,yn,zn,TempCoin,Ext)
        end do 
        Temp = Temp + TempCoin
        Ext = 0;
    end do 
    Temp = Temp/dble(N)
end subroutine Pbstudio

Subroutine MarcheSphere(N,x0,y0,z0,eps,Temp)
    implicit none
    Integer   :: N
    real(dp)  :: G(3) = 0._dp,Temp,Tempcoin =0._dp
    real(dp)  :: x0,y0,z0,xn,yn,zn,theta,r,eps
    do i = 1,N
        xn = x0; yn = y0; zn=z0;
        call rayon_temperature(x0,y0,z0,r,Tempcoin)
        do while (r>eps) 
            G = Gaussien3D(random(),random(),random(),random());
            xn = xn + r*G(1)/sqrt(G(1)**2+G(2)**2+G(3)**2);
            yn = yn + r*G(2)/sqrt(G(1)**2+G(2)**2+G(3)**2);
            zn = zn + r*G(3)/sqrt(G(1)**2+G(2)**2+G(3)**2);
            call rayon_temperature(xn,yn,zn,r,Tempcoin)
        end do 
        Temp = Temp + Tempcoin
    end do 
    Temp = Temp/dble(N)
end subroutine MarcheSphere



Subroutine rayon_temperature(x,y,z,r,temp)
 implicit none
 real(dp) :: x,y,z,temp,r,d(9)

 d = (/ abs(x),abs(x-5._dp),abs(y),abs(y-5._dp),abs(z), abs(z-2.5_dp), &
       sqrt(cercle(x-5._dp,y-5._dp))-1,sqrt(cercle(x,y))-1, 0._dp/)
    if (y >= 3._dp) then 
        if (x <= 4._dp) then
            d(9) = sqrt((4._dp - x)**2 + (3._dp - y)**2)
        else
            d(9) = abs(y - 3._dp)
        end if
    else if(y<=1) then
        if (x<= 4._dp)  then 
            d(9) = sqrt((4._dp - x)**2 + (3._dp - y)**2)
        else
            d(9) = abs(y - 1._dp)
        end if
    else
        d(9) = abs(x - 4._dp)
    end if
    r = minval(d,dim=1)

    ! Cercles
    if ((r == d(7)).or.(r==d(8))) then
        temp = 40._dp
    ! Sud/Nord
    else if (r == d(9)) then
        temp = 16._dp
    !Ouest
    else if(r == d(1)) then 
        temp = 10._dp
    ! Plafond 
    else if (r == d(6)) then
        temp = 14._dp 
    else ! le reste
        temp = 16._dp 
    end if
end subroutine rayon_temperature


Subroutine Exits(x,y,z,TempCoin,Ext)
    real(dp) :: x,y,z, TempCoin
    integer  :: Ext
    ! Face Est (9)
    if ((x<=0._dp).and.(y>=1._dp)) then 
        TempCoin = 16._dp; Ext = 1;
    ! Face ouest avant (4) + Faces Ouest sur les cotées du rectangle (3,5)
    else if((x>=4._dp).and.(1._dp<=y).and.(y<=3._dp)) then
        TempCoin = 16._dp; Ext = 1;
    ! Faces Ouest arriére (2,6)
    else if((x>=5._dp).and.(3._dp<=y).and.(y<=4._dp).and.(0._dp<=y).and.(y<=1._dp) ) then 
        TempCoin = 16._dp; Ext = 1;
    ! Cercle du bas (10) ou cercle du haut (7)
    else if((cercle(x,y)<=1._dp).or.(cercle(x-5._dp,y-5._dp) <= 1._dp)) then
        TempCoin = 40._dp; Ext = 1;
    ! Face sud (1)
    else if ((x>=1._dp).and.(y<=0._dp)) then
        TempCoin = 16._dp; Ext = 1;
    ! Face Nord (8)
    else if ((x<=4._dp).and.(y>=5._dp)) then
        TempCoin = 16._dp; Ext = 1;
    ! SOL (11)
    else if((z<=0._dp)) then
        TempCoin = 16._dp; Ext = 1;
    ! Plafond (12)
    else if (z>=2.5_dp) then
        TempCoin = 14._dp; Ext = 1;
    end if
    ! print *, x,y,z, Tempcoin
end subroutine exits

! Equation d'un cercle
real(dp) function cercle(x,y)
    real(dp) :: x,y
    cercle = x**2 + y**2
end function cercle

subroutine input(x0,y0,z0)
    real(dp) :: x0,y0,z0
    write(*,*) "Entrer la position initiale x0 y0 z0 dans la piece (0-5;0-5;0-2.5)"
    read(*,*) x0, y0, z0
end subroutine input

end program MN_3