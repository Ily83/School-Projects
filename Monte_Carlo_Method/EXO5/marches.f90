module marches
use, intrinsic :: iso_fortran_env    
use ogpf ! module Traçage 
use generateur_aleatoire
implicit none

contains

subroutine MvtBrownien_1(N,x0,y0,dt,Temp)
    integer               :: k,N,i, CL
    real (dp)             :: G(2),dt,Temp,x0,y0,xn,yn
    real(dp), parameter   :: PI = 4._dp*DATAN(1.d0)

    Temp=0.d0
    xn = x0; yn = y0;

    do i=1,N
        do while(dedans(xn,yn) == 1)
            G=Gaussien(random(),random())
            xn=xn+sqrt(dt)*G(1); yn=yn+sqrt(dt)*G(2)

            if(yn >= 1._dp) then ! Sortie Nord
                Temp = Temp + 100._dp
            endif

        enddo
        xn = x0; yn = y0;
    enddo
		
    Temp=Temp/dble(N)


end subroutine MvtBrownien_1

subroutine MvtBrownien_1_CL(N,x0,y0,dt,Temp)
    integer               :: k,N,i, CL
    real (dp)             :: G(2),dt,Temp,x0,y0,xn,yn
    real(dp), parameter   :: PI = 4._dp*DATAN(1.d0)

    Temp=0.d0
    xn = x0; yn = y0;

    do i=1,N
        do while(dedans(xn,yn) == 1)
            G=Gaussien(random(),random())
            xn=xn+sqrt(dt)*G(1); yn=yn+sqrt(dt)*G(2)
        enddo
        Temp = Temp + xn**2._dp - yn**2._dp
        xn = x0; yn = y0;
    enddo
		
    Temp=Temp/dble(N)


end subroutine MvtBrownien_1_CL

subroutine MvtBrownien_2(N,x0,y0,dt,Temp)
    integer               :: k,N,i, CL
    real(dp)             :: G(2),dt,Temp,x0,y0,xn,yn,distm,distp,u

    Temp=0.d0
    xn = x0; yn = y0;

    do i=1,N
        do while(dedans(xn,yn) == 1)
            G=Gaussien(random(),random())
            ! distance min avec les anciennes valeurs
            distm = min(abs(xn),abs(1._dp-xn),abs(yn),abs(1._dp-yn))
            ! Tirage des nouvelles valeurs
            xn=xn+sqrt(dt)*G(1); yn=yn+sqrt(dt)*G(2)
            ! distance min avec les nouvelles valeurs
            distp = min(abs(xn),abs(1._dp-xn),abs(yn),abs(1._dp-yn))
            u = random()

            if(yn >= 1.d0) then ! Sortie Nord
                Temp = Temp + 100._dp
            else if ( u < exp((-distm*distp)/(2*dt))) then 
                if (dp == abs(1._dp-yn)) Temp = Temp + 100._dp
            endif
        enddo
        xn = x0; yn = y0;
    enddo
		
    Temp=Temp/dble(N)


end subroutine MvtBrownien_2



subroutine MvtBrownien_2_CL(N,x0,y0,dt,Temp)
    integer               :: k,N,i, CL
    real(dp)             :: G(2),dt,Temp,x0,y0,xn,yn,distm,distp,u

    Temp=0.d0
    xn = x0; yn = y0;

    do i=1,N
        ! xn = random(); yn = random();
        do while(dedans(xn,yn) == 1)
            G=Gaussien(random(),random())
            ! distance min avec les anciennes valeurs
            distm = min(abs(xn),abs(1._dp-xn),abs(yn),abs(1._dp-yn))
            ! Tirage des nouvelles valeurs
            xn=xn+sqrt(dt)*G(1); yn=yn+sqrt(dt)*G(2)
            ! distance min avec les nouvelles valeurs
            distp = min(abs(xn),abs(1._dp-xn),abs(yn),abs(1._dp-yn))
            u = random()

            if(yn >= 1._dp) then ! Sortie Nord
                Temp = Temp + xn**2._dp - yn**2._dp
            else if ( u < exp((-distm*distp)/(2*dt))) then 
                if (dp == abs(1._dp-yn)) Temp = Temp + xn**2._dp - yn**2._dp
            endif
        enddo
        xn = x0; yn = y0;
    enddo
		
    Temp=Temp/dble(N)


end subroutine MvtBrownien_2_CL


integer function dedans(x,y)
    real(dp) :: x,y 
    if ((1.d0 >= x .and. x >= 0.d0).and.(1.d0 >= y .and. y >= 0.d0 )) then 
		dedans = 1
	else
		dedans = 0
	end if
end function dedans


subroutine MarcheCercle(N,x0,y0,eps,Temp)
    implicit none
    real(dp)             :: dt,Temp,x0,y0,xn,yn,eps,r,theta,d(4)
    integer               :: k,N,i,idx_bord
    real(dp), parameter   :: PI = 4._dp*DATAN(1.d0)

    Temp=0.d0

    Do i = 1,N
        xn = x0; yn = y0; r = 1.d0;
        do while(r>eps)
            theta = 2.d0*PI*random()
            ! Ouest , Sud, East, Nord 
            d = (/ abs(xn),abs(yn),abs(1.0_dp-xn),abs(1.0_dp-yn) /)
            r = minval(d,1)
            idx_bord = minloc(d,1)

            xn = xn + r*cos(theta)
            yn = yn + r*sin(theta) 
        end do 
        
        if (idx_bord == 4) Temp = Temp + 100
        
    end do 

    Temp=Temp/dble(N)

end subroutine MarcheCercle

subroutine MarcheCercleCL(N,x0,y0,eps,Temp)
    implicit none
    real (dp)             :: dt,Temp,x0,y0,xn,yn,eps,r,theta,d(4)
    integer               :: k,N,i,idx_bord
    real(dp), parameter   :: PI = 4._dp*DATAN(1.d0)
    
    Temp=0.d0

    Do i = 1,N
        xn = x0; yn = y0; r = 5.d-1;
        do while(r>eps)
            theta = 2*PI*random()
            ! Ouest , Sud, East, Nord 
            d = (/ abs(xn),abs(yn),abs(1._dp-xn),abs(1._dp-yn) /)
            r = minval(d,1)
            idx_bord = minloc(d,1)

            xn = xn + r*cos(theta)
            yn = yn + r*sin(theta) 
        end do 
        
        Temp = Temp + xn**2 - yn**2
        
    end do 

    Temp=Temp/dble(N)

end subroutine MarcheCercleCL

subroutine MvtBrownien_3(N,dt,Temp)
    integer               :: k,N,i, CL
    real(dp)             :: G(2),dt,Temp,xn,yn,distm,distp,u

    Temp=0.d0

    do i=1,N
        xn = random(); yn = random();
        do while(dedans(xn,yn) == 1)
            G=Gaussien(random(),random())
            ! distance min avec les anciennes valeurs
            distm = min(abs(xn),abs(1._dp-xn),abs(yn),abs(1._dp-yn))
            ! Tirage des nouvelles valeurs
            xn=xn+sqrt(dt)*G(1); yn=yn+sqrt(dt)*G(2)
            ! distance min avec les nouvelles valeurs
            distp = min(abs(xn),abs(1._dp-xn),abs(yn),abs(1._dp-yn))
            u = random()

            if(yn >= 1.d0) then ! Sortie Nord
                Temp = Temp + 100._dp
            else if ( u < exp((-distm*distp)/(2*dt))) then 
                if (dp == abs(1._dp-yn)) Temp = Temp + 100._dp
            endif
        enddo
    enddo
		
    Temp=Temp/dble(N)


end subroutine MvtBrownien_3

subroutine MvtBrownien_3_CL(N,dt,Temp)
    integer               :: k,N,i, CL
    real(dp)             :: G(2),dt,Temp,xn,yn,distm,distp,u

    Temp=0.d0

    do i=1,N
        xn = random(); yn = random();
        do while(dedans(xn,yn) == 1)
            G=Gaussien(random(),random())
            ! distance min avec les anciennes valeurs
            distm = min(abs(xn),abs(1._dp-xn),abs(yn),abs(1._dp-yn))
            ! Tirage des nouvelles valeurs
            xn=xn+sqrt(dt)*G(1); yn=yn+sqrt(dt)*G(2)
            ! distance min avec les nouvelles valeurs
            distp = min(abs(xn),abs(1._dp-xn),abs(yn),abs(1._dp-yn))
            u = random()

            if(yn >= 1._dp) then ! Sortie Nord
                Temp = Temp + xn**2._dp - yn**2._dp
            else if ( u < exp((-distm*distp)/(2*dt))) then 
                if (dp == abs(1._dp-yn)) Temp = Temp + xn**2._dp - yn**2._dp
            endif
        enddo
    enddo
		
    Temp=Temp/dble(N)


end subroutine MvtBrownien_3_CL

!--------------------------------------!
real(dp) function exacte(x0,y0)
	real(dp), parameter  :: PI = 4._dp*DATAN(1.d0)
	real(dp),intent(in)  :: x0,y0
    integer              :: i,N
	N =50
	exacte = 0._dp
    do concurrent(i=0:N)
        exacte=exacte+(400._dp/pi)*( (sin((2*i+1)*pi*x0)*sinh((2*i+1)*pi*y0))/((2*i+1)*sinh((2*i+1)*pi)))
    enddo
end function exacte

real(dp) function exacte2(x0,y0)
	real(dp) :: x0,y0
    integer              :: i
	exacte2 = x0**2 - y0**2
end function exacte2

end module marches