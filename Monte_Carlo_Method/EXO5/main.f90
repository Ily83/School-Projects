program MN_3
    use marches
    use, intrinsic :: iso_fortran_env    
    use ogpf ! module Traçage 
    use generateur_aleatoire
    implicit none
    type(gpf)  :: gp
    real(dp)   :: t1,t2
! ----------------------
    integer    :: N,Nit, i,plot,choix
	real(dp)   :: Temp,x0,y0,dt, d(10)
    ! real(wp), allocatable :: u(:,:)
    ! call input(choix,x0,y0)
    call cpu_time(t1)
    choix = 4
	! x0 = 0.5_dp; y0 = 1.00_dp;
    x0 = 0.5_dp; y0 = 0.5_dp;

	N = 100000;
	dt = 1.d-4;
    Select Case(choix)
    case(0)
    call MvtBrownien_1(N,x0,y0,dt,Temp)
    case(1)
    call MvtBrownien_1_CL(N,x0,y0,dt,Temp)
    case(2)
    call MvtBrownien_2(N,x0,y0,dt,Temp)
    case(3)
    call MvtBrownien_2_CL(N,x0,y0,dt,Temp)
    case(4)
    call MarcheCercle(N,x0,y0,1.d-5,Temp)
    case(5)
    call MarcheCercleCL(N,x0,y0,1.d-7,Temp)
    case(6)
    ! call MvtBrownien_3(N,dt,Temp)
    call MvtBrownien_3_CL(N,dt,Temp)
    end select
    Select Case(choix)
    case(0,2,4)
    print *, "Temperature:        ", Temp
    print *, "Temperature exacte: ", exacte(x0,y0)
    print *, "   Erreur absolu  : ",abs(Temp-exacte(x0,y0))
    case(1,3,5)
    print *, "Temperature:        ", Temp
    print *, "Temperature exacte: ", exacte2(x0,y0)
    print *, "   Erreur absolu  : ",abs(Temp-exacte2(x0,y0))
    case(6)
    print *, "Temperature Moyenne: ", Temp
    end select
    call cpu_time(t2)
    print *, "Temps d'éxecution: ",t2-t1
    contains




subroutine input(choix,x0,y0)
    integer  :: choix
    real(dp) :: x0,y0
    write(*,*) "Choix possibles (0-1): "  
    write(*,*) "(0) Mouvement Brownien Methode 1"
    write(*,*) "(1) Mouvement Brownien Methode 1 avec CL"
    write(*,*) "(2) Mouvement Brownien Methode 2"
    write(*,*) "(3) Mouvement Brownien Methode 2 avec CL"
    write(*,*) "(4) Marche Sur Le Cercle"
    write(*,*) "(5) Marche Sur Le Cercle avec CL"
    write(*,*) "(6) Mouvement Brownien Methode 2 avec position initiale aleatoire"
    write(*,*)
    write(*,'(a7)', advance='no') "Choix: "
    read(*,*)  choix
    write(*,*)
    write(*,*) "Entrer la position initiale  x0 y0 entre 0 et 1 (separe d'un espace): "
    read(*,*) x0, y0
end subroutine input

end program MN_3