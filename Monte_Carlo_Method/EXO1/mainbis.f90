program MN_1_bis
    !!!!!!!!!!!!!
    ! Ce programme est utilisé pour avoir l'aire en fonction du nombre d'itérations.
    !!!!!!!!!!!!!
    use airevolume
    use ogpf ! module Traçage 
    use, intrinsic :: iso_fortran_env    
    implicit none
! ----------------------
    integer :: N,i,plot,choix
    real(wp), allocatable :: Ntab(:), airetab(:),attendu(:)
    write(*,'(a44)', advance='no') "Choix: (Cercle = 1, Ellipse =2, sphere =3): "
    read(*,*)  choix
    call cpu_time(t1)
! ----- EXERCICE 1 -----   
    R = 1
    plot = 0 ! off
    ! plot = 1 ! on
    N = 10
    allocate(Ntab(1000),airetab(1000),attendu(1000))
    ntab = 0_wp; airetab=0_wp;
    ! call gp%multiplot(1,1)
    do i = 1,999
        N = N+100
        a = 2.0; b = 5.0; x0=0.0; y0=0.0;z0=0.0;
        Ntab(i) = N
        select case(choix)
        case(1)
        call cercle(N,x0,y0,R,aire,plot)
        airetab(i) = aire
        attendu(i) = aire_cercle(R)
        case(2)
        call ellipse(N,x0,y0,a,b,aire,plot)
        airetab(i) = aire
        attendu(i) = aire_ellipse(a,b)
        case(3)
        call sphere(N,x0,y0,z0,R,aire,plot)
        airetab(i) = aire
        attendu(i) = aire_sphere(R)
        end select
    end do
    call cpu_time(t2)
    print *,"Temps d'execution: ", t2-t1

    call gp%title("L'aire/volume en fonction du nombre d'iterations")
    call gp%xlabel('x')
    call gp%ylabel('y')
    call gp%options('set key off')
    !!!!!!!
    !! à changer selon le choix manuellement si les valeurs sont changées
    !!!!!!
    select case(choix)
    case(1)
    call gp%axis([0._wp,50000._wp,2.9_wp, 3.4_wp])
    case(2)
    call gp%axis([0._wp,50000._wp,29._wp, 34._wp])
    case(3)
    call gp%axis([0._wp,50000._wp,3.8_wp, 4.5_wp])
    end select
    call gp%plot(Ntab,airetab,'title "Aire en fonction de N" w lines lc "blue" lw 2 lt 1','',&
         Ntab,attendu,'title "pi" w lines lc "red" lw 2 lt 2 dt 2')
    deallocate(Ntab,airetab)
! ----- EXERCICE 2 -----   
end program MN_1_bis