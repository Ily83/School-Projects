program MN_1
    use airevolume
    use ogpf ! module Traçage 
    use, intrinsic :: iso_fortran_env    
    implicit none
! ----------------------
    integer :: N,i,plot,choix
    real(wp), allocatable :: Ntab(:), Errtab(:)
    write(*,'(a44)', advance='no') "Choix: (Cercle = 1, Ellipse =2, sphere =3): "
    read(*,*)  choix
    call cpu_time(t1)
! ----- EXERCICE 1 -----   
    R = 1.0
    ! plot = 0 ! off
    plot = 1 ! on
    N = 100 ! 10^3
    allocate(Ntab(4),Errtab(4))
    ntab = 0_wp; errtab=0_wp;
    a = 2._dp; b = 5._dp; x0=1._dp; y0=2._dp;z0=0._dp;
    do i = 1,4
        N = N*10
        Ntab(i) = N
        if (N>100000) plot = 0 ! Consome trop de RAM sinon!
        select case(choix)
        case(1)
        call cercle(N,x0,y0,R,aire,plot)
        Errtab(i) = abs(aire_cercle(R)-aire)
        print *, "Apres ", N, "iterations, l'aire du cercle calculee vaut: ", & 
        aire, ". L'erreur est de ", Errtab(i)
        case(2)
        call ellipse(N,x0,y0,a,b,aire,plot)
        Errtab(i) = abs(aire_ellipse(a,b)-aire)
        print *, "Apres ", N, "iterations, l'aire de l'ellipse calculee vaut: ", aire, & 
              ". L'erreur est de ", Errtab(i)
        case(3)
        call sphere(N,x0,y0,z0,R,aire,plot)
        errtab(i) = abs(aire_sphere(R)-aire)
        print *, "Apres ", N, "iterations, le volume de la sphere calculee vaut: ", aire, & 
              ". L'erreur est de ", Errtab(i)
        end select
    end do
    call cpu_time(t2)
    print *,"Temps d'execution: ", t2-t1
    call gp%title("Erreur en fonction du nombre d'iterations")
    call gp%options('set key off')
    call gp%options('set key noautotitle')
    call gp%axis([8.5_wp,17._wp,-9._wp, -1.5_wp])
    ! print*,Ntab
    ! print *,errtab
    call gp%plot(log(Ntab), log(errtab), 'title "Erreur" w lines lc "red" lw 2 lt 2 dt 2','',&
                log(Ntab),-0.5*log(Ntab),'title "1/sqrt(N)" w lines lc "blue" lw 3', '')
    deallocate(Ntab,errtab)
end program MN_1