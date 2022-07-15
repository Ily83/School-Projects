program MN_3
    use, intrinsic :: iso_fortran_env    
    use ogpf ! module Traçage 
    use generateur_aleatoire
    implicit none
    type(gpf)  :: gp
    real(dp) :: t1,t2
! ----------------------
    integer :: N,i,d,choix,plot
    real(wp), allocatable :: Ntab(:), integrale(:),attendu(:),vect(:,:)
    call input(choix)
    call cpu_time(t1)
! ----- EXERCICE 3 -----   
    d = 7
    ! d = 5
    ! d = 1000
!  --- Choix de tracer les points générés pour une N donnée  0 - Off , 1 - On  !! à changer pour tracer les points
    plot = 1
    if (plot == 1) d = 2  
! ----------------------------------------------------------------------------
    allocate(Ntab(d),integrale(d),attendu(d))
    Ntab = 0_wp; integrale=0_wp; attendu = 0_wp;
    N = 0
    if (plot == 1) N = 1  
    do i = 1,d
        print *,i,N
        N = 10**i
        if (plot == 1) N = N*10 ! 10 & 100 points
        Ntab(i) = N
        allocate(vect(2,N))

        select case(choix)
        case(0)
        call MC(N,vect,plot)
        case(1)
        call qmc(N,vect,plot)
        end select

        integrale(i) = sum(exp(vect(1,:)+vect(2,:)))/N
        deallocate(vect)

        attendu(i) = (exp(1._dp)-1._dp)**2._dp
        print *, "Apres ",N, "iterations, l'integrale calculee vaut: ",integrale(i),&
              ". L'erreur est de ", abs(integrale(i)-attendu(i));
    end do
    call cpu_time(t2)
    print *,"Temps d'execution: ", t2-t1
    call gp%title("Valeur de l'integrale en fonction du nombre d'iterations")
    call gp%options('set key top right')
    call gp%options('set autoscale fix')
    select case(choix)
    case(1) ! QuasiMC
    call gp%axis([0._wp, 1000001._wp, 2.952_wp, 2.954_wp])
    case(0) ! MC
    call gp%axis([0._wp, 1000001._wp, 2.94_wp, 3.0_wp])
    end select

    call gp%plot(Ntab, integrale, 'title "Int° Calcule" w lines lc "red" lw 2 lt 2 dt 2','',&
                Ntab,attendu,'title "Attendu" w lines lc "blue" lw 3', '')

                
    deallocate(Ntab,integrale,attendu)


contains


!--------------------------------------
! Methode Monte Carlo classique (vectoriel)
!--------------------------------------

subroutine MC(N,vect,plot)
    integer   :: N,i, plot 
    real(dp)  :: vect(2,N)

    call gp%axis([0._wp, 1._wp, 0._wp, 1._wp])
    call gp%options('set key off')
    call gp%options('set autoscale fix')
    call gp%options('set style line 1 lc rgb "red" lt 1 lw 2 pt 6 ps 1.5')

    vect(1,:) = (/ ( random(), i=1,N ) /)
    vect(2,:) = (/ ( random(), i=1,N ) /)

if (plot == 1) call gp%plot(vect(1,:),vect(2,:), 'title "Points Monte Carlo" with points pt 7 ps 1.2 lc rgb "#011582"')

end subroutine

!--------------------------------------
! Methode Quasi-Monte Carlo 
!--------------------------------------

subroutine QMC(N,vect,plot)
    integer   :: N, i, plot
    real(dp)  :: vect(2,N)

    call gp%axis([0._wp, 1._wp, 0._wp, 1._wp])
    call gp%options('set key off')
    call gp%options('set autoscale fix')
    call gp%options('set style line 1 lc rgb "red" lt 1 lw 2 pt 6 ps 1.5')

    vect(1,:) = (/( i*(2.d0**(0.5d0))-int(i*(2.d0**(0.5d0))), i =1,N )/)
    vect(2,:) = (/( i*(3.d0**(0.5d0))-int(i*(3.d0**(0.5d0))), i =1,N )/)

    if (plot == 1) call gp%plot(vect(1,:),vect(2,:),'title "Points Quasi-Monte Carlo" with points pt 7 ps 1.2 lc rgb "#011582"')
end subroutine QMC

subroutine input(choix)
    integer :: choix
    write(*,*) "Choix possibles (0-1): "  
    write(*,*) "(0) Monte Carlo Classique"
    write(*,*) "(1) Quasi-Monte Carlo"
    write(*,*)
    write(*,'(a7)', advance='no') "Choix: "
    read(*,*)  choix
end subroutine input


end program MN_3