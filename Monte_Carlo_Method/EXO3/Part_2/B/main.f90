program MN_3
    use, intrinsic :: iso_fortran_env    
    use ogpf ! module Traçage 
    use generateur_aleatoire
    implicit none
    type(gpf)  :: gp
    real(dp) :: t1,t2
! ----------------------
    integer :: N,i,d,choix,plot,Ne
    real(wp), allocatable :: Ntab(:), integrale(:),attendu(:),vect(:,:)
    call input(choix)
    call cpu_time(t1)
! ----- EXERCICE 3 -----   
    d = 5
!  --- Choix de tracer les points générés pour une N donnée  0 - Off , 1 - On  !! à changer pour tracer les points
    plot = 0
    if (plot == 1) d = 2  
! ----------------------------------------------------------------------------
    allocate(Ntab(d),integrale(d),attendu(d))
    Ntab = 0_wp; integrale=0_wp; attendu = 0_wp;

    N = 0  ! Nombre de points
    Ne = 20 ! Nombre d'étapes
    do i = 1,d
        N = 10**i
        Ntab(i) = N
        allocate(vect(2,N))

        select case(choix)
        case(0)
        call Kohonen(N,Ne,vect,plot,integrale(i))
        case(1)
        call kmeans(N,vect,Ne,plot,integrale(i))
        case(2)
        call kmeansquad(N,vect,Ne,plot,integrale(i))
        end select
        attendu(i) = exp(1._dp)
        deallocate(vect)

        select case(choix)
        case(0)
        print *, "Avec ",N, "points et ", Ne ," etapes, l'integrale calculee vaut: ",integrale(i),&
              ". L'erreur est de ", abs(integrale(i)-attendu(i));
        case(1,2)
        print *, "Avec ",N, "points, l'integrale calculee vaut: ",integrale(i),&
            ". L'erreur est de ", abs(integrale(i)-attendu(i));
        end select
    end do
    call cpu_time(t2)
    print *,"Temps d'execution: ", t2-t1
                
    deallocate(Ntab,integrale,attendu)


contains

! ---------------------------------
! Kouhonen
! ---------------------------------
subroutine Kohonen(N,Ne,vect,plot,integrale)
integer   :: N,i,plot,j0, Ne
real(dp)  :: ei, xa,ya,sei, seicarre,integrale
real(dp)  :: vect(2,N)

vect = 0_dp

! Placement de points aleatoires
vect(1,:) = (/ ( G1(), i=1,N ) /) ! abscisse x  
vect(2,:) = (/ ( G2(), i=1,N ) /) ! ordonnee y 
! print *,vect(1,:),vect(2,:)

! Parcours de proche en proche en Ne etape
sei = 0
seicarre = 0
do i = 1,Ne
    ei = 0.01 - ( dble(i)/dble(Ne) )*(0.01  - 0.0001 )
    sei = sei + ei 
    seicarre = seicarre + ei**2
    xa = G1(); ya = G2();
    ! Indice de la source la plus proche
    j0 = minloc(sqrt( (xa-vect(1,:))**2 + (ya-vect(2,:))**2 ),1)
    ! Rapprochement de cette indice
    vect(1,j0) = ( 1-ei )*vect(1,j0) + ei*xa
    vect(2,j0) = ( 1-ei )*vect(2,j0) + ei*ya
end do

print *, "somme epsilon: ", sei
print *, "somme epsilon carre: ", seicarre

if (plot == 1) then
    call gp%axis([-3._wp, 3._wp, -3._wp, 3._wp])
    call gp%options('set key off')
    call gp%options('set autoscale fix')
    call gp%plot(vect(1,:),vect(2,:), 'title "Points Kouhonen" with points pt 7 ps 1.5 lc rgb "#011582"')
end if

integrale = sum(exp(vect(1,:)+vect(2,:)))/N

end subroutine Kohonen


subroutine kmeans(N,vect,Ne,plot,integrale)
    integer   :: N,i,j,plot,j0,r,Ne
    real(dp)  :: xa,ya,sx,sy,cpt,integrale
    real(dp)  :: vect(2,N), temp(3,N)

    vect = 0_dp; temp = 0_dp; 
    call gp%axis([-3._wp, 3._wp, -3._wp, 3._wp])
    call gp%options('set key off')
    call gp%options('set autoscale fix')

    vect(1,:) = (/ ( G1() , i=1,N ) /) ! abscisse x  
    vect(2,:) = (/ ( G2() , i=1,N ) /) ! ordonnee y 
    
    do r = 1,Ne
        do i = 1,N
            xa = G1();ya = G2();
            j0 = minloc(sqrt( (vect(1,:)-xa)**2 + (vect(2,:)-ya)**2 ),1)
            temp(1,i) = xa
            temp(2,i) = ya
            temp(3,i) = j0
            cpt=0_dp;sx=0_dp;sy=0_dp;
            do concurrent (j = 1:N)
                if ( temp(3,j) == i) then 
                    sx = sx + temp(1,j)
                    sy = sy + temp(2,j)
                    cpt = cpt + 1_dp
                end if
            end do 
            if (cpt /= 0) then
                vect(1,i) = sx/cpt
                vect(2,i) = sy/cpt 
            end if 
        end do
    end do

if (plot == 1) call gp%plot(vect(1,:),vect(2,:), 'title "Points K-means" with points pt 7 ps 1.5 lc rgb "#011582"')


integrale = sum(exp(vect(1,:)+vect(2,:)))/N

end subroutine kmeans



subroutine kmeansquad(N,vect,Ne,plot,integrale)
    integer   :: N,i,j,plot,j0,r,Ne
    real(dp)  :: xa,ya,sx,sy,cpt,integrale
    real(dp)  :: vect(2,N), temp(4,N)

    vect = 0_dp; temp = 0_dp; 
    !initialisaiton
    vect(1,:) = (/ ( G1(), i=1,N ) /) ! abscisse x  
    vect(2,:) = (/ ( G2(), i=1,N ) /) ! ordonnee y 
    
    do r = 1,Ne
        do i = 1,N
            xa = G1();ya = G2();
            j0 = minloc(sqrt( (vect(1,:)-xa)**2 + (vect(2,:)-ya)**2 ),1)
            temp(1,i) = xa
            temp(2,i) = ya
            temp(3,i) = j0
            cpt=0_dp;sx=0_dp;sy=0_dp;
            do concurrent (j = 1:N)
                if ( temp(3,j) == i) then 
                    sx = sx + temp(1,j)
                    sy = sy + temp(2,j)
                    cpt = cpt + 1_dp
                end if
            end do 
            temp(4,i) = temp(i,4) + cpt
            if (cpt /= 0) then
                vect(1,i) = sx/cpt
                vect(2,i) = sy/cpt 
            end if 
        end do
    end do

if (plot == 1) call gp%plot(vect(1,:),vect(2,:), 'title "Points K-means quadratique" with points pt 7 ps 1.5 lc rgb "#011582"')


integrale = sum(temp(4,:)*exp(vect(1,:)+vect(2,:)))/sum(temp(4,:))

end subroutine kmeansquad


subroutine input(choix)
    integer :: choix
    write(*,*) "Choix possibles (0-1): "  
    write(*,*) "(0) Algorithme de Kohonen"
    write(*,*) "(1) Algorithme K-means"
    write(*,*) "(2) Algorithme K-means quadratique"
    write(*,*)
    write(*,'(a7)', advance='no') "Choix: "
    read(*,*)  choix
end subroutine input


end program MN_3
