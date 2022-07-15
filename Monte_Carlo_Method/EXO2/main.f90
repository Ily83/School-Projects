program MN_2
    use Integrales
    ! use ogpf ! module Traçage
    implicit none 
    ! type(gpf)  :: gp
    
! ----------------------
    integer               :: N,i,d,choix
    real(wp), allocatable :: Ntab(:), integrale(:),attendu(:),supe(:),infe(:)
    real(dp)              :: t1,t2
    call input(choix)
    call cpu_time(t1)
! ----- EXERCICE 2 -----   
    d = 5
    allocate(Ntab(d),integrale(d),attendu(d),supe(d),infe(d))
    Ntab = 0_wp; integrale=0_wp; attendu = 0_wp;
    N = 100
    do i = 1,d
        N = N*10
        Ntab(i) = N 
        ! call exemple1D(N) ! exemple tracé
        select case(choix)
        !! Somme 
        case(0)
        call integrals(N,integrale(i),ecart)
        attendu(i) = exacte("s")   
        !! Produit
        case(1)
        call integralp(N,integrale(i),ecart)
        attendu(i) = exacte("p") 
        !! produit controle 1
        case(2)
        call intcontrolef(N,integrale(i),var,ecart)
        attendu(i) = exacte("p") 
        ! produit controle 2
        case(3)
        call intcontroleG(N,integrale(i),var,ecart)
        attendu(i) = exacte("p")
        end select
        infe(i) = integrale(i) - ecart
        supe(i) = integrale(i) + ecart  
        !! Erreur
        print *, "Apres ",N, "iterations, l'integrale calculee vaut: ",integrale(i),"."
        print *, "L'erreur est de ", abs(integrale(i)-attendu(i));
        ! print *,"intervale de confiance: [",integrale(i)-ecart,",",integrale(i)+ecart,"]"
        print *,  " ========================================================================================  "
    end do
    call cpu_time(t2)
    print *,"Temps d'execution: ", t2-t1

    call gp%title("Resultat de l'integrale en fonction du nombre d'iterations N.")
    call gp%options('set key bottom right')
    call gp%options('set autoscale')
    select case(choix)
    case(0)
    call gp%axis([50._wp, 2000000._wp, 8.710_wp, 8.73_wp])
    case(1:3)
    call gp%axis([-250000._wp, 4500000._wp, 1.0689_wp, 1.0698_wp])
    end select
    ! print*,log(Ntab)
    ! print *,log(integrale(i)-attendu(i))
    ! print *,-0.5*log(Ntab)
    ! print *,integrale
    call gp%plot(Ntab, integrale, 'title "Int° Calcule" w lines lc "red" lw 2 lt 2 dt 2','',&
                Ntab,attendu,'title "Attendu" w lines lc "blue" lw 3','',&
                Ntab,infe,'title "infe" w lines lc "green" lw 2 lt 2 dt 2','',&
                Ntab,supe,'title "supe" w lines lc "brown" lw 2 lt 2 dt 2')
    deallocate(Ntab,integrale,attendu)
end program MN_2