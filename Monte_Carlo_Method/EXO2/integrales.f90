module Integrales
    use generateur_aleatoire   
    use ogpf ! module Traçage
    implicit none 
    type(gpf)  :: gp
    real(dp)   :: ecart,var

contains

! Exemple 1D
subroutine exemple1D(N)
    integer    :: N,i
    real(dp)   :: x,somme,xtab(N),ytab(N)
    Somme = 0._dp
    call gp%axis([0._wp, 1._wp,0._wp,1._wp])
    call gp%options('set key top left')
    call gp%options('set autoscale fix')
    call gp%options('set style line 1 lc rgb "blue" lt 1 lw 2 pt 6 ps 1.5')
    do i = 1,N
        ! Generation des nombres aléatoires
        x = random(); 
        Somme=Somme+sin(x)
        xtab(i) = x 
        ytab(i) = sin(x)
    end do
    print *, Somme/dble(N)
    call gp%plot(xtab,ytab,'title "x y" with points pt 7 ps 1.2 lc rgb "#ad2060"')
end subroutine exemple1D   

subroutine integrals(N,somme,ecart)
    integer    :: N,i
    real(dp)   :: x,y,z,t,nv, var, ecart,somme
    nv = 0._dp
    Somme=0._dp
    do i = 1,N
        ! Generation des nombres aléatoires
        x = random(); y = random()
        z = random(); t = random()

        Somme=Somme+fs(x,y,z,t)
        nv = nv + fs(x,y,z,t)*fs(x,y,z,t)
    end do
    Somme=Somme/real(N)
    var = sqrt(nv/real(N)-Somme*Somme)
    ecart = sqrt(var/real(N))*1.96
    ! print *,"intervale de confiance: [",somme-ecart,",",somme+ecart,"]"
end subroutine integrals

subroutine integralp(N,produit,ecart)
    integer    :: N,i
    real(dp)   :: x,y,z,t, nv,var,ecart,produit
    produit=0._dp
    nv = 0._dp
    do i = 1,n
        ! Generation des nombres aléatoires
        x = random(); y = random()
        z = random(); t = random()
        produit=produit+fp(x,y,z,t)
        nv = nv + fp(x,y,z,t)*fp(x,y,z,t)
    end do
    produit=produit/real(N)
    var = nv/real(N)-produit*produit
    print *, "Variance: ", var
    ecart = sqrt(var/real(N))*1.96
    ! print *,"intervale de confiance: [",produit-ecart,",",produit+ecart,"]"
end subroutine integralp

!---------------------
! Integrales avec Controle
!---------------------

subroutine intcontrolef(N,produit,var,ecart)
    implicit none
    integer        :: N,i
    real(dp)   :: x,y,z,t, nv,var,ecart,produit
    

    produit = 0._dp
    do i = 1,n
        ! Generation des nombres aléatoires
        x = random(); y = random()
        z = random(); t = random()
        
        produit = produit + fp(x,y,z,t) - h(x,y,z,t)    ! Controle avec l'approximation à l'ordre 1 xyzt
        nv = nv + (fp(x,y,z,t) - h(x,y,z,t))**2         ! Pour le calcul de la variance
    enddo

    ! var=(nv-produit)/(dble(N)-1._dp)
    var = nv/real(N)-(produit/N)**2
    produit=produit/real(N) + 1._dp/16._dp  ! moyenne
    print *, "Variance: ", var
    ecart = sqrt(var/real(N))*1.96
    ! print *, "ecart: ", ecart
    ! print *,"intervale de confiance: [",produit-ecart,",",produit+ecart,"]"
end subroutine  intcontrolef


subroutine intcontroleG(N,produit,var,ecart)
implicit none
integer        :: N,i
real(dp)   :: x,y,z,t, nv,var,ecart,produit

produit = 0._dp
do i = 1,n
    ! Generation des nombres aléatoires
    x = random(); y = random()
    z = random(); t = random()
    ! Controle avec l'approximation à l'ordre 2 xyzt + (xyzt)^2/2
    produit=produit+fp(x,y,z,t) - g(x,y,z,t) 
    ! Pour le calcul de la variance
    nv=nv+(fp(x,y,z,t)-g(x,y,z,t))**2         
enddo
var = nv/real(N)-(produit/N)**2
print *, "Variance: ", var
produit=produit/real(N) + 89._dp/1296._dp
ecart = sqrt(var/real(N))*1.96
! print *,"intervale de confiance: [",produit-ecart,",",produit+ecart,"]"
end subroutine  intcontroleG




!---------------------
! Fonctions
!---------------------

! Fonction exp produit
real(dp) function fp(x,y,z,t)
real(dp), intent(in) :: x,y,z,t
fp = exp(x*y*z*t)
end function fp

! Fonction exp somme
real(dp) function fs(x,y,z,t)
real(dp), intent(in) :: x,y,z,t
fs = exp(x+y+z+t)
end function fs

! Fonction de contrôle à l'ordre 1
real(dp) function h(x,y,z,t)
real(dp), intent(in) :: x,y,z,t
h = x*y*z*t
end function h

! Fonction de contrôle à l'ordre 2
real(dp) function g(x,y,z,t)
real(dp), intent(in) :: x,y,z,t
g = x*y*z*t+(x*y*z*t)**2._dp/2._dp
end function g

real(dp) function exacte(p)
character(1) :: p 
select case(p)
case("s")
exacte = (exp(1._dp)-1)**4
case("p")
exacte = 1.0693976088597708_dp  ! Calculé dans le programme python avec scipy
end select 
end function exacte

subroutine input(choix)
integer :: choix
write(*,*) "Choix possibles (0-3): "  
write(*,*) "(0) exp(x+y+z+t)"
write(*,*) "(1) exp(x*y*z*t)"
write(*,*) "(2) exp(x*y*z*t) controle au 1er ordre"
write(*,*) "(3) exp(x*y*z*t) controle au 2nd ordre"
write(*,*)
write(*,'(a7)', advance='no') "Choix: "
read(*,*)  choix
end subroutine input

end module Integrales