!module_params.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE  type_params
  IMPLICIT NONE

!! Initialisation des parametres
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Nombre de points en et y
  INTEGER                                   :: ntx, nty
  !Indices de debut et de fin dans chaque direction
  INTEGER                                   :: sx, ex, sy, ey
  !Nombre iterations maximum en temps
  INTEGER, PARAMETER                        :: it_max=100000

  integer, parameter                   :: dp = kind(1.d0)
  !Argument muet de la fonction F90 EPSILON
  REAL(kind=dp), PARAMETER                   :: eps=EPSILON(1._8)

END MODULE  type_params
