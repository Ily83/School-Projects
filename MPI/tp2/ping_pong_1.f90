!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ping_pong_1.f90 --- TP2 : Communications point ? point :
!!                           envoi d'un message du processus 0 au processus 1
!! 
!! Auteur          : Denis GIROU (CNRS/IDRIS - France) <Denis.Girou@idris.fr>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program ping_pong_1
  USE MPI
  implicit none

  integer, dimension(MPI_STATUS_SIZE)  :: statut
  integer, parameter                   :: nb_valeurs=1000,etiquette=99
  integer                              :: rang,code
  integer, parameter                   :: dp = kind(1.d0)
  real(kind=dp), dimension(nb_valeurs) :: valeurs

  ! Initialisation
  CALL MPI_INIT(code)

  ! Processus
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
  
!......................................................................
	if ( mod(rang,2) == 0 ) then  ! Cas pair (0)
      call random_number(valeurs)
      ! Envoi
	  call MPI_SEND(valeurs,nb_valeurs,MPI_REAL8,1,etiquette,MPI_COMM_WORLD,code)
      write(*,'(a7,I2,a4,I2)') "Ping ! ", rang, " -> ", rang+1
!......................................................................
    else if (rang == 1) then  ! cas impair (1)
	  ! Reception
      call MPI_RECV(valeurs,nb_valeurs,MPI_REAL8,0,etiquette,MPI_COMM_WORLD,statut,code)
      ! write(*,'(a7,I2,a4,I2)') "Pong ! ", rang, " <- ", rang-1
      write(*,*)
      print ('("Moi, processus 1, j''ai recu ",i4," valeurs (derniere = ", &
      & f4.2,") du processus 0.")'), nb_valeurs,valeurs(nb_valeurs)
!......................................................................
    end if
	
    ! Fin
    call MPI_FINALIZE(code)
end program ping_pong_1