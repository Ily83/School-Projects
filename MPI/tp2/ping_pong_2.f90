!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ping_pong_2.f90 --- TP2 : Communications point ? point : ping-pong
!! 
!! Auteur          : Denis GIROU (CNRS/IDRIS - France) <Denis.Girou@idris.fr>
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program ping_pong_2
  USE MPI
  implicit none

  integer, dimension(MPI_STATUS_SIZE)  :: statut
  integer, parameter                   :: nb_valeurs=1000,etiquette=99
  integer                              :: rang,code
  integer, parameter                   :: dp = kind(1.d0)
  real(kind=dp)                        :: temps_debut,temps_fin
  real(kind=dp), dimension(nb_valeurs) :: valeurs

  ! Initialisation
  CALL MPI_INIT(code)

  ! Processus
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)

!......................................................................

  if (rang == 0) then
    call random_number(valeurs)
    temps_debut=MPI_WTIME()
    call MPI_SEND(valeurs,nb_valeurs,MPI_REAL8,1,etiquette,MPI_COMM_WORLD,code)
    write(*,'(a7,I2,a4,I2)') "Ping ! ", rang, " -> ", rang+1
    call MPI_RECV(valeurs,nb_valeurs,MPI_REAL8,1,etiquette,MPI_COMM_WORLD,statut,code)
    temps_fin=MPI_WTIME()
    print ('("Moi, processus 0, j''ai envoye et recu ",i5, &
        & " valeurs (derniere = ",f4.2,") du processus 1", &
        & " en ",f8.6," secondes.")'), &
          nb_valeurs,valeurs(nb_valeurs),temps_fin-temps_debut
!......................................................................
  else if (rang == 1) then 
    call MPI_RECV(valeurs,nb_valeurs,MPI_REAL8,0,etiquette,MPI_COMM_WORLD,statut,code)
    temps_fin=MPI_WTIME()
    call MPI_SEND(valeurs,nb_valeurs,MPI_REAL8,0,etiquette,MPI_COMM_WORLD,code)
    write(*,'(a7I2,a4,I2)') "Pong ! ", rang, " <- ", rang-1
!......................................................................  
   end if

!......................................................................
   call MPI_FINALIZE(code)
   
end program ping_pong_2
