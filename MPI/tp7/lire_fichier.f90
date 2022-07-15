!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! lire_fichier.f90 --- T.P. 7 du cours MPI (à exécuter sur 4 processus)
!! 
!! Auteur          : Denis GIROU (CNRS/IDRIS - France) <Denis.Girou@idris.fr>
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program lire_fichier

  use mpi

  implicit none

  integer, parameter                  :: nb_valeurs=121
  INTEGER(KIND=MPI_OFFSET_KIND)       :: position_fichier
  integer, dimension(nb_valeurs)      :: valeurs
  integer, dimension(MPI_STATUS_SIZE) :: statut
  integer                             :: rang,nb_octets_entier,descripteur,code

  call MPI_INIT(code)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)

  ! Ouverture du fichier "donnees.dat" en lecture
  call MPI_FILE_OPEN(MPI_COMM_WORLD,"donnees.dat",MPI_MODE_RDONLY,MPI_INFO_NULL,descripteur,code)

  !Test pour verifier que le fichier a bien pu etre ouvert
   if (code == MPI_SUCCESS) write(*,*) "Le fichier a bien pu etre ouvert"; 

   ! Lecture de nb_valeurs valeurs sur chacun des processus
  valeurs(:)=0

  ! Lecture via des déplacements explicites, en mode individuel
  CALL MPI_TYPE_SIZE(MPI_INTEGER,nb_octets_entier,code)
  position_fichier = rang*nb_valeurs*nb_octets_entier
  CALL MPI_FILE_READ_AT(descripteur,position_fichier,valeurs,nb_valeurs,MPI_INTEGER,statut,code)
  
  open(unit=45,file="fichier_dei"//achar(48+rang)//".dat")
  write(unit=45,fmt='(I3)') valeurs(:)
  close(unit=45)

  ! Lecture via les pointeurs partagés, en mode collectif
  valeurs(:)=0
  CALL MPI_FILE_READ_ORDERED(descripteur,valeurs,nb_valeurs,MPI_INTEGER,statut,code)

  open(unit=45,file="fichier_ppc"//achar(48+rang)//".dat")
  write(unit=45,fmt='(I3)') valeurs(:)
  close(unit=45)

  ! Lecture via les pointeurs individuels, en mode individuel
  valeurs = 0
  position_fichier = rang*nb_valeurs*nb_octets_entier
  CALL MPI_FILE_SEEK(descripteur,position_fichier,MPI_SEEK_SET,code)
  CALL MPI_FILE_READ(descripteur,valeurs,nb_valeurs,MPI_INTEGER,statut,code)
  

  open(unit=45,file="fichier_pii"//achar(48+rang)//".dat")
  write(unit=45,fmt='(I3)') valeurs(:)
  close(unit=45)

  ! Lecture via les pointeurs partagés, en mode individuel
  ! (on doit tout d'abord repositionner le pointeur partagé au début du fichier)
  valeurs = 0
  position_fichier = 0
  CALL MPI_FILE_SEEK_SHARED(descripteur,position_fichier,MPI_SEEK_SET,code)
  CALL MPI_FILE_READ_SHARED(descripteur, valeurs, nb_valeurs, MPI_INTEGER,MPI_STATUS_IGNORE,code)

  open(unit=45,file="fichier_ppi"//achar(48+rang)//".dat")
  write(unit=45,fmt='(I3)') valeurs(:)
  close(unit=45)

  ! Fermeture du fichier
  CALL MPI_FILE_CLOSE(descripteur,code)

  ! Fin MPI
  call MPI_FINALIZE(code)

end program lire_fichier
