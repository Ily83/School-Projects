PROGRAM lire_tab

  USE mpi

  IMPLICIT NONE

  INTEGER, DIMENSION(MPI_STATUS_SIZE)            :: statut
  INTEGER                                        :: rang, code, descripteur, ntx, nty
  integer, parameter                             :: dp = kind(1.d0)
  REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)    :: u_lu
  INTEGER(KIND=MPI_OFFSET_KIND)                  :: taille_fichier
  INTEGER                                        :: taille_reel
  
  CALL MPI_INIT(code)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)

  OPEN(11, FILE='poisson.data', STATUS='OLD')
  READ(11, *) ntx
  READ(11, *) nty
  CLOSE(11)

  ALLOCATE(u_lu(ntx, nty))
  u_lu(:, :) = 0.d0

  !Ouverture du fichier "donnees.dat" en écriture
  CALL MPI_FILE_OPEN(MPI_COMM_WORLD, "donnees.dat", &
       MPI_MODE_RDONLY, &
       MPI_INFO_NULL, descripteur, code)

  !Test pour savoir si ouverture du fichier est correcte
  IF (code /= MPI_SUCCESS) THEN
     PRINT *, 'ATTENTION erreur lors ouverture du fichier'
     CALL MPI_ABORT(MPI_COMM_WORLD, 2, code)
     CALL MPI_FINALIZE(code)
  END IF

  CALL MPI_FILE_GET_SIZE(descripteur, taille_fichier, code)
  CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, taille_reel, code)
  if (taille_fichier /= ntx*nty*taille_reel) then
    print *, " ATTENTION Donnees.dat n'a pas la bonne longueur (",taille_fichier,",",&
            ntx*nty*taille_reel,")"
    write(11,*) 0
  else
    CALL MPI_FILE_READ(descripteur, u_lu, SIZE(u_lu), &
                       MPI_DOUBLE_PRECISION, statut, code)
    WRITE(11, 101)  u_lu
101 FORMAT (E19.12)
  end if

  CALL MPI_FILE_CLOSE(descripteur, code)
  
  CALL MPI_FINALIZE(code)

END PROGRAM lire_tab


