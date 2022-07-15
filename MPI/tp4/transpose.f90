!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! transpose.f90  --- Utilisation d'un type derive (type_transpose)
!!                    pour transposer une matrice.
!!
!!
!! Auteur          : Isabelle DUPAYS (CNRS/IDRIS - France)
!!                   <Isabelle.Dupays@idris.fr>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PROGRAM transpose
  USE MPI
  IMPLICIT NONE
  INTEGER, PARAMETER                     :: nb_lignes=5,nb_colonnes=4,&
                                            etiquette=1000
  INTEGER                                :: code,rang,type_ligne, nbr_procs,&
                                            type_transpose,taille_reel,i,j
  REAL, DIMENSION(nb_lignes,nb_colonnes) :: A
  REAL, DIMENSION(nb_colonnes,nb_lignes) :: AT
  INTEGER(kind=MPI_ADDRESS_KIND)         :: pas
  INTEGER, DIMENSION(MPI_STATUS_SIZE)    :: statut


  !Initialisation de MPI
  CALL MPI_INIT(code)

  !-- Savoir qui je suis
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nbr_procs,code)

  !-- Initialisation de la matrice AT
  AT(:,:) = 0.

  !Connaitre la taille du type de base MPI_REAL
  call MPI_TYPE_SIZE(MPI_REAL,taille_reel,code)
  call MPI_TYPE_VECTOR(nb_colonnes,1,nb_lignes,MPI_REAL,type_ligne,code)
  
  if(rang == 0) then
	write (* ,*) " un reel fait : " , taille_reel , " octets "
  end if
     !Construction du type derive type_transpose pour transposer la
  !matrice A composee de nb_lignes et de nb_colonnes   pas = 1*taille_real
   pas = 1*taille_reel

   call MPI_TYPE_CREATE_HVECTOR(nb_lignes,1,pas,type_ligne,type_transpose,code)

  !Validation du type cree type_transpose
  call MPI_TYPE_COMMIT(type_transpose,code)



  IF (rang == 0) THEN
    !Initialisation de la matrice A sur le processus 0
    A(:,:) = RESHAPE( (/ (i,i=1,nb_lignes*nb_colonnes) /), &
                      (/ nb_lignes,nb_colonnes /) )

    PRINT *,'Matrice A'
    DO i=1,nb_lignes
      PRINT *,A(i,:)
    END DO

    !Envoi de la matrice A au processus 1 avec le type type_transpose
	call MPI_SEND(A,1,type_transpose,1,100,MPI_COMM_WORLD,code)

  ELSE
    !Reception pour le processus 1 dans la matrice AT
	call MPI_RECV(AT,nb_colonnes*nb_lignes,MPI_REAL,0,100,MPI_COMM_WORLD,statut,code)

    PRINT *,'Matrice transposee AT'
    DO i=1,nb_colonnes
      PRINT *,AT(i,:)
    END DO

  END IF

  !Sortie de MPI
   call MPI_TYPE_FREE(type_transpose,code)
   call MPI_TYPE_FREE(type_ligne,code)
   call MPI_FINALIZE(code)

END PROGRAM transpose
