!module_parallel_mpi.f90
!!!!
!!subroutine initialisation_mpi
!!subroutine creation_topologie
!!subroutine domaine
!!subroutine voisinage
!!subroutine type_derive
!!subroutine communication
!!function erreur_globale
!!subroutine ecrire_mpi
!!subroutine finalisation_mpi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE parallel
  USE TYPE_PARAMS
  USE MPI
  IMPLICIT NONE

  !Rang du sous-domaine local
  INTEGER                                   :: rang
  !Nombre de processus
  INTEGER                                   :: nb_procs
  !Communicateur de la topologie cartesienne
  INTEGER                                   :: comm2d
  !Nombre de dimensions de la grille
  INTEGER, PARAMETER                        :: ndims = 2
  !Nombre de processus dans chaque dimension definissant la topologie
  INTEGER, DIMENSION(ndims)                 :: dims
  !Periodicite de  la topologie
  LOGICAL, DIMENSION(ndims)                 :: periodes
  ! Coordonnees du domaine locale dans la topologie cartesienne
  INTEGER, DIMENSION(ndims)                 :: coords
  ! Tableau definissant les voisins de chaque processus
  INTEGER, PARAMETER                        :: NB_VOISINS = 4
  INTEGER, PARAMETER                        :: N=1, E=2, S=3, O=4
  INTEGER, DIMENSION(NB_VOISINS)            :: voisin
  ! Types derives
  INTEGER                                   :: type_ligne, type_colonne
  !Constantes MPI
  INTEGER                                   :: code

CONTAINS

  SUBROUTINE initialisation_mpi
    !************
    !Initialisation pour chaque processus de son rang et du
    !nombre total de processus nb_procs
    !************

    !Initialisation de MPI
	Call MPI_INIT(code)

    !Savoir quel processus je suis et le nombre total de processus
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rang,code)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,code)

  END SUBROUTINE initialisation_mpi

  SUBROUTINE creation_topologie
    !************
    !Creation de la topologie cartesienne
    !************

    !Constantes MPI
    LOGICAL, PARAMETER                        :: reorganisation = .FALSE.

    ! Lecture du nombre de points ntx en x et nty en y
    OPEN(10, FILE='poisson.data', STATUS='OLD')
    READ(10, *) ntx
    READ(10, *) nty
    CLOSE(10)

    !Connaitre le nombre de processus selon x et le nombre de processus
    !selon y en fonction du nombre total de processus
     dims = 0
	 CALL MPI_DIMS_CREATE(nb_procs,ndims,dims,code)
    !Creation de la grille de processus 2D sans periodicite
	periodes = .FALSE.
	CALL MPI_CART_CREATE(MPI_COMM_WORLD,ndims,dims,periodes,reorganisation,comm2d,code)
	
	
    IF (rang == 0) THEN
      WRITE (*,'(A)') '-----------------------------------------'
      WRITE (*,'(A,i4,A)') 'Execution code poisson avec ', nb_procs, ' processus MPI'
      WRITE (*,'(A,i4,A,i4)') 'Taille du domaine : ntx=', ntx, ' nty=', nty
      WRITE (*,'(A,i4,A,i4,A)') 'Dimension de la topologie : ', &
                                 dims(1), ' suivant x, ', dims(2), ' suivant y'
      WRITE (*,'(A)') '-----------------------------------------'
    END IF
	
  END SUBROUTINE creation_topologie


  SUBROUTINE domaine
    !************
    !Calcul des coordonnées globales limites du sous domaine local
    !************

    ! Connaitre mes coordonnees dans la topologie
    CALL MPI_CART_COORDS(comm2d, rang, ndims, coords, code)


    !Calcul pour chaque processus de ses indices de debut et de fin suivant x
    sx = (coords(1)*ntx)/dims(1)+1
    ex = ((coords(1)+1)*ntx)/dims(1)

    !Calcul pour chaque processus de ses indices de debut et de fin suivant y
    sy = (coords(2)*nty)/dims(2)+1
    ey = ((coords(2)+1)*nty)/dims(2)

    WRITE (*,'(A,i4,A,i4,A,i4,A,i4,A,i4,A)') 'Rang dans la topologie : ', rang, &
         ' Indice des tableaux :', sx, ' a', ex, ' suivant x, ', &
         sy, ' a', ey, ' suivant y'

  END SUBROUTINE domaine

  SUBROUTINE voisinage
    !************
    !Calcul des processus voisins pour chaque processus
    !************
	 voisin = MPI_PROC_NULL
	! call MPI_CART_SHIFT(comm_nouveau, direction, pas, rang_precedent, rang_suivant, code)
    !Recherche des voisins Nord et Sud
	call MPI_CART_SHIFT(comm2d,0,1,voisin(N),voisin(S),code)
	
    !Recherche des voisins Ouest et Est
	call MPI_CART_SHIFT(comm2d,1,1,voisin(O),voisin(E),code)

    WRITE (*,'(A,i4,A,i4,A,i4,A,i4,A,i4)') "Processus ", rang, " a pour voisin : N", voisin(N), " E", voisin(E), &
            " S", voisin(S), " O", voisin(O)

  END SUBROUTINE voisinage


  SUBROUTINE type_derive
    !************
    !Creation des types derives type_ligne et type_colonne
    !************

    !Creation du type type_ligne pour echanger les points
    !au nord et au sud
	CALL MPI_TYPE_VECTOR(ey-sy+1,1,ex-sx+3,MPI_REAL8,type_ligne, code)
	CALL MPI_TYPE_COMMIT(type_ligne,code)


    !Creation du type type_colonne pour echanger
    !les points  a l'ouest et a l'est
	CALL MPI_TYPE_CONTIGUOUS(ex-sx+1,MPI_REAL8,type_colonne, code)
	CALL MPI_TYPE_COMMIT(type_colonne,code)

  END SUBROUTINE type_derive


  SUBROUTINE communication(u)
    !************
    !Echange des points aux interfaces
    !************

    REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :), INTENT(inout) :: u

    !Constantes MPI
    INTEGER, PARAMETER                   :: etiquette=100
    INTEGER, DIMENSION(MPI_STATUS_SIZE)  :: statut
	
	
	! MPI_SENDRECV(message_emis,longueur_message_emis,type_message_emis,rang_dest,etiq_message_emis, & 
	       ! message_recu, longueur_message_recu, type_message_recu, rang_source, etiq_message_recu, &
             ! comm, statut, code)

    !Envoi au voisin N et reception du voisin S
    CALL MPI_SENDRECV(u(sx,  sy), 1, type_ligne, voisin(N), etiquette, & 
				      u(ex+1,sy), 1, type_ligne, voisin(S), etiquette, &
					  comm2d, statut, code)

    !Envoi au voisin S et reception du voisin N
    CALL MPI_SENDRECV(u(ex,  sy), 1, type_ligne, voisin(S), etiquette, & 
				      u(sx-1,sy), 1, type_ligne, voisin(N), etiquette, &
					  comm2d, statut, code)

    !Envoi au voisin O et reception du voisin E 
    CALL MPI_SENDRECV(u(sx,  sy), 1, type_colonne, voisin(O), etiquette, & 
				      u(sx,ey+1), 1, type_colonne, voisin(E), etiquette, &
					  comm2d, statut, code)

    !Envoi au voisin E et reception du voisin O 
    CALL MPI_SENDRECV(u(sx,  ey), 1, type_colonne, voisin(E), etiquette, & 
				      u(sx, sy-1), 1, type_colonne, voisin(O), etiquette, &
					  comm2d, statut, code)

  END SUBROUTINE communication

  FUNCTION erreur_globale(u, u_nouveau)
    !************
    !Calcul de l'erreur globale (maximum des erreurs locales)
    !************
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :), INTENT(in) :: u, u_nouveau

    REAL(kind=dp)              :: erreur_globale, erreur_locale

    erreur_locale = MAXVAL(ABS(u(sx:ex, sy:ey) - u_nouveau(sx:ex, sy:ey)))

    !Calcul de l'erreur sur tous les sous-domaines
	! MPI MAX Recherche du maximum
    CALL MPI_AllREDUCE(erreur_locale,    &
	                   erreur_globale,   &
					   1,                &
					   MPI_REAL8,        &
					   MPI_MAX,          & 
					   comm2d,           &
					   code)


  END FUNCTION erreur_globale

  SUBROUTINE ecrire_mpi(u,nom_fichier)
    !********************
    ! Ecriture du tableau u a l'interieur d'un domaine pour chaque processus
    ! dans le fichier donnees.dat
    !********************
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :), INTENT(inout) :: u

    !Constantes MPI
    INTEGER, DIMENSION(MPI_STATUS_SIZE)   :: statut
    INTEGER                               :: descripteur
    INTEGER(kind = MPI_OFFSET_KIND)       :: deplacement_initial
    INTEGER, PARAMETER                    :: rang_tableau=2
    INTEGER, DIMENSION(rang_tableau)      :: profil_tab, profil_sous_tab, coord_debut
    INTEGER, DIMENSION(rang_tableau)      :: profil_tab_vue, profil_sous_tab_vue, coord_debut_vue
    INTEGER                               :: type_sous_tab, type_sous_tab_vue
	character(len = *), intent(in)        :: nom_fichier
    !Ouverture du fichier "donnees.dat" en écriture
	 CALL MPI_FILE_OPEN(comm2d,"donnees.dat",MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,descripteur,code)
    !Test pour verifier que le fichier a bien pu etre ouvert
    IF (code /= MPI_SUCCESS) THEN
      PRINT *, 'Erreur ouverture fichier'
      CALL MPI_ABORT(comm2d, 2, code)
    END IF

    !Définition de la vue sur le fichier a partir du debut
	deplacement_initial = 0
    CALL MPI_FILE_SET_VIEW(descripteur, deplacement_initial, MPI_REAL8, &
	                       type_sous_tab_vue, "native", MPI_INFO_NULL, code)

    !Creation du type derive type_sous_tab correspondant a la matrice u
    !sans les cellules fantomes
	
	profil_tab = size(u)
	profil_sous_tab = size(u)
	coord_debut = (/ 1 , 1/) ! Debut
	
	! call MPI_TYPE_CREATE_SUBARRAY (nb_dims,profil_tab,profil_sous_tab,coord_debut,
! ordre,ancien_type,nouveau_type,code)

	Call MPI_TYPE_CREATE_SUBARRAY(rang_tableau,profil_tab, profil_sous_tab, &
	     coord_debut,MPI_ORDER_FORTRAN, MPI_REAL8, type_sous_tab,code)
	CALL MPI_TYPE_COMMIT(type_sous_tab,code)


    !Ecriture du tableau u par tous les processus avec la vue
    CALL MPI_FILE_WRITE_ALL(descripteur, u, 1, type_sous_tab, statut, code)

    ! Fermeture du fichier
	CALL MPI_FILE_CLOSE(descripteur,code)

  END SUBROUTINE ecrire_mpi

  SUBROUTINE finalisation_mpi
    !************
    !Desactivation de l'environnement MPI
    !************
	
	Call MPI_TYPE_FREE(type_colonne,code)
	call MPI_TYPE_FREE(type_ligne,code)
	CALL MPI_COMM_FREE(comm2d,code)
    ! Desactivation de MPI
	CALL MPI_FINALIZE(code)

  END SUBROUTINE finalisation_mpi

END MODULE parallel

