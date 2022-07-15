!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! produit_matrices.f90 --- TP5 : produit de matrices
!
! Auteur          : Jalel Chergui (CNRS/IDRIS - France) <Jalel.Chergui@idris.fr>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Remarques :
! ---------
!
!   * On veut réaliser le produit de matrices C = A * B en parallèle.
!
!   * On suppose que ces matrices sont carrées et que leur ordre N
!     est divisible par le nombre Nprocs de processus.
!
!   * Le processus 0 initialise les matrices A et B qu'il distribue
!     ensuite aux autres processus.
!
!   * La distribution de A se fait par bandes horizontales.
!     La distribution de B se fait par bandes verticales.
!
!   * Chaque processus possède une bande des matrices A et B.
!
!   * Chaque processus calcule ainsi un bloc de la diagonale principale
!     de C avec les éléments qu'il possède. Le calcul des blocs
!     extra-diagonaux nécessite des communications avec les autres
!     processus.
!
!   * En fait, l'opération se reduit ici à un produit de matrices par bloc.

program produit_matrices

  use mpi
  implicit none
	
  integer, parameter                          :: etiquette=1000, Nprocs_max=16
  integer                                     :: i,j,rang, Nprocs,N, NL, code, k, type_temp, type_tranche
  integer                                     :: rang_suivant, rang_precedent, taille_type_reel,res_i,res_o
  integer, dimension(MPI_STATUS_SIZE) 	      :: statut
  real(kind=8)                                :: Emax
  real(kind=8), allocatable, dimension(:,:)   :: A, B, C, CC, E,A_T
  real(kind=8), allocatable, dimension(:,:)   :: AL, BL, CL
  integer(kind=MPI_ADDRESS_KIND)              :: borne_inferieure=0, taille_deplacement_type_tranche
  real(kind=8)                        		  :: temps_debut,temps_fin

   write(*,*)
   
  ! Initialisation de MPI
  call MPI_INIT(code)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, Nprocs, code)


  if (rang == 0) then
	print *, "ALGORITHME 2: "
    print *, 'Entrez l''ordre N global des matrices :'
    read *,N
  end if
  temps_debut=MPI_WTIME()
  ! Envoie à tous les processus le nombre "N" entré
  call MPI_BCAST(N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, code)
  
  !Determination du nombre de lignes à prendre 
  if (mod(N,Nprocs) == 0) then
	Nl = N/Nprocs
	if (rang==0) print *,"NL = ",Nl
  else
	write(*,*) "Le nombre N choisi doit etre divisible par le nombre de processeurs"
    call MPI_ABORT(MPI_COMM_WORLD, 1, code)
  end if

   ! Allocation dynamique de mémoire, entre autres, des matrices A, B et C
    allocate( A(N,N), B(N,N), C(N,N), CC(N,N),A_T(N,N) )
    A = 0.d0; B = 0.d0;  C = 0.d0; CC = 0.d0;
    ! Initialisation de A et B
    ! call RANDOM_NUMBER(A)
    ! call RANDOM_NUMBER(B)
	! A = floor(A*6)
	! B = floor(B*6)

	A = reshape((/ (i, i=1,N**2)  /), (/ N, N /))
	B = reshape((/ (i, i=1,N**2) /), (/ N, N /))
	if (rang == 0) then 
	! Call afficher(A)
	! call afficher(B)
		A_T = transpose(A)
	end if

  ! Allocation dynamique de mémoire des divers tableaux locaux
  allocate( AL(N,NL), BL(N,NL), CL(N,NL))
  AL = 0.d0; BL = 0.d0;  CL = 0.d0;

  ! Construction du type qui correspond à une colonne de N*NL valeurs
  call MPI_TYPE_CONTIGUOUS(N*NL, MPI_REAL8, type_tranche, code)  
  call MPI_TYPE_COMMIT(type_tranche, code)
  
  

  ! Le processus 0 distribue dans AL les tranches horizontales de la matrice A
    call MPI_SCATTER(A_T(1,1), 1, type_tranche, AL(1,1), 1, type_tranche, 0, MPI_COMM_WORLD, code)
  ! Le processus 0 distribue dans BL les tranches verticales de la matrice B
    call MPI_SCATTER(B(1,1), 1, type_tranche, BL(1,1), 1, type_tranche, 0, MPI_COMM_WORLD, code)  
	
! -------------------   Second algorithme ------------------------------------------
 rang_suivant = modulo(rang+1, Nprocs)
 rang_precedent = modulo(rang-1, Nprocs)
 res_i = rang*NL
 res_o = res_i
 do 
 
   ! Chaque processus calcule les blocs situés au-dessus
   ! et en dessous du bloc de la diagonale principale
   CALL  MATMUL2(AL, BL, CL(res_o+1:res_o+NL, :),rang)
   
   ! Passage au rang suivant 
   res_o = modulo(res_o+NL, N)
   
     ! Processus fini pour cette tranche
    if(res_o == res_i) exit
	
    ! Chaque processus ENVOIE sa tranche AL au processus précédent
    ! et REÇOIT la tranche AL du processus suivant (mais les contenus changent)
	call MPI_SENDRECV_REPLACE(AL(1,1),1,type_tranche,rang_precedent,etiquette,&
							  rang_suivant,etiquette,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
	end do
! -------------------   Fin Second algorithme ------------------------------------------

  ! Le processus 0 collecte les tranches CL de tous les processus
  ! pour former la matrice résultante C
  ! if (rang == 0) call afficher(CL)
  ! if (rang == 1) call afficher(Cl)
  call MPI_GATHER(CL, 1, type_tranche, C, 1, type_tranche, 0, MPI_COMM_WORLD, code)




  ! Les tableaux locaux sont désormais inutiles
  deallocate( AL, BL, CL)
  if (ALLOCATED(A_T)) deallocate(A_t)

  temps_fin=MPI_WTIME()
  
  ! Vérification des résultats
  if (rang == 0) then
    CC = matmul(A, B) ! Résultat monoprocesseur
  	! write(*,*)
	! write(*,*) "Multi= "
	! call afficher(C)
	! write(*,*) 
	! write(*,*) "Mono = "
	! Call afficher(CC)
    allocate( E(N,N) )
	! Calcul monoprocesseur du produit matriciel A*B
    E(:,:) = abs(C(:,:) - CC(:,:))
	! write(*,*) "ERREUR: "
	! call afficher(E)
    Emax   = maxval( E(:,:) ) / N**2
    deallocate( A, B, C, CC, E )

    if ( Emax <= epsilon(1.0) ) then
      print'(/,40X,"Bravo !",/,  &
            & 20X,"Le produit matriciel A*B calcule en parallele",/, &
            & 20X,"est bien egal a celui calcule en monoprocesseur")'
    else
      print'(/,33X,"Resultat incorrect !",/, &
            & 20X,"Le produit matriciel A*B calcule en parallele",/, &
            & 20X,"est different de celui calcule en monoprocesseur")'
    end if
	
	print ('(" Produit de matrice realise en ",f8.6," secondes.")'),temps_fin-temps_debut
  end if
  
  ! Fin MPI 
  call MPI_FINALIZE(code)
  
  contains
  
      subroutine MATMUL2(AL, BL, CL,rang)
        real(kind=8), dimension(:,:), intent(in) :: AL, BL
        real(kind=8), dimension(:,:), intent(out) :: CL
        integer :: i, j,rang
		forall (i = 1:size(CL,1),j=1:size(CL,2))
			    CL(i,j) = dot_product(AL(:,i), BL(:,j))
		end forall
    end subroutine MATMUL2
	
	subroutine Afficher(M)
	real(kind=8) :: M(:,:)
	integer       :: i,j
	write( * , "(*(g0))" ) ( (M(i,j)," ",j=1,size(M,2)), new_line("A"), i=1,size(M,1) )
	end subroutine
end program produit_matrices