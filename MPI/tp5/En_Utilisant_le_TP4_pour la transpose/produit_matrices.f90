!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! produit_matrices.f90 --- TP5 : produit de matrices
!
! Auteur          : Jalel Chergui (CNRS/IDRIS - France) <Jalel.Chergui@idris.fr>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Remarques :
! ---------
!
!   * On veut r�aliser le produit de matrices C = A * B en parall�le.
!
!   * On suppose que ces matrices sont carr�es et que leur ordre N
!     est divisible par le nombre Nprocs de processus.
!
!   * Le processus 0 initialise les matrices A et B qu'il distribue
!     ensuite aux autres processus.
!
!   * La distribution de A se fait par bandes horizontales.
!     La distribution de B se fait par bandes verticales.
!
!   * Chaque processus poss�de une bande des matrices A et B.
!
!   * Chaque processus calcule ainsi un bloc de la diagonale principale
!     de C avec les �l�ments qu'il poss�de. Le calcul des blocs
!     extra-diagonaux n�cessite des communications avec les autres
!     processus.
!
!   * En fait, l'op�ration se reduit ici � un produit de matrices par bloc.

program produit_matrices
  use transposempi
  use mpi
  implicit none
	
  integer, parameter                          :: etiquette=1000, Nprocs_max=16
  integer                                     :: i,j,rang, Nprocs,N, NL, code, k, type_temp, type_tranche
  integer                                     :: rang_suivant, rang_precedent, taille_type_reel,res_i,res_o
  integer, dimension(MPI_STATUS_SIZE) 	      :: statut
  real(kind=8)                                :: Emax
  real(kind=8), allocatable, dimension(:,:)   :: A, B, C, CC, E,A_T
  real(kind=8), allocatable, dimension(:,:)   :: AL, BL, CL, TEMP
  integer(kind=MPI_ADDRESS_KIND)              :: borne_inferieure=0, taille_deplacement_type_tranche
	   
  ! Initialisation de MPI
  call MPI_INIT(code)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, Nprocs, code)

  if (rang == 0) then
    print *, "ALGORITHME 1: "
    print *, 'Entrez l''ordre N global des matrices :'
    read *,N
  end if
  
  call MPI_BCAST(N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, code)
  
  if (mod(N,Nprocs) == 0) then
	Nl = N/Nprocs
  else
	write(*,*) "Le nombre N choisi doit etre divisible par le nombre de processeurs"
    call MPI_ABORT(MPI_COMM_WORLD, 1, code)
  end if
    
  ! Le processus 0 initialise les matrices A et B
    allocate( A(N,N), B(N,N), C(N,N), CC(N,N),A_T(N,N) )
    A = 0.d0; B = 0.d0;  C = 0.d0; CC = 0.d0; A_T = 0.d0;
    ! Initialisation de A et B
    ! call RANDOM_NUMBER(A)
    ! call RANDOM_NUMBER(B)
    ! Allocation dynamique de m�moire, entre autres, des matrices A, B et C


	
	! A = floor(A*6)
	! B = floor(B*6)
	A = reshape((/ (i, i=1,N**2)  /), (/ N, N /))
	B = reshape((/ (i, i=1,N**2) /), (/ N, N /))
	if (rang == 0) then 
    write(*,*) "A = "
	write( * , "(*(g0))" ) ( (A(i,j)," ",j=1,N), new_line("A"), i=1,N )
	write(*,*) "B = "
	write( * , "(*(g0))" ) ( (B(i,j)," ",j=1,N), new_line("A"), i=1,N )
	end if
  
  call transpose2(code,rang,Nprocs,A,A_T)
  ! if (rang==0) write( * , "(*(g0))" ) "A_T r0:", ( (A_T(i,j)," ",j=1,N), new_line("A"), i=1,N )
  ! if (rang==1) write( * , "(*(g0))" ) "A_T r1 :", ( (A_T(i,j)," ",j=1,N), new_line("A"), i=1,N )



  ! Allocation dynamique de m�moire des divers tableaux locaux
  allocate( AL(N,NL), BL(N,NL), CL(N,NL), TEMP(N,NL) )
  AL = 0.d0; BL = 0.d0;  CL = 0.d0; TEMP = 0.d0;

  ! call MPI_TYPE_SIZE(MPI_REAL8, taille_type_reel, code)

  ! Construction du type qui correspond a 1 bloc de NL lignes et N colonnes
  call MPI_TYPE_CONTIGUOUS(N*NL, MPI_REAL8, type_tranche, code)  
  call MPI_TYPE_CONTIGUOUS(N*NL, MPI_REAL8, type_temp, code)
  call MPI_TYPE_COMMIT(type_tranche, code)
  call MPI_TYPE_COMMIT(type_temp, code)
  
  

  ! Le processus 0 distribue dans AL les tranches horizontales de la matrice A
    call MPI_SCATTER(A_T(1,1), 1, type_tranche, AL(1,1), 1, type_tranche, 0, MPI_COMM_WORLD, code)
  ! Le processus 0 distribue dans BL les tranches verticales de la matrice B
    call MPI_SCATTER(B(1,1), 1, type_tranche, BL(1,1), 1, type_tranche, 0, MPI_COMM_WORLD, code)  
  ! write(*,*)
  ! write (* ,*)"AL: ", rang , " : " , AL
  ! Write (*, *) "BL: ", rang , " : " , BL
  ! write(*,*)

  ! Calcul des blocs diagonaux de la matrice resultante.
  call matmul2(AL,BL, CL(rang*NL+1:rang*NL+NL,:),rang)


  ! Premier algorithme (deux fois plus co�teux que le second)
  write(*,*) "Nprocs: ", Nprocs
  do k = 0, Nprocs-1
    ! Chaque processus ENVOIE sa tranche AL au processus k
    ! et RE�OIT dans TEMP la tranche AL du processus k
	call MPI_SENDRECV(AL(1,1),1,type_tranche,k,etiquette,&
		              TEMP(1,1),1,type_temp,k,etiquette,MPI_COMM_WORLD,statut,code)
	
    if (rang /= k) then


		if (rang == 0) then
			write(*,*) "k = ", k, "NL = ", NL
			! write(*,*) "Temp = ", Temp
			! write(*,*)  "BL =", BL
			! write(*,*) 
		end if
		! write(*,*)
		! if (rang == 1) then
			! write(*,*) "k = ", k, "NL = ", NL
			! write(*,*) "Temp = ", Temp
			! write(*,*)  "BL =", BL
			! ! write(*,*) 
		! end if
		
      ! Chaque processus calcule les blocs situ�s au-dessus
      ! et en dessous du bloc de la diagonale principale
      ! CL(k*NL+1:(k+1)*NL,:)=matmul(TEMP(:,:),BL(:,:))
	  call matmul2(TEMP(:,:),BL(:,:), CL(k*NL+1:(k+1)*NL,:),rang)
    end if
  end do


  ! Le processus 0 collecte les tranches CL de tous les processus
  ! pour former la matrice r�sultante C
  
  call MPI_GATHER(CL, 1, type_tranche, C, 1, type_tranche, 0, MPI_COMM_WORLD, code)




  ! Les tableaux locaux sont d�sormais inutiles
  deallocate( AL, BL, CL, TEMP)
  if (ALLOCATED(A_T)) deallocate(A_t)


  ! V�rification des r�sultats
  if (rang == 0) then
    CC = matmul(A, B)
  	write(*,*)
	write(*,*) "Multi= "
	write( * , "(*(g0))" ) ( (C(i,j)," ",j=1,n), new_line("A"), i=1,n )
	write(*,*) 
	write(*,*) "Mono = "
	write( * , "(*(g0))" ) ( (CC(i,j)," ",j=1,n), new_line("A"), i=1,n )
    allocate( E(N,N) )
	! Calcul monoprocesseur du produit matriciel A*B
    E(:,:) = abs(C(:,:) - CC(:,:))
	write(*,*) "ERREUR: "
	write( * , "(*(g0))" ) ( (E(i,j)," ",j=1,n), new_line("A"), i=1,n )
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

  end if

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

end program produit_matrices