module transposempi
  use mpi

  contains  
  subroutine transpose2(code,rang,nbr_procs,A,AT)
    real(kind=8),intent(in)             :: A(:,:)
    real(kind=8)                        :: AT(size(A,2),size(A,1))
    INTEGER                             :: code,rang, type_ligne, nbr_procs, type_transpose, taille_reel,i, j, M, N
    INTEGER(kind=MPI_ADDRESS_KIND)      :: pas
    M=size(A,1); N=size(A,2);
    print *, M,N
    AT(:,:) = 0.
    !Connaitre la taille du type de base MPI_REAL8
    call MPI_TYPE_SIZE(MPI_REAL8,taille_reel,code)
    ! Construction du type ligne
    call MPI_TYPE_VECTOR(N,1,M,MPI_REAL8,type_ligne,code)
    
    !Construction du type derive type_transpose pour transposer la
    !matrice A composee de M et de N   pas = 1*taille_real
     pas = 1.*taille_reel !  = 4
     call MPI_TYPE_CREATE_HVECTOR(M,1,pas,type_ligne,type_transpose,code)
  
    !Validation du type cree type_transpose
    call MPI_TYPE_COMMIT(type_transpose,code)

    IF (rang == 1) THEN
      !Envoi de la matrice A au processus 1 avec le type type_transpose
        call MPI_SEND(A,1,type_transpose,0,100,MPI_COMM_WORLD,code)
      ELSE
      !Reception pour le processus 1 dans la matrice AT
        call MPI_RECV(AT,N*M,MPI_REAL8,1,100,MPI_COMM_WORLD,MPI_STATUS_IGNORE ,code)
    END IF
    !Sortie de MPI
     call MPI_TYPE_FREE(type_transpose,code)
     call MPI_TYPE_FREE(type_ligne,code)
end subroutine
end module TransposeMPI
