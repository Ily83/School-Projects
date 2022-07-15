! Sur 1 processus !

program creation

  use mpi

  implicit none

  integer, parameter                  :: nb_valeurs=242
  integer, dimension(2*nb_valeurs)    :: valeurs
  integer, dimension(MPI_STATUS_SIZE) :: statut
  integer                             :: i,rang,descripteur,code

  open(unit=45,file="carte.dat")
  do i=1,nb_valeurs
     read(45,*) valeurs(i),valeurs(nb_valeurs+i)
  end do
  close(45)

  call MPI_INIT(code)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)

  call MPI_FILE_OPEN(MPI_COMM_WORLD,"donnees.dat",MPI_MODE_RDWR + MPI_MODE_CREATE, &
                     MPI_INFO_NULL,descripteur,code)

  call MPI_FILE_WRITE_ALL(descripteur,valeurs,2*nb_valeurs,MPI_INTEGER, &
                          statut,code)

  call MPI_FILE_CLOSE(descripteur,code)
  call MPI_FINALIZE(code)

end program creation
