program pair_impair
	use mpi
	implicit none
	integer 		   :: rang,code,nbprocs

	call MPI_INIT(code) !! Initialisation de l’environnement MPI
	call MPI_Comm_size(MPI_COMM_WORLD,nbprocs,code) ! Nombre de rang
	call MPI_COMM_rank(MPI_COMM_WORLD,rang,code) !! Nombre de processus
	
	! Pair ou Impair
	if (mod(rang,2) == 1) then ! Phase impair
        write(*,'(a35,I2)') "Je suis le processus impair de rang", rang
	else                       ! Phase pair
        write(*,'(a33,I2)') "Je suis le processus pair de rang", rang
	end if 

	! Desactivation de l’environnement MPI
	call MPI_FINALIZE(code)
end program pair_impair