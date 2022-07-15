program pi
  use mpi
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer, parameter :: li = selected_int_kind(15)
  real(dp), parameter :: pi_ex = 4_dp*atan(1.0)
  double precision ::   ecart, pi_f
  integer(kind=li):: nbbloc,i
  double precision :: t1 , t2 , temps
  real(kind=dp) :: largeur,somme,x
  integer :: code , rang , nbr_procs



 ! debut de la region parallele
 call MPI_INIT(code)
 call MPI_Comm_size(MPI_COMM_WORLD,nbr_procs,code)
 call MPI_Comm_rank(MPI_COMM_WORLD,rang,code)
 t1 = MPI_WTIME()


  ! Nombre d'intervalles
  nbbloc = 3*1000*1000_li*100
  ! largeur des intervalles
  largeur = 1._dp / real(nbbloc,dp)

  somme = 0._dp

 do i = rang +1 , nbbloc , nbr_procs
    ! Point au milieu de l'intervalle
    x = largeur*(real(i,dp)-0.5_dp)
    ! Calcul de l'aire
    somme = somme + (4._dp / (1._dp + x*x))
 end do
 somme = largeur * somme
 
 ! Reduction avec somme (MPI_sum)
 call MPI_Reduce(somme,pi_f,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,code);
 t2 = MPI_WTIME()
 temps = t2 - t1
 
 ecart = abs(pi_ex-pi_f)

 if(rang ==0) then
 write(* ,*) " pi_calculé           : " , pi_f
 write(* ,*) " Nombre de processus  : " , nbr_procs
 write(* ,*) " Nombre d intervalles : " , nbbloc
 write(* ,*) " |  pi_ex - pi_f  |   : " , ecart
 write(* ,*) " Temps reel           : " , temps
 endif

 ! fin de la region parallele
 call MPI_Finalize ( code )

end program pi