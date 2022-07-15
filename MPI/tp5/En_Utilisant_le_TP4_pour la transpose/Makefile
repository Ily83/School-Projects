############################ -*- Mode: Makefile -*- ###########################
## Makefile --- Cours MPI : TP5 : Communications collectives
##
## Auteur          : Jalel Chergui (CNRS/IDRIS)
## 
###############################################################################

# Compilateurs, options de compilation et d'�dition de liens
include ../../arch/make_inc

OBJS = produit_matrices.o

# R�gle implicite de compilation
.SUFFIXES: .o .f90
.f90.o:
	$(CF95) -c $(FFLAGS_TP5) $<

all : transposempi.mod transposempi.o produit_matrices

produit_matrices: $(OBJS)
	$(CF95) -o $@ transposempi.o $(LDFLAGS_TP5) $(OBJS)
	mpiexec -n 2 ./produit_matrices

transposempi.o: transposempi.f90
	mpif90 -c transposempi.f90 -o transposempi.o

transposempi.mod: transposempi.o transposempi.f90
	mpif90 -c transposempi.f90

clean:
	rm -f $(OBJS) produit_matrices *.o *.mod core
