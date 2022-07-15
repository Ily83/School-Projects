############################ -*- Mode: Makefile -*- ###########################
## Makefile --- Cours MPI : TP5 : Communications collectives
##
## Auteur          : Jalel Chergui (CNRS/IDRIS)
## 
###############################################################################

# Compilateurs, options de compilation et d'édition de liens
include ../arch/make_inc

OBJS1 = produit_matrices_1.o
OBJS2 = produit_matrices_2.o

# Règle implicite de compilation
.SUFFIXES: .o .f90
.f90.o:
	$(CF95) -c $(FFLAGS_TP5) $< 

all : produit_matrices_1 produit_matrices_2

produit_matrices_1: $(OBJS1)
	$(CF95) -o $@ $(LDFLAGS_TP5) $(OBJS1)
	$(MPIEXEC_TP5) ./produit_matrices_1
	
produit_matrices_2: $(OBJS2)
	$(CF95) -o $@ $(LDFLAGS_TP5) $(OBJS2) 
	$(MPIEXEC_TP5) ./produit_matrices_2

clean:
	rm -f $(OBJS1) $(OBJS2) produit_matrices_1 produit_matrices_2 core
