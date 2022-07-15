############################ -*- Mode: Makefile -*- ###########################
## Makefile --- Cours MPI : TP6 : Subdiviser une grille 2D avec MPI_COMM_SPLIT
##
## Auteur          : Jalel Chergui (CNRS/IDRIS)
## 
###############################################################################

# Compilateurs, options de compilation et d'édition de liens
include ../arch/make_inc

OBJS = commsplit.o

# Règle implicite de compilation
.SUFFIXES: .o .f90
.f90.o:
	$(CF95) -c $(FFLAGS_TP6) $<

all : commsplit

commsplit: $(OBJS)
	$(CF95) -o $@ $(LDFLAGS_TP6) $(OBJS)
	$(MPIEXEC_TP6) ./commsplit

clean :
	rm -f $(OBJS) commsplit core
