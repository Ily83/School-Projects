############################ -*- Mode: Makefile -*- ###########################
## Makefile --- Cours MPI : TP1 : Environnement
##
## Auteur          : Denis Girou (CNRS/IDRIS) <Denis.Girou@idris.fr>
###############################################################################

# Compilateurs, options de compilation et d'�dition de liens
include ../arch/make_inc

OBJS = pair_impair.o

# R�gle implicite de compilation
.SUFFIXES: .o .f95
.f95.o:
	$(CF95) -g -c $(FFLAGS_TP1) $<

all: pair_impair

pair_impair:$(OBJS)
	$(CF95) -o $@ $(LDFLAGS_TP1) $(OBJS)
	$(MPIEXEC_TP1) ./pair_impair

clean:
	rm -f $(OBJS) pair_impair core
