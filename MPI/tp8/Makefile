############################ -*- Mode: Makefile -*- ###########################
## Makefile --- Cours MPI : TP8 : Equation de poisson
##
## Auteur          : Isabelle Dupays (CNRS/IDRIS) <Isabelle.Dupays@idris.fr>
## 
###############################################################################
SHELL := /bin/bash
# Compilateur et options de compilation et de link
include ../arch/make_inc

OBJS = module_params.o module_parallel_mpi.o module_calcul.o poisson.o
OBJS1 = read.o
OBJS2 = calcul_exact.o
OBJS3 = module_params.o module_parallel_mpi.o module_calcul.o poisson_non_bloquant.o

# R�gle implicite de compilation
.SUFFIXES: .o .f90
.f90.o:
	$(CF95) -c  $(FFLAGS_TP8) $<

default: poisson

all: poisson read calcul_exact exec

poisson: $(OBJS)
	$(CF95) $(LDFLAGS_TP8) -o $@ $(OBJS)

poisson_non_bloquant: $(OBJS3)
	$(CF95) $(LDFLAGS_TP8) -o $@ $(OBJS3)

read: $(OBJS1)
	$(CF95) $(LDFLAGS_TP8) -o $@ $(OBJS1)

calcul_exact: $(OBJS2)
	$(CF95) $(LDFLAGS_TP8) -o $@ $(OBJS2)

verification: fort.11 calcul_exact
	$(EXEC_TP8_VERIF) ./calcul_exact

fort.11: read exec
	$(EXEC_TP8_VERIF) ./read

exec: poisson poisson.data
	rm -f donnees.dat
	$(MPIEXEC_TP8) ./poisson

clean:
	rm -f $(OBJS) $(OBJS1) $(OBJS2) $(OBJS3) poisson poisson_non_bloquant calcul_exact read *.mod core donnees.dat fort.10 fort.11 *.mod
