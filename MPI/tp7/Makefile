############################ -*- Mode: Makefile -*- ###########################
## Makefile --- T.P. 7 
## 
## Auteur          : Denis GIROU (CNRS/IDRIS - France) <Denis.Girou@idris.fr>
## 
###############################################################################

# Compilateur et options de compilation et de link
include ../arch/make_inc

OBJS    = lire_fichier.o

# Règle implicite de compilation
.SUFFIXES:
.SUFFIXES: .o .f90
.f90.o: 
	$(CF95) -c $(FFLAGS_TP7) $<
	
all: donnee lire_fichier verification

lire_fichier: donnee $(OBJS)
	$(CF95) $(LDFLAGS_TP7) -o $@ $(OBJS)
	$(MPIEXEC_TP7) ./lire_fichier
	rm -f $(OBJS) core

donnee: 
	cd creation;make clean;make

verification: 
	$(CC) verification.c $(LDCFLAGS) -o verification
	$(EXEC_TP7_VERIF) ./verification
	rm -f verification


.PHONY : clean cleanall exec
clean:
	rm -f *.o lire_fichier fichier*.dat verification core*

cleanall:
	make clean
	rm -f donnees.dat
