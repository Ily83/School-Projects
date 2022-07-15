# Makefile
CC = gfortran
CFLAGS1 = -c
CFLAGS2 = 
# CFLAGS1 = -Wall -g -c
# CFLAGS2 = -g

MN: ogpf.o generateur_aleatoire.o main.o
	${CC} ${CFLAGS2} -o MN generateur_aleatoire.o ogpf.o main.o -static
	./MN

ogpf.mod: ogpf.o ogpf.f90
	${CC} ${CFLAGS1} ogpf.f90 -static

ogpf.o: ogpf.f90
	${CC} ${CFLAGS1} ogpf.f90 -static

main.o:	generateur_aleatoire.mod ogpf.mod main.f90
	${CC} ${CFLAGS1} main.f90 -static

generateur_aleatoire_aleatoire.mod: generateur_aleatoire.o generateur_aleatoire.f90
	${CC} ${CFLAGS1} ogpf.f90 -static

generateur_aleatoire.o: generateur_aleatoire.f90
	${CC} ${CFLAGS1} generateur_aleatoire.f90 -static

clean:
	rm *.o *.mod *.exe *.plt *.gp
