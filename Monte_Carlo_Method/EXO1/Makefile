# Makefile
CC = gfortran
CFLAGS1 = -c
CFLAGS2 = 
# CFLAGS1 = -Wall -g -c
# CFLAGS2 = -g

MN: ogpf.o main.o
	${CC} ${CFLAGS2} -o MN airevolume.o ogpf.o main.o -static
	./MN

main.o:	airevolume.mod ogpf.mod main.f90
	${CC} ${CFLAGS1} main.f90 -static

# MN: ogpf.o mainbis.o
# 	${CC} ${CFLAGS2} -o MN airevolume.o ogpf.o mainbis.o -static

# mainbis.o:	airevolume.mod ogpf.mod main.f90
# 	${CC} ${CFLAGS1} mainbis.f90 -static

ogpf.mod: ogpf.o ogpf.f90
	${CC} ${CFLAGS1} ogpf.f90 -static

ogpf.o: ogpf.f90
	${CC} ${CFLAGS1} ogpf.f90 -static

airevolume.mod: airevolume.o airevolume.f90
	${CC} ${CFLAGS1} airevolume.f90 -static

airevolume.o: airevolume.f90
	${CC} ${CFLAGS1} airevolume.f90 -static

clean:
	rm *.o *.mod *.exe *.plt *.gp
