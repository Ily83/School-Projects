#!/bin/bash
scorep mpif90 -c produit_matrices_1.f90;
scorep mpif90 -o produit_matrices_1 produit_matrices_1.o;
mpiexec --oversubscribe -np 8 -x SCOREP_ENABLE_TRACING=1 -x SCOREP_TIMER=clock_gettime -x SCOREP_ENABLE_PROFILING=0 produit_matrices_1;