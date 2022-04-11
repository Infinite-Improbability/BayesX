# Parallel with Intel Compilers

#FC = mpif90
#CC = mpicc
#CXX = mpiCC
#FFLAGS = -w -O3 -DMPI -extend_source
#CFLAGS = -I. -O3
#CXXFLAGS = -I. -O3
#LIBS = -llapack -lcfitsio


# Serial with Intel Compilers

#FC = ifort
#CC = icc
#CXX = icpc
#FFLAGS = -w -O3 -extend_source
#CFLAGS = -I. -O3
#CXXFLAGS = -I. -O3
#LIBS = -llapack -lcfitsio


# Parallel with GNU Compilers

FC = mpif90
CC = mpicc
CXX = mpic++
FFLAGS = -w -O3 -ffree-line-length-none -DMPI
CFLAGS = -I. -O3
CXXFLAGS = -I. -O3
LIBS = -llapack -lcfitsio


# Serial with GNU Compilers

#FC = gfortran
#CC = gcc
#CXX = g++
#FFLAGS = -w -O3 -ffree-line-length-none
#CFLAGS = -I. -O3
#CXXFLAGS = -I. -O3
#LIBS = -llapack -lcfitsio



export FC FFLAGS CC CFLAGS CXX CXXFLAGS LIBS
 
.PHONY: multinest BayesX
 
BINDIR = bin
 
all: multinest BayesX

multinest:
	gmake -C multinest all
      
BayesX:
	gmake -C src BayesX

clean:
	gmake -C src clean 
	gmake -C multinest clean
	-rm BayesX

