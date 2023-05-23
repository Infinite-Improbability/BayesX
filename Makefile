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

FC = $(PREP) mpifort
CC = $(PREP) mpicc
CXX = $(PREP) mpic++
FFLAGS = -w -O3 -ffree-line-length-none -DMPI -fallow-argument-mismatch -fopenmp
CFLAGS = -I. -O3
CXXFLAGS = -I. -O3
LIBS = -llapack -lcfitsio


# Serial with GNU Compilers

# FC = gfortran
# CC = gcc
# CXX = g++
# FFLAGS = -w -O3 -ffree-line-length-none
# CFLAGS = -I. -O3
# CXXFLAGS = -I. -O3
# LIBS = -llapack -lcfitsio


# Serial with GNU Compilers debug flags

# FC = gfortran
# CC = gcc
# CXX = g++
# FFLAGS = -ggdb3 -ffree-line-length-none
# CFLAGS = -ggdb3 -I. -Og
# CXXFLAGS = -I. -Og
# LIBS = -llapack -lcfitsio

export FC FFLAGS CC CFLAGS CXX CXXFLAGS LIBS

.PHONY: multinest BayesX test

BINDIR = bin

all: multinest BayesX

multinest: | multinestDirs
	$(MAKE) -C multinest all

multinestDirs:
	mkdir -p lib

BayesX: multinest | BayesXdirs
	$(MAKE) -C src BayesX

BayesXdirs:
	mkdir -p $(BINDIR) chains

test: multinest | BayesXdirs
	$(MAKE) -C src libbayesx.a
	$(MAKE) -C tests tests

clean:
	$(MAKE) -C src clean
	$(MAKE) -C multinest clean
	-rm BayesX
