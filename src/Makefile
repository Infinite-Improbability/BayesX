PROGS = BayesX

INCLUDEDIRS = -I../multinest
 
FFLAGS += $(INCLUDEDIRS)
 
AR = ar r  
LINKLIB = ld -shared  
LIBDIR = ../lib
BINDIR = ../bin
 
MCLIBS = -L$(LIBDIR) -lnest3 $(LIBS)
 
MOBJECTS = params.o ReadInputs.o constants.o random.o rewind.o matrix_utils.o utilities.o cosmology.o massfunction.o MassModels.o ReadWrite.o GasModels.o Priors.o CheckPars.o like.o nestwrap.o BayesX.o
 
%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $^
	
%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $^
 
 
all:  BayesX
 
BayesX: $(MOBJECTS)
	$(FC) $(FFLAGS) -o BayesX $(MOBJECTS) $(MCLIBS)
	-mv $@ $(BINDIR)
	-rm *.o

libbayesx.a: $(MOBJECTS) 
	$(AR) $@ $^ 
	-mv $@ $(LIBDIR)
	-rm *.o
 
clean: 
	-rm *.o *.mod
	-rm $(BINDIR)/BayesX

 

