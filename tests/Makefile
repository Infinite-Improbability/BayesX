PROGS = BayesX-tests

INCLUDEDIRS = -I../src
 
override FFLAGS += $(INCLUDEDIRS)
 
AR = ar r  
LINKLIB = ld -shared  
LIBDIR = ../lib
BINDIR = ../bin
 
MCLIBS = -L$(LIBDIR) -lnest3 -lbayesx $(LIBS)
 
MOBJECTS = testdrive.o testdrive_version.o test_utilities.o main.o
 
%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $^
	
%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $^
 
 
all:  tests
 
tests: $(MOBJECTS)
	$(FC) $(FFLAGS) -o BayesX-tests $(MOBJECTS) $(MCLIBS)
	-mv BayesX-tests $(BINDIR)
	-rm *.o
 
clean: 
	-rm *.o *.mod
	-rm $(BINDIR)/BayesX-tests
