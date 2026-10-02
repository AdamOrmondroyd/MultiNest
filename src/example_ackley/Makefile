NESTLIBDIR = ../
LIBS = -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB)
 
OBJFILES = params.o like.o nestwrap.o main.o

all: ackley 

%.o: %.f90
	$(FC) $(FFLAGS) -I$(NESTLIBDIR) -c $*.f90
 
ackley: $(OBJFILES)
	$(FC) -o ../../bin/ackley $(OBJFILES) \
	$(FFLAGS) $(LIBS)

clean:
	rm -f *.o *.mod ../../bin/ackley
