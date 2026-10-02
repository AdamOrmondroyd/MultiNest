NESTLIBDIR = ../
LIBS = -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB)
 
OBJFILES = params.o like.o nestwrap.o main.o

all: loggamma 

%.o: %.f90
	$(FC) $(FFLAGS) -I$(NESTLIBDIR) -c $*.f90
 
loggamma: $(OBJFILES)
	$(FC) -o ../../bin/loggamma $(OBJFILES) \
	$(FFLAGS) $(LIBS)

clean:
	rm -f *.o *.mod ../../bin/loggamma
