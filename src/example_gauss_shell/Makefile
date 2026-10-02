NESTLIBDIR = ../
LIBS = -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB)
 
OBJFILES = params.o like.o nestwrap.o main.o

all: gauss_shell 

%.o: %.f90
	$(FC) $(FFLAGS) -I$(NESTLIBDIR) -c $*.f90
 
gauss_shell: $(OBJFILES)
	$(FC) -o ../../bin/gauss_shell $(OBJFILES) \
	$(FFLAGS) $(LIBS)

clean:
	rm -f *.o *.mod ../../bin/gauss_shell
