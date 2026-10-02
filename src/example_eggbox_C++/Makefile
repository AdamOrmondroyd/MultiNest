NESTLIBDIR = ../
LIBS =  -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB) -lstdc++
HEADERDIR = ../../include
 
OBJFILES = eggbox.o

all: eggboxC++

%.o: %.cc
	$(CXX) $(CFLAGS) -I$(HEADERDIR) -c $*.cc
 
eggboxC++: $(OBJFILES)
	$(FC) $(FFLAGS) $(NOFORMAIN) -o ../../bin/eggboxC++ $(OBJFILES) $(LIBS)

clean:
	rm -f *.o *.mod ../../bin/eggboxC++
