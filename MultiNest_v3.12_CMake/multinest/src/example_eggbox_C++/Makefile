NESTLIBDIR = ../
LIBS =  -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB) -lstdc++
HEADERDIR = ../../include
 
OBJFILES = eggbox.o

all: eggboxC++

%.o: %.cc
	$(CXX) $(CFLAGS) -I$(HEADERDIR) -c $*.cc
 
eggbox: $(OBJFILES)
	$(FC) $(FFLAGS) $(NOFORMAIN) -o ../eggboxC++ $(OBJFILES) $(LIBS)

clean:
	rm -f *.o *.mod ../eggboxC++
