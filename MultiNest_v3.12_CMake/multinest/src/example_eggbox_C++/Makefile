NESTLIBDIR = ../
LIBS =  -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB) -lstdc++
HEADERDIR = ../../include
 
OBJFILES = eggbox.o

all: eggbox 

%.o: %.cc
	$(CXX) $(CFLAGS) -I$(HEADERDIR) -c $*.cc
 
eggbox: $(OBJFILES)
	$(FC) $(FFLAGS) -nofor_main -o ../eggboxC++ $(OBJFILES) $(LIBS)

clean:
	rm -f *.o *.mod ../eggboxC++
