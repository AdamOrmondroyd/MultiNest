NESTLIBDIR = ../
LIBS =  -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB)
HEADERDIR = ../../include
 
OBJFILES = eggbox.o

all: eggboxC

%.o: %.c
	$(CC) $(CFLAGS) -I$(HEADERDIR) -c $*.c
 
eggboxC: $(OBJFILES)
	$(FC) $(FFLAGS) $(NOFORMAIN) -o ../../bin/eggboxC $(OBJFILES) $(LIBS)

clean:
	rm -f *.o *.mod ../../bin/eggboxC
