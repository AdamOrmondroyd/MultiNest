NESTLIBDIR = ../
LIBS =  -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB)
HEADERDIR = ../../include
 
OBJFILES = eggbox.o

all: eggbox 

%.o: %.c
	$(CC) $(CFLAGS) -I$(HEADERDIR) -c $*.c
 
eggbox: $(OBJFILES)
	$(FC) $(FFLAGS) $(NOFORMAIN) -o ../eggboxC $(OBJFILES) $(LIBS)

clean:
	rm -f *.o *.mod ../eggboxC
