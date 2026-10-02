NESTLIBDIR = ../
LIBS =  -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB)
HEADERDIR = ../../include
 
OBJFILES = eggbox.o

all: eggbox 

%.o: %.c
	$(CC) $(CFLAGS) -I$(HEADERDIR) -c $*.c
 
eggbox: $(OBJFILES)
	$(FC) $(FFLAGS) -nofor_main -o ../eggboxC $(OBJFILES) $(LIBS)

clean:
	rm -f *.o *.mod ../eggboxC
