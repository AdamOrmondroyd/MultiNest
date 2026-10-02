#FC = ifort -lmpi
#CC = icc -lmpi
#CXX = icpc -lmpi
#FFLAGS += -O3 -DMPI
#CFLAGS += -O3 -DMPI
## newer intel compilers spell this -nofor-main
#NOFORMAIN = -nofor_main

FC = mpif90
CC = gcc
CXX = g++
FFLAGS += -ffree-line-length-none -O3 -DMPI
CFLAGS += -O3 -DMPI

# gfortran >= 10 errors on argument rank mismatches in MPI calls
FFLAGS += $(shell [ $$($(FC) -dumpversion) -ge 10 ] && echo -fallow-argument-mismatch)

LAPACKLIB = -llapack

NESTLIBDIR = ./

export FC CC CXX FFLAGS CFLAGS LAPACKLIB NOFORMAIN

 
AR = ar r  
LINKLIB = ld -shared
 
NSOBJECTS = utils.o utils1.o priors.o kmeans_clstr.o xmeans_clstr.o posterior.o nested.o

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $^ 

%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $^ 

.PHONY: default all clean cleanall \
      obj_detect gaussian gaussian_weinberg loggamma rosenbrock ackley \
      himmelblau eggbox gauss_shell eggboxC eggboxC++ \
      clean_obj_detect clean_gaussian clean_gaussian_weinberg clean_loggamma \
      clean_rosenbrock clean_ackley clean_himmelblau clean_eggbox \
      clean_gauss_shell clean_eggbox_C clean_eggbox_C++

default: libnest3.a

all: libnest3.a obj_detect eggbox eggboxC eggboxC++ gaussian gaussian_weinberg \
loggamma gauss_shell rosenbrock himmelblau ackley
 
libnest3.so: $(NSOBJECTS) 
	$(LINKLIB) -o $(LIBS) $@ $^ 
 
libnest3.a: $(NSOBJECTS) 
	$(AR) $@ $^ 
 
obj_detect:
	make -C example_obj_detect
 
gaussian:
	make -C example_gaussian
 
gaussian_weinberg:
	make -C example_gaussian_weinberg
 
loggamma:
	make -C example_loggamma
 
rosenbrock:
	make -C example_rosenbrock
 
ackley:
	make -C example_ackley
 
himmelblau:
	make -C example_himmelblau
 
eggbox:
	make -C example_eggbox
 
gauss_shell:
	make -C example_gauss_shell
	
eggboxC:
	make -C example_eggbox_C
	
eggboxC++:
	make -C example_eggbox_C++

clean: 
	-rm $(NESTLIBDIR)/libnest3.*  *.o *.mod
	
cleanall: clean clean_obj_detect clean_gaussian clean_gaussian_weinberg \
clean_loggamma clean_gauss_shell clean_eggbox clean_eggbox_C clean_eggbox_C++ \
clean_rosenbrock clean_himmelblau \
clean_ackley

clean_obj_detect:
	make -C example_obj_detect clean
	
clean_gaussian:
	make -C example_gaussian clean
	
clean_gaussian_weinberg:
	make -C example_gaussian_weinberg clean
	
clean_loggamma:
	make -C example_loggamma clean
	
clean_rosenbrock:
	make -C example_rosenbrock clean
	
clean_ackley:
	make -C example_ackley clean
	
clean_himmelblau:
	make -C example_himmelblau clean
	
clean_eggbox:
	make -C example_eggbox clean
	
clean_gauss_shell:
	make -C example_gauss_shell clean
	
clean_eggbox_C:
	make -C example_eggbox_C clean
	
clean_eggbox_C++:
	make -C example_eggbox_C++ clean
