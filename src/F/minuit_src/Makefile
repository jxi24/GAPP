
FC=gfortran
#FCOPT= -fbacktrace
FCOPT= -g -w  -C  -O
PATHLIB=../../../lib/

all: lib

lib: 
	$(FC) $(FCOPT) -c *.F
	ar rcs libminuit.a *.o

install: all
	cp libminuit.a  $(PATHLIB)

clean: 
	rm *.o 
	rm *.mod
	rm *.a
