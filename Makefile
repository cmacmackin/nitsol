ALL: bratu cavity pormed

# Macro definitions.

FC 		= gfortran
FFLAGS 		= -g
FLINKER 	= gfortran

NITSOL = ./Nitsol/libnitsol.a

LAPACK = -llapack

BLAS = -lblas

# Default compilation rules.

.f.o:
	$(FC) -c $(FFLAGS) $*.f

# Rules to build libraries.

nitsol_lib:
	cd Nitsol; $(MAKE) "FFLAGS=$(FFLAGS)" "FC=$(FC)"

lapack_lib:
	cd Lapack; $(MAKE) "FFLAGS=$(FFLAGS)" "FC=$(FC)"

blas_lib:
	cd Blas; $(MAKE) "FFLAGS=$(FFLAGS)" "FC=$(FC)"

# Rules to build test programs.

bratu_obj:
	cd Bratu; $(MAKE) "FFLAGS=$(FFLAGS)" "FC=$(FC)"

bratu:  nitsol_lib lapack_lib blas_lib bratu_obj
	$(FLINKER) -o bratu Bratu/*.o $(NITSOL) $(LAPACK) $(BLAS)

cavity_obj:
	cd Cavity; $(MAKE) "FFLAGS=$(FFLAGS)" "FC=$(FC)"

cavity:  nitsol_lib lapack_lib blas_lib cavity_obj
	$(FLINKER) -o cavity Cavity/*.o $(NITSOL) $(LAPACK) $(BLAS)

pormed_obj:
	cd Pormed; $(MAKE) "FFLAGS=$(FFLAGS)" "FC=$(FC)"

pormed:  nitsol_lib lapack_lib blas_lib pormed_obj
	$(FLINKER) -o pormed Pormed/*.o $(NITSOL) $(LAPACK) $(BLAS)

clean:
	- /bin/rm -f bratu cavity pormed core
	cd Bratu; $(MAKE) clean
	cd Cavity; $(MAKE) clean
	cd Pormed; $(MAKE) clean

veryclean:
	make clean
	cd Nitsol; $(MAKE) clean
	cd Lapack; $(MAKE) clean
	cd Blas; $(MAKE) clean

remake:
	make clean
	make

help:
	cat make.help
