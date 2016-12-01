ALL: libblas.a

BLAS = daxpy.o dcopy.o ddot.o dnrm2.o dscal.o dswap.o idamax.o lsame.o

# Dependencies.

daxpy.o:	daxpy.f
dcopy.o:	dcopy.f
ddot.o:		ddot.f 
dnrm2.o:	dnrm2.f
dscal.o:	dscal.f
dswap.o:	dswap.f
idamax.o:	idamax.f
lsame.o:	lsame.f

# Default compilation rules.

.f.o:
	$(FC) -c $(FFLAGS) $*.f

# Rules to build library.

libblas.a:  $(BLAS)
	ar rs $@ $?

clean:
	- /bin/rm $(BLAS) libblas.a

remake:
	make clean
	make
