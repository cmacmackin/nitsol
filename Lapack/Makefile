ALL: liblapack.a

LAPACK = dlaic1.o dlamch.o

# Dependencies.

dlaic1.o:	dlaic1.f
dlamch.o:	dlamch.f

# Default compilation rules.

.f.o:
	$(FC) -c $(FFLAGS) $*.f

# Rules to build library.

liblapack.a:  $(LAPACK)
	ar rs $@ $?

clean:
	- /bin/rm $(LAPACK) liblapack.a

remake:
	make clean
	make
