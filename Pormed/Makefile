
# Macro definitions.

PORMED = pormed.o fpm.o jacvpm.o getrow.o
ILUTF = ilutf.o ilutf_aux.o

# Default target.

ALL:  $(PORMED) $(ILUTF)

# Default compilation rules.

.f.o:
	$(FC) -c $(FFLAGS) $*.f

# Dependencies.

fpm.o:		fpm.f
getrow.o:	getrow.f scratch.h
ilutf.o:	ilutf.f
ilutf_aux.o:	ilutf_aux.f
jacvpm.o:	jacvpm.f ilut_pormed.h scratch.h ../Nitsol/nitinfo.h 
pormed.o:	pormed.f ilut_pormed.h ../Nitsol/nitprint.h \
                ../Nitsol/nitparam.h

clean:
	- /bin/rm -f $(PORMED) $(ILUTF)

remake:
	make clean
	make
