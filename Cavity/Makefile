
# Macro definitions.

CAVITY = cavity.o fcav.o jacvcav.o init.o init2.o psol.o bihsl.o resid.o \
         wload.o 
LINPACK = dspfa.o dspsl.o dppfa.o dppsl.o
FBISLV = dbihar.o $(LINPACK)

# Default target.

ALL:  $(CAVITY) $(FBISLV)

# Default compilation rules.

.f.o:
	$(FC) -c $(FFLAGS) $*.f

# Dependencies.

cavity.o:	cavity.f ../Nitsol/nitprint.h ../Nitsol/nitparam.h
fcav.o:		fcav.f
jacvcav.o:	jacvcav.f
init.o:		init.f
init2.o:	init2.f
psol.o:		psol.f
resid.o:	resid.f
wload.o:	wload.f

dbihar.o:	dbihar.f

dppfa.o:	dppfa.f
dppsl.o:	dppsl.f
dspfa.o:	dspfa.f
dspsl.o:	dspsl.f

clean:
	- /bin/rm -f $(CAVITY) $(FBISLV)

remake:
	make clean
	make
