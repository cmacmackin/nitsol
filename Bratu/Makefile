
# Macro definitions.

BRATU = bratu.o fbratu.o jacvbratu.o psolb.o
FPSLV = hwscrtd.o

# Default target.

ALL:  $(BRATU) $(FPSLV)

# Default compilation rules.

.f.o:
	$(FC) -c $(FFLAGS) $*.f

# Dependencies.

bratu.o:	bratu.f ../Nitsol/nitprint.h ../Nitsol/nitparam.h
fbratu.o:	fbratu.f
javbratu.o:	jacvbratu.f
psolb.o:	psolb.f
hwscrtd.o:	hwscrtd.f

clean:
	- /bin/rm -f $(BRATU) $(FPSLV)

remake:
	make clean
	make
