.SUFFIX:
.SUFFIX: .f90 .o

.PHONY: dftd3 lib testlib

all: lib dftd3

include make.arch

lib:
	$(MAKE) -C lib FC="$(FC)" FCFLAGS="$(FCFLAGS)" LN="$(LN)" \
	LNFLAGS="$(LNFLAGS)"

dftd3: lib
	$(MAKE) -C prg FC="$(FC)" FCFLAGS="$(FCFLAGS)" LN="$(LN)" \
	LNFLAGS="$(LNFLAGS)" dftd3


.PHONY: clean distclean
clean:
	make -C lib clean

distclean: 
	make -C lib distclean
