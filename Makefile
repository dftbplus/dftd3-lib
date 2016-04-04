.SUFFIX:
.SUFFIX: .f90 .o

PROG = dftd3 
LIBD3 = libdftd3.a
TEST_LIBD3 = test_libdftd3


.PHONY: all lib

all: $(PROG) $(TEST_LIBD3)

include make.arch

lib:
	$(MAKE) -f Makefile.lib FC="$(FC)" FCFLAGS="$(FCFLAGS)" LN="$(LN)" \
	LNFLAGS="$(LNFLAGS)" $(LIBD3)


%.o: %.f90
	$(FC) $(FCFLAGS) -c $< -o $@

$(PROG): main.o $(LIBD3) 
	$(LN) $^ -o $@ $(LNFLAGS)

$(TEST_LIBD3): test_libdftd3.o $(LIBD3)
	$(LN) $^ -o $@ $(LNFLAGS)


.PHONY: clean distclean
clean:
	make -f Makefile.lib clean

distclean: clean
	make -f Makefile.lib distclean
	rm -f $(PROG)


test_libdftd3.o: lib
main.o: lib
