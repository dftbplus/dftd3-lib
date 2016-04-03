PROG = dftd3 
LIBD3 = libdftd3.a

OBJS = main.o


.PHONY: all lib

all: $(PROG)


.SUFFIX:
.SUFFIX: .f90 .o

include make.arch

lib:
	$(MAKE) -f Makefile.lib FC="$(FC)" FCFLAGS="$(FFLAGS)" LN="$(LINKER)" \
	LNFLAGS="$(LFLAGS)"


%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(PROG): $(OBJS) $(LIBD3) 
	$(LINKER) $^ -o $@ $(LFLAGS)


.PHONY: clean
clean:
	make -f Makefile.lib clean
	rm -f $(PROG)


main.o: lib

