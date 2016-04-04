################################################################################
#
# Program makefile
#
# Needs following variables
#
# FC  Fortran compiler
# FCFLAGS  Fortran compiler flags
# LN  Linker
# LNFLAGS  Linker flags
#
################################################################################
.SUFFIX:
.SUFFIX: .f90 .o

.PHONY: clean distclean

DFTD3 = dftd3
LIBDIR = ../lib

all: $(DFTD3)

OBJS_DFTD3  = extras.o dftd3.o

$(DFTD3): $(OBJS_DFTD3)
	$(LN) $(LNFLAGS) -o $@ $^ -L$(LIBDIR) -ldftd3

clean:
	rm -f *.o *.mod

distclean: clean
	rm -f $(DFTD3)


%.o: %.f90
	$(FC) $(FCFLAGS) -I$(LIBDIR) -c $< -o $@


# Dependencies
dftd3.o: extras.o
