################################################################################
#
# Library makefile
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

LIBD3 = libdftd3.a

OBJS  = api.o common.o core.o pars.o sizes.o

all: $(LIBD3)

$(LIBD3): $(OBJS)
	ar r $@ $^

clean:
	rm -f *.o

distclean: clean
	rm -f *.mod $(LIBD3)


%.o: %.f90
	$(FC) $(FCFLAGS) -c $< -o $@


# Dependencies
api.o: common.o core.o sizes.o
common.o:
core.o: common.o pars.o sizes.o
pars.o: common.o sizes.o
sizes.o:
