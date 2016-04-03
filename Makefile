
PROG = ./dftd3 

#--------------------------------------------------------------------------
 OSTYPE=LINUXL
# source /usr/qc/lf95/csh_setup nicht vergessen (Lahey compiler)
# source /usr/qc/intel/compiler70/ia32/bin/ifcvars.csh (Intel compiler)
#--------------------------------------------------------------------------

OBJS= main.o dftd3.o copyc6.o param.o pars.o common.o

#--------------------------------------------------------------------------

ifeq ($(OSTYPE),LINUXL)
  #FC = lf95 
  FC = ifort
  #FC = gfortran
  #LINKER = lf95
  LINKER = ifort -static
  #LINKER = gfortran
  FFLAGS= -O -C -traceback -g
  #FFLAGS = -O -openmp -I$(MKLROOT)/include -mkl=parallel
  #LFLAGS = -openmp -I$(MKLROOT)/include -mkl=parallel
  #FFLAGS = -O  --chk a,e,s,u
  #FFLAGS = -O0 -g 
  #LFLAGS = -O0 -g
endif

ifeq ($(OSTYPE),LINUXI)
  FC = ifort 
  LINKER = ifort   
  FFLAGS = -w90 -O
  #LFLAGS =
endif                     

# diese ziele gibts:
.PHONY: all
.PHONY: clean

# dieses ist das erste auftretende,
# wird also beim aufruf von make erzeugt (default)
all: $(PROG)


#--------------------------------------------------------------------------
.SUFFIX:
.SUFFIX: .f90 .o

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(PROG): $(OBJS) 
	$(LINKER) $(OBJS) $(LIBS) -o $(PROG) $(LFLAGS)


#aufraeumen
clean:
	rm -f *.o $(PROG) 


# Abhaengigkeiten
main.o: dftd3.o
dftd3.o: param.o copyc6.o common.o
copyc6.o: pars.o common.o
