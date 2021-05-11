F90 = gfortran
FOPT = -O2
LIB = -llapack -lblas

OBJ = \
calculate_one_integral.o \
Fgamma.o \
pair_index.o \
construct_JK_matrices.o \
calculate_two_integral.o \


.PHONY:	h2
h2: $(OBJ)
	$(F90) $(FOPT) -o h2 hydrogen_molecule.f90 $(OBJ) $(LIB)

.SUFFIXES: .o .f90

.f90.o:
	$(F90) $(FOPT) -c $<

.PHONY:	clean
clean:
	rm -f *.o h2*
