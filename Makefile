# Makefile for Fortran F90 program

# Compiler
FC = gfortran

# Flags
FFLAGS = -O3 -g

# Source files
SRC = config.f90 panel3d.f90 libgeom.f90 main.f90 

# Library files

LIB = liblapack.a librefblas.a 

# Object files
OBJ = $(SRC:.f90=.o)

# Executable name
EXEC = solver

# Rule to build the executable
$(EXEC): $(OBJ) $(LIB)
	$(FC) $(FFLAGS) -o $@ $^

# Rule to build object files
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# Clean rule
clean:
	del /Q $(OBJ) $(EXEC)  
