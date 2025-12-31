# Compiler and flags
FC = gfortran
FFLAGS = -fcheck=all -g
LDFLAGS = -llapack -lblas

# link object file and output the program program_main
runccfull: main.o initialization.o Interaction.o CoulombWave.o Numerov.o DiscreteBasis.o ModifiedDiscreteBasis.o Rmatrix.o
	$(FC) $(FFLAGS) -o runccfull main.o initialization.o Interaction.o CoulombWave.o Numerov.o DiscreteBasis.o ModifiedDiscreteBasis.o Rmatrix.o $(LDFLAGS)

main.o: main.f90 ccfull_initialization_mod.mod
	$(FC) $(FFLAGS) -c main.f90

initialization.o ccfull_initialization_mod.mod: initialization.f90
	$(FC) $(FFLAGS) -c initialization.f90

Interaction.o: Interaction.f90
	$(FC) $(FFLAGS) -c Interaction.f90

CoulombWave.o: CoulombWave.f90
	$(FC) $(FFLAGS) -c CoulombWave.f90

Numerov.o: Numerov.f90
	$(FC) $(FFLAGS) -c Numerov.f90

DiscreteBasis.o: DiscreteBasis.f90
	$(FC) $(FFLAGS) -c DiscreteBasis.f90

ModifiedDiscreteBasis.o: ModifiedDiscreteBasis.f90
	$(FC) $(FFLAGS) -c ModifiedDiscreteBasis.f90

Rmatrix.o: Rmatrix.f90
	$(FC) $(FFLAGS) -c Rmatrix.f90
# delete object and mod file
clean:
	rm -f *.o *.mod 

.PHONY: clean