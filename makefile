COMPILER= gfortran

FLAGS = 

EXEC = nuclear-reactor

SRC = $(wildcard *.f90) 

OBJ = $(SRC:.f90=.o)

$(EXEC): $(OBJ)
	$(COMPILER) $(FLAGS) -o $@ $^

types.o: types.f90
	$(COMPILER) $(FLAGS) -c $<

quadrature.o: quadrature.f90 types.o
	$(COMPILER) $(FLAGS) -c $<

neutron_flux.o: neutron_flux.f90 types.o quadrature.o
	$(COMPILER) $(FLAGS) -c $<

read_write.o: read_write.f90 types.o neutron_flux.o
	$(COMPILER) $(FLAGS) -c $<

main.o: main.f90 types.o read_write.o
	$(COMPILER) $(FLAGS) -c $<

clean:
	rm -rf *.o *.mod

mrproper: clean
	rm -rf $(EXEC)
