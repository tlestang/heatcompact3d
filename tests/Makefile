test_field:  module_thomas.o module_differentiate.o \
             module_field.o test_field.o
	gfortran -o test_field module_thomas.o module_differentiate.o \
                               module_field.o test_field.o

test_differentiate: module_thomas.o module_differentiate.o \
                    test_differentiate.o
	gfortran -o test_differentiate module_thomas.o \
                    module_differentiate.o test_differentiate.o

module_%.o : ../src/module_%.f90
	gfortran -c $< -o $@
test_%.o : test_%.f90
	gfortran -c $< -o $@

.PHONY: clean
clean:
	rm *.o *.mod
