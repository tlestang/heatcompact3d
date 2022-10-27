test_thomas:  module_thomas.o test_thomas.o
	gfortran -g -o test_thomas module_thomas.o test_thomas.o

test_field:  module_thomas.o module_differentiate.o \
             module_field.o test_field.o
	gfortran -g -o test_field module_thomas.o module_differentiate.o \
                               module_field.o test_field.o

test_differentiate: module_thomas.o module_differentiate.o \
                    test_differentiate.o
	gfortran -g -o test_differentiate module_thomas.o \
                    module_differentiate.o test_differentiate.o

test_time_integration: module_thomas.o module_differentiate.o \
                       module_field.o module_time_integration.o \
                       test_time_integration.o
	gfortran -g -o test_time_integration module_thomas.o \
                    module_differentiate.o module_field.o \
                    module_time_integration.o test_time_integration.o

module_%.o : ../src/module_%.f90
	gfortran -g -c $< -o $@
test_%.o : test_%.f90
	gfortran -g -c $< -o $@

.PHONY: clean test
clean:
	rm *.o *.mod
test: test_thomas test_differentiate test_field test_time_integration
	echo "-*-*-*- TESTING THOMAS SOLVER MODULE -*-*-*-"
	./test_thomas
	echo "-*-*-*- TESTING SPATIAL DIFFERENTIATION MODULE -*-*-*-"
	./test_differentiate
	echo "-*-*-*- TESTING FIELD MODULE -*-*-*-"
	./test_field
	echo "-*-*-*- TESTING TIME INTEGRATION MODULE -*-*-*-"
	./test_time_integration
