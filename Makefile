FC = gfortran
BINFLD = ./bin/
SRCFLD = ./src/
TSTFLD = $(SRCFLD)test/


about:
	@echo The fortran compiler to be used is $(FC)
	@echo The folder where binaries are created is $(BINFLD) 
	@echo The folder where the source files are is $(SRCFLD)
	@echo The folder where the test files are is $(TSTFLD)
