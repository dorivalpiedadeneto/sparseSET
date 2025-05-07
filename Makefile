FC = gfortran
BINFLD = ./bin/
SRCFLD = ./src/
TSTFLD = $(SRCFLD)test/


all:
	$(MAKE) $(BINFLD)sparseset.mod

about:
	@echo The fortran compiler to be used is $(FC)
	@echo The folder where binaries are created is $(BINFLD) 
	@echo The folder where the source files are is $(SRCFLD)
	@echo The folder where the test files are is $(TSTFLD)

$(BINFLD)sparseset.mod: $(SRCFLD)sparseSET.f90
	@echo The sparse_set file is $< 
	@echo The outputfile is $@
	$(FC) -c $< -o $@
