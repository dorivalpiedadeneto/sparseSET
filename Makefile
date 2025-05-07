FC = gfortran
BINFLD = ./bin/
SRCFLD = ./src/
TSTFLD = $(SRCFLD)test/


all:
	$(MAKE) $(BINFLD)sparseset.mod
	$(MAKE) $(BINFLD)dev_tests.mod

about:
	@echo The fortran compiler to be used is $(FC)
	@echo The folder where binaries are created is $(BINFLD) 
	@echo The folder where the source files are is $(SRCFLD)
	@echo The folder where the test files are is $(TSTFLD)

$(BINFLD)sparseset.mod: $(SRCFLD)sparseSET.f90
	@echo The sparse_set file is $< 
	@echo The outputfile is $@
	$(FC) -c $< -o $@

$(BINFLD)dev_tests.mod: $(TSTFLD)dev_tests.f90 $(BINFLD)sparseset.mod
	$(FC) -J$(BINFLD) -c $< -o $@

clean-sparseset:
	rm -rf $(BINFLD)sparseset.mod

clean-dev_tests:
	rm -rf $(BINFLD)dev_tests.mod

clean:
	$(MAKE) clean-sparseset
	$(MAKE) clean-dev_tests
