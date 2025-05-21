FC = gfortran
BINFLD = ./bin
SRCFLD = ./src
TSTFLD = $(SRCFLD)/test
OPTFLAG = -O2

all:
	$(MAKE) $(BINFLD)/tests

about:
	@echo The fortran compiler to be used is $(FC)
	@echo The folder where binaries are created is $(BINFLD) 
	@echo The folder where the source files are is $(SRCFLD)
	@echo The folder where the test files are is $(TSTFLD)

$(BINFLD)/sparseset.o: $(SRCFLD)/sparseSET.f90
	@echo The sparse_set file is $< 
	@echo The outputfile is $@
	$(FC) -J$(BINFLD) -c $< -o $@

$(BINFLD)/dev_tests.o: $(TSTFLD)/dev_tests.f90 $(BINFLD)/sparseset.o
	$(FC) -J$(BINFLD) -c $< -o $@

$(BINFLD)/tests: $(TSTFLD)/tests.f90 $(BINFLD)/dev_tests.o $(BINFLD)/sparseset.o
	$(FC) -I$(BINFLD)  $^ -o $@

$(BINFLD)/ptests: $(TSTFLD)/ptests.f90
	$(FC) $(OPTFLAG) $< -o $@

run-ptest: $(BINFLD)/ptests
	@$<

run-test: $(BINFLD)/tests
	@$<

clean:
	rm -rf $(BINFLD)/*.mod
	rm -rf $(BINFLD)/*.o
	rm -rf $(BINFLD)/tests
	rm -rf $(BINFLD)/ptests
