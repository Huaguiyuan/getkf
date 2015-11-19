.SUFFIXES: .f90
.f90.o:
	$(F90) $(FOPTS) -c $<


FOPTS=-fpp -I$(MKLROOT)/include -openmp
LIBS=-L$(MKLROOT)/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core

test.x: main.f90 constants.o banddata.o wanndata.o solve_kf.o
	$(F90) $(FOPTS) -o test.x main.f90 banddata.o constants.o wanndata.o solve_kf.o $(LIBS) 

clean:
	rm -f *.x *.o *.mod
