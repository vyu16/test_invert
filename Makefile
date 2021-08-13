#FC = ifort
#LDFLAGS = -L/opt/intel/oneapi/mkl/2021.3.0/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

#FC = nvfortran
#FCFLAGS = -cuda -gpu=cc70 -Minfo
#LDFLAGS = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -cudalib=cusolver,cublas

SRC = test_lapack.f90 test_cusolver.f90 test_cublas.f90
EXE = $(SRC:.f90=.x)

all: $(EXE)

%.x: %.f90
	$(FC) $(FCFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f *.[ox]
