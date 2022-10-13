#FC = ifort -g -traceback -check bounds -check uninit -check pointers -fpe0
#LDFLAGS = -L/opt/intel/oneapi/mkl/2021.3.0/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

FC = nvfortran
FCFLAGS = -cuda -gpu=cc80,cuda11.7 -Minfo -cudalib=cublas,cusolver -acc

all: test_real.x test_cmplx.x

%.x: %.o test_cublas.o test_custom_gpu_newapi.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ -c $<

clean:
	rm -f *.[ox]
