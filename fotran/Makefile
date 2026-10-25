FC     := gfortran
FFLAGS := -O3 -march=native -fopenmp -Wall -Wextra

all: riemann_tol riemann_traj riemann_dn

riemann_lib.o: riemann_lib.f90
	$(FC) $(FFLAGS) -c $<

riemann_tol: riemann_tol.f90 riemann_lib.o
	$(FC) $(FFLAGS) -o $@ $^

riemann_traj: riemann_traj.f90 riemann_lib.o
	$(FC) $(FFLAGS) -o $@ $^

riemann_dn: riemann_dn.f90 riemann_lib.o
	$(FC) $(FFLAGS) -o $@ $^ -llapack -lblas

clean:
	rm -f *.o *.mod riemann_tol riemann_traj riemann_dn
