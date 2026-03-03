gfortran riemann.f90 -O2 $(pkg-config --cflags fortran_stdlib) $(pkg-config --libs fortran_stdlib) -o riemann
./riemann