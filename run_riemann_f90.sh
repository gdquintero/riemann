gfortran riemann.f90 -O2 \
  -I$HOME/.local/include -L$HOME/.local/lib -lquadpack \
  $(pkg-config --cflags fortran_stdlib) $(pkg-config --libs fortran_stdlib) \
  -o riemann
./riemann