export PKG_CONFIG_PATH := /home/gustavo/.local/lib/pkgconfig:/home/gustavo/.local/lib64/pkgconfig:$(PKG_CONFIG_PATH)

STDLIB_CFLAGS := $(shell pkg-config --cflags fortran_stdlib)
STDLIB_LIBS   := $(shell pkg-config --libs   fortran_stdlib)

FC := gfortran
FFLAGS := -O2 -Wall

EXE := riemann
SRC := riemann.f90

all: $(EXE)

$(EXE): $(SRC:.f90=.o)
	$(FC) -o $@ $^ $(FFLAGS) $(STDLIB_LIBS)

%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(STDLIB_CFLAGS)

clean:
	rm -f *.o *.mod $(EXE)
