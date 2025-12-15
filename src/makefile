INCLUDES_MC = flux.h observables.h sampling.h \
	lattice_hyperbolic.h io/io.h rng/pcg_basic.h monte_carlo.h \
	replica.h
SOURCES_MC = main.c flux.c observables.c sampling.c \
  lattice_hyperbolic.c io/io.c rng/pcg_basic.c monte_carlo.c \
	replica.c
INCLUDES_ANA = lattice_hyperbolic.h io/io.h analysis.h observables.h
SOURCES_ANA = ana.c lattice_hyperbolic.c io/io.c analysis.c observables.c

FLAGS = -O2 -fopenmp
DEBUG_FLAGS = -Wall -Wshadow
LIBS = -lgsl -lcblas -lm -lgfortran

ifeq ($(OS),Windows_NT)
  $(Windows is not supported. Use a UNIX OS.)
else
  UNAME_S := $(shell uname -s)
  ifeq ($(UNAME_S),Linux)
    CC = gcc
    FORT = gfortran
  endif
  ifeq ($(UNAME_S),Darwin)
    CC = gcc-15
    FORT = gfortran
  endif
endif

all: main ana

main: read_parameters $(SOURCES_MC) $(INCLUDES_MC)
	$(CC) -o main $(FLAGS) $(DEBUG_FLAGS) io/read_parameters.o io/read_parameters_analysis.o io/read_parameters_replica.o $(SOURCES_MC) $(LIBS)

ana: read_parameters $(SOURCES_ANA) $(INCLUDES_ANA) 
	$(CC) -o ana $(FLAGS) $(DEBUG_FLAGS) io/read_parameters.o io/read_parameters_analysis.o io/read_parameters_replica.o $(SOURCES_ANA) $(LIBS)

read_parameters: io/read_parameters.f90 io/read_parameters_analysis.f90
	$(FORT) -o io/read_parameters.o -c io/read_parameters.f90
	$(FORT) -o io/read_parameters_analysis.o -c io/read_parameters_analysis.f90
	$(FORT) -o io/read_parameters_replica.o -c io/read_parameters_replica.f90

clean: 
	rm io/read_parameters.o 
	rm io/read_parameters_analysis.o 
	rm io/read_parameters_replica.o 
	rm main
	rm ana
