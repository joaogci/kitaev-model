#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <stdbool.h>

#include "monte_carlo.h"
#include "flux.h"
#include "rng/pcg_basic.h"
#include "sampling.h"
#include "observables.h"
#include "lattice.h"
#include "io/io.h"

#ifndef TESTING
#define SEED (u_int64_t) time(NULL)
#else
#define SEED (u_int64_t) 2
#endif

// Lattice
int Lx, Ly;
Lattice latt;

// Simulation
Sim_info sim;
double beta;

// Hamiltonian
Flux flux_conf, new_flux_conf;
int *flux;
double K[3];

// Observables
Obs_scalar* obs_scal;
Obs_latt* obs_eq;
int n_scal = 2;
int n_eq = 0;


int main(int argc, char **argv)
{
  int n, i, b;
  long t;
  char filename[BUFFER_SIZE];
  pcg32_random_t rng;
  clock_t start_clock, end_clock;
  time_t time_;

  // Read and process inputs
  if (argc < 1) {
    // printf("Please provide the number of threads and the path to the main directory. \n");
    // printf("Usage: %s MC_DIR \n", argv[0]);
    exit(1);
  }

  read_parameters(&Lx, &Ly, &beta, &(sim.n_bins), &(sim.mc_sweeps), &(K[0]), &(K[1]), &(K[2]));

  // Lattice
  make_lattice(Lx, Ly, &latt);

  // Hamiltonian
  init_flux(K, &latt, &flux_conf);
  init_flux(K, &latt, &new_flux_conf);
  flux = (int *) malloc(latt.Nb * sizeof(double));

  // Rng
  pcg32_srandom_r(&rng, SEED ^ (intptr_t)&rng, SEED);

  // Observables
  obs_scal = (Obs_scalar*) malloc(n_scal * sizeof(Obs_scalar));
  obs_eq = (Obs_latt*) malloc(n_eq * sizeof(Obs_latt));
  set_observables(obs_scal, n_scal, obs_eq, n_eq, &latt);

  time_ = time(NULL);
  printf("Starting SSE simulation\n");
  printf("Simulation started at: %s \n", ctime(&time_));
  fflush(stdout);

  // Read confin file
  sprintf(filename, "confin");
  if (access(filename, F_OK) == 0) {
    read_configuration(flux);
    printf("Configration file read. \n");
  } else {
    for (n = 0; n < latt.Nb; n++) {
      flux[n] = 1;
    }
  }
  set_flux(flux, &flux_conf);
  diagonalise_flux(&flux_conf);

  // Monte Carlo Sampling
  start_clock = clock();

  for (n = 0; n < sim.n_bins; n++) {
    reset_observables(obs_scal, n_scal, obs_eq, n_eq);

    for (t = 0; t < sim.mc_sweeps; t++) {
      for (i = 0; i < latt.N; i++) {
        b = pcg32_random_r(&rng) % latt.Nb;
        monte_carlo_step(b, beta, pcg32_double_r(&rng), &flux_conf, &new_flux_conf);
      }

      if (n_scal > 0 || n_eq > 0) {
        sample(obs_scal, n_scal, obs_eq, n_eq, beta, &flux_conf);
      }
    }

    // Write bin to file
    if (n_scal > 0 || n_eq > 0) {
      write_observables(obs_scal, n_scal, obs_eq, n_eq);
    }
    write_configuration(&flux_conf);
  }

  end_clock = clock();
  sim.wall_time = ((double) (end_clock - start_clock)) / (CLOCKS_PER_SEC);

  write_sim_info(sim);

  printf("Simulation finished in %.5lfs \n", sim.wall_time);
  fflush(stdout);

  // FREE THE VARIABLES
  free_lattice(&latt);
  free_flux(&flux_conf);
  free(obs_eq);
  free(obs_scal);

  exit(0);
}
