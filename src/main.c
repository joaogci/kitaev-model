#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <stdbool.h>
#include <omp.h>

#include "monte_carlo.h"
#include "flux.h"
#include "replica.h"
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
Flux flux_conf;
int *flux;
double K[3];

// Observables
Obs_scalar* obs_scal;
Obs_latt* obs_eq;
int n_scal = 2;
int n_eq = 0;

// Replica Exchange 
Replica *replicas;
Flux *flux_replicas;
bool replica_exchange;
int n_exchange_steps;
int n_exchange_frequency;
int n_replica_walkers;
double beta_init = 10.0;


int main(int argc, char **argv)
{
  int n, i, b;
  long t;
  char filename[BUFFER_SIZE];
  pcg32_random_t rng;
  clock_t start_clock, end_clock;
  time_t time_;
  Flux new_flux_conf;

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

  // Replica Exchange
  replica_exchange = false;
  n_replica_walkers = 1;

  if (beta >= 75.0) {
    read_parameters_replica(&n_exchange_steps, &n_exchange_frequency, &n_replica_walkers);
    replica_exchange = true;

    replicas = (Replica *) malloc(n_replica_walkers * sizeof(Replica));
    flux_replicas = (Flux *) malloc(n_replica_walkers * sizeof(Flux));
    for (i = 0; i < n_replica_walkers; i++) {
      init_flux(K, &latt, &(flux_replicas[i]));
      init_replica(&(flux_replicas[i]), beta_init + i * (beta - beta_init) / (n_replica_walkers - 1), &(replicas[i]));
    }
  }
  
  omp_set_num_threads(n_replica_walkers);

  time_ = time(NULL);
  printf("Starting Finite T Kitaev simulation\n");
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
  if (!replica_exchange) {
    set_flux(flux, &flux_conf);
    diagonalise_flux(&flux_conf);
  } else {
    for (i = 0; i < n_replica_walkers; i++) {
      set_flux(flux, replicas[i].flux);
      diagonalise_flux(replicas[i].flux);
    }
  }

  // Monte Carlo Sampling
  start_clock = clock();

  if (!replica_exchange) {
    // Normal Monte Carlo
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
  } else {
    // Replica Exchange Monte Carlo for low temperatures
    #pragma omp parallel default(shared) private(obs_scal, obs_eq, rng, new_flux_conf, n, t, i, b)
    {
      int n_ex, i_ex;
      int t_id;

      t_id = omp_get_thread_num();
      pcg32_srandom_r(&rng, SEED * (t_id + 1) ^ (intptr_t)&rng, SEED * (t_id + 1));

      init_flux(K, &latt, &new_flux_conf);

      obs_scal = (Obs_scalar*) malloc(n_scal * sizeof(Obs_scalar));
      obs_eq = (Obs_latt*) malloc(n_eq * sizeof(Obs_latt));
      set_observables(obs_scal, n_scal, obs_eq, n_eq, &latt);

      for (n = 0; n < sim.n_bins; n++) {
        if (t_id == n_replica_walkers - 1) {
          reset_observables(obs_scal, n_scal, obs_eq, n_eq);
        }

        for (t = 0; t < sim.mc_sweeps; t++) {
          for (i = 0; i < latt.N; i++) {
            b = pcg32_random_r(&rng) % latt.Nb;
            monte_carlo_step(b, replicas[t_id].beta, pcg32_double_r(&rng), replicas[t_id].flux, &new_flux_conf);
          }

          if (t_id == n_replica_walkers - 1 && (n_scal > 0 || n_eq > 0)) {
            sample(obs_scal, n_scal, obs_eq, n_eq, replicas[t_id].beta, replicas[t_id].flux);
          }

          // Exchange configurations
          if (t % n_exchange_frequency == 0) {
            #pragma omp barrier
            #pragma omp master 
            {
              for (n_ex = 0; n_ex < n_exchange_steps; n_ex++) {
                i_ex = pcg32_random_r(&rng) % (n_replica_walkers - 1);
                exchange(pcg32_double_r(&rng), &(replicas[i_ex]), &(replicas[i_ex + 1]));
              }
            }
            #pragma omp barrier
          }
          
        }

        // Write bin to file
        if (t_id == n_replica_walkers - 1 && (n_scal > 0 || n_eq > 0)) {
          write_observables(obs_scal, n_scal, obs_eq, n_eq);
        }
        if (t_id == n_replica_walkers - 1) {
          write_configuration(&flux_conf);
        }
      }

      free(obs_scal);
      free(obs_eq);
      free_flux(&new_flux_conf);
    }
  }

  end_clock = clock();
  sim.wall_time = ((double) (end_clock - start_clock)) / (CLOCKS_PER_SEC * n_replica_walkers);

  write_sim_info(sim);

  printf("Simulation finished in %.5lfs \n", sim.wall_time);
  fflush(stdout);

  // FREE THE VARIABLES
  free_lattice(&latt);
  free_flux(&flux_conf);
  free(obs_eq);
  free(obs_scal);
  free(flux);

  if (replica_exchange) {
    free(replicas);
    free(flux_replicas);
  }

  exit(0);
}
