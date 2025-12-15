#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../sampling.h"
#include "../flux.h"

#define BUFFER_SIZE 256

typedef struct Sim_info {
  long mc_sweeps;
  int n_bins;
  double wall_time;
} Sim_info;

typedef struct Ana_info {
  int n_skip;
  int n_rebin;
} Ana_info;

/*
 *  reads input file for the simulation - Fortran function 
 */
void read_parameters(int *N, double* beta, int* n_bins, long* mc_sweeps, double* alpha);

/* 
 * reads parameters for analysis
 */
void read_parameters_analysis(int *n_skip, int* n_rebin);

 /* 
 * reads parameters for replica exchange
 */
void read_parameters_replica(int *n_exchange_steps, int* n_exchange_frequency, int *n_replica_walkers);

/*
 * Reads hyperbolic lattice
 */
void read_hyperbolic_lattice(int N, int *Nb, int*** adj, char* sse_path);

/*
 * writes observables
 */
void write_observables(Obs_scalar* obs_scal, int n_scal, Obs_latt* obs_eq, int n_eq, Obs_spectral *obs_spec, int n_spec);

// /*
//  * Writes transport observables
//  */
// void write_transport_obeservables(Obs_transport* obs, int n_transp);

/*
 * writes simulation statistics to file
 */
void write_sim_info(Sim_info sim);

/*
 * Writes SSE configuration
 */
void write_configuration(Flux* flux_conf);

/*
 * Reads SSE configuration
 */
void read_configuration(int* flux);

/*
 * Returns the number of lines of a file
 */
int num_lines(FILE* fp);

#endif // IO_H
