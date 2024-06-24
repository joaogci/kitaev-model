#ifndef SAMPLING_H
#define SAMPLING_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "flux.h"
#include "observables.h"
#include "lattice.h"

/*
 * Sets the obervables to sample
 */
void set_observables(Obs_scalar *obs_scalar, int n_scal, Obs_latt *obs_eq, int n_eq, Obs_spectral *obs_spec, int n_spec, Lattice *latt);

/*
 * Reset observables
 */
void reset_observables(Obs_scalar *obs_scalar, int n_scal, Obs_latt *obs_eq, int n_eq, Obs_spectral *obs_spec, int n_spec);

/*
 * Samples the set obsevables
 */
void sample(Obs_scalar *obs_scal, int n_scal, Obs_latt *obs_eq, int n_eq, Obs_spectral *obs_spec, int n_spec, double beta, Flux *flux_conf);

/*
 * Sample scalar observables
 */ 
void sample_obs_scalar(Obs_scalar *obs, int n_scal, double beta, Flux *flux_conf);

/*
  *Sample equal time observables
 */
void sample_obs_eq(Obs_latt *obs, int n_eq, double beta, Flux *flux_conf);

/*
  *Sample equal time observables
 */
void sample_obs_spec(Obs_spectral *obs, int n_spec, double beta, Flux *flux_conf);

/*
 * Free the observable structs
 */
void free_observables(Obs_scalar *obs_scalar, int n_scal, Obs_latt *obs_eq, int n_eq, Obs_spectral *obs_spec, int n_spec);

double _Complex delta(double om, double eta, double dE);

#endif // SAMPLING_H
