#ifndef FLUX_H
#define FLUX_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "lattice.h"

/*
 * Flux struct
 */
typedef struct Flux
{
  double K;

  int *flux;
  double _Complex* H_majorana;

  double _Complex* U;
  double *E;

  double *fermi_func;

  Lattice *latt;
} Flux;

/*
 * Initialises flux struct
 */
void init_flux(double K, Lattice *latt, Flux *flux_conf);

/*
 * Sets flux
 */
void set_flux(int *flux, Flux *flux_conf);

/*
 * Diagonalise Hamiltonian
 */
void diagonalise_flux(Flux *flux_conf);

/*
 * Returns flux weight for inverse temperature beta
 */
double weight(double beta, Flux *flux_conf);

/*
 * Returns the logarithm of the flux weight for inverse temperature beta
 */
double log_weight(double beta, Flux *flux_conf);

/*
 * Fermi function for single particle eigenstate k at inverse temperature beta
 */
double fermi_function(double beta, double Ek);

/*
 * Frees flux struct
 */
void free_flux(Flux *flux_conf);

#endif // FLUX_H
