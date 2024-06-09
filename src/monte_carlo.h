#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "flux.h"

/*
 * Makes one Monte Carlo step
 */
void monte_carlo_step(int b, double beta, double r, Flux *flux_conf, Flux *new_flux_conf);

/*
 * Calculates the Monte Carlo ratio between F and F^\prime
 */
double ratio(double beta, Flux *flux_conf, Flux *new_flux_conf);

#endif // MONTE_CARLO_H
