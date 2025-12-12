#include "monte_carlo.h"

void monte_carlo_step(int b, double beta, double r, Flux *flux_conf, Flux *new_flux_conf)
{
  int flux[flux_conf->latt->Nb];

  memcpy(flux, flux_conf->flux, flux_conf->latt->Nb * sizeof(int));
  
  // Flip flux on bond b
  flux[b] = - flux[b];
  set_flux(flux, new_flux_conf);
  diagonalise_flux(new_flux_conf);

  if (r <= ratio(beta, flux_conf, new_flux_conf)) {
    set_flux(flux, flux_conf);
    memcpy(flux_conf->U, new_flux_conf->U, flux_conf->latt->Ns * flux_conf->latt->N * sizeof(double _Complex));
    memcpy(flux_conf->E, new_flux_conf->E, flux_conf->latt->N * sizeof(double));
  }
}

double ratio(double beta, Flux *flux_conf, Flux *new_flux_conf)
{
  return exp(log_weight(beta, new_flux_conf) - log_weight(beta, flux_conf));
}
