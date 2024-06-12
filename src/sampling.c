#include "sampling.h"

void set_observables(Obs_scalar *obs_scalar, int n_scal, Obs_latt *obs_eq, int n_eq, Lattice *latt)
{
  int n;

  for (n = 0; n < n_scal; n++) {
    switch (n)
    {
    case 0:
      init_obs_scalar("Sz", &(obs_scalar[n]));
      break;
    case 1:
      init_obs_scalar("W", &(obs_scalar[n]));
      break;
    case 2:
      init_obs_scalar("E", &(obs_scalar[n]));
      break;
    case 3:
      init_obs_scalar("E2", &(obs_scalar[n]));
      break;
    case 4:
      init_obs_scalar("dE_dbeta", &(obs_scalar[n]));
      break;
    default:
      printf("Observable not found. \n");
      exit(1);
      break;
    }
  }

  for (n = 0; n < n_eq; n++) {
    switch (n)
    {
    case 0:
      init_obs_latt("Sz", latt, &(obs_eq[n]));
      break;
    default:
      printf("Observable not found. \n");
      exit(1);
      break;
    }
  }
}

void reset_observables(Obs_scalar *obs_scalar, int n_scal, Obs_latt *obs_eq, int n_eq)
{
  int n;

  for (n = 0; n < n_scal; n++) {
    reset_obs_scalar(&(obs_scalar[n]));
  }

  for (n = 0; n < n_eq; n++) {
    reset_obs_latt(&(obs_eq[n]));
  }
}

void sample(Obs_scalar *obs_scal, int n_scal, Obs_latt *obs_eq, int n_eq, double beta, Flux *flux_conf)
{
  int i; 
  for (i = 0; i < flux_conf->latt->N * flux_conf->latt->N; i++) {
    flux_conf->A[i] = flux_conf->XT[i] + flux_conf->YT[i];
    flux_conf->B[i] = flux_conf->XT[i] - flux_conf->YT[i];
  }

  if (n_scal > 0) {
    sample_obs_scalar(obs_scal, n_scal, beta, flux_conf);
  }

  if (n_eq > 0) {
    sample_obs_eq(obs_eq, n_eq, beta, flux_conf);
  }
}

void sample_obs_scalar(Obs_scalar *obs, int n_scal, double beta, Flux *flux_conf)
{
  int i, j, b, m, p, n;
  int Fij;
  double res, tmp;
  double E, dEdb;

  for (i = 0; i < n_scal; i++) {
    obs[i].N++;
  }

  res = 0.0;
  for (p = 0; p < flux_conf->latt->N; p++) {
    b = flux_conf->latt->z_bonds[p];
    i = flux_conf->latt->r_sites[flux_conf->latt->bond_list[b][0]];
    j = flux_conf->latt->r_sites[flux_conf->latt->bond_list[b][1]];
    Fij = flux_conf->flux[b];
    for (m = 0; m < flux_conf->latt->N; m++) {
      res += Fij * flux_conf->A[i * flux_conf->latt->N + m] * flux_conf->B[j * flux_conf->latt->N + m] * (1.0 - 2.0 * fermi_function(beta, flux_conf->E[m]));
    }
  }
  obs[0].obs_vec += res / flux_conf->latt->N;

  res = 0.0;
  for (p = 0; p < flux_conf->latt->Np; p++) {
    tmp = 1.0;
    for (b = 0; b < flux_conf->latt->size_p; b++) {
      tmp *= flux_conf->flux[flux_conf->latt->plaquette_bonds[p][b]];
    }
    res += tmp;
  }
  obs[1].obs_vec += res / flux_conf->latt->N;

  E = 0.0;
  dEdb = 0.0;
  for (n = 0; n < flux_conf->latt->N; n++) {
    // E += flux_conf->E[n] * fermi_function(beta, flux_conf->E[n]);
    // dEdb += - flux_conf->E[n] * flux_conf->E[n] * fermi_function(beta, flux_conf->E[n]) * (1.0 - fermi_function(beta, flux_conf->E[n]));
    E += - flux_conf->E[n] * tanh(beta * flux_conf->E[n] * 0.5) * 0.5;
    dEdb += - flux_conf->E[n] * flux_conf->E[n] * 0.25 / (cosh(beta * flux_conf->E[n] * 0.5) * cosh(beta * flux_conf->E[n] * 0.5));
  }
  obs[2].obs_vec += E;
  obs[3].obs_vec += E * E;
  obs[4].obs_vec += dEdb;
}

void sample_obs_eq(Obs_latt *obs, int n_eq, double beta, Flux *flux_conf)
{
  // int i, j; // , a, b, p;

  // for (i = 0; i < n_eq; i++) {
  //   obs[i].N++;
  // }
}

void free_observables(Obs_scalar *obs_scalar, int n_scal, Obs_latt *obs_eq, int n_eq)
{
  int n;

  for (n = 0; n < n_eq; n++) {
    free_obs_latt(&(obs_eq[n]));
  }
}
