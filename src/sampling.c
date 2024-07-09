#include "sampling.h"

void set_observables(Obs_scalar *obs_scalar, int n_scal, Obs_latt *obs_eq, int n_eq, Obs_spectral *obs_spec, int n_spec, Lattice *latt)
{
  int n;

  for (n = 0; n < n_scal; n++) {
    switch (n)
    {
    case 0:
      init_obs_scalar("SpinZZ", &(obs_scalar[n]));
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
    case 5:
      init_obs_scalar("SpinYY", &(obs_scalar[n]));
      break;
    case 6:
      init_obs_scalar("SpinXX", &(obs_scalar[n]));
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
      init_obs_latt("DimerXX", latt, &(obs_eq[n]));
      break;
    case 1:
      init_obs_latt("DimerXY", latt, &(obs_eq[n]));
      break;
    case 2:
      init_obs_latt("DimerXZ", latt, &(obs_eq[n]));
      break;
    case 3:
      init_obs_latt("DimerYY", latt, &(obs_eq[n]));
      break;
    case 4:
      init_obs_latt("DimerYZ", latt, &(obs_eq[n]));
      break;
    case 5:
      init_obs_latt("DimerZZ", latt, &(obs_eq[n]));
      break;
    default:
      printf("Observable not found. \n");
      exit(1);
      break;
    }
  }

  for (n = 0; n < n_spec; n++) {
    switch (n)
    {
    case 0:
      init_obs_spectral("DimerXX", 200, 4.0, latt, &(obs_spec[n]));
      break;
    case 1:
      init_obs_spectral("DimerXY", 200, 4.0, latt, &(obs_spec[n]));
      break;
    case 2:
      init_obs_spectral("DimerXZ", 200, 4.0, latt, &(obs_spec[n]));
      break;
    case 3:
      init_obs_spectral("DimerYY", 200, 4.0, latt, &(obs_spec[n]));
      break;
    case 4:
      init_obs_spectral("DimerYZ", 200, 4.0, latt, &(obs_spec[n]));
      break;
    case 5:
      init_obs_spectral("DimerZZ", 200, 4.0, latt, &(obs_spec[n]));
      break;
    default:
      printf("Observable not found. \n");
      exit(1);
      break;
    }
  }
}

void reset_observables(Obs_scalar *obs_scalar, int n_scal, Obs_latt *obs_eq, int n_eq, Obs_spectral *obs_spec, int n_spec)
{
  int n;

  for (n = 0; n < n_scal; n++) {
    reset_obs_scalar(&(obs_scalar[n]));
  }

  for (n = 0; n < n_eq; n++) {
    reset_obs_latt(&(obs_eq[n]));
  }

  for (n = 0; n < n_spec; n++) {
    reset_obs_spectral(&(obs_spec[n]));
  }
}

void sample(Obs_scalar *obs_scal, int n_scal, Obs_latt *obs_eq, int n_eq, Obs_spectral *obs_spec, int n_spec, double beta, Flux *flux_conf)
{
  int i; 
  for (i = 0; i < flux_conf->latt->N * flux_conf->latt->N; i++) {
    flux_conf->A[i] = flux_conf->XT[i] + flux_conf->YT[i];
    flux_conf->B[i] = flux_conf->XT[i] - flux_conf->YT[i];
  }
  for (i = 0; i < flux_conf->latt->N; i++) {
    flux_conf->fermi_func[i] = fermi_function(beta, flux_conf->E[i]);
  }

  if (n_scal > 0) {
    sample_obs_scalar(obs_scal, n_scal, beta, flux_conf);
  }

  if (n_eq > 0) {
    sample_obs_eq(obs_eq, n_eq, beta, flux_conf);
  }

  if (n_spec > 0) {
    sample_obs_spec(obs_spec, n_spec, beta, flux_conf);
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
      res += Fij * flux_conf->A[i * flux_conf->latt->N + m] * flux_conf->B[j * flux_conf->latt->N + m] * (1.0 - 2.0 * flux_conf->fermi_func[m]);
    }
  }
  obs[0].obs_vec += res / flux_conf->latt->N;

  res = 0.0;
  for (p = 0; p < flux_conf->latt->N; p++) {
    b = flux_conf->latt->y_bonds[p];
    i = flux_conf->latt->r_sites[flux_conf->latt->bond_list[b][0]];
    j = flux_conf->latt->r_sites[flux_conf->latt->bond_list[b][1]];
    Fij = flux_conf->flux[b];
    for (m = 0; m < flux_conf->latt->N; m++) {
      res += Fij * flux_conf->A[i * flux_conf->latt->N + m] * flux_conf->B[j * flux_conf->latt->N + m] * (1.0 - 2.0 * flux_conf->fermi_func[m]);
    }
  }
  obs[5].obs_vec += res / flux_conf->latt->N;

  res = 0.0;
  for (p = 0; p < flux_conf->latt->N; p++) {
    b = flux_conf->latt->x_bonds[p];
    i = flux_conf->latt->r_sites[flux_conf->latt->bond_list[b][0]];
    j = flux_conf->latt->r_sites[flux_conf->latt->bond_list[b][1]];
    Fij = flux_conf->flux[b];
    for (m = 0; m < flux_conf->latt->N; m++) {
      res += Fij * flux_conf->A[i * flux_conf->latt->N + m] * flux_conf->B[j * flux_conf->latt->N + m] * (1.0 - 2.0 * flux_conf->fermi_func[m]);
    }
  }
  obs[6].obs_vec += res / flux_conf->latt->N;

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
    // dEdb += - flux_conf->E[n] * flux_conf->E[n] * (fermi_function(beta, flux_conf->E[n]) * fermi_function(beta, flux_conf->E[n]));
    E += flux_conf->E[n] * 0.5 * tanh(beta * flux_conf->E[n] * 0.5);
    dEdb += - flux_conf->E[n] * flux_conf->E[n] * 0.25 / (cosh(beta * flux_conf->E[n] * 0.5) * cosh(beta * flux_conf->E[n] * 0.5));
  }
  obs[2].obs_vec += E;
  obs[3].obs_vec += E * E;
  obs[4].obs_vec += dEdb;
}

void sample_obs_eq(Obs_latt *obs, int n_eq, double beta, Flux *flux_conf)
{
  int a_, b_, bonds_a[flux_conf->latt->N], bonds_b[flux_conf->latt->N], n_;
  int i, j, ip, jp, b, bp, k, q, p, pp;
  int Fij, Fijp;
  double res_PP, res_PH;

  for (i = 0; i < n_eq; i++) {
    obs[i].N++;
  }

  n_ = 0;
  for (a_ = 0; a_ < 3; a_++) {
    if (a_ == 0) {
      memcpy(bonds_a, flux_conf->latt->x_bonds, flux_conf->latt->N * sizeof(int));
    } else if (a_ == 1) {
      memcpy(bonds_a, flux_conf->latt->y_bonds, flux_conf->latt->N * sizeof(int));
    } else if (a_ == 2) {
      memcpy(bonds_a, flux_conf->latt->z_bonds, flux_conf->latt->N * sizeof(int));
    }
    for (b_ = a_; b_ < 3; b_++) {
      if (b_ == 0) {
        memcpy(bonds_b, flux_conf->latt->x_bonds, flux_conf->latt->N * sizeof(int));
      } else if (b_ == 1) {
        memcpy(bonds_b, flux_conf->latt->y_bonds, flux_conf->latt->N * sizeof(int));
      } else if (b_ == 2) {
        memcpy(bonds_b, flux_conf->latt->z_bonds, flux_conf->latt->N * sizeof(int));
      }

      for (p = 0; p < 1; p++) { // 
        b = bonds_a[p];
        i = flux_conf->latt->r_sites[flux_conf->latt->bond_list[b][0]];
        j = flux_conf->latt->r_sites[flux_conf->latt->bond_list[b][1]];
        Fij = flux_conf->flux[b];

        for (pp = 0; pp < flux_conf->latt->N; pp++) { // flux_conf->latt->N
          bp = bonds_b[pp];
          ip = flux_conf->latt->r_sites[flux_conf->latt->bond_list[bp][0]];
          jp = flux_conf->latt->r_sites[flux_conf->latt->bond_list[bp][1]];
          Fijp = flux_conf->flux[bp];

          res_PP = 0.0;
          res_PH = 0.0;

          for (k = 0; k < flux_conf->latt->N; k++) {
            for (q = 0; q < flux_conf->latt->N; q++) {
              res_PP += flux_conf->A[i * flux_conf->latt->N + k] * flux_conf->B[j * flux_conf->latt->N + q] * (flux_conf->A[ip * flux_conf->latt->N + k] * flux_conf->B[jp * flux_conf->latt->N + q] - flux_conf->A[ip * flux_conf->latt->N + q] * flux_conf->B[jp * flux_conf->latt->N + k]) * (flux_conf->fermi_func[k]*flux_conf->fermi_func[q] + (1.0-flux_conf->fermi_func[k])*(1.0-flux_conf->fermi_func[q]));
              res_PH += flux_conf->A[i * flux_conf->latt->N + k] * flux_conf->B[j * flux_conf->latt->N + q] * (flux_conf->A[ip * flux_conf->latt->N + k] * flux_conf->B[jp * flux_conf->latt->N + q] + flux_conf->A[ip * flux_conf->latt->N + q] * flux_conf->B[jp * flux_conf->latt->N + k]) * (flux_conf->fermi_func[k]*(1.0-flux_conf->fermi_func[q]) + (1.0-flux_conf->fermi_func[k])*flux_conf->fermi_func[q]);
            }
          }
          obs[n_].obs_latt[p][pp] += Fij * Fijp * (res_PP + res_PH);
        }    
      }
      n_++;

    }
  }
}

void sample_obs_spec(Obs_spectral *obs, int n_spec, double beta, Flux *flux_conf)
{
  int a_, b_, bonds_a[flux_conf->latt->N], bonds_b[flux_conf->latt->N], n_;
  int i, j, ip, jp, b, bp, k, q, p, pp, n;
  int Fij, Fijp;
  double _Complex res_PP, res_PH;

  for (i = 0; i < n_spec; i++) {
    obs[i].N++;
  }

  n_ = 0;
  for (a_ = 0; a_ < 3; a_++) {
    if (a_ == 0) {
      memcpy(bonds_a, flux_conf->latt->x_bonds, flux_conf->latt->N * sizeof(int));
    } else if (a_ == 1) {
      memcpy(bonds_a, flux_conf->latt->y_bonds, flux_conf->latt->N * sizeof(int));
    } else if (a_ == 2) {
      memcpy(bonds_a, flux_conf->latt->z_bonds, flux_conf->latt->N * sizeof(int));
    }
    for (b_ = a_; b_ < 3; b_++) {
      if (b_ == 0) {
        memcpy(bonds_b, flux_conf->latt->x_bonds, flux_conf->latt->N * sizeof(int));
      } else if (b_ == 1) {
        memcpy(bonds_b, flux_conf->latt->y_bonds, flux_conf->latt->N * sizeof(int));
      } else if (b_ == 2) {
        memcpy(bonds_b, flux_conf->latt->z_bonds, flux_conf->latt->N * sizeof(int));
      }

      for (n = 0; n < obs[0].n_om; n++) {
        for (p = 0; p < 1; p++) { // flux_conf->latt->N
          b = bonds_a[p];
          i = flux_conf->latt->r_sites[flux_conf->latt->bond_list[b][0]];
          j = flux_conf->latt->r_sites[flux_conf->latt->bond_list[b][1]];
          Fij = flux_conf->flux[b];

          for (pp = 0; pp < flux_conf->latt->N; pp++) { // 
            bp = bonds_b[pp];
            ip = flux_conf->latt->r_sites[flux_conf->latt->bond_list[bp][0]];
            jp = flux_conf->latt->r_sites[flux_conf->latt->bond_list[bp][1]];
            Fijp = flux_conf->flux[bp];

            res_PP = 0.0;
            res_PH = 0.0;

            for (k = 0; k < flux_conf->latt->N; k++) {
              for (q = 0; q < flux_conf->latt->N; q++) {
                res_PP += flux_conf->A[i * flux_conf->latt->N + k] * flux_conf->B[j * flux_conf->latt->N + q] * 
                  (flux_conf->A[ip * flux_conf->latt->N + k] * flux_conf->B[jp * flux_conf->latt->N + q] - flux_conf->A[ip * flux_conf->latt->N + q] * flux_conf->B[jp * flux_conf->latt->N + k]) * 
                  (flux_conf->fermi_func[k]*flux_conf->fermi_func[q] * delta(obs[0].omega[n], obs[0].eta, - (flux_conf->E[q] + flux_conf->E[k])) +
                  (1.0-flux_conf->fermi_func[k])*(1.0-flux_conf->fermi_func[q]) * delta(obs[0].omega[n], obs[0].eta, (flux_conf->E[q] + flux_conf->E[k]))); 

                res_PH += flux_conf->A[i * flux_conf->latt->N + k] * flux_conf->B[j * flux_conf->latt->N + q] * 
                    (flux_conf->A[ip * flux_conf->latt->N + k] * flux_conf->B[jp * flux_conf->latt->N + q] + flux_conf->A[ip * flux_conf->latt->N + q] * flux_conf->B[jp * flux_conf->latt->N + k]) * 
                    (flux_conf->fermi_func[k]*(1.0-flux_conf->fermi_func[q]) * delta(obs[0].omega[n], obs[0].eta, (flux_conf->E[k] - flux_conf->E[q])) + 
                    (1.0-flux_conf->fermi_func[k])*flux_conf->fermi_func[q] * delta(obs[0].omega[n], obs[0].eta, - (flux_conf->E[k] - flux_conf->E[q])));
              }
            }
            obs[n_].obs_latt[n][p][pp] += - Fij * Fijp * cimag(res_PP + res_PH);
          }    
        }
      }
      n_++;

    }
  }
}

void free_observables(Obs_scalar *obs_scalar, int n_scal, Obs_latt *obs_eq, int n_eq, Obs_spectral *obs_spec, int n_spec)
{
  int n;

  for (n = 0; n < n_eq; n++) {
    free_obs_latt(&(obs_eq[n]));
  }

  for (n = 0; n < n_spec; n++) {
    free_obs_spectral(&(obs_spec[n]));
  }
}

double _Complex delta(double om, double eta, double dE)
{
  return  - (om - dE + I * eta) / ((om - dE) * (om - dE) + eta * eta);
}
