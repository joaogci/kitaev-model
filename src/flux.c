#include "flux.h"

void init_flux(double K[3], Lattice *latt, Flux *flux_conf)
{
  flux_conf->K[0] = K[0];
  flux_conf->K[1] = K[1];
  flux_conf->K[2] = K[2];

  flux_conf->latt = latt;

  flux_conf->flux = (int *) malloc(flux_conf->latt->Nb * sizeof(int));

  flux_conf->H_matter = (double *) malloc(flux_conf->latt->Ns * flux_conf->latt->Ns * sizeof(double));

  flux_conf->E = (double *) malloc(flux_conf->latt->N * sizeof(double));
  flux_conf->F = (double *) malloc(flux_conf->latt->N * flux_conf->latt->N * sizeof(double));
  flux_conf->XT = (double *) malloc(flux_conf->latt->N * flux_conf->latt->N * sizeof(double));
  flux_conf->YT = (double *) malloc(flux_conf->latt->N * flux_conf->latt->N * sizeof(double));
  flux_conf->A = (double *) malloc(flux_conf->latt->N * flux_conf->latt->N * sizeof(double));
  flux_conf->B = (double *) malloc(flux_conf->latt->N * flux_conf->latt->N * sizeof(double));
  flux_conf->fermi_func = (double *) malloc(flux_conf->latt->N * sizeof(double));
}

void set_flux(int *flux, Flux *flux_conf)
{
  int b, i, j, i_, j_;

  memcpy(flux_conf->flux, flux, flux_conf->latt->Nb * sizeof(int));
  memset(flux_conf->F, 0.0, flux_conf->latt->N * flux_conf->latt->N * sizeof(double));

  /* Set flux matrix */
  for (b = 0; b < flux_conf->latt->Nb; b++) {
    i = flux_conf->latt->r_sites[flux_conf->latt->bond_list[b][0]]; 
    j = flux_conf->latt->r_sites[flux_conf->latt->bond_list[b][1]];

    flux_conf->F[i * flux_conf->latt->N + j] = flux_conf->K[flux_conf->latt->bond_type[b]] * flux_conf->flux[b];
  }    

  /* Construct the block Hamiltonian */
  memset(flux_conf->H_matter, 0.0, flux_conf->latt->Ns * flux_conf->latt->Ns * sizeof(double));

  for (i = 0; i < flux_conf->latt->N; i++) {
    i_ = i;

    for (j = 0; j < flux_conf->latt->N; j++){
      j_ = j;
      flux_conf->H_matter[i * flux_conf->latt->Ns + j] = flux_conf->F[j_ * flux_conf->latt->N + i_] + flux_conf->F[i_ * flux_conf->latt->N + j_];
    }

    for (j = flux_conf->latt->N; j < flux_conf->latt->Ns; j++) {
      j_ = j % flux_conf->latt->N;
      flux_conf->H_matter[i * flux_conf->latt->Ns + j] = flux_conf->F[j_ * flux_conf->latt->N + i_] - flux_conf->F[i_ * flux_conf->latt->N + j_];
    }
  }

  for (i = flux_conf->latt->N; i < flux_conf->latt->Ns; i++) {
    i_ = i % flux_conf->latt->N;

    for (j = 0; j < flux_conf->latt->N; j++){
      j_ = j;
      flux_conf->H_matter[i * flux_conf->latt->Ns + j] = flux_conf->F[i_ * flux_conf->latt->N + j_] - flux_conf->F[j_ * flux_conf->latt->N + i_];
    }

    for (j = flux_conf->latt->N; j < flux_conf->latt->Ns; j++) {
      j_ = j % flux_conf->latt->N;
      flux_conf->H_matter[i * flux_conf->latt->Ns + j] = - (flux_conf->F[j_ * flux_conf->latt->N + i_] + flux_conf->F[i_ * flux_conf->latt->N + j_]);
    }
  }
}

void diagonalise_flux(Flux *flux_conf)
{
  int i, j, i_, j_;
  gsl_matrix_view H_matter;
  gsl_vector *eval;
  gsl_matrix *evec;
  gsl_eigen_symmv_workspace *ws;

  ws = gsl_eigen_symmv_alloc(flux_conf->latt->Ns);

  H_matter = gsl_matrix_view_array(flux_conf->H_matter, flux_conf->latt->Ns, flux_conf->latt->Ns);
  eval = gsl_vector_alloc(flux_conf->latt->Ns);
  evec = gsl_matrix_alloc(flux_conf->latt->Ns, flux_conf->latt->Ns);

  /* Diagonalise matrix */
  gsl_eigen_symmv(&H_matter.matrix, eval, evec, ws);

  gsl_eigen_symmv_free(ws);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);

  /* Get BdG form */
  for (j = flux_conf->latt->N; j < flux_conf->latt->Ns; j++) {
    j_ = j % flux_conf->latt->N;

    flux_conf->E[j_] = gsl_vector_get(eval, j);

    for (i = 0; i < flux_conf->latt->N; i++) {
      i_ = i % flux_conf->latt->N;
      flux_conf->XT[i_ * flux_conf->latt->N + j_] = gsl_matrix_get(evec, i, j);
    }
    for (i = flux_conf->latt->N; i < flux_conf->latt->Ns; i++) {
      i_ = i % flux_conf->latt->N;
      flux_conf->YT[i_ * flux_conf->latt->N + j_] = gsl_matrix_get(evec, i, j);
    }
  }

  gsl_vector_free(eval);
  gsl_matrix_free(evec);
}

double weight(double beta, Flux *flux_conf)
{
  int i;
  double weight;

  weight = 1.0;
  for (i = 0; i < flux_conf->latt->N; i++) {
    weight *= 2.0 * cosh(0.5 * beta * flux_conf->E[i]);
  }

  return weight;
}

double log_weight(double beta, Flux *flux_conf)
{
  int i;
  double log_weight;

  log_weight = 0.0;
  for (i = 0; i < flux_conf->latt->N; i++) {
    if (beta * flux_conf->E[i] <= 100.0) {
      log_weight += log(exp(beta * flux_conf->E[i]) + 1.0) - 0.5 * beta * flux_conf->E[i];
    } else {
      log_weight += 0.5 * beta * flux_conf->E[i];
    }
  }

  return log_weight;
}

double fermi_function(double beta, double Ek)
{
  if (beta * Ek <= 100.0) {
    return 1.0 / (1.0 + exp(beta * Ek));
  } else {
    return 0.0;
  }
}

void free_flux(Flux *flux_conf)
{
  free(flux_conf->H_matter);
  free(flux_conf->F);
  free(flux_conf->XT);
  free(flux_conf->YT);
  free(flux_conf->A);
  free(flux_conf->B);
  free(flux_conf->fermi_func);
  free(flux_conf->E);
  free(flux_conf->flux);
}
