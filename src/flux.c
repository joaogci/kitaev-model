#include "flux.h"

void init_flux(double K, LatticeHyperbolic *latt, Flux *flux_conf)
{
  flux_conf->K = K;
  flux_conf->latt = latt;

  flux_conf->flux = (int *) malloc(flux_conf->latt->Nb * sizeof(int));

  flux_conf->H_majorana = (double _Complex*) malloc(flux_conf->latt->Ns * flux_conf->latt->Ns * sizeof(double));

  flux_conf->E = (double *) malloc(flux_conf->latt->N * sizeof(double));
  flux_conf->U = (double _Complex*) malloc(flux_conf->latt->Ns * flux_conf->latt->N * sizeof(double));
  flux_conf->fermi_func = (double *) malloc(flux_conf->latt->N * sizeof(double));
}

void set_flux(int *flux, Flux *flux_conf)
{
  int b, i, j;

  memcpy(flux_conf->flux, flux, flux_conf->latt->Nb * sizeof(int));
  memset(flux_conf->H_majorana, 0.0, flux_conf->latt->Ns * flux_conf->latt->Ns * sizeof(double _Complex));
  
  /* Build the Majorana Hamiltonian */
  for (b = 0; b < flux_conf->latt->Nb; b++) {
    i = flux_conf->latt->bond_list[b][0];
    j = flux_conf->latt->bond_list[b][1];

    flux_conf->H_majorana[i * flux_conf->latt->Ns + j] = 2.0 * flux_conf->K * I * flux_conf->flux[b];
    flux_conf->H_majorana[j * flux_conf->latt->Ns + i] = - 2.0 * flux_conf->K * I * flux_conf->flux[b];
  }
}

void diagonalise_flux(Flux *flux_conf)
{
  int i, j, j_;
  gsl_matrix_complex_view H_matter;
  gsl_vector *eval;
  gsl_matrix_complex *evec;
  gsl_eigen_hermv_workspace *ws;
  gsl_complex z;

  ws = gsl_eigen_hermv_alloc(flux_conf->latt->Ns);

  H_matter = gsl_matrix_complex_view_array((double *) flux_conf->H_majorana, flux_conf->latt->Ns, flux_conf->latt->Ns);
  eval = gsl_vector_alloc(flux_conf->latt->Ns);
  evec = gsl_matrix_complex_alloc(flux_conf->latt->Ns, flux_conf->latt->Ns);

  /* Diagonalise matrix */
  gsl_eigen_hermv(&H_matter.matrix, eval, evec, ws);

  gsl_eigen_hermv_free(ws);
  gsl_eigen_hermv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);

  /* Get the positive eigenvalues and eigenvectors */
  for (j = flux_conf->latt->N; j < flux_conf->latt->Ns; j++) {
    j_ = j - flux_conf->latt->N;

    flux_conf->E[j_] = gsl_vector_get(eval, j);

    for (i = 0; i < flux_conf->latt->Ns; i++) {
      z = gsl_matrix_complex_get(evec, i, j);
      flux_conf->U[i * flux_conf->latt->N + j_] = GSL_REAL(z) + I * GSL_IMAG(z);
    }
  }

  gsl_vector_free(eval);
  gsl_matrix_complex_free(evec);
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
  free(flux_conf->H_majorana);
  free(flux_conf->U);
  free(flux_conf->fermi_func);
  free(flux_conf->E);
  free(flux_conf->flux);
}
