#include "observables.h"

void init_obs_scalar(char* filename, Obs_scalar* obs)
{
  strcpy(obs->filename, filename);
}

void init_obs_latt(char* filename, Lattice* latt, Obs_latt* obs)
{
  int i;

  strcpy(obs->filename, filename);
  obs->latt = latt;

  obs->obs_latt = (double _Complex**) malloc(latt->N * sizeof(double _Complex*));
  for (i = 0; i < latt->N; i++) {
    obs->obs_latt[i] = (double _Complex*) malloc(latt->N * sizeof(double _Complex));
  }
  obs->obs_latt0 = (double _Complex*) malloc(latt->N * sizeof(double _Complex));

  obs->obs_i = (double _Complex*) malloc(latt->N * sizeof(double _Complex));
  obs->obs_k = (double _Complex*) malloc(latt->N * sizeof(double _Complex));
}

void init_obs_spectral(char* filename, int n_om, double om_f, Lattice *latt, Obs_spectral* obs)
{
  int n, i;

  strcpy(obs->filename, filename);
  obs->n_om = n_om;
  obs->om_f = om_f;
  obs->eta = om_f / (n_om - 1);
  obs->latt = latt;

  obs->omega = (double*) malloc(n_om * sizeof(double));
  for (n = 0; n < n_om; n++) {
    obs->omega[n] = n * om_f / (n_om - 1);
  }

  obs->obs_latt = (double _Complex***) malloc(obs->n_om * sizeof(double _Complex**));
  obs->obs_i = (double _Complex**) malloc(obs->n_om * sizeof(double _Complex*));
  obs->obs_k = (double _Complex**) malloc(obs->n_om * sizeof(double _Complex*));
  for (n = 0; n < obs->n_om; n++) {
    obs->obs_latt[n] = (double _Complex**) malloc(latt->N * sizeof(double _Complex*));
    obs->obs_k[n] = (double _Complex*) malloc(latt->N * sizeof(double _Complex));
    obs->obs_i[n] = (double _Complex*) malloc(latt->N * sizeof(double _Complex));

    for (i = 0; i < latt->N; i++) {
      obs->obs_latt[n][i] = (double _Complex*) malloc(latt->N * sizeof(double _Complex));
    }
  }
}

void reset_obs_scalar(Obs_scalar* obs)
{
  obs->N = 0;
  obs->obs_vec = 0.0;
}

void reset_obs_latt(Obs_latt* obs)
{
  int i, j;

  obs->N = 0;
  for (i = 0; i < obs->latt->N; i++) {
    for (j = 0; j < obs->latt->N; j++) {
      obs->obs_latt[i][j] = 0.0;
    }
    obs->obs_latt0[i] = 0.0;
    obs->obs_i[i] = 0.0;
    obs->obs_k[i] = 0.0;
  }
}

void reset_obs_spectral(Obs_spectral* obs)
{
  int n, i, j;

  obs->N = 0;
  for (n = 0; n < obs->n_om; n++) {
    for (i = 0; i < obs->latt->N; i++) {
      obs->obs_k[n][i] = 0.0;
      obs->obs_i[n][i] = 0.0;

      for (j = 0; j < obs->latt->N; j++) {
        obs->obs_latt[n][i][j] = 0.0;
      }
    }
  }
}

void write_obs_scalar(FILE* out, Obs_scalar* obs)
{ 
  obs->obs_vec = obs->obs_vec / obs->N;
  fprintf(out, "(%lf, %lf) \n", creal(obs->obs_vec), cimag(obs->obs_vec));
}

void write_obs_scalar_info(FILE* info, Obs_scalar* obs)
{
  fprintf(info, "-- Analysis Mode --\n");
  fprintf(info, "identity\n");
}

void write_obs_latt(FILE* out_i, FILE* out_k, Obs_latt* obs)
{
  int i, j;

  for (i = 0; i < obs->latt->N; i++) {
    for (j = 0; j < obs->latt->N; j++) {
      obs->obs_latt[i][j] = obs->obs_latt[i][j] / obs->N;
    }
    obs->obs_latt0[i] = obs->obs_latt0[i] / obs->N;
  }

  fourier_trans(obs);
  inv_fourier_trans(obs);

  for (i = 0; i < obs->latt->N; i++) {
    fprintf(out_i, "(%lf, %lf) (%lf, %lf) \n", (double) (obs->latt->r[i][0] * obs->latt->a_1[0] + obs->latt->r[i][1] * obs->latt->a_2[0]), 
                                               (double) (obs->latt->r[i][0] * obs->latt->a_1[1] + obs->latt->r[i][1] * obs->latt->a_2[1]), 
                                               creal(obs->obs_i[i]), 
                                               cimag(obs->obs_i[i]));
  }

  for (i = 0; i < obs->latt->N; i++) {
    fprintf(out_k, "(%lf, %lf) (%lf, %lf) \n", (double) (obs->latt->k[i][0] * obs->latt->b_1[0] + obs->latt->k[i][1] * obs->latt->b_2[0]), 
                                               (double) (obs->latt->k[i][0] * obs->latt->b_1[1] + obs->latt->k[i][1] * obs->latt->b_2[1]), 
                                               creal(obs->obs_k[i]), 
                                               cimag(obs->obs_k[i]));
  }
}

void write_obs_eq_info(FILE* info, Obs_latt* obs)
{
  fprintf(info, "-- Analysis Mode --\n");
  fprintf(info, "equal time\n");
  fprintf(info, "-- Lattice --\n");
  fprintf(info, "Lx: %d\n", obs->latt->Lx);
  fprintf(info, "Ly: %d\n", obs->latt->Ly);
}

void write_obs_spectral(FILE* out_i, FILE *out_k, Obs_spectral* obs)
{
  int n, i, j; 

  for (n = 0; n < obs->n_om; n++) {
    for (i = 0; i < obs->latt->N; i++) {
      for (j = 0; j < obs->latt->N; j++) {
        obs->obs_latt[n][i][j] = obs->obs_latt[n][i][j] / obs->N;
      }
    }
  }

  spec_fourier_trans(obs);
  spec_inv_fourier_trans(obs);

  for (n = 0; n < obs->n_om; n++) {
    for (i = 0; i < obs->latt->N; i++) {
      fprintf(out_i, "%lf (%lf, %lf) (%lf, %lf) \n", creal(obs->omega[n]), (double) (obs->latt->r[i][0] * obs->latt->a_1[0] + obs->latt->r[i][1] * obs->latt->a_2[0]), 
                                                (double) (obs->latt->r[i][0] * obs->latt->a_1[1] + obs->latt->r[i][1] * obs->latt->a_2[1]), 
                                                creal(obs->obs_i[n][i]), 
                                                cimag(obs->obs_i[n][i]));
    }

    for (i = 0; i < obs->latt->N; i++) {
      fprintf(out_k, "%lf (%lf, %lf) (%lf, %lf) \n", creal(obs->omega[n]), (double) (obs->latt->k[i][0] * obs->latt->b_1[0] + obs->latt->k[i][1] * obs->latt->b_2[0]), 
                                                (double) (obs->latt->k[i][0] * obs->latt->b_1[1] + obs->latt->k[i][1] * obs->latt->b_2[1]), 
                                                creal(obs->obs_k[n][i]), 
                                                cimag(obs->obs_k[n][i]));
    }
  }
}

void write_obs_spectral_info(FILE* info, Obs_spectral* obs)
{
  fprintf(info, "-- Analysis Mode --\n");
  fprintf(info, "spectral\n");
  fprintf(info, "-- Spectral --\n");
  fprintf(info, "Lx: %d\n", obs->latt->Lx);
  fprintf(info, "Ly: %d\n", obs->latt->Ly);
  fprintf(info, "n_om: %d\n", obs->n_om);
  fprintf(info, "om_f: %lf\n", obs->om_f);
}

void free_obs_latt(Obs_latt* obs)
{
  int i;
  
  for (i = 0; i < obs->latt->N; i++) {
    free(obs->obs_latt[i]);
  }
  free(obs->obs_latt);
  free(obs->obs_latt0);
  free(obs->obs_i);
  free(obs->obs_k);
}

void free_obs_spectral(Obs_spectral* obs)
{
  int n, i;
  free(obs->omega);

  for (n = 0; n < obs->n_om; n++) {
    free(obs->obs_k[n]);
    free(obs->obs_i[n]);

    for (i = 0; i < obs->latt->N; i++) {
      free(obs->obs_latt[n][i]);
    }
  }
}

void fourier_trans(Obs_latt* obs) 
{
  int i, j, n;
  double a[2], b[2];

  for (n = 0; n < obs->latt->N; n++) {
    for (i = 0; i < obs->latt->N; i++) {
      for (j = 0; j < obs->latt->N; j++) {
        a[0] = obs->latt->k[n][0] * obs->latt->b_1[0] + obs->latt->k[n][1] * obs->latt->b_2[0];
        a[1] = obs->latt->k[n][0] * obs->latt->b_1[1] + obs->latt->k[n][1] * obs->latt->b_2[1];
        b[0] = obs->latt->r_ij[i][j][0] * obs->latt->a_1[0] + obs->latt->r_ij[i][j][1] * obs->latt->a_2[0];
        b[1] = obs->latt->r_ij[i][j][0] * obs->latt->a_1[1] + obs->latt->r_ij[i][j][1] * obs->latt->a_2[1];

        obs->obs_k[n] += cexp(- I * (a[0] * b[0] + a[1] * b[1])) * (obs->obs_latt[i][j] - obs->obs_latt0[i] * obs->obs_latt0[j]) / (obs->latt->N);
      }
    }
  }
}

void inv_fourier_trans(Obs_latt* obs)
{
  int i, n;
  double a[2], b[2];

  for (i = 0; i < obs->latt->N; i++) {
    for (n = 0; n < obs->latt->N; n++) {
      a[0] = obs->latt->k[n][0] * obs->latt->b_1[0] + obs->latt->k[n][1] * obs->latt->b_2[0];
      a[1] = obs->latt->k[n][0] * obs->latt->b_1[1] + obs->latt->k[n][1] * obs->latt->b_2[1];
      b[0] = obs->latt->r[i][0] * obs->latt->a_1[0] + obs->latt->r[i][1] * obs->latt->a_2[0];
      b[1] = obs->latt->r[i][0] * obs->latt->a_1[1] + obs->latt->r[i][1] * obs->latt->a_2[1];

      obs->obs_i[i] += cexp(I * (a[0] * b[0] + a[1] * b[1])) * obs->obs_k[n] / (obs->latt->N);
    }
  }
}

void spec_fourier_trans(Obs_spectral* obs) 
{
  int i, j, n, om;
  double a[2], b[2];

  for (om = 0; om < obs->n_om; om++) {
    for (n = 0; n < obs->latt->N; n++) {
      for (i = 0; i < obs->latt->N; i++) {
        for (j = 0; j < obs->latt->N; j++) {
          a[0] = obs->latt->k[n][0] * obs->latt->b_1[0] + obs->latt->k[n][1] * obs->latt->b_2[0];
          a[1] = obs->latt->k[n][0] * obs->latt->b_1[1] + obs->latt->k[n][1] * obs->latt->b_2[1];
          b[0] = obs->latt->r_ij[i][j][0] * obs->latt->a_1[0] + obs->latt->r_ij[i][j][1] * obs->latt->a_2[0];
          b[1] = obs->latt->r_ij[i][j][0] * obs->latt->a_1[1] + obs->latt->r_ij[i][j][1] * obs->latt->a_2[1];

          obs->obs_k[om][n] += cexp(- I * (a[0] * b[0] + a[1] * b[1])) * obs->obs_latt[om][i][j] / (obs->latt->N);
        }
      }
    }
  }
}

void spec_inv_fourier_trans(Obs_spectral* obs)
{
  int i, n, om;
  double a[2], b[2];

  for (om = 0; om < obs->n_om; om++) {
    for (i = 0; i < obs->latt->N; i++) {
      for (n = 0; n < obs->latt->N; n++) {
        a[0] = obs->latt->k[n][0] * obs->latt->b_1[0] + obs->latt->k[n][1] * obs->latt->b_2[0];
        a[1] = obs->latt->k[n][0] * obs->latt->b_1[1] + obs->latt->k[n][1] * obs->latt->b_2[1];
        b[0] = obs->latt->r[i][0] * obs->latt->a_1[0] + obs->latt->r[i][1] * obs->latt->a_2[0];
        b[1] = obs->latt->r[i][0] * obs->latt->a_1[1] + obs->latt->r[i][1] * obs->latt->a_2[1];

        obs->obs_i[om][i] += cexp(I * (a[0] * b[0] + a[1] * b[1])) * obs->obs_k[om][n] / (obs->latt->N);
      }
    }
  }
}


