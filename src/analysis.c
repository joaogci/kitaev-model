#include "analysis.h"

void analyse_scal(FILE* fp, Ana_scalar* obs, int n_bins, Ana_info *ana_info)
{
  int i, j;
  double real, imag;
  double _Complex* measurements;
  double _Complex* measurements_rebin;

  n_bins = n_bins - ana_info->n_skip;
  measurements = (double _Complex*) malloc(n_bins * sizeof(double _Complex));
  for (i = 0; i < n_bins + ana_info->n_skip; i++) {
    fscanf(fp, "(%lf, %lf) ", &real, &imag);
    if (i >= ana_info->n_skip) {
      measurements[i - ana_info->n_skip] = real + I * imag;
    }
  }

  measurements_rebin = (double _Complex*) malloc(n_bins / ana_info->n_rebin * sizeof(double _Complex));
  for (i = 0; i < n_bins / ana_info->n_rebin; i++) {
    measurements_rebin[i] = 0.0;

    for (j = 0; j < ana_info->n_rebin; j++) {
      measurements_rebin[i] += measurements[i * ana_info->n_rebin + j];
    }
    measurements_rebin[i] = measurements_rebin[i] / ana_info->n_rebin;
  }

  for (i = 0; i < n_bins/ana_info->n_rebin; i++) {
    obs->obs_mean += measurements_rebin[i];
  }  
  obs->obs_mean = obs->obs_mean / (n_bins/ana_info->n_rebin);

  for (i = 0; i < n_bins/ana_info->n_rebin; i++) {
    obs->obs_std += cpow(measurements_rebin[i] - obs->obs_mean, 2.0);
  }
  obs->obs_std = csqrt(obs->obs_std / (n_bins/ana_info->n_rebin));

  free(measurements);
  free(measurements_rebin);
}

void analyse_latt(FILE* fp, Ana_latt* obs, LatticeHyperbolic *latt, int n_bins, Ana_info *ana_info)
{
  int i, j, n;
  double real, imag, x, y;
  double _Complex** measurements;
  double _Complex** measurements_rebin;

  n_bins = n_bins - ana_info->n_skip;
  measurements = (double _Complex**) malloc(n_bins * sizeof(double _Complex*));
  for (i = 0; i < n_bins + ana_info->n_skip; i++) {
    if (i >= ana_info->n_skip) {
      measurements[i - ana_info->n_skip] = (double _Complex*) malloc(latt->N * sizeof(double _Complex));
    }

    for (n = 0; n < latt->N; n++) {
      fscanf(fp, "(%lf, %lf) (%lf, %lf) \n", &x, &y, &real, &imag);
      if (i >= ana_info->n_skip) {
        measurements[i - ana_info->n_skip][n] = real + I * imag;
      }
    }
  }

  measurements_rebin = (double _Complex**) malloc(n_bins / ana_info->n_rebin * sizeof(double _Complex*));
  for (i = 0; i < n_bins / ana_info->n_rebin; i++) {
    measurements_rebin[i] = (double _Complex*) malloc(latt->N * sizeof(double _Complex));
    
    for (n = 0; n < latt->N; n++) {
      for (j = 0; j < ana_info->n_rebin; j++) {
        measurements_rebin[i][n] += measurements[i * ana_info->n_rebin + j][n];
      }
      measurements_rebin[i][n] = measurements_rebin[i][n] / ana_info->n_rebin;
    }
  }

  for (n = 0; n < latt->N; n++) {
    for (i = 0; i < n_bins / ana_info->n_rebin; i++) {
      obs->obs_mean[n] += measurements_rebin[i][n];
    }
    obs->obs_mean[n] = obs->obs_mean[n] / (n_bins / ana_info->n_rebin);
  }  

  for (n = 0; n < latt->N; n++) {
    for (i = 0; i < n_bins / ana_info->n_rebin; i++) {
      obs->obs_std[n] += cpow(measurements_rebin[i][n] - obs->obs_mean[n], 2.0);
    }
    obs->obs_std[n] = csqrt(obs->obs_std[n] / (n_bins / ana_info->n_rebin));
  }

  for (i = 0; i < n_bins; i++) {
    free(measurements[i]);
  }
  free(measurements);
  for (i = 0; i < n_bins/ana_info->n_rebin; i++) {
    free(measurements_rebin[i]);
  }
  free(measurements_rebin);
}

void analyse_spectral(FILE* fp, Ana_spectral* obs, LatticeHyperbolic *latt, int n_om, double om_f, int n_bins, Ana_info *ana_info)
{
  int i, j, n, b;
  double real, imag, x, y, om;
  double _Complex*** measurements;
  double _Complex*** measurements_rebin;

  n_bins = n_bins - ana_info->n_skip;
  measurements = (double _Complex***) malloc(n_bins * sizeof(double _Complex**));
  for (b = 0; b < n_bins + ana_info->n_skip; b++) {
    if (b >= ana_info->n_skip) {
      measurements[b - ana_info->n_skip] = (double _Complex**) malloc(n_om * sizeof(double _Complex*));

      for (n = 0; n < n_om; n++) {
        measurements[b - ana_info->n_skip][n] = (double _Complex*) malloc(latt->N * sizeof(double _Complex));
        for (i = 0; i < latt->N; i++) {
          fscanf(fp, "%lf (%lf, %lf) (%lf, %lf) \n", &om, &x, &y, &real, &imag);
          measurements[b - ana_info->n_skip][n][i] = real + I * imag;
        }
      }
    }
  }

  measurements_rebin = (double _Complex***) malloc(n_bins / ana_info->n_rebin * sizeof(double _Complex**));
  for (b = 0; b < n_bins / ana_info->n_rebin; b++) {
    measurements_rebin[b] = (double _Complex**) malloc(n_om * sizeof(double _Complex*));
    
    for (n = 0; n < n_om; n++) {
      measurements_rebin[b][n] = (double _Complex*) malloc(latt->N * sizeof(double _Complex));
      for (i = 0; i < latt->N; i++) {
        for (j = 0; j < ana_info->n_rebin; j++) {
          measurements_rebin[b][n][i] += measurements[b * ana_info->n_rebin + j][n][i];
        }
        measurements_rebin[b][n][i] = measurements_rebin[b][n][i] / ana_info->n_rebin;
      }
    }
  }

  for (n = 0; n < n_om; n++) {
    for (i = 0; i < latt->N; i++) {
      for (b = 0; b < n_bins / ana_info->n_rebin; b++) {
        obs->obs_mean[n][i] += measurements_rebin[b][n][i];
      }
      obs->obs_mean[n][i] = obs->obs_mean[n][i] / (n_bins / ana_info->n_rebin);
    }
  }  

  for (n = 0; n < n_om; n++) {
    for (i = 0; i < latt->N; i++) {
      for (b = 0; b < n_bins / ana_info->n_rebin; b++) {
        obs->obs_std[n][i] += cpow(measurements_rebin[b][n][i] - obs->obs_mean[n][i], 2.0);
      }
      obs->obs_std[n][i] = csqrt(obs->obs_std[n][i] / (n_bins / ana_info->n_rebin));
    }
  }

  for (i = 0; i < n_bins; i++) {
    for (n = 0; n < n_om; n++) {
      free(measurements[i][n]);
    }
    free(measurements[i]);
  }
  free(measurements);
  for (i = 0; i < n_bins/ana_info->n_rebin; i++) {
    for (n = 0; n < n_om; n++) {
      free(measurements_rebin[i][n]);
    }
    free(measurements_rebin[i]);
  }
  free(measurements_rebin);
}

void write_scal(FILE* fp, Ana_scalar* obs)
{
  fprintf(fp, "%lf %lf %lf %lf \n", creal(obs->obs_mean), creal(obs->obs_std), cimag(obs->obs_mean), cimag(obs->obs_std));
}

void write_latt_r(FILE* fp, Ana_latt* obs, LatticeHyperbolic* latt)
{
  // int i;

  // for (i = 0; i < latt->N; i++) {
  //   fprintf(fp, "%lf %lf \n", (double) (latt->r[i][0] * latt->a_1[0] + latt->r[i][1] * latt->a_2[0]), (double) (latt->r[i][0] * latt->a_1[1] + latt->r[i][1] * latt->a_2[1]));
  //   fprintf(fp, "%lf %lf %lf %lf \n", creal(obs->obs_mean[i]), creal(obs->obs_std[i]), cimag(obs->obs_mean[i]), cimag(obs->obs_std[i]));
  // }
}

void write_latt_k(FILE* fp, Ana_latt* obs, LatticeHyperbolic* latt)
{
  // int i;

  // for (i = 0; i < latt->N; i++) {
  //   fprintf(fp, "%lf %lf \n", (double) (latt->k[i][0] * latt->b_1[0] + latt->k[i][1] * latt->b_2[0]), (double) (latt->k[i][0] * latt->b_1[1] + latt->k[i][1] * latt->b_2[1]));
  //   fprintf(fp, "%lf %lf %lf %lf \n", creal(obs->obs_mean[i]), creal(obs->obs_std[i]), cimag(obs->obs_mean[i]), cimag(obs->obs_std[i]));
  // }
}

void write_spectral_r(FILE* fp, Ana_spectral* obs, LatticeHyperbolic *latt, int n_om, double om_f)
{
  // int n, i;

  // for (n = 0; n < n_om; n++) {
  //   fprintf(fp, "%lf \n", n * om_f / (n_om - 1));
  //   for (i = 0; i < latt->N; i++) {
  //     fprintf(fp, "%lf %lf \n", (double) (latt->r[i][0] * latt->a_1[0] + latt->r[i][1] * latt->a_2[0]), (double) (latt->r[i][0] * latt->a_1[1] + latt->r[i][1] * latt->a_2[1]));
  //     fprintf(fp, "%lf %lf %lf %lf \n", creal(obs->obs_mean[n][i]), creal(obs->obs_std[n][i]), cimag(obs->obs_mean[n][i]), cimag(obs->obs_std[n][i]));
  //   }
  // }
}

void write_spectral_k(FILE* fp, Ana_spectral* obs, LatticeHyperbolic *latt, int n_om, double om_f)
{
  // int i, n;

  // for (n = 0; n < n_om; n++) {
  //   fprintf(fp, "%lf \n", n * om_f / (n_om - 1));
  //   for (i = 0; i < latt->N; i++) {
  //     fprintf(fp, "%lf %lf \n", (double) (latt->k[i][0] * latt->b_1[0] + latt->k[i][1] * latt->b_2[0]), (double) (latt->k[i][0] * latt->b_1[1] + latt->k[i][1] * latt->b_2[1]));
  //     fprintf(fp, "%lf %lf %lf %lf \n", creal(obs->obs_mean[n][i]), creal(obs->obs_std[n][i]), cimag(obs->obs_mean[n][i]), cimag(obs->obs_std[n][i]));
  //   }
  // }
}
