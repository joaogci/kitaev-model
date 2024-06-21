#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "observables.h"
#include "lattice.h"
#include "io/io.h"
#include "analysis.h"

// Lattice
int Lx, Ly;
Lattice latt;

// Spectral
int n_om;
double om_f;

// Observables
Ana_scalar obs_scal;
Ana_latt obs_eq;
Ana_spectral obs_spec;

Ana_info ana_info;

int main(int argc, char** argv)
{
  int i, n, j, n_bins;
  size_t len;
  char filename[BUFFER];
  char* line;
  FILE* out;
  FILE* info,* res;

  if (argc <= 1) {
    printf("Please provide the files to be analysed. \n");
    printf("Usage: %s file_1 file_2 ... \n ", argv[0]);
    exit(1);
  }

  read_parameters_analysis(&(ana_info.n_skip), &(ana_info.n_rebin));

  // Analysis
  for (i = 1; i < argc; i++) {
    if (strstr(argv[i], "_info") == NULL && (strstr(argv[i], "_scal") != NULL || strstr(argv[i], "_eq") != NULL || strstr(argv[i], "_spec") != NULL)) {
    // strcmp(argv[i], "vertices_info") != 0 && strcmp(argv[i], "parameters") != 0 && strcmp(argv[i], "info") != 0) {
      if (access(argv[i], F_OK) != 0) {
        printf("%s file does not exist. \n", argv[i]);
        continue;
      }
      
      if (strstr(argv[i], "J") != NULL || strstr(argv[i], "_eqJR") != NULL || strstr(argv[i], "_eqKJ") != NULL || strstr(argv[i], "_specJ") != NULL) {
        continue;
      }

      strcpy(filename, argv[i]);
      strcat(filename, "_info");
      info = fopen(filename, "r");

      line = NULL;
      getline(&line, &len, info);
      getline(&line, &len, info);
      line[strlen(line) - 1] = '\0';

      // strcpy(filename, argv[i]);
      // filename[strlen(filename) - 5] = '\0';
      
      switch (line[0]) 
      {
        case 'i':
          obs_scal.obs_mean = 0.0;
          obs_scal.obs_std = 0.0;

          strcpy(filename, argv[i]);
          res = fopen(filename, "r");
          len = num_lines(res);
          n_bins = len;
          fclose(res);
          
          printf("Analysing %s \n", filename);
          printf("Number of bins: %d \n", n_bins - ana_info.n_skip);
          printf("Number of effective bins: %d \n", (n_bins - ana_info.n_skip) / ana_info.n_rebin);

          res = fopen(filename, "r");
          analyse_scal(res, &obs_scal, n_bins, &ana_info);
          fclose(res);

          strcat(filename, "J");
          out = fopen(filename, "w");
          write_scal(out, &obs_scal);
          fclose(out);
          
          break;
        case 'e':
          getline(&line, &len, info);
          fscanf(info, "Lx: %d\n", &Lx);
          fscanf(info, "Ly: %d\n", &Ly);

          make_lattice(Lx, Ly, &latt);
          
          obs_eq.obs_mean = (double _Complex*) malloc(latt.N * sizeof(double _Complex));
          obs_eq.obs_std = (double _Complex*) malloc(latt.N * sizeof(double _Complex));
          for (n = 0; n < latt.N; n++) {
            obs_eq.obs_mean[n] = 0.0;
            obs_eq.obs_std[n] = 0.0;
          }

          strcpy(filename, argv[i]);
          // strcat(filename, "R");
          res = fopen(filename, "r");
          len = num_lines(res);
          n_bins = len / latt.N;
          fclose(res);

          printf("Analysing %s \n", filename);
          printf("Number of bins: %d \n", n_bins - ana_info.n_skip);
          printf("Number of effective bins: %d \n", (n_bins - ana_info.n_skip) / ana_info.n_rebin);

          res = fopen(filename, "r");
          analyse_latt(res, &obs_eq, &latt, n_bins, &ana_info);
          fclose(res);

          filename[strlen(filename) - 1] = '\0';
          strcat(filename, "JR");
          out = fopen(filename, "w");
          write_latt_r(out, &obs_eq, &latt);
          fclose(out);

          for (n = 0; n < latt.N; n++) {
            obs_eq.obs_mean[n] = 0.0;
            obs_eq.obs_std[n] = 0.0;
          }

          filename[strlen(filename) - 2] = '\0';
          strcat(filename, "K");
          res = fopen(filename, "r");
          len = num_lines(res);
          n_bins = len / latt.N;
          fclose(res);

          printf("Analysing %s \n", filename);
          printf("Number of bins: %d \n", n_bins - ana_info.n_skip);
          printf("Number of effective bins: %d \n", (n_bins - ana_info.n_skip) / ana_info.n_rebin);

          res = fopen(filename, "r");
          analyse_latt(res, &obs_eq, &latt, n_bins, &ana_info);
          fclose(res);


          filename[strlen(filename) - 1] = '\0';
          strcat(filename, "JK");
          out = fopen(filename, "w");
          write_latt_k(out, &obs_eq, &latt);
          fclose(out);

          free(obs_eq.obs_mean);
          free(obs_eq.obs_std);

          free_lattice(&latt);
          break;
        case 's':
          getline(&line, &len, info);
          fscanf(info, "Lx: %d\n", &Lx);
          fscanf(info, "Ly: %d\n", &Ly);
          fscanf(info, "n_om: %d\n", &n_om);
          fscanf(info, "om_f: %lf\n", &om_f);

          make_lattice(Lx, Ly, &latt);

          obs_spec.obs_mean = (double _Complex**) malloc(n_om * sizeof(double _Complex*));
          obs_spec.obs_std = (double _Complex**) malloc(n_om * sizeof(double _Complex*));
          for (n = 0; n < n_om; n++) {
            obs_spec.obs_mean[n] = (double _Complex*) malloc(latt.N * sizeof(double _Complex));
            obs_spec.obs_std[n] = (double _Complex*) malloc(latt.N * sizeof(double _Complex));
            for (j = 0; j < latt.N; j++) {
              obs_spec.obs_mean[n][j] = 0.0;
              obs_spec.obs_std[n][j] = 0.0;
            }
          }

          strcpy(filename, argv[i]);
          // strcat(filename, "R");
          res = fopen(filename, "r");
          len = num_lines(res);
          n_bins = len / (latt.N * n_om);
          fclose(res);

          printf("Analysing %s \n", filename);
          printf("Number of bins: %d \n", n_bins - ana_info.n_skip);
          printf("Number of effective bins: %d \n", (n_bins - ana_info.n_skip) / ana_info.n_rebin);

          res = fopen(filename, "r");
          analyse_spectral(res, &obs_spec, &latt, n_om, om_f, n_bins, &ana_info);
          fclose(res);

          filename[strlen(filename) - 1] = '\0';
          strcat(filename, "JR");
          out = fopen(filename, "w");
          write_spectral_r(out, &obs_spec, &latt, n_om, om_f);
          fclose(out);

          for (n = 0; n < n_om; n++) {
            for (j = 0; j < latt.N; j++) {
              obs_spec.obs_mean[n][j] = 0.0;
              obs_spec.obs_std[n][j] = 0.0;
            }
          } 

          filename[strlen(filename) - 2] = '\0';
          strcat(filename, "K");
          res = fopen(filename, "r");
          len = num_lines(res);
          n_bins = len / (latt.N * n_om);
          fclose(res);

          printf("Analysing %s \n", filename);
          printf("Number of bins: %d \n", n_bins - ana_info.n_skip);
          printf("Number of effective bins: %d \n", (n_bins - ana_info.n_skip) / ana_info.n_rebin);

          res = fopen(filename, "r");
          analyse_spectral(res, &obs_spec, &latt, n_om, om_f, n_bins, &ana_info);
          fclose(res);


          filename[strlen(filename) - 1] = '\0';
          strcat(filename, "JK");
          out = fopen(filename, "w");
          write_spectral_k(out, &obs_spec, &latt, n_om, om_f);
          fclose(out);

          for (n = 0; n < n_om; n++) {
            free(obs_spec.obs_mean[n]);
            free(obs_spec.obs_std[n]);
          }
          free(obs_spec.obs_mean);
          free(obs_spec.obs_std);
          free_lattice(&latt);
          break;
      }
      
      fclose(info);
    }
  }

  return 0;
}
