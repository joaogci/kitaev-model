#include "io.h"

void write_observables(Obs_scalar* obs_scal, int n_scal, Obs_latt* obs_eq, int n_eq, Obs_spectral *obs_spec, int n_spec)
{
  int i;
  FILE* out,* out_i,* out_k,* info;
  char filename[BUFFER_SIZE];

  for (i = 0; i < n_scal; i++) {
    strcpy(filename, obs_scal[i].filename);
    strcat(filename, "_scal");

    out = fopen(filename, "a");
    write_obs_scalar(out, &(obs_scal[i]));
    fclose(out);

    strcat(filename, "_info");
    info = fopen(filename, "w");
    write_obs_scalar_info(info, &(obs_scal[i]));
    fclose(info);
  }

  for (i = 0; i < n_eq; i++) {
    strcpy(filename, obs_eq[i].filename);
    strcat(filename, "_eqR");
    out_i = fopen(filename, "a");

    strcpy(filename, obs_eq[i].filename);
    strcat(filename, "_eqK");
    out_k = fopen(filename, "a");

    write_obs_latt(out_i, out_k, &(obs_eq[i]));

    fclose(out_i);
    fclose(out_k);

    strcpy(filename, obs_eq[i].filename);
    strcat(filename, "_eqR_info");
    info = fopen(filename, "w");
    write_obs_eq_info(info, &(obs_eq[i]));
    fclose(info);

    strcpy(filename, obs_eq[i].filename);
    strcat(filename, "_eqK_info");
    info = fopen(filename, "w");
    write_obs_eq_info(info, &(obs_eq[i]));
    fclose(info);
  }

  for (i = 0; i < n_spec; i++) {
    strcpy(filename, obs_spec[i].filename);
    strcat(filename, "_specR");
    out_i = fopen(filename, "a");

    strcpy(filename, obs_spec[i].filename);
    strcat(filename, "_specK");
    out_k = fopen(filename, "a");

    write_obs_spectral(out_i, out_k, &(obs_spec[i]));

    fclose(out_i);
    fclose(out_k);

    strcpy(filename, obs_spec[i].filename);
    strcat(filename, "_specR_info");
    info = fopen(filename, "w");
    write_obs_spectral_info(info, &(obs_spec[i]));
    fclose(info);

    strcpy(filename, obs_spec[i].filename);
    strcat(filename, "_specK_info");
    info = fopen(filename, "w");
    write_obs_spectral_info(info, &(obs_spec[i]));
    fclose(info);
  }
}

// void write_transport_obeservables(Obs_transport* obs, int n_transp)
// {
//   int i;
//   FILE* out,* info;
//   char filename[BUFFER_SIZE], tmp[BUFFER_SIZE];

//   for (i = 0; i < n_transp; i++) {
//     strcpy(filename, obs[i].filename);
    
//     strcat(filename, "_");
//     sprintf(tmp, "%d", obs[i].x);
//     strcat(filename, tmp);
    
//     strcat(filename, "_");
//     sprintf(tmp, "%d", obs[i].y);
//     strcat(filename, tmp);
//     strcat(filename, "_transp");

//     out = fopen(filename, "a");
//     write_obs_transport(out, &(obs[i]));
//     fclose(out);

//     strcat(filename, "_info");
//     info = fopen(filename, "w");
//     write_obs_transport_info(info, &(obs[i]));
//     fclose(info);
//   }
// }

void write_sim_info(Sim_info sim)
{
  FILE* out;

  out = fopen("info", "a");
  fprintf(out, "  -------------------------------  \n");
  fprintf(out, "Number of bins, sweeps:  %d, %ld\n", sim.n_bins, sim.mc_sweeps);
  // fprintf(out, "Number of threads:       %d \n", sim.n_threads);
  fprintf(out, "Wall Time:               %lfs \n", sim.wall_time);
  fclose(out);
}

void write_configuration(Flux* flux_conf)
{
  FILE* out;
  int n;
  char filename[BUFFER_SIZE];

  sprintf(filename, "confout");
  out = fopen(filename, "w");
  
  for (n = 0; n < flux_conf->latt->Nb; n++) {
    fprintf(out, "%d \n", flux_conf->flux[n]);
  }

  fclose(out);
}

void read_configuration(int *flux)
{
  FILE* out;
  int n;
  char filename[BUFFER_SIZE];

  sprintf(filename, "confin");
  out = fopen(filename, "r");

  if (out != NULL) {
    for (n = 0; n < num_lines(out); n++) {
      fscanf(out, "%d \n", &(flux[n]));
    }
  } else {
    printf("Error opening the %s file. \n", filename);
    exit(1);
  } 

  fclose(out);
}

int num_lines(FILE* fp)
{
  int ch, lines;

  lines = 0;
  while(!feof(fp))
  {
    ch = fgetc(fp);
    if(ch == '\n')
    {
      lines++;
    }
  }

  return lines;
}

void read_hyperbolic_lattice(int N, int *Nb, int*** adj_mat, char* sse_path)
{
  int i, j;
  char filename[BUFFER_SIZE];
  FILE* file;

  sprintf(filename, "%s/src/lattices/C%d", sse_path, N);
  file = fopen(filename, "r");

  if (file != NULL) {
    fscanf(file, "%d %d\n", &N, Nb);
    (*adj_mat) = (int**) malloc((N) * sizeof(int*));
    
    for (i = 0; i < N; i++) {
      (*adj_mat)[i] = (int*) malloc(N * sizeof(int));

      for (j = 0; j < N; j++) {
        fscanf(file, "%d.", &((*adj_mat)[i][j]));
      }
    }
  } else {
    printf("Error opening the adjacency matrix file. \n");
    exit(1);
  }

  fclose(file); 
}

