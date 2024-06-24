#ifndef OBSERVABLES_H 
#define OBSERVABLES_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <stdint.h>

#include "lattice.h"

#define BUFFER 32

/*
 * Struct for scalar observables
 */
typedef struct Obs_scalar
{
  u_int64_t N;
  char filename[BUFFER];

  double _Complex obs_vec;
} Obs_scalar;

/*
 * Struct for equal time correlation functions
 * or static susceptibilities 
 */
typedef struct Obs_latt
{
  u_int64_t N;
  char filename[BUFFER];

  double _Complex** obs_latt;
  double _Complex* obs_latt0;

  double _Complex* obs_k;
  double _Complex* obs_i;

  Lattice* latt;
} Obs_latt;

/* 
 * Struct for spectral functions
 */
typedef struct Obs_spectral
{
  u_int64_t N;
  char filename[BUFFER];
  
  int n_om;
  double om_f, eta;
  double* omega;

  double _Complex*** obs_latt;
  double _Complex **obs_k, **obs_i;
  Lattice *latt;
} Obs_spectral;


void init_obs_scalar(char* filename, Obs_scalar* obs);
void init_obs_latt(char* filename, Lattice* latt, Obs_latt* obs);
void init_obs_spectral(char* filename, int n_om, double om_f, Lattice *latt, Obs_spectral* obs);

void reset_obs_scalar(Obs_scalar* obs);
void reset_obs_latt(Obs_latt* obs);
void reset_obs_spectral(Obs_spectral* obs);

void write_obs_scalar(FILE* out, Obs_scalar* obs);
void write_obs_latt(FILE* out_i, FILE* out_k, Obs_latt* obs);
void write_obs_spectral(FILE* out_i, FILE* out_k, Obs_spectral* obs);

void write_obs_scalar_info(FILE* info, Obs_scalar* obs);
void write_obs_eq_info(FILE* info, Obs_latt* obs);
void write_obs_spectral_info(FILE* out, Obs_spectral* obs);

void free_obs_latt(Obs_latt* obs);
void free_obs_spectral(Obs_spectral* obs);

void fourier_trans(Obs_latt* obs);
void inv_fourier_trans(Obs_latt* obs);

void spec_fourier_trans(Obs_spectral* obs);
void spec_inv_fourier_trans(Obs_spectral* obs);

#endif // OBSRVABLES_H
