#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "io/io.h"
#include "flux.h"
#include "observables.h"
#include "lattice_hyperbolic.h"

#define BUFFER 32

/*
 * Struct for scalar observables
 */
typedef struct Ana_scalar
{
  double _Complex obs_mean;
  double _Complex obs_std;
} Ana_scalar;

/*
 * Struct for LatticeHyperbolic observables
 */
typedef struct Ana_latt
{
  double _Complex* obs_mean;
  double _Complex* obs_std;
} Ana_latt;

/*
 * Struct for spectral observables
 */
typedef struct Ana_spectral
{
  double _Complex** obs_mean;
  double _Complex** obs_std;
} Ana_spectral;


/*
 * Analyses scalar observable
 */
void analyse_scal(FILE* fp, Ana_scalar* obs, int n_bins, Ana_info *ana_info);

/*
 * Analyses LatticeHyperbolic observable
 */
void analyse_latt(FILE* fp, Ana_latt* obs, LatticeHyperbolic *latt, int n_bins, Ana_info *ana_info);

/*
 * Analyses spectral observable
 */
void analyse_spectral(FILE* fp, Ana_spectral* obs, LatticeHyperbolic *latt, int n_om, double om_f, int n_bins, Ana_info *ana_info);

/*
 * Writes the results for scalar observable
 */
void write_scal(FILE* fp, Ana_scalar* obs);

/*
 * Writes the results for LatticeHyperbolic observable
 */
void write_latt_r(FILE* fp, Ana_latt* obs, LatticeHyperbolic* latt);
void write_latt_k(FILE* fp, Ana_latt* obs, LatticeHyperbolic* latt);

/*
 * Wrties the results for transport observable
 */
void write_spectral_r(FILE* fp, Ana_spectral* obs, LatticeHyperbolic *latt, int n_om, double om_f);
void write_spectral_k(FILE* fp, Ana_spectral* obs, LatticeHyperbolic *latt, int n_om, double om_f);

#endif // ANALYSIS_H
