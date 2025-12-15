#ifndef LATTICE_HYPERBOLIC_H
#define LATTICE_HYPERBOLIC_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

/*
 * Hyperbolic lattice struct
 */
typedef struct LatticeHyperbolic
{
  int Ns;
  int Nb;
  int N;

  int** bond_list;  
} LatticeHyperbolic;

/*
 * makes hyperbolic lattice
 */
void make_lattice_hyperbolic(int N, int Nb, int*** adj_mat, LatticeHyperbolic* latt);

/*
 * frees hyperbolic lattice
 */
void free_lattice_hyperbolic(LatticeHyperbolic* latt);

#endif // LATTICE_HYPERBOLIC_H
