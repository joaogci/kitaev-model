#ifndef LATTICE_H 
#define LATTICE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * Honeycomb lattice struct
 */
typedef struct Lattice
{
  int Lx;
  int Ly;
  int z;

  int N;
  int Ns;
  int Nb;
  
  double a_1[2];
  double a_2[2];
  double b_1[2];
  double b_2[2];
  double p_a[2];
  double p_b[2];

  int **r;
  int **inv_r;
  int ***r_ij;

  int **k;
  int **inv_k;
  int ***k_ij;

  int *r_sites;
  int *sublattice;
  int **bond_list;
  int *bond_type;
  int *z_bonds;
  int *y_bonds;
  int *x_bonds;

  int Np;
  int size_p;
  int **plaquette_bonds;
  int **plaquette_sites;
} Lattice;

/* 
 * Initializes Lattice 
 */
void make_lattice(int Lx, int Ly, Lattice *lattice);

/* 
 * Frees the lattice 
 */
void free_lattice(Lattice *lattice);

#endif // LATTICE_H
