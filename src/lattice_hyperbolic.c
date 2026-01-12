#include "lattice_hyperbolic.h"

void make_lattice_hyperbolic(int Ns, int Nb, int Np, int Nvp, int*** adj_mat, int*** plaquette_list, LatticeHyperbolic* latt)
{
  int i, j, b;

  latt->Ns = Ns;
  latt->N = Ns/2;
  latt->Nb = Nb;
  latt->Np = Np;
  latt->Nvp = Nvp;
  
  /* bond list */
  latt->bond_list = (int**) malloc(latt->Nb * sizeof(int*));
  for (i = 0; i < latt->Nb; i++) {
    latt->bond_list[i] = (int*) malloc(2 * sizeof(int));
  }

  b = 0;
  for (i = 0; i < latt->Ns; i++) {
    for (j = 0; j <= i; j++) {
      if ((*adj_mat)[i][j] != 0) {
        latt->bond_list[b][0] = i;
        latt->bond_list[b][1] = j;
        b++;
      }
    }
  }

  /* plaquette list */
  latt->plaqutte_list = (int**) malloc(latt->Np * sizeof(int*));
  for (i = 0; i < latt->Np; i++) {
    latt->plaqutte_list[i] = (int*) malloc(latt->Nvp * sizeof(int));
  }

  for (i = 0; i < latt->Np; i++) {
    for (j = 0; j < latt->Nvp; j++) {
      latt->plaqutte_list[i][j] = (*plaquette_list)[i][j];
    }
  }

  for (i = 0; i < latt->Ns; i++) {
    free((*adj_mat)[i]);
  }
  for (i = 0; i < latt->Np; i++) {
    free((*plaquette_list)[i]);
  }
  free((*adj_mat));
  free((*plaquette_list));
}

void free_lattice_hyperbolic(LatticeHyperbolic* latt)
{
  int i;
  for (i = 0; i < latt->Nb; i++) {
    free(latt->bond_list[i]);
  }
  for (i = 0; i < latt->Np; i++) {
    free(latt->plaqutte_list[i]);
  }
  free(latt->bond_list);
  free(latt->plaqutte_list);
}


