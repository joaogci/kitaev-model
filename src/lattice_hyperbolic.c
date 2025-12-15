#include "lattice_hyperbolic.h"

void make_lattice_hyperbolic(int Ns, int Nb, int*** adj_mat, LatticeHyperbolic* latt)
{
  int i, j, b;

  latt->Ns = Ns;
  latt->N = Ns/2;
  latt->Nb = Nb;
  
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

  for (i = 0; i < latt->Ns; i++) {
    free((*adj_mat)[i]);
  }
  free((*adj_mat));
}

void free_lattice_hyperbolic(LatticeHyperbolic* latt)
{
  int i;
  for (i = 0; i < latt->Nb; i++) {
    free(latt->bond_list[i]);
  }
  free(latt->bond_list);
}


