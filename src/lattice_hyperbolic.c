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
  latt->inv_bond_list = (int**) malloc(latt->Ns * sizeof(int*));
  for (i = 0; i < latt->Ns; i++) {
    latt->inv_bond_list[i] = (int*) malloc(latt->Ns * sizeof(int));
    for (j = 0; j < latt->Ns; j++) {
      latt->inv_bond_list[i][j] = -1;
    }
  }

  b = 0;
  for (i = 0; i < latt->Ns; i++) {
    for (j = 0; j <= i; j++) {
      if ((*adj_mat)[i][j] != 0) {
        latt->bond_list[b][0] = i;
        latt->bond_list[b][1] = j;
        latt->inv_bond_list[i][j] = b;
        b++;
      }
    }
  }

  /* plaquette list */
  latt->plaquette_list = (int**) malloc(latt->Np * sizeof(int*));
  latt->plaquette_bond_list = (int**) malloc(latt->Np * sizeof(int*));
  latt->bond_sign_plaquette = (int**) malloc(latt->Np * sizeof(int*));
  for (i = 0; i < latt->Np; i++) {
    latt->plaquette_list[i] = (int*) malloc(latt->Nvp * sizeof(int));
    latt->plaquette_bond_list[i] = (int*) malloc(latt->Nvp * sizeof(int));
    latt->bond_sign_plaquette[i] = (int*) malloc(latt->Nvp * sizeof(int));
  }

  for (i = 0; i < latt->Np; i++) {
    for (j = 0; j < latt->Nvp; j++) {
      latt->plaquette_list[i][j] = (*plaquette_list)[i][j];
    }
  }

  for (i = 0; i < latt->Np; i++) {
    for (j = 0; j < latt->Nvp; j++) {
      b = latt->inv_bond_list[latt->plaquette_list[i][j]][latt->plaquette_list[i][(j+1)%latt->Nvp]];
      if (b > -1) { 
        latt->bond_sign_plaquette[i][j] = 1;
        latt->plaquette_bond_list[i][j] = b;
      }
      b = latt->inv_bond_list[latt->plaquette_list[i][(j+1)%latt->Nvp]][latt->plaquette_list[i][j]];
      if (b > -1) {
        latt->bond_sign_plaquette[i][j] = -1;
        latt->plaquette_bond_list[i][j] = b;
      }
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
    free(latt->plaquette_list[i]);
    free(latt->plaquette_bond_list[i]);
    free(latt->bond_sign_plaquette[i]);
  }
  for (i = 0; i < latt->Ns; i++) {
    free(latt->inv_bond_list[i]);
  }
  free(latt->bond_list);
  free(latt->plaquette_list);
  free(latt->inv_bond_list);
  free(latt->plaquette_bond_list);
  free(latt->bond_sign_plaquette);
}


