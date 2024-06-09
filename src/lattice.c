#include "lattice.h"

void make_lattice(int Lx, int Ly, Lattice *latt) 
{
  int i, j, n, n_, b, p;

  latt->Lx = Lx;
  latt->Ly = Ly;
  latt->N = Lx * Ly;
  latt->Ns = 2 * latt->N;
  latt->z = 3;
  latt->Nb = latt->z * latt->N;

  latt->a_1[0] = 0.5; 
  latt->a_1[1] = 0.5 * sqrt(3.0);
  latt->a_2[0] = - 0.5;
  latt->a_2[1] = 0.5 * sqrt(3.0);

  latt->b_1[0] = 2.0 * M_PI * latt->a_2[1] / (latt->Lx * (latt->a_1[0] * latt->a_2[1] - latt->a_1[1] * latt->a_2[0])); 
  latt->b_1[1] = - 2.0 * M_PI * latt->a_2[0] / (latt->Lx * (latt->a_1[0] * latt->a_2[1] - latt->a_1[1] * latt->a_2[0])); 
  latt->b_2[0] = - 2.0 * M_PI * latt->a_1[1] / (latt->Ly * (latt->a_1[0] * latt->a_2[1] - latt->a_1[1] * latt->a_2[0])); 
  latt->b_2[1] = 2.0 * M_PI * latt->a_1[0] / (latt->Ly * (latt->a_1[0] * latt->a_2[1] - latt->a_1[1] * latt->a_2[0])); 

  latt->p_a[0] = 0.0;
  latt->p_a[1] = - 0.5;
  latt->p_b[0] = 0.0;
  latt->p_b[1] = 0.5;

  /* Real-space lattice */
  latt->r = (int **) malloc(latt->N * sizeof(int *));
  for (i = 0; i < latt->N; i++) latt->r[i] = (int *) malloc(2 * sizeof(int));

  latt->inv_r = (int **) malloc(latt->Lx * sizeof(int *));
  for (i = 0; i < latt->Lx; i++) latt->inv_r[i] = (int *) malloc(latt->Ly * sizeof(int));

  latt->r_ij = (int ***) malloc(latt->N * sizeof(int **));
  for (i = 0; i < latt->N; i++) { 
    latt->r_ij[i] = (int **) malloc(latt->N * sizeof(int *));
    for (j = 0; j < latt->N; j++) {
      latt->r_ij[i][j] = (int *) malloc(2 * sizeof(int));
    }
  }

  n = 0;
  for (j = 0; j < Ly; j++) {
    for (i = 0; i < Lx; i++) {
      latt->r[n][0] = i;
      latt->r[n][1] = j;
      latt->inv_r[i][j] = n;
      n++;
    }
  }

  for (i = 0; i < latt->N; i++) {
    for (j = 0; j < latt->N; j++) {
      latt->r_ij[i][j][0] = latt->r[i][0] - latt->r[j][0];
      latt->r_ij[i][j][1] = latt->r[i][1] - latt->r[j][1];
    }
  }

  /* Reciprocal lattice */
  latt->k = (int **) malloc(latt->N * sizeof(int *));
  for (i = 0; i < latt->N; i++) latt->k[i] = (int *) malloc(2 * sizeof(int));

  latt->inv_k = (int **) malloc(latt->Lx * sizeof(int *));
  for (i = 0; i < latt->Lx; i++) latt->inv_k[i] = (int *) malloc(latt->Ly * sizeof(int));

  latt->k_ij = (int ***) malloc(latt->N * sizeof(int **));
  for (i = 0; i < latt->N; i++) { 
    latt->k_ij[i] = (int **) malloc(latt->N * sizeof(int *));
    for (j = 0; j < latt->N; j++) {
      latt->k_ij[i][j] = (int *) malloc(2 * sizeof(int));
    }
  }

  n = 0;
  for (j = 0; j < Ly; j++) {
    for (i = 0; i < Lx; i++) {
      latt->k[n][0] = i;
      latt->k[n][1] = j;
      latt->inv_k[i][j] = n;
      n++;
    }
  }

  for (i = 0; i < latt->N; i++) {
    for (j = 0; j < latt->N; j++) {
      latt->k_ij[i][j][0] = latt->k[i][0] - latt->k[j][0];
      latt->k_ij[i][j][1] = latt->k[i][1] - latt->k[j][1];
    }
  }

  /* Sites for lattice and bond list */
  latt->r_sites = (int *) malloc(latt->Ns * sizeof(int));
  latt->sublattice = (int *) malloc(latt->Ns * sizeof(int));

  n = 0;
  for (i = 0; i < latt->Lx; i++) {
    for (j = 0; j < latt->Ly; j++) {
      n_ = latt->inv_r[i][j];

      latt->r_sites[n] = n_;
      latt->sublattice[n] = 0;
      n++;

      latt->r_sites[n] = n_;
      latt->sublattice[n] = 1;
      n++;
    }
  }

  latt->bond_list = (int **) malloc(latt->Nb * sizeof(int *));
  for (i = 0; i < latt->Nb; i++) {
    latt->bond_list[i] = (int *) malloc(2 * sizeof(int));
  }

  latt->bond_type = (int *) malloc(latt->Nb * sizeof(int));
  latt->z_bonds = (int *) malloc(latt->N * sizeof(int));
  latt->y_bonds = (int *) malloc(latt->N * sizeof(int));
  latt->x_bonds = (int *) malloc(latt->N * sizeof(int));

  b = 0;
  p = 0;
  for (i = 0; i < latt->Lx; i++) {
    for (j = 0; j < latt->Ly; j++) {
      n = 2 * latt->inv_r[i][j];

      n_ = 2 * latt->inv_r[i][j] + 1;
      latt->bond_list[b][0] = n;
      latt->bond_list[b][1] = n_;
      latt->bond_type[b] = 2;
      latt->z_bonds[p] = b;
      b++;

      n_ = 2 * latt->inv_r[i][(j + 1) % latt->Ly] + 1;
      latt->bond_list[b][0] = n;
      latt->bond_list[b][1] = n_;
      latt->bond_type[b] = 1;
      latt->y_bonds[p] = b;
      b++;

      n_ = 2 * latt->inv_r[(i + 1) % latt->Lx][j] + 1;
      latt->bond_list[b][0] = n;
      latt->bond_list[b][1] = n_;
      latt->bond_type[b] = 0;
      latt->x_bonds[p] = b;
      b++;
      p++;
    }
  }

  /* Plaquettes */  
  latt->Np = latt->N;
  latt->size_p = 6;
  latt->plaquette_bonds = (int **) malloc(latt->Np * sizeof(int *));
  latt->plaquette_sites = (int **) malloc(latt->Np * sizeof(int *));
  for (p = 0; p < latt->Np; p++) {
    latt->plaquette_bonds[p] = (int *) malloc(latt->size_p * sizeof(int));
    latt->plaquette_sites[p] = (int *) malloc(latt->size_p * sizeof(int));
  }

  for (p = 0; p < latt->Np; p++) {
    latt->plaquette_sites[p][0] = latt->bond_list[latt->z_bonds[p]][0];
    latt->plaquette_sites[p][1] = latt->bond_list[latt->z_bonds[p]][1];
    latt->plaquette_bonds[p][0] = latt->z_bonds[p];

    for (i = 0; i < latt->N; i++) {
      if (latt->bond_list[latt->x_bonds[i]][1] == latt->plaquette_sites[p][1]) {
        latt->plaquette_sites[p][2] = latt->bond_list[latt->x_bonds[i]][0];
        latt->plaquette_bonds[p][1] = latt->x_bonds[i];
        break;
      }
    }

    for (i = 0; i < latt->N; i++) {
      if (latt->bond_list[latt->y_bonds[i]][0] == latt->plaquette_sites[p][2]) {
        latt->plaquette_sites[p][3] = latt->bond_list[latt->y_bonds[i]][1];
        latt->plaquette_bonds[p][2] = latt->y_bonds[i];
        break;
      }
    }

    for (i = 0; i < latt->N; i++) {
      if (latt->bond_list[latt->z_bonds[i]][1] == latt->plaquette_sites[p][3]) {
        latt->plaquette_sites[p][4] = latt->bond_list[latt->z_bonds[i]][0];
        latt->plaquette_bonds[p][3] = latt->z_bonds[i];
        break;
      }
    }

    for (i = 0; i < latt->N; i++) {
      if (latt->bond_list[latt->x_bonds[i]][0] == latt->plaquette_sites[p][4]) {
        latt->plaquette_sites[p][5] = latt->bond_list[latt->x_bonds[i]][1];
        latt->plaquette_bonds[p][4] = latt->x_bonds[i];
        break;
      }
    }

    for (i = 0; i < latt->N; i++) {
      if (latt->bond_list[latt->y_bonds[i]][1] == latt->plaquette_sites[p][5]) {
        latt->plaquette_bonds[p][5] = latt->y_bonds[i];
        break;
      }
    }
  }
}

void free_lattice(Lattice* latt) 
{
  int i, j;

  for (i = 0; i < latt->N; i++) {
    for (j = 0; j < latt->N; j++) {
      free(latt->r_ij[i][j]);
      free(latt->k_ij[i][j]);
    }
    free(latt->r_ij[i]);
    free(latt->k_ij[i]);
    free(latt->r[i]);
    free(latt->k[i]);
  }   
  free(latt->r);
  free(latt->r_ij);
  free(latt->k);

  for (i = 0; i < latt->Lx; i++) {
    free(latt->inv_r[i]);
    free(latt->inv_k[i]);
  }
  free(latt->inv_r);
  free(latt->inv_k);

  for (i = 0; i < latt->Nb; i++) {
    free(latt->bond_list[i]);
  }
  free(latt->bond_list);
  free(latt->r_sites);
  free(latt->sublattice);
  free(latt->bond_type);
  free(latt->x_bonds);
  free(latt->y_bonds);
  free(latt->z_bonds);

  for (i = 0; i < latt->Np; i++) {
    free(latt->plaquette_bonds[i]);
    free(latt->plaquette_sites[i]);
  }
  free(latt->plaquette_bonds);
  free(latt->plaquette_sites);
}
