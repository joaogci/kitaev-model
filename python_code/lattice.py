import numpy as np

class Lattice():
  # Real space unit cell vectors
  n1 = 0.5 * np.array([  1, np.sqrt(3)])
  n2 = 0.5 * np.array([- 1, np.sqrt(3)])
  # Reciprocal space unit cell vectors
  b1 = 2.0 * np.pi * np.array([n2[1], - n2[0]]) / (n1[0] * n2[1] - n1[1] * n2[0])
  b2 = 2.0 * np.pi * np.array([- n1[1], n1[0]]) / (n1[0] * n2[1] - n1[1] * n2[0])
  # Position of the lattice sites inside the unit cell
  __p1 = np.array([0.0, - 0.5])
  __p2 = np.array([0.0, 0.5])

  def __init__(self, Lx: int, Ly: int):
    self.Lx = Lx
    self.Ly = Ly
    self.N = self.Lx * self.Ly
    self.z = 3
    self.Ns = 2 * self.N
    self.Nb = self.z * self.Ns // 2
    
    # Real space and k space unit cell coordinates
    self.r_unit_cell = np.zeros((self.N, 2))
    self.inv_r_unit_cell = np.zeros((self.Lx, self.Ly), dtype=int)
    self.k_unit_cell = np.zeros((self.N, 2))
    self.inv_k_unit_cell = np.zeros((self.Lx, self.Ly), dtype=int)
   
    n = 0 
    for j in range(self.Ly):
      for i in range(self.Lx):
        self.r_unit_cell[n, :] = i * self.n1 + j * self.n2
        self.inv_r_unit_cell[i, j] = n
        self.k_unit_cell[n, :] = i * self.b1 / self.Lx + j * self.b2 / self.Ly
        self.inv_k_unit_cell[i, j] = n
        
        n += 1
   
    # Lattice site coordinates and bond matrix
    self.r_sites = np.zeros((self.Ns, 2))
    self.unit_cell = np.zeros((self.Ns), dtype=int)
    self.sublattice = np.zeros((self.Ns), dtype=int)
    self.bond_list = np.zeros((self.Nb, 2), dtype=int)
    self.bond_type = np.zeros((self.Nb), dtype=int)
    
    n = 0
    for i in range(self.Lx):
      for j in range(self.Ly):
        n_ = self.inv_r_unit_cell[i, j]
        
        self.r_sites[n, :] = self.r_unit_cell[n_, :] + self.__p1
        self.unit_cell[n] = n_
        self.sublattice[n] = 1
        n += 1
        
        self.r_sites[n, :] = self.r_unit_cell[n_, :] + self.__p2
        self.unit_cell[n] = n_
        self.sublattice[n] = 2
        n += 1

    nb = 0
    for i in range(self.Lx):
      for j in range(self.Ly):
        n = 2 * self.inv_r_unit_cell[i, j]
        
        # nn in z dir
        n1 = 2 * self.inv_r_unit_cell[i, j] + 1
        self.bond_list[nb, 0] = n
        self.bond_list[nb, 1] = n1
        self.bond_type[nb] = 3
        nb += 1
        
        # nn in y dir
        n1 = 2 * self.inv_r_unit_cell[i, (j + 1)%self.Ly] + 1
        self.bond_list[nb, 0] = n
        self.bond_list[nb, 1] = n1
        self.bond_type[nb] = 2
        nb += 1
        
        # nn in x dir
        n1 = 2 * self.inv_r_unit_cell[(i + 1)%self.Lx, j] + 1
        self.bond_list[nb, 0] = n
        self.bond_list[nb, 1] = n1
        self.bond_type[nb] = 1
        nb += 1
    
    # Brillouin zone for plotting
    gamma = 0.0 * self.b1 + 0.0 * self.b2
    gamma_prime = 1.0 * self.b1 + 0.0 * self.b2
    m1 = 0.5 * self.b1 + 0.5 * self.b2
    m3 = 0.5 * self.b1 + 0 * self.b2
    x = -0.5 * self.b1 + 0.5 * self.b2
    
    self.q_points = list()
 
    self.gamma_idx = 0 
    self.q_points.append(gamma)
    for n in range(self.N):
      if self.k_unit_cell[n, 0] == gamma[0] and self.k_unit_cell[n, 1] < m1[1] and self.k_unit_cell[n, 1] > gamma[1]:
        self.q_points.append(self.k_unit_cell[n, :])
    self.q_points.append(m1)
    self.m1_idx = len(self.q_points) - 1
    for n in range(self.N):
      if self.k_unit_cell[n, 1] == m1[1] and self.k_unit_cell[n, 0] > m1[0] and self.k_unit_cell[n, 0] < gamma_prime[0]:
        self.q_points.append(self.k_unit_cell[n, :])
    self.q_points.append(gamma_prime)
    self.gamma_prime_idx = len(self.q_points) - 1
    for i in range(self.Lx-1,-1,-1):
      n = self.inv_k_unit_cell[i, 0]
      if self.k_unit_cell[n, 1] > m3[1] and self.k_unit_cell[n, 0] > m3[0]:
        self.q_points.append(self.k_unit_cell[n, :])
    self.q_points.append(m3)
    self.m3_idx = len(self.q_points) - 1
    for i in range(self.Lx-1,-1,-1):
      n = self.inv_k_unit_cell[i, 0]
      if self.k_unit_cell[n, 1] < m3[1] and self.k_unit_cell[n, 0] < m3[0]:
        self.q_points.append(self.k_unit_cell[n, :])

    self.q_points = np.array(self.q_points)

    # Plaquettes 
    self.Np = self.N
    self.size_p = 6
    self.plaquettes_sites = np.zeros((self.Np, self.size_p), dtype=int)
    self.plaquettes_bonds = np.zeros((self.Np, self.size_p), dtype=int)
    
    self.z_bonds = list()
    self.y_bonds = list()
    self.x_bonds = list()
    p = 0
    for b in range(self.Nb):
      if self.bond_type[b] == 3:
        self.z_bonds.append(b)
      if self.bond_type[b] == 2:
        self.y_bonds.append(b)
      if self.bond_type[b] == 1:
        self.x_bonds.append(b)
    self.x_bonds = np.array(self.x_bonds)
    self.y_bonds = np.array(self.y_bonds)
    self.z_bonds = np.array(self.z_bonds)
    
    for p in range(self.Np):
      self.plaquettes_sites[p, 0] = self.bond_list[self.z_bonds[p], 0]
      self.plaquettes_sites[p, 1] = self.bond_list[self.z_bonds[p], 1]
      self.plaquettes_bonds[p, 0] = self.z_bonds[p]
      
      for i in range(self.x_bonds.shape[0]):
        if self.bond_list[self.x_bonds[i], 1] == self.plaquettes_sites[p, 1]:
          self.plaquettes_sites[p, 2] = self.bond_list[self.x_bonds[i], 0]
          self.plaquettes_bonds[p, 1] = self.x_bonds[i]
          break

      for i in range(self.y_bonds.shape[0]):
        if self.bond_list[self.y_bonds[i], 0] == self.plaquettes_sites[p, 2]:
          self.plaquettes_sites[p, 3] = self.bond_list[self.y_bonds[i], 1]
          self.plaquettes_bonds[p, 2] = self.y_bonds[i]
          break
        
      for i in range(self.z_bonds.shape[0]):
        if self.bond_list[self.z_bonds[i], 1] == self.plaquettes_sites[p, 3]:
          self.plaquettes_sites[p, 4] = self.bond_list[self.z_bonds[i], 0]
          self.plaquettes_bonds[p, 3] = self.z_bonds[i]
          break
        
      for i in range(self.x_bonds.shape[0]):
        if self.bond_list[self.x_bonds[i], 0] == self.plaquettes_sites[p, 4]:
          self.plaquettes_sites[p, 5] = self.bond_list[self.x_bonds[i], 1]
          self.plaquettes_bonds[p, 4] = self.x_bonds[i]
          break
        
      for i in range(self.y_bonds.shape[0]):
        if self.bond_list[self.y_bonds[i], 1] == self.plaquettes_sites[p, 5]:
          self.plaquettes_bonds[p, 5] = self.y_bonds[i]
          break
    

