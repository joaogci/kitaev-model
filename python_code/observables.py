import numpy as np
from lattice import Lattice

class Obs():
  def __init__(self, name: str):
    self.name = name
    self.N = 0
  
  def write_to_file(self):
    pass
  
  def reset(self):
    pass

class ObsScalar(Obs):
  def __init__(self, name: str):
    super().__init__(name)
    self.obs_vec = complex(0.0)
     
  def write_to_file(self):
    self.obs_vec = self.obs_vec / self.N
    
    with open(self.name + "_scal", "a") as f:
      f.write(f"{self.obs_vec.real} {self.obs_vec.imag}\n")
      
    with open(self.name + "_scal_info", "w") as f:
      f.write("-- Analysis Mode --\n")
      f.write("identity\n")
  
  def reset(self):
    self.N = 0
    self.obs_vec = complex(0.0)

class ObsLatt(Obs):
  def __init__(self, name: str, latt: Lattice):
    super().__init__(name)
    self.latt = latt
    self.obs_latt = np.zeros((self.latt.N, self.latt.N), dtype=complex)
    
    self.obs_k = np.zeros(self.latt.N, dtype=complex)
    self.obs_r = np.zeros(self.latt.N, dtype=complex)
  
  def write_to_file(self):
    for i in range(self.latt.N):
      for j in range(self.latt.N):
        self.obs_latt[i, j] = self.obs_latt[i, j] / self.N
    
    self.fourier_transform() 
    self.inv_fourier_transform()
    
    with open(self.name + "_eqR", "a") as f:
      for i in range(self.latt.N):
        f.write(f"{self.latt.r_unit_cell[i, 0]} {self.latt.r_unit_cell[i, 1]} {np.real(self.obs_r[i])} {np.imag(self.obs_r[i])}\n")

    with open(self.name + "_eqK", "a") as f:
      for i in range(self.latt.N):
        f.write(f"{self.latt.k_unit_cell[i, 0]} {self.latt.k_unit_cell[i, 1]} {np.real(self.obs_k[i])} {np.imag(self.obs_k[i])}\n")
  
  def reset(self):
    self.N = 0
    self.obs_latt = np.zeros((self.latt.N, self.latt.N), dtype=complex)

  def fourier_transform(self):
    for n in range(self.latt.N):
      for i in range(self.latt.N):
        for j in range(self.latt.N):
          self.obs_k[n] += np.exp(- 1.0j * self.latt.k_unit_cell[n,:]@(self.latt.r_unit_cell[i,:] - self.latt.r_unit_cell[j,:])) * self.obs_latt[i,j] / self.latt.N

  def inv_fourier_transform(self):
    for i in range(self.latt.N):
      for n in range(self.latt.N):
        self.obs_r[i] += np.exp(1.0j * self.latt.k_unit_cell[n,:]@self.latt.r_unit_cell[i,:]) * self.obs_k[n] / self.latt.N
