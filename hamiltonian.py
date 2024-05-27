import numpy as np
from lattice import Lattice

class Hamiltonian():
  def __init__(self, K: list[int], latt: Lattice):
    self.__latt = latt
    self.K = K

  def set_gs_flux(self):
    self.F = np.zeros((self.__latt.N, self.__latt.N), dtype=int)
    for nb in range(self.__latt.Nb):
      self.F[self.__latt.unit_cell[self.__latt.bond_list[nb, 0]], self.__latt.unit_cell[self.__latt.bond_list[nb, 1]]] = 1
        
    self.__set_H()
    
  def set_random_flux(self):
    ...

  def set_flux(self, flux: list[int]):
    self.F = np.zeros((self.__latt.N, self.__latt.N), dtype=int)
    for nb in range(self.__latt.Nb):
      self.F[self.__latt.unit_cell[self.__latt.bond_list[nb, 0]], self.__latt.unit_cell[self.__latt.bond_list[nb, 1]]] = flux[nb]
      
    self.__set_H()
  
  def diagonalise(self):
    self.E_full, self.U = np.linalg.eigh(self.H_matter)
  
    self.E = self.E_full[self.__latt.N:]
    
    V = self.U[:, self.__latt.N:]
    self.X = V[:self.__latt.N, :]
    self.Y = V[self.__latt.N:, :]
    self.X = self.X.T
    self.Y = self.Y.T
    
    self.T = np.block([[self.X.conj(), self.Y.conj()], [self.Y, self.X]])
    
  def obs_gs_energy(self):
    return - 0.5 * np.sum(self.E) / self.__latt.N
  
  def obs_gs_spinZ_eq(self):
    i = self.__latt.unit_cell[self.__latt.bond_list[0, 0]]
    j = self.__latt.unit_cell[self.__latt.bond_list[0, 1]]
    Bij = self.F[i, j]
    
    A = self.X.T + self.Y.T
    B = self.X.T - self.Y.T
    
    res = 0.0
    for m in range(len(self.E)):
      res += Bij * A[i, m] * B[j, m]

    return res

  def ground_state_energy(self):
    self.set_gs_flux()
    self.diagonalise()
    
    return - 0.5 * np.sum(self.E) / self.__latt.N

  def ana_ground_state_dispersion(self, q: np.ndarray):
    return self.K[2] + self.K[1] * np.exp(q @ self.__latt.n2 * 1j) + self.K[0] * np.exp(q @ self.__latt.n1 * 1j)
  
  def ana_gound_state_energy(self):
    res = 0.0
    for n in range(self.__latt.N):
      res += np.abs(self.ana_ground_state_dispersion(self.__latt.k_unit_cell[n, :]))
      
    return - res / self.__latt.N

  def __set_H(self):
    self.H_majoranas = 2.0j * np.block([[np.zeros(self.F.shape), self.F], [- self.F.T, np.zeros(self.F.shape)]])
    
    h = self.F + self.F.T
    d = self.F.T - self.F
    self.H_matter = np.block([[h, d], [d.conj().T, - h.T]])

  def __fermi_function(self, m: int, T: float):
    beta = self.__T_to_beta(T)
    if T == 0.0:
      return 1.0
    return 1.0 / (np.exp(beta * self.E[m]) + 1.0)

  def __T_to_beta(self, T: float):
    return 1.0 / T if (T != 0.0) else 1e16
