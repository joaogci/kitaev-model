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
    self.E_full, self.U_full = np.linalg.eigh(self.H_majoranas)
    
    self.E = np.zeros(self.__latt.N)
    self.U = np.zeros((self.__latt.Ns, self.__latt.N), dtype=complex)
    n = 0
    for i in range(self.__latt.Ns):
      if self.E_full[i] > 0.0:
        self.E[n] = self.E_full[i]
        self.U[:, n] = self.U_full[:, i]
        n += 1

  def obs_energy(self):
    return - 0.5 * np.sum(self.E) / self.__latt.N
  
  def obs_spinZ_eq(self, T: float):
    i = self.__latt.bond_list[0, 0]
    j = self.__latt.bond_list[0, 1]
    
    res = 0.0
    for m in range(self.__latt.N):
      res += np.conj(self.U[i, m]) * self.U[j, m] * self.F[self.__latt.unit_cell[i], self.__latt.unit_cell[j]] * (2.0 * self.__fermi_function(m, T) - 1.0)
  
    return - 1.0j * res

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
    self.H_matter = np.block([[h, d], [np.conj(d), - h.T]])

  def __fermi_function(self, m: int, T: float):
    beta = self.__T_to_beta(T)
    if T == 0.0:
      return 1.0
    return 1.0 / (np.exp(beta * self.E[m]) + 1.0)

  def __T_to_beta(self, T: float):
    return 1.0 / T if (T != 0.0) else 1e16
