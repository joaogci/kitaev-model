import numpy as np
from .lattice import Lattice

class Flux():
  """
      Flux class
    Stores the flux configuration 

    Diagonalises the Matter Hamiltonian for a given flux configuration
  """
  
  def __init__(self, K: list[float], latt: Lattice):
    self.K = K
    self.__latt = latt
    
  def set_flux(self, flux: list[int]):
    self.flux = flux
    
    self.__F = np.zeros((self.__latt.N, self.__latt.N), dtype=int)
    for nb in range(self.__latt.Nb):
      i = self.__latt.unit_cell[self.__latt.bond_list[nb, 0]]
      j = self.__latt.unit_cell[self.__latt.bond_list[nb, 1]]
      self.__F[i, j] = flux[nb]
      
    self.__set_H()
  
  def diagonalise_flux(self):
    E_tmp, U = np.linalg.eigh(self.__H_matter)
    
    self.X = U[:self.__latt.N, self.__latt.N:].T
    self.Y = U[self.__latt.N:, self.__latt.N:].T
    self.E = E_tmp[self.__latt.N:]
    self.T = np.block([[self.X.conj(), self.Y.conj()], [self.Y, self.X]])
    
  def weight(self, beta: float):
    """
      Returns W(F) = prod_n 2 cosh(beta E_n / 2)
    """
    return np.prod(2.0 * np.cosh(0.5 * beta * self.E))

  def ln_weight(self, beta: float):
    """
      Returns ln(W(F)) = sum_n ln(2) + ln(cosh(beta E_n / 2))
    """
    return np.sum(np.log(2.0) + np.log(np.cosh(0.5 * beta * self.E)))

  def fermi_function(self, m: int, beta: float):
    return 1.0 / (np.exp(beta * self.E[m]) + 1.0)
  
  def __set_H(self):
    #self.H_majoranas = 2.0j * np.block([[np.zeros(self.F.shape), self.F], [- self.F.T, np.zeros(self.F.shape)]])
    
    h = self.__F + self.__F.T
    d = self.__F.T - self.__F
    self.__H_matter = np.block([[h, d], [d.conj().T, - h.T]])
