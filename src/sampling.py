import numpy as np
from .observables import Obs, ObsScalar, ObsLatt
from .flux import Flux
from .lattice import Lattice

class Sampling():
  def __init__(self, n_scal: int, n_eq: int, latt: Lattice):
    self.n_scal = n_scal
    self.n_eq = n_eq
    self.__latt = latt

    self.obs_scal = list()
    for i in range(self.n_scal):
      match i:
        case 0:
          self.obs_scal.append(ObsScalar("Sz_0A0B"))
        case 1:
          ...
        case _:
          print("Observable not found. ")
          exit(1) 
    
    self.obs_eq = list()
    for i in range(self.n_eq):
      match i:
        case 0:
          ...
        case _:
          print("Observable not found. ")
          exit(1) 

  def sample(self, flux: Flux, beta: float):
    for i in range(self.n_scal):
      self.obs_scal[i].N += 1
    
    A = flux.X.T + flux.Y.T
    B = flux.X.T - flux.Y.T
    
    i = self.__latt.unit_cell[self.__latt.bond_list[0, 0]]
    j = self.__latt.unit_cell[self.__latt.bond_list[0, 1]]
    Fij = flux.flux[0]
        
    res = complex(0.0)
    for m in range(self.__latt.N):
      res += Fij * A[i, m] * B[j, m] * (1.0 - 2.0 * flux.fermi_function(m, beta))
    
    self.obs_scal[0].obs_vec += res
  
  def write_to_file(self):
    for i in range(self.n_scal):
      self.obs_scal[i].write_to_file()
  
  def reset(self):
    for i in range(self.n_scal):
      self.obs_scal[i].reset()



  # def obs_spinZ_eq(self, T: float):
  #     Equal time spin-spin correlations 
  #     S_{ij}^{zz} = \langle \sigma^z_{iA} \sigma^z_{jB} \rangle
    
  #   i = self.__latt.unit_cell[self.__latt.bond_list[0, 0]]
  #   j = self.__latt.unit_cell[self.__latt.bond_list[0, 1]]
  #   Bij = self.F[i, j]
        
  #   res = 0.0
  #   for m in range(len(self.E)):
  #     res += Bij * self.A[i, m] * self.B[j, m] * (1.0 - 2.0 * self.__fermi_function(m, T))

  #   return res
