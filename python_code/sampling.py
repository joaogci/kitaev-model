import numpy as np
from observables import Obs, ObsScalar, ObsLatt
from flux import Flux
from lattice import Lattice

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
          self.obs_scal.append(ObsScalar("Sz"))
        case 2:
          self.obs_scal.append(ObsScalar("W"))
        case 3:
          ...
        case _:
          print("Observable not found. ")
          exit(1) 
    
    self.obs_eq = list()
    for i in range(self.n_eq):
      match i:
        case 0:
          self.obs_eq.append(ObsLatt("Dimer", latt))
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
    
    res = complex(0.0)
    for b in self.__latt.z_bonds:
      i = self.__latt.unit_cell[self.__latt.bond_list[b, 0]]
      j = self.__latt.unit_cell[self.__latt.bond_list[b, 1]]
      Fij = flux.flux[b]
      for m in range(self.__latt.N):
        res += Fij * A[i, m] * B[j, m] * (1.0 - 2.0 * flux.fermi_function(m, beta))
    res = res / self.__latt.N

    self.obs_scal[1].obs_vec += res
    
    res = complex(0.0)
    for p in range(self.__latt.Np):
      tmp = complex(1.0)
      for b in range(self.__latt.size_p):
        tmp *= flux.flux[self.__latt.plaquettes_bonds[p, b]]
      res += tmp
    res = np.abs(res) / self.__latt.N
    
    self.obs_scal[2].obs_vec += res
    
  def sample_exact(self, flux: Flux, beta: float):
    for i in range(self.n_scal):
      self.obs_scal[i].N += np.exp(flux.ln_weight(beta))
    for i in range(self.n_eq):
      self.obs_eq[i].N += np.exp(flux.ln_weight(beta))
    
    A = flux.X.T + flux.Y.T
    B = flux.X.T - flux.Y.T
    
    # i = self.__latt.unit_cell[self.__latt.bond_list[0, 0]]
    # j = self.__latt.unit_cell[self.__latt.bond_list[0, 1]]
    # Fij = flux.flux[0]
    res = complex(0.0)
    for m in range(self.__latt.N):
      res += - flux.E[m] * np.tanh(0.5 * beta * flux.E[m]) * 0.5
    
    self.obs_scal[0].obs_vec += np.exp(flux.ln_weight(beta)) * res
    
    res = complex(0.0)
    for b in self.__latt.z_bonds:
      i = self.__latt.unit_cell[self.__latt.bond_list[b, 0]]
      j = self.__latt.unit_cell[self.__latt.bond_list[b, 1]]
      Fij = flux.flux[b]
      for m in range(self.__latt.N):
        res += Fij * A[i, m] * B[j, m] * (1.0 - 2.0 * flux.fermi_function(m, beta))
    res = res / self.__latt.N

    self.obs_scal[1].obs_vec += np.exp(flux.ln_weight(beta)) * res
    
    res = complex(0.0)
    for p in range(self.__latt.Np):
      tmp = complex(1.0)
      for b in range(self.__latt.size_p):
        tmp *= flux.flux[self.__latt.plaquettes_bonds[p, b]]
      res += tmp
    res = np.abs(res) / self.__latt.N
    
    self.obs_scal[2].obs_vec += np.exp(flux.ln_weight(beta)) * res

    for p, b in enumerate(self.__latt.z_bonds):
      i = self.__latt.unit_cell[self.__latt.bond_list[b, 0]]
      j = self.__latt.unit_cell[self.__latt.bond_list[b, 1]]
      Fij = flux.flux[b]
      
      for pp, bp in enumerate(self.__latt.z_bonds):
        ip = self.__latt.unit_cell[self.__latt.bond_list[bp, 0]]
        jp = self.__latt.unit_cell[self.__latt.bond_list[bp, 1]]
        Fijp = flux.flux[bp]
        
        res_PP = complex(0.0)
        res_PH = complex(0.0)
        
        for k in range(self.__latt.N):
          for q in range(self.__latt.N):
            res_PP += A[i, k] * B[j, q] * (A[ip, k] * B[jp, q] - A[ip, q] * B[jp, k]) * (flux.fermi_function(k, beta) * flux.fermi_function(q, beta) + (1.0 - flux.fermi_function(k, beta)) * (1.0 - flux.fermi_function(q, beta)))
            res_PH += A[i, k] * B[j, q] * (A[ip, k] * B[jp, q] + A[ip, q] * B[jp, k]) * ((1.0 - flux.fermi_function(k, beta)) * flux.fermi_function(q, beta) + (1.0 - flux.fermi_function(q, beta)) * flux.fermi_function(k, beta))

        self.obs_eq[0].obs_latt[p][pp] += np.exp(flux.ln_weight(beta)) * Fij * Fijp * (res_PP + res_PH)

  def write_to_file(self):
    for i in range(self.n_scal):
      self.obs_scal[i].write_to_file()
    for i in range(self.n_eq):
      self.obs_eq[i].write_to_file()
  
  def reset(self):
    for i in range(self.n_scal):
      self.obs_scal[i].reset()
    for i in range(self.n_eq):
      self.obs_eq[i].reset()

