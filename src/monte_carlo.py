import numpy as np
from .flux import Flux
from .lattice import Lattice
from .sampling import Sampling
import os
from tqdm import tqdm
import time

class MonteCarlo():
  def __init__(self, n_sweeps: int, n_bins: int, beta: float):
    self.n_sweeps = n_sweeps
    self.n_bins = n_bins
    self.beta = beta
    self.wall_time = 0.0

  def set_flux(self, flux: Flux):
    self.__flux = flux
  
  def init_observables(self, n_scal: int, n_eq: int, latt: Lattice):
    self.__sampling_obs = Sampling(n_scal, n_eq, latt)
    self.__latt = latt
  
  def simulation(self):
    self.accepted = 0

    if os.path.exists("confout"):
      self.__read_configuration()
    else:
      flux = np.ones(self.__latt.Nb)
      self.__flux.set_flux(flux)
    self.__flux.diagonalise_flux()
    
    start = time.time()
     
    for b in tqdm(range(self.n_bins)):
      self.__sampling_obs.reset()
      
      for t in range(self.n_sweeps):
        for i in range(self.__latt.N):
          b_ = np.random.randint(self.__latt.Nb)
          self.__mc_step(b_)
          
        self.__sampling_obs.sample(self.__flux, self.beta)
      
      self.__sampling_obs.write_to_file()
      self.__write_configuration()
    
    self.wall_time = time.time() - start
    self.__write_info()
 
  def __mc_step(self, b: int):
    flux_new = Flux(self.__flux.K, self.__latt)
    flux = self.__flux.flux
    new_flux = flux.copy()
    new_flux[b] = - new_flux[b]
    flux_new.set_flux(new_flux)
    
    flux_new.set_flux(new_flux)
    flux_new.diagonalise_flux()
    
    ratio = self.__ratio(flux_new)
    
    if np.random.rand() <= ratio:
      self.__flux.set_flux(new_flux)
      self.__flux.X = flux_new.X.copy()
      self.__flux.Y = flux_new.Y.copy()
      self.__flux.T = flux_new.T.copy()
      self.__flux.E = flux_new.E.copy()
      self.accepted += 1
  
  def __ratio(self, flux_new: Flux):
    # return flux_new.weight(self.beta) / self.__flux.weight(self.beta)
    # print("new: ", flux_new.ln_weight(self.beta), "old: ", self.__flux.ln_weight(self.beta), "diff: ", flux_new.ln_weight(self.beta) - self.__flux.ln_weight(self.beta))
    return np.exp(flux_new.ln_weight(self.beta) - self.__flux.ln_weight(self.beta))
 
  def __write_info(self):
    with open("info", "a") as f:
      f.write("Simulation Finished \n")
      f.write(f"Lx, Ly: {self.__latt.Lx}, {self.__latt.Ly} \n")
      f.write(f"Kx, Ky, Kz: {self.__flux.K[0]}, {self.__flux.K[1]}, {self.__flux.K[2]} \n")
      f.write(f"beta: {self.beta} \n")
      f.write(f"n_bins, mc_sweeps: {self.n_bins}, {self.n_sweeps} \n")
      f.write(f"Accepted ration: {self.accepted/(self.n_bins * self.n_sweeps * self.__latt.N)} \n")
      f.write(f"Time: {self.wall_time}s \n")
      f.write("----------------------- \n")
  
  def __write_configuration(self):
    np.savetxt("confout", self.__flux.flux)
  
  def __read_configuration(self):
    flux = np.loadtxt("confout")
    self.__flux.set_flux(flux)