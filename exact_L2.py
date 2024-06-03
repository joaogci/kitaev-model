import numpy as np
import matplotlib.pyplot as plt
import os
from itertools import product

from src.sampling import Sampling
from src.lattice import Lattice
from src.flux import Flux
from src.analysis import Analysis

import scienceplots
plt.style.use("science")

SEED = 2024
np.random.seed(SEED)

Lx, Ly = 2, 2
K = [1.0, 1.0, 1.0]
beta = 1.0

n_sweeps = 100
n_bins = 25
n_rebin = 1

n_sims = 1

for i in range(n_sims): 
  latt = Lattice(Lx, Ly)
  flux = Flux(K, latt)
  sampling_obs = Sampling(1, 0, latt)

  all_flux_configurations = product([-1.0, 1.0], repeat=latt.Nb)

  sim_name = f"Lx{Lx}_Ly{Ly}_Kx{K[0]}_Ky{K[1]}_Kz{K[2]}_beta{int(beta)}_exact"
  if not os.path.exists(sim_name):
    os.mkdir(sim_name)
  os.chdir(sim_name)
  
  print(f"Starting: {sim_name}")  
  
  sampling_obs.reset()
  for flux_conf in all_flux_configurations:
    flux.set_flux(flux_conf)
    flux.diagonalise_flux()
    
    sampling_obs.sample_exact(flux, beta)
  
  sampling_obs.write_to_file()

  print(f"Analysing: {sim_name} ")

  ana = Analysis(n_rebin)
  ana.analyse()

  os.chdir("../")
