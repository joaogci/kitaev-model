import numpy as np
import matplotlib.pyplot as plt
import os
from itertools import product

from sampling import Sampling
from lattice import Lattice
from flux import Flux
from analysis import Analysis

import scienceplots
plt.style.use("science")

SEED = 2024
np.random.seed(SEED)

Lx, Ly = 2, 2
K = [1.0/3.0, 1.0/3.0, 1.0/3.0]

# T_vals = np.arange(0.01, 4.0 + 0.05, 0.05)
T_vals = np.logspace(-2.5, 1, 400, base=10)
beta_vals = 1.0 / T_vals

n_sims = len(T_vals)
n_sims = 4
beta_vals[0] = 100
beta_vals[1] = 50
beta_vals[2] = 10
beta_vals[3] = 1
T_vals = 1.0 / beta_vals

for i in range(n_sims): 
  latt = Lattice(Lx, Ly)
  flux = Flux(K, latt)
  sampling_obs = Sampling(3, 1, latt)
  
  # all_flux_configurations = product([-1.0, 1.0], repeat=latt.Nb)

  sim_name = f"Lx{Lx}_Ly{Ly}_alpha1.0_beta{beta_vals[i]:.0f}_gs_flux"
  if not os.path.exists(sim_name):
    os.mkdir(sim_name)
  os.chdir(sim_name)
  
  print(f"Starting: {sim_name}")  
  
  sampling_obs.reset()
#  for flux_conf in all_flux_configurations:
  flux_conf = np.ones(latt.Nb)
  flux.set_flux(flux_conf)
  flux.diagonalise_flux()

  sampling_obs.sample_exact(flux, beta_vals[i])
  
  sampling_obs.write_to_file()

  print(f"Analysing: {sim_name} ")

  ana = Analysis(1, 0)
  ana.analyse()

  os.chdir("../")
