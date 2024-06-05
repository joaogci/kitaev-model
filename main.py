import numpy as np
import matplotlib.pyplot as plt
import os

from src.monte_carlo import MonteCarlo
from src.lattice import Lattice
from src.flux import Flux
from src.analysis import Analysis

import scienceplots
plt.style.use("science")

SEED = 2024
np.random.seed(SEED)

Lx, Ly = 12, 12
K = [1.0/3.0, 1.0/3.0, 1.0/3.0]

T_vals = np.logspace(-2.5, 1, 50, base=10)
# T_vals = np.arange(0.01, 4.0 + 0.05, 0.05)
beta_vals = 1.0 / T_vals

n_sims = len(T_vals)

n_sweeps = 100
n_bins = 500
n_skip = 100
n_rebin = 1

for i in range(n_sims): 
  latt = Lattice(Lx, Ly)
  flux = Flux(K, latt)

  mc = MonteCarlo(n_sweeps, n_bins, beta_vals[i])
  mc.init_observables(3, 0, latt)
  mc.set_flux(flux)

  sim_name = f"Lx{Lx}_Ly{Ly}_Kx{K[0]:.1f}_Ky{K[1]:.1f}_Kz{K[2]:.1f}_logT{np.log10(T_vals[i]):.3f}"
  if not os.path.exists(sim_name):
    os.mkdir(sim_name)
  os.chdir(sim_name)
  
  print(f"Starting: {sim_name}")  
  mc.simulation()
  print(f"Analysing: {sim_name} ")

  ana = Analysis(n_rebin, n_skip)
  ana.analyse()

  os.chdir("../")
