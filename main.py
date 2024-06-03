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

Lx, Ly = 4, 4
K = [1.0, 1.0, 1.0]
beta = 1.0

n_sweeps = 100
n_bins = 10
n_rebin = 1

n_sims = 1

for i in range(n_sims): 
  latt = Lattice(Lx, Ly)
  flux = Flux(K, latt)

  mc = MonteCarlo(n_sweeps, n_bins, beta)
  mc.init_observables(1, 0, latt)
  mc.set_flux(flux)

  sim_name = f"Lx{Lx}_Ly{Ly}_Kx{K[0]}_Ky{K[1]}_Kz{K[2]}_beta{int(beta)}"
  if not os.path.exists(sim_name):
    os.mkdir(sim_name)
  os.chdir(sim_name)
  
  print(f"Starting: {sim_name}")  
  #mc.simulation()
  print(f"Analysing: {sim_name} ")

  ana = Analysis(n_rebin)
  ana.analyse()

  os.chdir("../")
