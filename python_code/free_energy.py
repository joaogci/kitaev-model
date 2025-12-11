import numpy as np
import matplotlib.pyplot as plt
import os
from itertools import product
from tqdm import tqdm

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
T_vals = np.logspace(-1.5, -0.5, 100, base=10)
beta_vals = 1.0 / T_vals
n_sims = len(T_vals)

beta_vals[0] = 1
beta_vals[1] = 10
beta_vals[3] = 100
n_sims = 3

F = np.zeros(n_sims)
Z = np.zeros(n_sims)

for i in tqdm(range(n_sims)):
  latt = Lattice(Lx, Ly)
  flux = Flux(K, latt)
  
  all_flux_configurations = product([-1.0, 1.0], repeat=latt.Nb)
  
  weights = list()
  
  for flux_conf in all_flux_configurations:
    flux.set_flux(flux_conf)
    flux.diagonalise_flux()
    
    weights.append(flux.weight(beta_vals[i]))
    
  F[i] = - np.log(np.sum(weights)) / beta_vals[i]
  Z[i] = np.sum(weights)

print(Z, F, F / latt.N, F / latt.Ns)

exit()

for i in tqdm(range(n_sims)):
  latt = Lattice(Lx, Ly)
  flux = Flux(K, latt)
  
  weights = list()
  
  for _ in range(1000):
    flux_conf = np.random.randint(2, size=latt.Nb, dtype=int)
    flux_conf[flux_conf == 0] = -1
    
    flux.set_flux(flux_conf)
    flux.diagonalise_flux()
    
    weights.append(flux.weight(beta_vals[i]))
  
  F_random[i] = - np.log(np.sum(weights) * np.power(2, latt.Nb) / len(weights)) / beta_vals[i]

F_random = F_random / latt.Ns

plt.figure(1)
# plt.plot(np.log10(T_vals), F_exact, label="F exact")
plt.plot(np.log10(T_vals), F_random, label="F random")
plt.legend()

# plt.figure(2)
# plt.plot(np.log10(T_vals), np.abs(F_exact - F_random) / F_exact, label="Relative Error of F random")
# plt.legend()

plt.figure(3)

# S = - np.gradient(F_exact, T_vals)
# C = T_vals * np.gradient(S, T_vals)
# plt.plot(np.log10(T_vals), S - S[0], label="S exact")
# plt.plot(np.log10(T_vals), C, label="C exact")

S = - np.gradient(F_random, T_vals)
C = T_vals * np.gradient(S, T_vals)
plt.plot(np.log10(T_vals), S - S[0], label="S random")
plt.plot(np.log10(T_vals), C, label="C random")


plt.legend()

plt.show()
