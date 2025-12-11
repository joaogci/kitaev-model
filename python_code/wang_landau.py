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

# Kitaev parameters
Lx, Ly = 2, 2
K = [1.0/3.0, 1.0/3.0, 1.0/3.0]
beta = 10

# WL parameters
f = np.exp(1)
p = 0.9
p_warmup = 0.75

def find_F_bin(e):
  for i in range(nF):
    if e >= F_vals[i] and e <= F_vals[i + 1]:
      return i

  print(f"Error: IDX NOT FOUND! {e} {i} {nF}")
  exit()

latt = Lattice(Lx, Ly)

# Warm up
# Find the energies which g(F) > 0
dF = 0.05
F_min = 0.0 * beta * latt.N
F_max = 1.0 * beta * latt.N
nF = int((F_max - F_min) / dF)
F_vals = np.linspace(F_min, F_max, nF + 1)

flux = Flux(K, latt)
flux.set_flux(np.ones(latt.Nb))
flux.diagonalise_flux()

F_old = flux.ln_weight(beta)
idx_old = find_F_bin(F_old)

n_warmup = 1000
mc_sweeps = 0
ln_gF = np.zeros(nF)
hist = np.zeros(nF)
visited_energies = list()

while True:
  for i in range(latt.N):
    b = np.random.randint(latt.Nb)

    flux.flux[b] = - flux.flux[b]
    flux.set_flux(flux.flux)
    flux.diagonalise_flux()

    F_new = flux.ln_weight(beta)
    idx_new = find_F_bin(F_new)

    ratio = np.exp(ln_gF[idx_old] - ln_gF[idx_new])

    if np.random.rand() <= ratio:
      idx_old = idx_new
      F_old = F_new
    else:
      flux.flux[b] = - flux.flux[b]

    hist[idx_old] += 1
    ln_gF[idx_old] += np.log(f)
    if np.round(F_old, 8) not in visited_energies:
      visited_energies.append(np.round(F_old, 8))

  mc_sweeps += 1

  if mc_sweeps % n_warmup == 0:
    if np.min(hist[hist > 0]) > np.average(hist[hist > 0]) * p_warmup:
      print(f"warm up done; found {len(np.where(hist != 0)[0])} energies")
      break


# Wang-Landau simulation

F_vals = np.array(sorted(visited_energies))
nF = F_vals.shape[0]

flux = Flux(K, latt)
flux.set_flux(np.ones(latt.Nb))
flux.diagonalise_flux()

F_old = np.round(flux.ln_weight(beta), 8)
tmp = np.where(F_vals == F_old)[0]
if len(tmp) != 1:
  print(f"ERROR: More than one index found for energy {F_old} {tmp}")
  print(f"Known energies: {F_vals}")
  exit()
idx_old = tmp[0]

mc_sweeps = 0
ln_gF = np.zeros(nF)
hist = np.zeros(nF)

while f > 1 + 1e-8:
  for i in range(latt.N):
    b = np.random.randint(latt.Nb)
    
    flux.flux[b] = - flux.flux[b]
    flux.set_flux(flux.flux)
    flux.diagonalise_flux()

    F_new = np.round(flux.ln_weight(beta), 8)
    tmp = np.where(F_vals == F_new)[0]
    if len(tmp) != 1:
      print(f"ERROR: More than one index found for energy {F_old} {tmp}")
      exit()
    idx_new = tmp[0]
    
    ratio = np.exp(ln_gF[idx_old] - ln_gF[idx_new])
    
    if np.random.rand() <= ratio:
      idx_old = idx_new
    else:
      flux.flux[b] = - flux.flux[b]
      
    hist[idx_old] += 1
    ln_gF[idx_old] += np.log(f)
  
  mc_sweeps += 1
  
  if mc_sweeps % 1000 == 0:
    if np.min(hist) > np.average(hist) * p:
      hist = np.zeros(nF)
      print(f"after {mc_sweeps} the histogram is flat. f: {f:.8f}")
      f = np.sqrt(f)
      mc_sweeps = 0

Z = 0.0
ln_gFn = np.zeros(nF)
s = 0
for i in range(1, nF):
  s += np.exp(ln_gF[i] - ln_gF[0])
norm = ln_gF[0] + np.log(1 + s)

for i in range(nF):
  ln_gFn[i] = ln_gF[i] + latt.Nb * np.log(2) - norm
  Z += np.exp(F_vals[i]) * np.exp(ln_gFn[i])

for i in range(nF):
  print(f"g({F_vals[i]}) = {np.exp(ln_gFn[i])}")

F = - np.log(Z) / beta

print(f"Z: {Z}")
print(f"F / N: {F / latt.N}")

