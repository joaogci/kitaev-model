import numpy as np
import matplotlib.pyplot as plt
from itertools import product

from lattice import Lattice
from hamiltonian import Hamiltonian

SEED = 2024
np.random.seed(SEED)

Lx, Ly = 32, 32
K = [1.0, 1.0, 1.0]

latt = Lattice(Lx, Ly)
ham = Hamiltonian(K, latt)
ham.set_gs_flux()
ham.diagonalise()

print(f"E_gs = {ham.obs_gs_energy()}")
print(f"S_ij(T = 0) = {ham.obs_spinZ_eq(0.0)}")
print(f"S_ij(T = 0.01) = {ham.obs_spinZ_eq(0.01)}")
print(f"S_ij(T = 0.1) = {ham.obs_spinZ_eq(0.1)}")
print(f"S_ij(T = 1.0) = {ham.obs_spinZ_eq(1.0)}")

plt.figure(1)
plt.title("Eigenvalues for ground state energy")
plt.plot(ham.E_full, ".")
plt.xlabel(r"$i$")
plt.ylabel(r"$E_i$")
plt.show()


# Generate all possible flux configurations
# all_flux_configurations = list(product([1, -1], repeat=latt.Nb))

# energies = list()
# sum_flux = list()

# for flux in all_flux_configurations:
#   ham.set_flux(flux)
#   ham.diagonalise()
#   erg = np.round(ham.obs_gs_energy(), 8)
  
#   if np.sum(flux) not in sum_flux:
#     sum_flux.append(np.sum(flux))
#   if erg not in energies:
#     energies.append(erg)

# energies.sort()
# sum_flux.sort()

# histogram = np.zeros((len(sum_flux), len(energies)))

# for flux in all_flux_configurations:
#   ham.set_flux(flux)
#   ham.diagonalise()
#   erg = np.round(ham.obs_gs_energy(), 8)
  
#   idx_erg = energies.index(erg)
#   idx_flux = sum_flux.index(np.sum(flux))
#   histogram[idx_flux, idx_erg] += 1

# import seaborn as sns

# plt.figure(2)
# plt.title("Histogram of ground state energy for all flux sectors")
# sns.heatmap(histogram, annot=True)
# plt.xlabel(r"$E$")
# plt.ylabel(r"Sum of Flux")
# plt.show() 

