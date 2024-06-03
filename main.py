import numpy as np
import matplotlib.pyplot as plt
from itertools import product

from src.monte_carlo import MonteCarlo
from src.lattice import Lattice
from src.flux import Flux
from src.sampling import Sampling

SEED = 2024
np.random.seed(SEED)

Lx, Ly = 4, 4
K = [0.5, 0.5, 0.5]
beta = 1.0

n_sweeps = 1_00
n_bins = 10

latt = Lattice(Lx, Ly)
flux = Flux(K, latt)

mc = MonteCarlo(n_sweeps, n_bins, beta)
mc.init_observables(1, 0, latt)
mc.set_flux(flux)

mc.simulation()

