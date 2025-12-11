import numpy as np
import sys 

def n_conf(L):
  Nb = 3 * L * L
  return np.power(2, Nb)
  

n = n_conf(int(sys.argv[1]))
print(n, np.log10(n))

