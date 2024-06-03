import numpy as np
from .lattice import Lattice

class Obs():
  def __init__(self, name: str):
    self.name = name
    self.N = 0
  
  def write_to_file(self):
    pass
  
  def reset(self):
    pass

class ObsScalar(Obs):
  def __init__(self, name: str):
    super().__init__(name)
    self.obs_vec = complex(0.0)
     
  def write_to_file(self):
    self.obs_vec = self.obs_vec / self.N
    
    with open(self.name + "_scal", "a") as f:
      f.write(f"{self.obs_vec.real} {self.obs_vec.imag}\n")
      
    with open(self.name + "_scal_info", "w") as f:
      f.write("-- Analysis Mode --\n")
      f.write("identity\n")
  
  def reset(self):
    self.N = 0
    self.obs_vec = complex(0.0)

class ObsLatt(Obs):
  pass

