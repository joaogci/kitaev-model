import numpy as np
import os
import sys

class Analysis():
  def __init__(self, n_rebin: int):
    self.n_rebin = n_rebin
    
  def analyse(self):
    files = os.listdir(".")
    
    for file in files: 
      if file.endswith("_info"):
        
        if file.find("scal") != -1:
          self.__analyse_scal(file[:-5])
        
        if file.find("eq") != -1: 
          self.__analyse_eq(file[:-5])
  
  def __analyse_scal(self, filename: str):
    raw_data = np.loadtxt(filename)
    
    data = np.zeros(raw_data.shape[0] // self.n_rebin, dtype=complex)
    for i in range(data.shape[0]):
      for j in range(self.n_rebin):
        data[i] += complex(raw_data[i * self.n_rebin + j, 0], raw_data[i * self.n_rebin + j, 1])
      data[i] = data[i] / self.n_rebin
    
    mean = np.mean(data)
    std = np.std(data)
    
    with open(filename + "J", "w") as f:
      f.write(f"{mean.real} {std.real} {mean.imag} {std.imag} \n")
    
  def __analyse_eq(self, filename: str):
    ...

