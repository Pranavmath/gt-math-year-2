import scipy.integrate
import matplotlib.pyplot as plt
import scipy
import numpy as np

# palmitoyl_coa concentraion (over time) array
palmitoyl_coa = np.load("palmitoyl-coa.npy")

# time array
t = np.load("time.npy")

P_total = 100
kcat_pal = 1.4
Km_pal = 23

kcat_depal = 1.4
Km_depal = 23

E_pal = 1
E_depal = 1

Km_C = 10
E_depal = 1


print(palmitoyl_coa)
print(t)