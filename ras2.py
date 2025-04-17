import scipy.integrate
import matplotlib.pyplot as plt
import scipy
import numpy as np
from time import time

# long-term rate/slope of palmitic acid
PALMITIC_RATE = 1.5
# initial time at which palmitic reachs linear behaviour and we start this simulation 
INITIAL_TIME = 20
# intital palmitic acid concentration at 20 s
INITIAL_PALMITIC = 327.06
# inital Rc concentration
INITIAL_Rc = 100

k_palm_max = 1
k_depal_max = 90
Km_palm = 23
Km_depal = 89

# concentrations values match the order in the variables list
variables = ["PalmiticAcid", "R_c", "R_p", "R_r"]
intital = [INITIAL_PALMITIC, INITIAL_Rc, 0, 0]

def k_palm(concentrations):
    return k_palm_max * concentrations[0]/(Km_palm + concentrations[0])

def k_depal(concentrations):
    return k_depal_max * concentrations[0]/(Km_depal + concentrations[0])


def rate(t, concentrations):
    return [
        PALMITIC_RATE + (-k_palm_max * concentrations[0] * concentrations[1] / (Km_palm + concentrations[0])),
        -k_palm(concentrations) * concentrations[1] + k_depal(concentrations) * concentrations[2],
        k_palm(concentrations) * concentrations[1] - k_depal(concentrations) * concentrations[2],
        0 # placeholder for now
    ]


FINAL_TIME = 50
num_iter = FINAL_TIME-INITIAL_TIME
t = np.linspace(INITIAL_TIME, FINAL_TIME, num_iter * 10)


sol = scipy.integrate.solve_ivp(rate, [INITIAL_TIME, FINAL_TIME], intital, method="LSODA", dense_output=True)
z = sol.sol(t)

for idx, metabolite in enumerate(variables):
    datapoints = z[idx]
    
    # Plot the iterative data with markers and a solid line
    plt.plot(t, datapoints, color="purple", linestyle='-', linewidth=1.5)
    
    # Add title and axis labels with LaTeX for µM symbol
    plt.title(f"Comparison for {metabolite}", fontsize=16, fontweight='bold')
    plt.xlabel("Iteration", fontsize=14)
    plt.ylabel(r'Concentration ($\mu$M)', fontsize=14)

    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(f"rasgraphs/{metabolite}", dpi=600)
    plt.clf()

