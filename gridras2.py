import matplotlib.pyplot as plt
import numpy as np

# long-term rate/slope of palmitic acid
PALMITIC_RATE = 1.5
# intital palmitic acid concentration at 20 s
INITIAL_PALMITIC = 327.06


DELTA_T = 0.1
D_p = 0.1
KOFF = 
AREA_SQUARE = 0.01
SIGMA = 3
KON_BASE = 10**6
HILL = 1.5
N = 5
k_depal_max = 90
Km_depal = 89

# initial time at which palmitic reachs linear behaviour and we start this simulation 
INITIAL_TIME = 20
FINAL_TIME = 50
times = np.arange(INITIAL_TIME, FINAL_TIME, DELTA_T)

# grid side size
GRID_SIZE = 120

"""
Grid code now
"""

# colors
WHITE = [1.0, 1.0, 1.0]
LIGHT_RED = [0.9, 0.4, 0.4]
LIGHT_BLUE = [0.3, 0.5, 0.8]

# grid is white initially
colors = np.tile(WHITE, (GRID_SIZE, GRID_SIZE, 1))

center_grid = GRID_SIZE//2
# 1x1 golgi zone
golgi = (center_grid, center_grid)
colors[golgi] = LIGHT_RED

# indices other than golgi zone
outside_indices = [(i, j) for i in range(GRID_SIZE) for j in range(GRID_SIZE) if not (i == center_grid and j == center_grid)]

# choose a random # for lipid rafts (30-1800ish depending on density)
num_rafts = 1200

# list of (xi, yi)
rafts = np.random.choice(len(outside_indices), num_rafts, replace=False)

for idx in rafts:
    i, j = outside_indices[idx]
    colors[i, j] = LIGHT_BLUE


"""
These are the variables we are tracking over time
"""
# 2d array
Rp = np.zeros((GRID_SIZE, GRID_SIZE))
# rafts R
Rrs = np.zeros(num_rafts)
# Rc
Rc = 100


def laplacian_rp(u, v):
    
    total = -4 * Rp[v, u]

    for (i, j) in [(u+1, v), (u-1, v), (u, v+1), (u, v-1)]:
        if (0 <= i < GRID_SIZE) and (0 <= j < GRID_SIZE):
            total += Rp[j, i]
    
    return total / AREA_SQUARE

def kon(raft_idx):
    xi, yi = rafts[raft_idx]
    dist_squared = xi**2 + yi**2
    return KON_BASE * np.exp(-dist_squared/(SIGMA**2))


def k_depal(palmitic_acid):
    return k_depal_max * palmitic_acid/(Km_depal + palmitic_acid)


print("Simulation Started")

for t in times:
    # write injection code here??


    palmitic_acid = (t - INITIAL_TIME) * PALMITIC_RATE + INITIAL_PALMITIC

    Rp_change = np.array([
        [(D_p * laplacian_rp(u, v) - k_depal(palmitic_acid) * Rp[v, u]) for v in range(GRID_SIZE)]
        for u in range(GRID_SIZE)
    ])

    for raft_idx, (xi, yi) in enumerate(rafts):
        Rp_change += (-kon(raft_idx) * Rp[yi, xi]**HILL * (1 - Rrs[raft_idx]/N))
        Rp_change += (KOFF * Rrs[raft_idx])

    Rr_change = np.array([kon(raft_idx) * Rp[yi, xi]**HILL * (1 - (Rrs[raft_idx]/N)) - KOFF * Rrs[raft_idx] for raft_idx, (xi, yi) in enumerate(rafts)])

    Rc_change = k_depal(palmitic_acid) * np.sum(Rp) * AREA_SQUARE

    Rp += Rp_change * DELTA_T
    Rrs += Rr_change * DELTA_T
    Rc += Rc_change * DELTA_T


