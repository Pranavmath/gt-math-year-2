import matplotlib.pyplot as plt
import numpy as np
import cv2
import imageio
from tqdm import tqdm

# long-term rate/slope of palmitic acid
PALMITIC_RATE = 1.5
# intital palmitic acid concentration at 20 s
INITIAL_PALMITIC = 327.06


DELTA_T = 0.1
D_p = 1
AREA_SQUARE = 1
k_depal_max = 0.015
Km_depal = 89

# these don't matter when num_rafts=0
KON_BASE = 100
KOFF = 10.75
SIGMA = 3
HILL = 1.5
N = 5

# initial time at which palmitic reachs linear behaviour and we start this simulation 
INITIAL_TIME = 0
FINAL_TIME = 2000
times = np.arange(INITIAL_TIME, FINAL_TIME, DELTA_T)

# grid side size
GRID_SIZE = 200

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
num_rafts = 1000

# list of (xi, yi)
rafts = np.random.choice(len(outside_indices), num_rafts, replace=False)

for idx in rafts:
    i, j = outside_indices[idx]
    colors[i, j] = LIGHT_BLUE

rafts = np.array([outside_indices[idx] for idx in rafts])

"""
These are the variables we are tracking over time
"""
# 2d array
Rp = np.zeros((GRID_SIZE, GRID_SIZE))
# rafts R
Rrs = np.zeros(num_rafts)
# Rc
Rc = 100


def compute_laplacian(Rp):
    laplacian = np.zeros_like(Rp)
    laplacian[1:-1, 1:-1] = (
        Rp[1:-1, 2:] + Rp[1:-1, :-2] +
        Rp[2:, 1:-1] + Rp[:-2, 1:-1] -
        4 * Rp[1:-1, 1:-1]
    )
    return laplacian / AREA_SQUARE


def kon(raft_idx):
    xi, yi = rafts[raft_idx]
    dist_squared = xi**2 + yi**2
    return KON_BASE * np.exp(-dist_squared/(SIGMA**2))


def k_depal(palmitic_acid):
    return k_depal_max * palmitic_acid/(Km_depal + palmitic_acid)


print("Simulation Started")
frames = []

for t in tqdm(times):
    Rp_normalized = np.uint8(255 * (Rp - np.min(Rp)) / (np.max(Rp) - np.min(Rp) + 1e-8))
    colored = cv2.applyColorMap(Rp_normalized, cv2.COLORMAP_INFERNO)
    colored = cv2.cvtColor(colored, cv2.COLOR_BGR2RGB)
    frames.append(colored)

    # ----------------------------------------------------

    if t < 100:
        Rp[center_grid, center_grid] = 1000

    palmitic_acid = INITIAL_PALMITIC #+ (t - INITIAL_TIME)/60 * PALMITIC_RATE 
    
    Rp_change = D_p * compute_laplacian(Rp) - k_depal(palmitic_acid) * Rp


    Rr_change = np.zeros(num_rafts)

    for raft_idx, (xi, yi) in enumerate(rafts):
        total = (-kon(raft_idx) * Rp[yi, xi]**HILL * (1 - Rrs[raft_idx]/N)) + (KOFF * Rrs[raft_idx])
        # clip so that its not too positive that it makes Rr negative and too negative that it makes Rp negative
        total = np.clip(total, -Rp[yi, xi]/DELTA_T, Rrs[raft_idx]/DELTA_T)

        Rp_change[yi, xi] += total
        Rr_change[raft_idx] = -total

    # clipping rp for conservation of matter
    #"""
    previous_sum = np.sum(Rp_change)
    Rp_change = np.maximum(Rp_change, -Rp/DELTA_T)
    new_sum = np.sum(Rp_change)
    assert new_sum >= previous_sum
    #"""

    Rc_change = k_depal(palmitic_acid) * np.sum(Rp) * AREA_SQUARE - (new_sum-previous_sum)

    Rp += Rp_change * DELTA_T
    Rrs += Rr_change * DELTA_T
    Rc += Rc_change * DELTA_T

    #Rp, Rrs, Rc = np.clip(Rp, 0, None), np.clip(Rrs, 0, None), min(0, Rc)

print(len(frames))

# 100x speed
out = cv2.VideoWriter("heatmap.mp4", cv2.VideoWriter_fourcc(*'mp4v'), 100/DELTA_T, (GRID_SIZE, GRID_SIZE))
for frame in frames:
    out.write(cv2.cvtColor(frame, cv2.COLOR_RGB2BGR))
out.release()
