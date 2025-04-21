import matplotlib.pyplot as plt
import numpy as np
import cv2
import imageio
from tqdm import tqdm
from scipy.ndimage import convolve

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
FINAL_TIME = 4000
times = np.arange(INITIAL_TIME, FINAL_TIME, DELTA_T)

# grid side size
GRID_SIZE = 100

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
num_rafts = 1300

# list of (xi, yi)
rafts = np.random.choice(len(outside_indices), num_rafts, replace=False)

for idx in rafts:
    i, j = outside_indices[idx]
    colors[i, j] = LIGHT_BLUE

# get raft coords and kons for each raft
rafts = np.array([outside_indices[idx] for idx in rafts])
xi_vals, yi_vals = rafts[:, 0], rafts[:, 1]
dist_squared =  100 + (xi_vals-center_grid)**2 + (yi_vals-center_grid)**2
kon_vals = KON_BASE * np.exp(-dist_squared / (SIGMA**2))

# mask where we don't do (- k_depal(palmitic_acid) * Rp)
no_depal_mask = np.zeros((GRID_SIZE, GRID_SIZE))
no_depal_mask[yi_vals, xi_vals] = 1
# dont do at golgi
no_depal_mask[center_grid, center_grid] = 1


"""
These are the variables we are tracking over time
"""
# 2d array
Rp = np.zeros((GRID_SIZE, GRID_SIZE))
# rafts R
Rrs = np.zeros(num_rafts)
# Rc
Rc = 100

Rcs = []
Rr_sums = []
Rr_heatmap = np.zeros((GRID_SIZE, GRID_SIZE))

def get_random_kernel():
    kernel = np.random.normal(loc=0, scale=2, size=(3, 3))
    kernel[1, 1] = 0
    # make sure its always going outwards
    kernel = np.sign(np.sum(kernel)) * kernel
    assert np.sum(kernel) >= 0
    kernel[1, 1] = -np.sum(kernel)
    return kernel

def compute_laplacian(Rp, laplacian_kernel):
    laplacian = convolve(Rp, laplacian_kernel, mode='constant', cval=0.0)
    return laplacian / AREA_SQUARE


def k_depal(palmitic_acid):
    return k_depal_max * palmitic_acid/(Km_depal + palmitic_acid)


print("Simulation Started")
frames = []
rr_frames = []

start_time = 0

change_sum_before = []
change_sum_after = []

def get_colored(arr):
    arr_normalized = np.uint8(255 * (arr - np.min(arr)) / (np.max(arr) - np.min(arr) + 1e-8))
    colored = cv2.applyColorMap(arr_normalized, cv2.COLORMAP_INFERNO)
    colored = cv2.cvtColor(colored, cv2.COLOR_BGR2RGB)
    return colored

Rp_total = []

for t in tqdm(times):

    #"""
    frames.append(get_colored(Rp))
    rr_frames.append(get_colored(Rr_heatmap))
    #"""

    # ----------------------------------------------------
    Rcs.append(Rc)
    Rr_sums.append(np.sum(Rrs))
    Rp_total.append(np.sum(Rp))
    
    """
    if (Rp[center_grid, center_grid] < 1):
        start_time = t
    
    if (start_time <= t <= start_time + 5):
        Rp[center_grid, center_grid] += 100
    """
    if (abs(t % 200) <= 0.5):
        Rp[center_grid, center_grid] = 100
    

    palmitic_acid = INITIAL_PALMITIC #+ (t - INITIAL_TIME)/60 * PALMITIC_RATE 

    #"""
    kernel = np.array([
        [0,  0.1, 0],
        [0.1, -0.4, 0.1],
        [0,  0.1, 0]
    ])
    #"""
    #kernel = get_random_kernel()
    #print(kernel)
    
    laplacian_rp = compute_laplacian(Rp, kernel)
    Rp_change = D_p * laplacian_rp  - (1-no_depal_mask) * k_depal(palmitic_acid) * Rp * 0.05

    """
    if this is commented then sinusodial if not then not sinusodial
    """
    Rp_change[center_grid, center_grid] = 0

    change_sum_before.append(np.sum(Rp_change))

    # calculate the changes for Rp and Rr
    Rp_vals = Rp[yi_vals, xi_vals]

    # can go to nan if Rp_vals is negative => complex
    binding = -kon_vals * Rp_vals**HILL * (1 - Rrs / N)
    unbinding = KOFF * Rrs
    totals = binding + unbinding
    """
    #clip for conservation of matter
    max_neg = -Rp_vals / DELTA_T
    max_pos = Rrs / DELTA_T
    totals = np.clip(totals, max_neg, max_pos)
    """
    # Apply changes
    np.add.at(Rp_change, (yi_vals, xi_vals), totals)
    Rr_change = -totals


    change_sum_after.append(np.sum(Rp_change))

    """
    # clipping rp_change for conservation of matter
    previous_sum = np.sum(Rp_change)
    Rp_change = np.maximum(Rp_change, -Rp/DELTA_T)
    new_sum = np.sum(Rp_change)
    assert new_sum >= previous_sum
    """

    Rc_change = k_depal(palmitic_acid) * np.sum(Rp * (1-no_depal_mask)) * AREA_SQUARE #- (new_sum-previous_sum)
    

    Rp += Rp_change * DELTA_T
    Rrs += Rr_change * DELTA_T
    Rc += Rc_change * DELTA_T

    Rp = np.clip(Rp, 0, a_max=None)
    #print(np.min(Rrs), np.max(Rrs))
    Rr_heatmap[yi_vals, xi_vals] = np.where(Rrs > 0, np.log(np.clip(Rrs, 1e-10, None)), 0)




def save_graph(a, b, c):
    plt.plot(times, a, color="purple", linestyle='-', linewidth=1.5, label="Rc Total")
    plt.plot(times, b, color="blue", linestyle='-', linewidth=1.5, label = "Rp Total")
    plt.plot(times, c, color="black", linestyle='-', linewidth=1.5, label="Rr Total")

    
    # Add title and axis labels with LaTeX for µM symbol
    plt.title(f"Comparison for Rp, Rr, and Rc", fontsize=16, fontweight='bold')
    plt.xlabel("Iteration", fontsize=14)
    plt.ylabel(r'Log Concentration ($\mu$M)', fontsize=14)
    plt.legend()

    # Enable grid for better readability of the plot
    plt.grid(True, linestyle='--', alpha=0.7)

    # Improve the layout to avoid overlap
    plt.tight_layout()

    # Save the plot with higher resolution (300 dpi)
    plt.savefig(f"gridrasgraphs/combined.jpg", dpi=600)
    
    # Clear figure for the next plot
    plt.clf()


#plt.plot(times, change_sum_before, label="before")
#plt.plot(times, change_sum_after, label="after")
#plt.legend()
#plt.show()


print(len(frames))

# 1000x speed
out = cv2.VideoWriter("heatmap.mp4", cv2.VideoWriter_fourcc(*'mp4v'), 1000/DELTA_T, (GRID_SIZE, GRID_SIZE))
for frame in frames:
    out.write(cv2.cvtColor(frame, cv2.COLOR_RGB2BGR))
out.release()

out = cv2.VideoWriter("rrheatmap.mp4", cv2.VideoWriter_fourcc(*'mp4v'), 10, (GRID_SIZE, GRID_SIZE))
for frame in rr_frames[:1000]:
    out.write(cv2.cvtColor(frame, cv2.COLOR_RGB2BGR))
out.release()


save_graph(np.log(Rcs), np.log(Rp_total), np.log(Rr_sums))