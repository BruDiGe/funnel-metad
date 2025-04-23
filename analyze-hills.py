import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import subprocess as sp
import glob
import shutil
import scipy.stats
from collections import deque
from math import exp, log

# Constants
kT = 0.593  # kcal/mol at 298 K
funnel_correction = 0.0

# Directory with renamed HILLS files (adapt if needed)
hills_dir = os.getcwd()
#hills_dir = '/scratch/project_465001666/elbrunis/ptp/all_HILLS'

# Explicit list of HILLS files
hills_files = [
    "HILLS-4erc-D65A-1", "HILLS-4erc-D65A-2", "HILLS-4erc-D65A-3",
    "HILLS-4erc-D65N-1", "HILLS-4erc-D65N-2", "HILLS-4erc-D65N-3",
    "HILLS-4erc-E134A-1", "HILLS-4erc-E134A-2", "HILLS-4erc-E134A-3",
    "HILLS-4erc-E134Q-1", "HILLS-4erc-E134Q-2", "HILLS-4erc-E134Q-3",
    "HILLS-4erc-P64A-1", "HILLS-4erc-P64A-2", "HILLS-4erc-P64A-3",
    "HILLS-4erc-P68V-1", "HILLS-4erc-P68V-2", "HILLS-4erc-P68V-3",
    "HILLS-4erc-P69A-1", "HILLS-4erc-P69A-2", "HILLS-4erc-P69A-3",
    "HILLS-4erc-Q138A-1", "HILLS-4erc-Q138A-2", "HILLS-4erc-Q138A-3"
]

# Output directory
fes_output = os.path.join(hills_dir, 'fes_outputs')
os.makedirs(fes_output, exist_ok=True)

# Function to calculate relative FES
def fes_rel(fes_list, fes_all=None, kT=0.593):
    results = []
    for i in fes_list:
        temp = sum(exp(-j / kT) for j in i)
        results.append(temp)
    total = sum(fes_all if fes_all else results)
    results_2 = [-kT * log(i / total) if i > 0 else np.inf for i in results]
    min_res = min(results_2)
    return [i - min_res for i in results_2]

# Process each HILLS file
for hills_file in hills_files:
    print(f"Processing {hills_file}...")

    hills_path = os.path.join(hills_dir, hills_file)
    temp_hills = os.path.join(hills_dir, 'HILLS_temp')
    shutil.copyfile(hills_path, temp_hills)

    out_1d = os.path.join(fes_output, f'{hills_file}_fes_1d.dat')
    out_2d = os.path.join(fes_output, f'{hills_file}_fes_2d.dat')

    # Generate 1D FES
    sp.call(f'plumed sum_hills --hills {temp_hills} --idw pp.proj --kt {kT} --outfile {out_1d} --stride 5000', shell=True)
    # Generate 2D FES
    sp.call(f'plumed sum_hills --hills {temp_hills} --kt {kT} --outfile {out_2d} --stride 5000', shell=True)

    # Plot pp.proj vs time
    try:
        df = pd.read_csv(temp_hills, delim_whitespace=True, comment='#')
        df.columns = ['time', 'pp.proj', 'pp.ext', 'sigma_proj', 'sigma_ext', 'height', 'biasf']
        df = df.apply(pd.to_numeric)

        plt.figure()
        plt.plot(df['time'] / 1000, df['pp.proj'])
        plt.xlabel('Time (ns)')
        plt.ylabel('pp.proj (Å)')
        plt.title(hills_file)
        plt.savefig(os.path.join(fes_output, f'{hills_file}_proj_vs_time.png'))
        plt.close()
    except Exception as e:
        print(f"Warning: could not plot pp.proj for {hills_file} - {e}")

    os.remove(temp_hills)

print("All HILLS files processed. FES data saved in:", fes_output)
