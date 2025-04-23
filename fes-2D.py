import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import glob
import re

# Directory and systems
input_dir = '/scratch/project_465001666/elbrunis/ptp/analysis/fes_outputs'
systems = ['4erc-D65A', '4erc-D65N', '4erc-E134A', '4erc-E134Q',
           '4erc-P64A', '4erc-P68V', '4erc-P69A', '4erc-Q138A']
reps = [1, 2, 3]

# Plot setup: 8 systems × 3 replicas = 24 subplots in 8×3 grid
fig, axes = plt.subplots(nrows=8, ncols=3, figsize=(18, 28), constrained_layout=True)
cmap = 'RdBu'

# Helper to extract number from .datXX.dat
def extract_dat_number(f):
    match = re.search(r'(\d+)\.dat$', f)
    return int(match.group(1)) if match else -1

# Loop over systems and replicas
for i_sys, system in enumerate(systems):
    for i_rep, rep in enumerate(reps):
        ax = axes[i_sys, i_rep]

        # Find all matching files like HILLS-4erc-D65A-1_fes_2d.dat30.dat
        pattern = os.path.join(input_dir, f'HILLS-{system}-{rep}_fes_2d.dat*.dat')
        matching_files = sorted(glob.glob(pattern), key=extract_dat_number)

        if not matching_files:
            ax.set_visible(False)
            continue

        try:
            file_path = matching_files[-1]
            fes = np.loadtxt(file_path)
            x, y, z = fes[:, 0], fes[:, 1], fes[:, 2] / 4.168  # Convert to kcal/mol

            N = 200
            xi = np.linspace(x.min(), x.max(), N)
            yi = np.linspace(y.min(), y.max(), N)
            zi = scipy.interpolate.griddata((x, y), (z - z.min()), (xi[None, :], yi[:, None]), method='cubic')

            cf = ax.contourf(xi, yi, zi, levels=np.linspace(0, 10, 21), cmap=cmap, vmin=0, vmax=10)
            cbar = fig.colorbar(cf, ax=ax)
            cbar.set_label('Free Energy (kcal/mol)')

            ax.set_title(f'{system} Rep {rep}')
            ax.set_xlabel('pp.proj (nm)')
            ax.set_ylabel('pp.ext (nm)')
        except Exception as e:
            print(f"Could not process {file_path}: {e}")
            ax.set_visible(False)

# Save plot
plt.savefig("fes_2d_all_replicas_individual_cbars.png")

