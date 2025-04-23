import os
import re
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from collections import deque

# Function to compute relative free energy from 1D FES
def fes_rel(fes_data, kT=2.49):
    integrals = np.array([np.exp(-val / kT) for val in fes_data[:, 1]])
    total = np.sum(integrals)
    deltaGs = -kT * np.log(integrals / total)
    return deltaGs - np.min(deltaGs)

# Define input directory and systems
input_dir = '/scratch/project_465001666/elbrunis/ptp/analysis/fes_outputs'
systems = ['4erc-D65A', '4erc-D65N', '4erc-E134A', '4erc-E134Q',
           '4erc-P64A', '4erc-P68V', '4erc-P69A', '4erc-Q138A']
reps = [1, 2, 3]
stride = 10  # Each FES file corresponds to 10 ns
funnel_correction = 2.45
start_from_idx = 5

# Experimental ΔG values in kcal/mol
exp_dict = {
    '4erc-D65A': -4.99,
    '4erc-D65N': -2.81,
    '4erc-P64A': -2.02,
    '4erc-P69A': -2.06,
    '4erc-P68V': -3.04,
    '4erc-E134Q': -2.88,
    '4erc-E134A': -2.72,
    '4erc-Q138A': -2.57,
}

# Plotting setup
# Plotting setup: 2 rows × 4 columns for 8 systems
fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(20, 10))
axes = axes.flatten()  # Make it easier to index in a loop

for plot_index, system in enumerate(systems):
    ax = axes[plot_index]
    ax.set_title(system)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('ΔG (kcal/mol)')
    ax.set_ylim(-20, 0)  # Set fixed y-axis range
    estimates_lists = deque()

    for rep in reps:
        rep_tag = f"HILLS-{system}-{rep}_fes_1d.dat"
        rep_files = sorted(
            [f for f in os.listdir(input_dir) if f.startswith(rep_tag)],
            key=lambda x: int(re.findall(r'\\d+', x)[-1]) if re.findall(r'\\d+', x) else 0
        )

        ns = []
        estimates = []
        oscillations = []

        for idx, file in enumerate(rep_files):
            fes_data = np.loadtxt(os.path.join(input_dir, file))
            fes_1d = fes_rel(fes_data)
            estimate = -(np.median(fes_1d[-65:-15]) - funnel_correction)

            ns.append((idx + 1) * stride)
            estimates.append(round(estimate, 2))

            if idx < start_from_idx:
                oscillations.append(0)
            else:
                if idx < 10:
                    osc = scipy.stats.sem(estimates)
                else:
                    osc = scipy.stats.sem(estimates[idx - 10:])
                oscillations.append(osc * 2)

        ax.errorbar(ns, estimates, yerr=oscillations, label=f'Rep {rep}')
        estimates_lists.append(estimates)

    # Mean estimate across replicas
    shortest = min(len(est) for est in estimates_lists)
    mean_estimates = [round(np.mean([est[i] for est in estimates_lists]), 2) for i in range(shortest)]

    ax.plot(ns[:shortest], mean_estimates, color='black', linewidth=2, label='Mean')
    ax.axhline(y=exp_dict[system], color='red', linestyle='--', label='Experimental ΔG')
    ax.legend()

plt.tight_layout()
plt.savefig("dg_vs_time_plot.png")
