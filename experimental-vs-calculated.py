import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import glob
import re

# Function to compute relative free energy from 1D FES
def fes_rel(fes_list, fes_all=None, kT=2.49):
    results = []
    for i in fes_list:
        temp = sum(np.exp(-j / kT) for j in i)
        results.append(temp)
    total = sum(fes_all) if fes_all else sum(results)
    deltaGs = [-kT * np.log(i / total) if i > 0 else np.inf for i in results]
    min_res = min(deltaGs)
    return [i - min_res for i in deltaGs]

# Extract numerical suffix from .dat files
def extract_dat_number(fname):
    match = re.search(r'dat(\d+)\.dat$', fname)
    return int(match.group(1)) if match else -1

# Experimental ΔG values (kcal/mol)
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

systems = list(exp_dict.keys())
n_replicas = 3
base_dir = '/scratch/project_465001666/elbrunis/ptp/analysis/fes_outputs'
funnel_correction = 2.45

# Plot setup
plt.figure(figsize=(8, 6))
exp_list = []
est_list = []
stdiv_list = []
error_list = []

for system in systems:
    estimates = []
    for rep in range(1, n_replicas + 1):
        pattern = os.path.join(base_dir, f'HILLS-{system}-{rep}_fes_1d.dat*.dat')
        files = sorted(glob.glob(pattern), key=extract_dat_number)
        if not files:
            continue

        fes_file = files[-1]  # Use the most recent datXX.dat file (e.g., dat30.dat)
        fes_data = np.loadtxt(fes_file)
        fes_1d = fes_rel(fes_data)
        estimate = -(np.median(fes_1d[-65:-15]) - funnel_correction)
        estimates.append(estimate)

    if estimates:
        mean_est = np.mean(estimates)
        stdev = np.std(estimates)
        exp_val = exp_dict[system]

        exp_list.append(exp_val)
        est_list.append(mean_est)
        stdiv_list.append(stdev)
        error_list.append(abs(exp_val - mean_est))

        print(f'{system}: ΔG = {mean_est:.2f} ± {stdev:.2f} kcal/mol | Error: {abs(exp_val - mean_est):.2f} kcal/mol')
    else:
        print(f'[WARNING] Missing FES data for {system}, skipping.')

# Plotting
plt.plot([-17, 1], [-17, 1], 'k--', label='y = x')
# ±1 kcal/mol band (light grey with higher alpha)
plt.fill_between([-17, 1], [-18, 0], [-16, 2],
                 facecolor='lightgrey', edgecolor=None, alpha=0.8, label='±1 kcal/mol')

# ±2 kcal/mol band (darker grey with lower alpha)
plt.fill_between([-17, 1], [-19, -1], [-15, 3],
                 facecolor='darkgrey', edgecolor=None, alpha=0.4, label='±2 kcal/mol')


plt.errorbar(exp_list, est_list, yerr=stdiv_list, fmt='o', color='black', capsize=5)
plt.xlabel('Experimental ΔG (kcal/mol)')
plt.ylabel('Estimated ΔG (kcal/mol)')
plt.title('Experimental vs Estimated ΔG')
plt.xlim(-17, 1)
plt.ylim(-17, 1)

# Pearson, Kendall, MUE
plt.text(-15, -5, f'Pearson r = {stats.pearsonr(exp_list, est_list)[0]:.2f}')
plt.text(-15, -6, f'Kendall τ = {stats.kendalltau(exp_list, est_list)[0]:.2f}')
plt.text(-15, -7, f'MUE = {np.mean(error_list):.2f} kcal/mol')
plt.legend()

plt.tight_layout()
plt.savefig("exp_vs_est_dG_fes_latest.png")

