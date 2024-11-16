import h5py
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit

# Define the exponential function for fitting
def exponential(x, a, b, c, d):
    return a * np.exp(-(x - b) / c) + d

def optimal_bins(data):
    n = len(data)
    if n == 0:
        return {}
    # Sturges' Rule
    sturges_bins = int(np.ceil(np.log2(n) + 1))
    return sturges_bins

# Constants and Initialization
t1_vals = {i: [] for i in range(6)}
t1_errs = {i: [] for i in range(6)}
qubit_for_this_index = []
rounds = []
reps = []
file_names = []
dates = {i: [] for i in range(6)}
mean_values = {}
dates_run = ['2024-11-14']
show_legends = False

# Grab all of the data
for date in dates_run:
    folder_path = f"/data/QICK_data/6transmon_run4a/{date}/T1_ge/"

    for file_name in os.listdir(folder_path):
        if file_name.endswith(".h5"):
            file_path = os.path.join(folder_path, file_name)

            with h5py.File(file_path, 'r') as h5_file:

                for qubit_index in range(1, 7):  # Adjusting for 0-based index
                    group_name = f"Q{qubit_index}"
                    if group_name in h5_file:
                        # Extracting the datasets
                        t1_data = h5_file[f"{group_name}/T1"][()]
                        t1_err_data = h5_file[f"{group_name}/Errors"][()]
                        dates_data = h5_file[f"{group_name}/Dates"][()]

                        # Appending values to the appropriate lists
                        t1_vals[qubit_index - 1].extend(t1_data)  # Store T1 values
                        t1_errs[qubit_index - 1].extend(t1_err_data)  # Store T1 error values
                        dates[qubit_index - 1].extend([d.decode('utf-8') for d in dates_data])  # Decode bytes to string

                        # You can also append qubit indices if needed
                        qubit_for_this_index.extend([qubit_index] * len(t1_data))

fig, axes = plt.subplots(2, 3, figsize=(12, 8))
axes = axes.flatten()
font = 14
titles = [f"Qubit {i+1}" for i in range(6)]
gaussian_xvals =  {i: [] for i in range(0, 6)}
gaussian_yvals =  {i: [] for i in range(0, 6)}
gaussian_colors = {i: [] for i in range(0, 6)}
gaussian_dates = {i: [] for i in range(0, 6)}
colors = ['orange','blue','purple','green','brown','pink']
for i, ax in enumerate(axes):
    if len(dates[i])>1:
        date_label = dates[i][0]
    else:
        date_label = ''


    if len(t1_vals[i]) >1:
        optimal_bin_num = optimal_bins(t1_vals[i])

        # Fit a Gaussian to the raw data instead of the histogram
        # get the mean and standard deviation of the data
        mu_1, std_1 = norm.fit(t1_vals[i])
        mean_values[f"Qubit {i + 1}"] = mu_1  # Store the mean value for each qubit

        # Generate x values for plotting a gaussian based on this mean and standard deviation
        x_1 = np.linspace(min(t1_vals[i]), max(t1_vals[i]), optimal_bin_num)
        p_1 = norm.pdf(x_1, mu_1, std_1)

        # Calculate histogram data for t1_vals[i]
        hist_data_1, bins_1 = np.histogram(t1_vals[i], bins=optimal_bin_num)
        bin_centers_1 = (bins_1[:-1] + bins_1[1:]) / 2

        # Scale the Gaussian curve to match the histogram
        # the gaussian height natrually doesnt match the bin heights in the histograms
        # np.diff(bins_1)  calculates the width of each bin by taking the difference between bin edges
        # the total counts are in hist_data_1.sum()
        # to scale, multiply data gaussian by bin width to convert the probability density to probability within each bin
        # then multiply by the total count to scale the probability to match the overall number of datapoints
        # https://mathematica.stackexchange.com/questions/262314/fit-function-to-histogram
        # https://stackoverflow.com/questions/23447262/fitting-a-gaussian-to-a-histogram-with-matplotlib-and-numpy-wrong-y-scaling
        ax.plot(x_1, p_1 * (np.diff(bins_1) * hist_data_1.sum()), 'b--', linewidth=2, color=colors[i])

        # Plot histogram and Gaussian fit for t1_vals[i]
        ax.hist(t1_vals[i], bins=optimal_bin_num, alpha=0.7,color=colors[i], edgecolor='black', label=date_label)

        #make a fuller gaussian to make smoother lotting for cumulative plot
        x_1_full = np.linspace(min(t1_vals[i]), max(t1_vals[i]), 2000)
        p_1_full = norm.pdf(x_1_full, mu_1, std_1)

        gaussian_xvals[i].append(x_1_full)
        gaussian_yvals[i].append(p_1_full )
        gaussian_colors[i].append(colors[i])
        gaussian_dates[i].append(date_label)

        #rough start at errors:
        #counts, bin_edges, _ = ax.hist(t1_vals[i], bins=20, alpha=0.7, color='blue', edgecolor='black')
        #bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
        #bin_errors = [np.sqrt(np.sum(t1_errs[i])) for _ in range(len(bin_centers))]  #using error propogation through the sum of valuse in each bin
        #ax.errorbar(bin_centers, counts, yerr=bin_errors, fmt='o', color='red', ecolor='black', capsize=3, linestyle='None')
        if show_legends:
            ax.legend()
        ax.set_title(titles[i] + f" $\mu$: {mu_1:.2f} $\sigma$:{std_1:.2f}",fontsize = font)
        ax.set_xlabel('T1 (µs)',fontsize = font)
        ax.set_ylabel('Frequency',fontsize = font)
        ax.tick_params(axis='both', which='major', labelsize=font)

plt.tight_layout()
plt.savefig('plots/hists.png', transparent=True, dpi=500)

fig, ax = plt.subplots(1, 1, figsize=(12, 8))
plt.title('Cumulative Distribution',fontsize = font)
for i in range(0, len(t1_vals)):
    if len(dates[i])>1:
        date_label = dates[i][0]
    else:
        date_label = ''

    if len(t1_vals[i]) > 1:
        t1_vals_sorted = np.sort(t1_vals[i])
        len_samples = len(t1_vals_sorted)
        var = np.linspace(1,len_samples,len_samples)/ len_samples

        cumulative_gaussian = np.cumsum(gaussian_yvals[i][0]) / np.sum(gaussian_yvals[i][0])
        ax.scatter(t1_vals_sorted,var,color = colors[i], label = f'Q{i+1}', s = 5)
        ax.plot(gaussian_xvals[i][0], cumulative_gaussian, color=colors[i], label='Gauss Fit ' + f'Q{i + 1}',
                linestyle='--')
        ax.tick_params(axis='both', which='major', labelsize=font)
#ax.set_title('')
ax.set_xlabel('T1 (us)',fontsize = font)
ax.set_ylabel('Cumulative Distribution',fontsize = font)
ax.loglog()
ax.legend(edgecolor='black')
ax.set_xlim(10**0, 10**3)
ax.set_ylim(10 ** -7, 10 ** 0) #to compare to johns plot, need to adjust a little
plt.tight_layout()
plt.savefig('plots/cumulative.png', transparent=True, dpi=500)



fig, axes = plt.subplots(2, 3, figsize=(12, 8))
plt.title('Fit Error vs T1 Time',fontsize = font)
axes = axes.flatten()
titles = [f"Qubit {i + 1}" for i in range(6)]
for i, ax in enumerate(axes):

    if len(dates[i])>1:
        date_label = dates[i][0]
    else:
        date_label = ''
    ax.set_title(titles[i], fontsize = font)
    ax.scatter(t1_vals[i],t1_errs[i], label = date_label, color = colors[i])
    if show_legends:
        ax.legend(edgecolor='black')
    ax.set_xlabel('T1 (us)', fontsize = font)
    ax.set_ylabel('Fit error (us)', fontsize = font)
    ax.tick_params(axis='both', which='major', labelsize=font)
plt.tight_layout()
plt.savefig('plots/errs.png', transparent=True, dpi=500)

plt.show()
