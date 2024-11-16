import datetime

import h5py
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.stats import norm
from scipy import optimize
from sklearn import preprocessing

def t2_fit(x_data, y_data, verbose=False, guess=None, plot=False):
    # Normalizing the vectors
    xn = preprocessing.normalize([x_data], return_norm=True)
    yn = preprocessing.normalize([y_data], return_norm=True)
    x = xn[0][0]
    y = yn[0][0]
    x_normal = xn[1][0]
    y_normal = yn[1][0]

    # Compute the FFT for guessing the frequency
    fft = np.fft.fft(y)
    f = np.fft.fftfreq(len(x))
    # Take the positive part only
    fft = fft[1: len(f) // 2]
    f = f[1: len(f) // 2]
    # Remove the DC peak if there is one
    if (np.abs(fft)[1:] - np.abs(fft)[:-1] > 0).any():
        first_read_data_ind = np.where(np.abs(fft)[1:] - np.abs(fft)[:-1] > 0)[0][0]  # away from the DC peak
        fft = fft[first_read_data_ind:]
        f = f[first_read_data_ind:]

    # Finding a guess for the frequency
    out_freq = f[np.argmax(np.abs(fft))]
    guess_freq = out_freq / (x[1] - x[0])

    # The period is 1 / guess_freq --> number of oscillations --> peaks decay to get guess_T2
    period = int(np.ceil(1 / out_freq))
    peaks = (
            np.array([np.std(y[i * period: (i + 1) * period]) for i in range(round(len(y) / period))]) * np.sqrt(
        2) * 2
    )

    # Finding a guess for the decay (slope of log(peaks))
    if len(peaks) > 1:
        guess_T2 = -1 / ((np.log(peaks)[-1] - np.log(peaks)[0]) / (period * (len(peaks) - 1))) * (x[1] - x[0])
    else:
        guess_T2 = 100 / x_normal

    # Finding a guess for the offsets
    initial_offset = np.mean(y[:period])
    final_offset = np.mean(y[-period:])

    # Finding a guess for the phase
    guess_phase = np.angle(fft[np.argmax(np.abs(fft))]) - guess_freq * 2 * np.pi * x[0]

    # Check user guess
    if guess is not None:
        for key in guess.keys():
            if key == "f":
                guess_freq = float(guess[key]) * x_normal
            elif key == "phase":
                guess_phase = float(guess[key])
            elif key == "T2":
                guess_T2 = float(guess[key]) * x_normal
            elif key == "amp":
                peaks[0] = float(guess[key]) / y_normal
            elif key == "initial_offset":
                initial_offset = float(guess[key]) / y_normal
            elif key == "final_offset":
                final_offset = float(guess[key]) / y_normal
            else:
                raise Exception(
                    f"The key '{key}' specified in 'guess' does not match a fitting parameters for this function."
                )

    # Print the initial guess if verbose=True
    if verbose:
        print(
            f"Initial guess:\n"
            f" f = {guess_freq / x_normal:.3f}, \n"
            f" phase = {guess_phase:.3f}, \n"
            f" T2 = {guess_T2 * x_normal:.3f}, \n"
            f" amp = {peaks[0] * y_normal:.3f}, \n"
            f" initial offset = {initial_offset * y_normal:.3f}, \n"
            f" final_offset = {final_offset * y_normal:.3f}"
        )

    # Fitting function
    def func(x_var, a0, a1, a2, a3, a4, a5):
        return final_offset * a4 * (1 - np.exp(-x_var / (guess_T2 * a1))) + peaks[0] / 2 * a2 * (
                np.exp(-x_var / (guess_T2 * a1))
                * (a5 * initial_offset / peaks[0] * 2 + np.cos(2 * np.pi * a0 * guess_freq * x + a3))
        )

    def fit_type(x_var, a):
        return func(x_var, a[0], a[1], a[2], a[3], a[4], a[5])

    popt, pcov = optimize.curve_fit(
        func,
        x,
        y,
        p0=[1, 1, 1, guess_phase, 1, 1],
    )

    perr = np.sqrt(np.diag(pcov))

    # Output the fitting function and its parameters
    out = {
        "fit_func": lambda x_var: fit_type(x_var / x_normal, popt) * y_normal,
        "f": [popt[0] * guess_freq / x_normal, perr[0] * guess_freq / x_normal],
        "phase": [popt[3] % (2 * np.pi), perr[3] % (2 * np.pi)],
        "T2": [(guess_T2 * popt[1]) * x_normal, perr[1] * guess_T2 * x_normal],
        "amp": [peaks[0] * popt[2] * y_normal, perr[2] * peaks[0] * y_normal],
        "initial_offset": [
            popt[5] * initial_offset * y_normal,
            perr[5] * initial_offset * y_normal,
        ],
        "final_offset": [
            final_offset * popt[4] * y_normal,
            perr[4] * final_offset * y_normal,
        ],
    }
    # Print the fitting results if verbose=True
    if verbose:
        print(
            f"Fitting results:\n"
            f" f = {out['f'][0] * 1000:.3f} +/- {out['f'][1] * 1000:.3f} MHz, \n"
            f" phase = {out['phase'][0]:.3f} +/- {out['phase'][1]:.3f} rad, \n"
            f" T2 = {out['T2'][0]:.2f} +/- {out['T2'][1]:.3f} ns, \n"
            f" amp = {out['amp'][0]:.2f} +/- {out['amp'][1]:.3f} a.u., \n"
            f" initial offset = {out['initial_offset'][0]:.2f} +/- {out['initial_offset'][1]:.3f}, \n"
            f" final_offset = {out['final_offset'][0]:.2f} +/- {out['final_offset'][1]:.3f} a.u."
        )
    # Plot the data and the fitting function if plot=True
    if plot:
        plt.plot(x_data, fit_type(x, popt) * y_normal)
        plt.plot(
            x_data,
            y_data,
            ".",
            label=f"T2  = {out['T2'][0]:.1f} +/- {out['T2'][1]:.1f}ns \n f = {out['f'][0] * 1000:.3f} +/- {out['f'][1] * 1000:.3f} MHz",
        )
        plt.legend(loc="upper right")
    t2r_est = out['T2'][0]  # in ns
    t2r_err = out['T2'][1]  # in ns
    return fit_type(x, popt) * y_normal, t2r_est, t2r_err

def optimal_bins(data):
    n = len(data)
    if n == 0:
        return {}
    # Sturges' Rule
    sturges_bins = int(np.ceil(np.log2(n) + 1))
    return sturges_bins

def plot_results(I, Q, delay_times, QubitIndex, fit, t2r_est, t2r_err, signal, save_figs):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    plt.rcParams.update({'font.size': 18})

    if 'I' in signal:
        signal = I
        plot_sig='I'
    elif 'Q' in signal:
        signal = Q
        plot_sig = 'Q'



    # I subplot

    ax1.plot(delay_times[0], I[0], label="Gain (a.u.)", linewidth=2)
    ax1.set_ylabel("I Amplitude (a.u.)", fontsize=20)
    ax1.tick_params(axis='both', which='major', labelsize=16)
    # ax1.axvline(freq_q, color='orange', linestyle='--', linewidth=2)

    # Q subplot
    ax2.plot(delay_times[0], Q[0], label="Q", linewidth=2)
    ax2.set_xlabel("Delay time (us)", fontsize=20)
    ax2.set_ylabel("Q Amplitude (a.u.)", fontsize=20)
    ax2.tick_params(axis='both', which='major', labelsize=16)
    # ax2.axvline(freq_q, color='orange', linestyle='--', linewidth=2)

    if 'I' in plot_sig:
        ax1.plot(delay_times[0], fit, '-', color='red', linewidth=3, label="Fit")
    else:
        ax2.plot(delay_times[0], fit, '-', color='red', linewidth=3, label="Fit")

    # Adjust spacing
    plt.tight_layout()

    # Calculate the middle of the plot area
    plot_middle = (ax1.get_position().x0 + ax1.get_position().x1) / 2

    # Add title, centered on the plot area
    fig.text(plot_middle, 0.98,
             f"T2 Q{QubitIndex + 1}, pi gain %.2f" + f" T2 = {t2r_est:.3f} ± {t2r_err:.3f} us",
             fontsize=24, ha='center', va='top')

    # Adjust the top margin to make room for the title
    plt.subplots_adjust(top=0.93)
    outerFolder_expt = "/data/QICK_data/6transmon_run4a/" + str(datetime.date.today())  + "/Replotting_T2_from_h5" + "/"
    create_folder_if_not_exists(outerFolder_expt)
    formatted_datetime = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    file_name = outerFolder_expt + f"Q_{QubitIndex+1}_" + f"{formatted_datetime}_" +  f"_q{QubitIndex + 1}.png"
    if save_figs:
        fig.savefig(file_name, dpi=300, bbox_inches='tight')  # , facecolor='white'
    plt.close(fig)


def create_folder_if_not_exists(folder_path):
    import os
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

# Constants and Initialization
t1_vals = {i: [] for i in range(6)}
t1_errs = {i: [] for i in range(6)}
I_vals = {i: [] for i in range(6)}
Q_vals = {i: [] for i in range(6)}
round_nums = {i: [] for i in range(6)}
batch_nums = {i: [] for i in range(6)}
delay_time_vals = {i: [] for i in range(6)}
qubit_for_this_index = []
rounds = []
reps = []
file_names = []
dates = {i: [] for i in range(6)}
mean_values = {}
dates_run = ['2024-11-11']
show_legends = False

# Grab all of the data
for date in dates_run:
    folder_path = f"/data/QICK_data/6transmon_run4a/{date}/T2_ge/"

    for file_name in os.listdir(folder_path):
        if file_name.endswith(".h5"):
            file_path = os.path.join(folder_path, file_name)

            with h5py.File(file_path, 'r') as h5_file:

                for qubit_index in range(1, 7):  # Adjusting for 0-based index
                    group_name = f"Q{qubit_index}"
                    if group_name in h5_file:
                        #print all of the groups in this file
                        #h5_file.visititems(lambda name, obj: print(name))
                        #print("-" * 40)

                        # Extracting the datasets
                        t1_data = h5_file[f"{group_name}/T2"][()]
                        t1_err_data = h5_file[f"{group_name}/Errors"][()]
                        dates_data = h5_file[f"{group_name}/Dates"][()]
                        delay_times = h5_file[f"{group_name}/Delay Times"][()]
                        I = h5_file[f"{group_name}/I"][()]
                        Q = h5_file[f"{group_name}/Q"][()]
                        round_num = h5_file[f"{group_name}/Round Num"][()]
                        batch_num = h5_file[f"{group_name}/Batch Num"][()]

                        # Appending values to the appropriate lists
                        Q_vals[qubit_index - 1].extend(Q)
                        I_vals[qubit_index - 1].extend(I)
                        round_nums[qubit_index - 1].extend(round_num)
                        batch_nums[qubit_index - 1].extend(batch_num)
                        delay_time_vals[qubit_index - 1].extend(delay_times)
                        t1_vals[qubit_index - 1].extend(t1_data)  # Store T1 values
                        t1_errs[qubit_index - 1].extend(t1_err_data)  # Store T1 error values
                        dates[qubit_index - 1].extend([d.decode('utf-8') for d in dates_data])  # Decode bytes to string

                        # You can also append qubit indices if needed
                        qubit_for_this_index.extend([qubit_index] * len(t1_data))

#recreate plots and save them for T2
for QubitIndex in range(0, 6):
    for i in range(0, len(I_vals[QubitIndex])-1):
        fit, t2r_est, t2r_err = t2_fit(delay_times[QubitIndex], I_vals[QubitIndex][i])
        plot_results(I, Q, delay_times, QubitIndex, fit, t2r_est, t2r_err, "I", True)

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
        ax.set_xlabel('T2 (µs)',fontsize = font)
        ax.set_ylabel('Frequency',fontsize = font)
        ax.tick_params(axis='both', which='major', labelsize=font)
plt.tight_layout()
plt.savefig('hists.png', transparent=True, dpi=500)

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
ax.set_xlabel('T2 (us)',fontsize = font)
ax.set_ylabel('Cumulative Distribution',fontsize = font)
ax.loglog()
ax.legend(edgecolor='black')
ax.set_xlim(10**0, 10**3)
ax.set_ylim(10 ** -7, 10 ** 0) #to compare to johns plot, need to adjust a little
plt.tight_layout()
plt.savefig('cumulative.png', transparent=True, dpi=500)



fig, axes = plt.subplots(2, 3, figsize=(12, 8))
plt.title('Fit Error vs T2 Time',fontsize = font)
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
    ax.set_xlabel('T2 (us)', fontsize = font)
    ax.set_ylabel('Fit error (us)', fontsize = font)
    ax.tick_params(axis='both', which='major', labelsize=font)
plt.tight_layout()
plt.savefig('errs.png', transparent=True, dpi=500)

plt.show()
