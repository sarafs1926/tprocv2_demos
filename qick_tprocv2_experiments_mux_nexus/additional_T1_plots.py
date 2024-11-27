import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

# Define the exponential function for fitting
def exponential(x, a, b, c, d):
    return a * np.exp(-(x - b) / c) + d

t1_vals = {i: [] for i in range(0, 6)}
t1_errs = {i: [] for i in range(0, 6)}
# Counter for files that required recalculation
recalculation_count = 0

# Define colors to match previous color scheme
colors = ['orange', 'blue', 'purple', 'green', 'brown', 'pink']

# Dictionary to store T1 values for each qubit
t1_vals = {i: [] for i in range(6)}

# Specify the dates to retrieve data from
dates_run = ['2024-10-23', '2024-10-24', '2024-10-25', '2024-10-26']
q6_files_with_high_t1_err = []

# Iterate over each date and file to collect T1 values
for date in dates_run:
    folder_path = f"/data/QICK_data/6transmon_run4a/{date}/T1_ge/"
    if '23' in date or '24' in date:
        updated_q2_folder = os.path.join(folder_path, "Updated_Q2_Data")
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".h5"):
            # Determine the file path based on qubit-specific subfolders
            file_path = os.path.join(updated_q2_folder, file_name) if "q2" in file_name else os.path.join(folder_path, file_name)
            with h5py.File(file_path, 'r') as h5_file:
                qubit_num = h5_file['Qubit_num'][0]  # Get qubit number
                t1_val = float(h5_file['T1_est'][0])  # Retrieve T1 estimate
                #t1_vals[qubit_num - 1].append(t1_val)  # Add T1 to the list for this qubit

                """"
                # Check T1_err for Qubit 6
                if qubit_num == 6:
                    t1_err = float(h5_file['T1_err'][0])
                    if t1_err > 10:
                        q6_files_with_high_t1_err.append(file_name)  # Store file name if T1_err > 10
                        print('T1 val: ',h5_file['T1_est'][0])
                        print('T1 err: ',t1_err)
                """
                T1_est = float(h5_file['T1_est'][0])
                T1_err = float(h5_file['T1_err'][0])

                if T1_err > 10:
                    print('old T1 val: ', T1_est)
                    print('old T1 err: ', T1_err)

                    signal = h5_file["Q"][:]
                    delay_times = h5_file["delay_times"][:]
                    # Increment recalculation count
                    recalculation_count += 1

                    # Initial guesses for the fit parameters
                    q1_a_guess = np.max(signal) - np.min(signal)  # Amplitude
                    q1_b_guess = 0  # Time shift
                    q1_c_guess = (delay_times[-1] - delay_times[0]) / 5  # Decay constant T1
                    q1_d_guess = np.min(signal)  # Baseline

                    # Perform the fit with bounds
                    q1_guess = [q1_a_guess, q1_b_guess, q1_c_guess, q1_d_guess]
                    lower_bounds = [-np.inf, -np.inf, 0, -np.inf]
                    upper_bounds = [np.inf, np.inf, np.inf, np.inf]

                    q1_popt, q1_pcov = curve_fit(exponential, delay_times, signal, p0=q1_guess,
                                                 bounds=(lower_bounds, upper_bounds), method='trf', maxfev=10000)

                    # Update T1 and T1_err with the recalculated values
                    T1_est = q1_popt[2]
                    T1_err = np.sqrt(q1_pcov[2][2]) if q1_pcov[2][2] >= 0 else float('inf')

                    print('new T1 val: ', T1_est)
                    print('new T1 err: ', T1_err)
                    print('Number of recalculated T1s: ', recalculation_count)

                    # Replace old values with recalculated T1 and T1_err
                    #t1_vals[qubit_num - 1][-1] = T1_est
                    #t1_errs[qubit_num - 1][-1] = T1_err

                t1_vals[qubit_num-1].append(T1_est)


# Print names of files with T1_err > 10 for Qubit 6
#print("Files with T1_err > 10 for Qubit 6:", q6_files_with_high_t1_err)
print("Qubit 6 T1 values over 40 µs:", [val for val in t1_vals[5] if val > 40])

# Plot T1 vs iteration number for each qubit
fig, axes = plt.subplots(2, 3, figsize=(12, 8))
axes = axes.flatten()
for i, ax in enumerate(axes):
    if t1_vals[i]:  # Check if there is data for this qubit
        iterations = range(1, len(t1_vals[i]) + 1)  # Define iteration numbers
        ax.scatter(iterations, t1_vals[i], color=colors[i], label=f"Qubit {i+1}")
        ax.set_title(f"Qubit {i+1}")
        ax.set_xlabel("Iteration Number")
        ax.set_ylabel("T1 (µs)")
        ax.legend(edgecolor='black')

        # Set a wider y-axis range
        y_min =  min(t1_vals[i]) * 0.3
        y_max = max(t1_vals[i]) * 1.3
        ax.set_ylim(y_min, y_max)

plt.tight_layout()
plt.show()

# Plot all T1 values on a single scatter plot
plt.figure(figsize=(10, 6))
for i in range(6):
    if t1_vals[i]:  # Check if there is data for this qubit
        iterations = range(1, len(t1_vals[i]) + 1)  # Define iteration numbers
        plt.scatter(iterations, t1_vals[i], color=colors[i], label=f"Qubit {i+1}")
# Customize the plot
plt.title("T1 vs Iteration Number for All Qubits")
plt.xlabel("Iteration Number")
plt.ylabel("T1 (µs)")
plt.legend(edgecolor='black')
plt.tight_layout()
plt.show()

#print("T1 values for Qubit 6:", t1_vals[5])
###########################################################################
def exponential(x, a, b, c, d):
    return a * np.exp(-(x - b) / c) + d

# Specify the target file name
target_file = "2024-10-26_21-47-50_T1_ge_q6.h5"

# Specify the folder path where your file is located
folder_path = "/data/QICK_data/6transmon_run4a/2024-10-26/T1_ge/"

# Check if the target file exists in the directory and load data if it does
if target_file in os.listdir(folder_path):
    file_path = os.path.join(folder_path, target_file)

    # Open the HDF5 file and load data
    with h5py.File(file_path, 'r') as f:
        I = f["I"][:]
        Q = f["Q"][:]
        delay_times = f["delay_times"][:]
        pi_amp = f['pi_amp'][0]
        sigma = f['sigma'][0]
        reps = f['reps'][0]
        rounds = f['rounds'][0]
        T1_est = f['T1_est'][0]
        T1_err = f['T1_err'][0]
        plot_middle = f['plot_middle'][0]
        Qubit_num = f['Qubit_num'][0]

    # Initial guesses for the fit parameters
    signal = Q
    q1_a_guess = np.max(signal) - np.min(signal)  # Amplitude
    q1_b_guess = 0  # Time shift
    q1_c_guess = (delay_times[-1] - delay_times[0]) / 5  # Decay constant T1
    q1_d_guess = np.min(signal)  # Baseline

    # Perform the fit with bounds
    q1_guess = [q1_a_guess, q1_b_guess, q1_c_guess, q1_d_guess]
    lower_bounds = [-np.inf, -np.inf, 0, -np.inf]
    upper_bounds = [np.inf, np.inf, np.inf, np.inf]

    q1_popt, q1_pcov = curve_fit(exponential, delay_times, signal, p0=q1_guess,
                                 bounds=(lower_bounds, upper_bounds), method='trf', maxfev=10000)

    # Generate the fitted exponential curve
    q1_fit_exponential = exponential(delay_times, *q1_popt)

    # Extract T1 and its error
    T1_est = q1_popt[2]
    T1_err = np.sqrt(q1_pcov[2][2]) if q1_pcov[2][2] >= 0 else float('inf')

    # Plot the data
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    plt.rcParams.update({'font.size': 18})

    # I subplot
    ax1.plot(delay_times, I, label="Gain (a.u.)", linewidth=2)
    ax1.set_ylabel("I Amplitude (a.u.)", fontsize=20)
    ax1.tick_params(axis='both', which='major', labelsize=16)
    #ax1.plot(delay_times, q1_fit_exponential, '-', color='red', linewidth=3, label="Fit")

    # Q subplot
    ax2.plot(delay_times, Q, label="Q", linewidth=2)
    ax2.plot(delay_times, q1_fit_exponential, '-', color='red', linewidth=3, label="Fit")
    ax2.set_xlabel("Delay time (us)", fontsize=20)
    ax2.set_ylabel("Q Amplitude (a.u.)", fontsize=20)
    ax2.tick_params(axis='both', which='major', labelsize=16)

    # Add text annotation with the fit parameters
    fig.text(plot_middle, 0.99,
             f"T1 Q{Qubit_num}, pi gain {pi_amp:.2f}, {sigma * 1000:.1f} ns sigma, "
             f"{reps}*{rounds} avgs, T1 = {T1_est:.3f} ± {T1_err:.3f} µs",
             fontsize=16, ha='center', va='top')

    plt.tight_layout()
    plt.show()
else:
    print(f"File {target_file} not found in the directory.")