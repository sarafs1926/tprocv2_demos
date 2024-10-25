import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os

T1_est_list = []
delay_times_list = []
I_list = []
Q_list = []

folder_path2 = "/data/QICK_data/6transmon_run4a/2024-10-24/T1_ge/"

for file_name in os.listdir(folder_path2):
    if file_name.endswith(".h5") and "q2" in file_name:
        file_path = os.path.join(folder_path2, file_name)
        with h5py.File(file_path, 'r') as f:
            I = f["I"][:]
            I_list.append(I)
            Q = f["Q"][:]
            Q_list.append(Q)
            delay_times = f["delay_times"][:]
            delay_times_list.append(delay_times)
            q1_fit_exponential = f["q1_fit_exponential"][:]

            # Load the parameters needed for plotting the text
            pi_amp = f['pi_amp'][0]
            sigma = f['sigma'][0]
            reps = f['reps'][0]
            rounds = f['rounds'][0]
            T1_est = f['T1_est'][0]
            T1_est_list.append(T1_est)
            T1_err = f['T1_err'][0]
            # Load plot_middle
            plot_middle = f['plot_middle'][0]
            # Load Qubit_num
            Qubit_num = f['Qubit_num'][0]


            def exponential(x, a, b, c, d):
                return a * np.exp(- (x - b) / c) + d

            signal = Q
            #average_last_nums = np.mean(signal[-20:])

            #if signal[-1] > signal[0]:
                #q1_a_guess = -( np.max(signal) - np.min(signal))
            #else:
                #q1_a_guess = (np.max(signal) - np.min(signal))
            q1_a_guess = (np.max(signal) - np.min(signal))
            #if "inf" in str(T1_err):
            q1_a_guess = np.max(signal) - np.min(signal)  # Initial guess for amplitude (a)
            q1_b_guess = 0  # Initial guess for time shift (b)
            q1_c_guess = (delay_times[-1] - delay_times[0]) / 5  # Initial guess for decay constant (T1)
            q1_d_guess = np.min(signal)  # Initial guess for baseline (d)

            # Form the guess array
            q1_guess = [q1_a_guess, q1_b_guess, q1_c_guess, q1_d_guess]

            # Define bounds to constrain T1 (c) to be positive, but allow amplitude (a) to be negative
            lower_bounds = [-np.inf, -np.inf, 0, -np.inf]  # Amplitude (a) can be negative/positive, but T1 (c) > 0
            upper_bounds = [np.inf, np.inf, np.inf, np.inf]  # No upper bound on parameters

            # Perform the fit using the 'trf' method with bounds
            q1_popt, q1_pcov = curve_fit(exponential, delay_times, signal,
                                             p0=q1_guess, bounds=(lower_bounds, upper_bounds),
                                             method='trf', maxfev=10000)

            # Generate the fitted exponential curve
            q1_fit_exponential = exponential(delay_times, *q1_popt)

            # Extract T1 and its error
            T1_est = q1_popt[2]  # Decay constant T1
            T1_err = np.sqrt(q1_pcov[2][2]) if q1_pcov[2][2] >= 0 else float('inf')  # Ensure error is valid

            """
            #divide_val = delay_times[-1] / (T1_est/5)
            q1_b_guess = 0
            q1_c_guess = delay_times[0]
            q1_d_guess = np.min(signal)

            q1_guess = [q1_a_guess, q1_b_guess, q1_c_guess, q1_d_guess]
            q1_popt, q1_pcov = curve_fit(exponential, delay_times, signal, maxfev=500000, p0=q1_guess)
            q1_fit_exponential = exponential(delay_times, *q1_popt)

            T1_est = q1_popt[2]
            T1_err = np.sqrt(q1_pcov[2][2])
            """
            # Plot the data
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
            plt.rcParams.update({'font.size': 18})

            # I subplot
            ax1.plot(delay_times, I, label="Gain (a.u.)", linewidth=2)
            ax1.set_ylabel("I Amplitude (a.u.)", fontsize=20)
            ax1.tick_params(axis='both', which='major', labelsize=16)

            # Q subplot
            ax2.plot(delay_times, Q, label="Q", linewidth=2)
            ax2.plot(delay_times, q1_fit_exponential, '-', color='red', linewidth=3, label="Fit")
            ax2.set_xlabel("Delay time (us)", fontsize=20)
            ax2.set_ylabel("Q Amplitude (a.u.)", fontsize=20)
            ax2.tick_params(axis='both', which='major', labelsize=16)

            #plt.subplots_adjust(top=0.75)

            fig.text(plot_middle, 0.99,
                     f"T1 Q{Qubit_num}, pi gain %.2f" % pi_amp +
                     f", {sigma * 1000} ns sigma" + f", {reps}*{rounds} avgs," +
                    f" T1 = {T1_est:.3f} ± {T1_err:.3f} µs", fontsize=16, ha='center', va='top')

            plt.tight_layout()
            save_folder = '/data/QICK_data/6transmon_run4a/2024-10-24/T1_ge/Updated_Q2_Data/'
            plt.savefig(save_folder+ file_name.split('.')[0] + '.png')

            # Save all data to a new .h5 file with the updated parameters
            h5_save_path = os.path.join(save_folder, file_name)  # Keep the same file name
            with h5py.File(file_path, 'r') as f_src, h5py.File(h5_save_path, 'w') as f_dest:
                # Copy all datasets from original file to new file
                for key in f_src.keys():
                    f_src.copy(key, f_dest)

                # Update modified parameters
                f_dest["T1_est"][...] = T1_est
                f_dest["T1_err"][...] = T1_err
                if "q1_fit_exponential" in f_dest:
                    del f_dest["q1_fit_exponential"]  # Remove old dataset if it exists
                f_dest.create_dataset("q1_fit_exponential", data=q1_fit_exponential)



