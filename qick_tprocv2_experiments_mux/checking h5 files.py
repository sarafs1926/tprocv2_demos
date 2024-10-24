import h5py
import matplotlib.pyplot as plt

#

# Open the HDF5 file and load the I, Q, and delay_times data
with h5py.File("/data/QICK_data/6transmon_run4a/2024-10-23/T1_ge/2024-10-23_17-23-24_T1_ge_q2.h5", "r") as f:
    I = f["I"][:]
    Q = f["Q"][:]
    delay_times = f["delay_times"][:]
    q1_fit_exponential = f["q1_fit_exponential"][:]

    # Load the parameters needed for plotting the text
    pi_amp = f['pi_amp'][0]
    sigma = f['sigma'][0]
    reps = f['reps'][0]
    rounds = f['rounds'][0]
    T1_est = f['T1_est'][0]
    T1_err = f['T1_err'][0]
    # Load plot_middle
    plot_middle = f['plot_middle'][0]
    # Load Qubit_num
    Qubit_num = f['Qubit_num'][0]

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
plt.show()
