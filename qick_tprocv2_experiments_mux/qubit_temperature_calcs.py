import numpy as np
import h5py
import os
from section_005_single_shot_ge import SingleShot
from system_config import QICK_experiment
import datetime
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt

def calculate_qubit_temperature(frequency_mhz, ground_state_population, excited_state_population):
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    h = 6.62607015e-34  # Planck's constant in J·s
    frequency_hz = frequency_mhz * 1e6
    T = (h * frequency_hz) / (k_B * np.log(ground_state_population / excited_state_population))
    return T


def fit_double_gaussian(iq_data):
    # Initialize a Gaussian Mixture Model with 2 components (for ground and excited states)
    gmm = GaussianMixture(n_components=2)

    # Fit the Gaussian Mixture Model to the IQ data, reshaped to a 2D array (required by GMM)
    gmm.fit(iq_data.reshape(-1, 1))

    # Extract the means of the two Gaussian components as a 1D array
    # This helps identify which Gaussian represents the ground state (lower mean)
    means = gmm.means_.flatten()  # Get the means of the Gaussians

    # Predict labels for each data point in IQ data, assigning each point to one of the two Gaussians
    # The result is an array of labels (0 or 1) indicating the Gaussian each point belongs to
    labels = gmm.predict(iq_data.reshape(-1, 1))

    # Determine which label corresponds to the ground state by identifying the Gaussian with the lower mean
    # This assigns 'ground_label' as the index of the Gaussian with the smallest mean
    ground_label = np.argmin(means)

    # Calculate the population in the ground state by computing the proportion of data points
    # assigned to the Gaussian identified as the ground state
    ground_state_population = np.mean(labels == ground_label)

    # Calculate the population in the excited state as the complement of the ground state population
    excited_state_population = np.mean(labels != ground_label)

    print(f"Ground state population: {ground_state_population}")
    print(f"Excited state population: {excited_state_population}")

    # Return the calculated ground and excited state populations
    return ground_state_population, excited_state_population

# Output folder configuration
date_str = str(datetime.date.today())
output_folder = f"/data/QICK_data/6transmon_run4a/{date_str}/Qubit_Temps/"
os.makedirs(output_folder, exist_ok=True)

# Optimized parameters for qubits 3 and 4
qubits = {
    3: {"length": 3.50, "gain": 0.8462, "freq_offset": -0.2000, "reference_freq": 6292.261},
    4: {"length": 4.75, "gain": 0.9615, "freq_offset": -0.2000, "reference_freq": 6405.79}
}

# Loop through specified qubits and apply settings
for qubit_index, params in qubits.items():
    length = params["length"]
    gain = params["gain"]
    frequency = params["reference_freq"] + params["freq_offset"]

    print(f"Processing Qubit {qubit_index} with Length {length}, Gain {gain}, and Frequency {frequency}")

    # Initialize experiment
    experiment = QICK_experiment(output_folder)
    experiment.readout_cfg['res_length'] = length
    res_gains = experiment.set_gain_filter_ge(qubit_index - 1, gain) # in this script qubit_index starts at 1 not zero, so we must subtract 1
    experiment.readout_cfg['res_gain_ge'] = res_gains

    # Run the single-shot experiment
    ss = SingleShot(qubit_index - 1, output_folder, experiment, round_num=0, save_figs=False) # in this script qubit_index starts at 1 not zero, so we must subtract 1
    fid, angle, iq_list_g, iq_list_e = ss.run(experiment.soccfg, experiment.soc)

    I_g = iq_list_g[qubit_index-1][0].T[0]
    Q_g = iq_list_g[qubit_index-1][0].T[1]
    I_e = iq_list_e[qubit_index-1][0].T[0]
    Q_e = iq_list_e[qubit_index-1][0].T[1]

    # Call hist_ssf to get rotated I data
    fid, threshold, rotation_angle, ig_new, ie_new = ss.hist_ssf(
        data=[I_g, Q_g, I_e, Q_e], cfg=ss.config, plot=True)

    # Use the rotated I data for Gaussian fitting
    combined_I_data = np.concatenate([ig_new, ie_new]) # array that includes all the I-values from both states. The function fit_double_gaussian will fit two Gaussian distributions over the combined data.
    ground_state_population, excited_state_population = fit_double_gaussian(combined_I_data) #gets ground state and excited state populations

    # Calculate temperature in Kelvin and millikelvin
    temperature_k = calculate_qubit_temperature(frequency, ground_state_population, excited_state_population)
    temperature_mk = temperature_k * 1e3

    # Save results in HDF5
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    h5_filename = os.path.join(output_folder, f"qubit_{qubit_index}_temperatureSSdata_{timestamp}.h5")
    with h5py.File(h5_filename, 'w') as f:
        f.create_dataset("fidelity", data=fid)
        f.create_dataset("gain", data=gain)
        f.create_dataset("frequency", data=frequency)
        f.create_dataset("ground_iq_data", data=iq_list_g)
        f.create_dataset("excited_iq_data", data=iq_list_e)
        f.create_dataset("temperature_kelvin", data=temperature_k)
        f.create_dataset("temperature_millikelvin", data=temperature_mk)
        f.attrs['ground_state_population'] = ground_state_population
        f.attrs['excited_state_population'] = excited_state_population
        f.attrs['rotation_angle'] = rotation_angle

    print(f"Data and temperature saved for Qubit {qubit_index} to {h5_filename}")
    print(f"Temperature for Qubit {qubit_index}: {temperature_mk:.2f} mK")
    del experiment
