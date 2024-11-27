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

def fit_double_gaussian_with_overlap_exclusion(iq_data):
    # Initialize and fit a Gaussian Mixture Model with 2 components
    gmm = GaussianMixture(n_components=2)
    gmm.fit(iq_data.reshape(-1, 1))

    # Get the means and covariances to find the "true" ground and excited (leakage) states
    means = gmm.means_.flatten()
    covariances = np.sqrt(gmm.covariances_).flatten()

    # Identify the ground state Gaussian as the one with the lower mean
    ground_gaussian = np.argmin(means)
    excited_gaussian = 1 - ground_gaussian  # Other Gaussian is considered the excited state overlap. Since there are only two Gaussians, the excited Gaussian will have the other index (either 0 or 1)

    # Calculate the threshold (midpoint) between the two Gaussian means
    midpoint = (means[ground_gaussian] + means[excited_gaussian]) / 2

    # Label each data point according to which Gaussian it belongs to
    labels = gmm.predict(iq_data.reshape(-1, 1))

    # Exclude data points within the overlapping region defined by the midpoint
    ground_data = iq_data[(labels == ground_gaussian) & (iq_data < midpoint)]
    excited_data = iq_data[(labels == excited_gaussian) & (iq_data > midpoint)]

    # Calculate populations as the fraction of data points in each Gaussian, excluding the shared region
    ground_state_population = len(ground_data) / len(iq_data)
    excited_state_population_overlap = len(excited_data) / len(iq_data)

    print(f"Ground state population (true): {ground_state_population}")
    print(f"Excited state population (overlap): {excited_state_population_overlap}")

    return ground_state_population, excited_state_population_overlap

# Output folder configuration
date_str = str(datetime.date.today())
output_folder = f"/data/QICK_data/6transmon_run4a/{date_str}/Qubit_Temps/"
os.makedirs(output_folder, exist_ok=True)

# Optimized parameters for qubits 3 and 4
qubits = {
    3: {"length": 3.50, "gain": 0.8462, "freq_offset": -0.2000, "reference_freq": 6292.261, "qubit_freq": 4156.57},
    4: {"length": 4.75, "gain": 0.9615, "freq_offset": -0.2000, "reference_freq": 6405.79, "qubit_freq": 4459.19}
}

# Loop through specified qubits and apply settings
for qubit_index, params in qubits.items():
    length = params["length"]
    gain = params["gain"]
    frequency = params["reference_freq"] + params["freq_offset"] #resonator offset frequency
    qubit_frequency = params["qubit_freq"]

    print(f"Processing Qubit {qubit_index} with Res_Length {length}, Res_Gain {gain}, Res Freq {frequency}, and Qubit Frequency {qubit_frequency}")

    # Initialize experiment
    QubitIndex = qubit_index - 1
    experiment = QICK_experiment(output_folder)
    experiment.readout_cfg['res_length'] = length
    experiment.readout_cfg['res_freq_ge'][QubitIndex] = frequency
    res_gains = experiment.mask_gain_res(QubitIndex, gain)
    experiment.readout_cfg['res_gain_ge'] = res_gains

    # Run the single-shot experiment
    ss = SingleShot(QubitIndex, output_folder, experiment, round_num=0, save_figs=False)
    fid, angle, iq_list_g, iq_list_e = ss.run(experiment.soccfg, experiment.soc)

    I_g = iq_list_g[QubitIndex][0].T[0]  # Ground-state I data
    Q_g = iq_list_g[QubitIndex][0].T[1]
    I_e = iq_list_e[QubitIndex][0].T[0]
    Q_e = iq_list_e[QubitIndex][0].T[1]

    # Call hist_ssf to get rotated I data
    fid, threshold, rotation_angle, ig_new, ie_new = ss.hist_ssf(
        data=[I_g, Q_g, I_e, Q_e], cfg=ss.config, plot=True)

    # Fit the rotated ground-state I data only to two Gaussians and compute overlap with excited Gaussian
    ground_state_population, excited_state_population_overlap = fit_double_gaussian_with_overlap_exclusion(ig_new)

    # Calculate temperature in Kelvin and millikelvin
    temperature_k = calculate_qubit_temperature(qubit_frequency, ground_state_population, excited_state_population_overlap)
    temperature_mk = temperature_k * 1e3

    # Save results in HDF5
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    h5_filename = os.path.join(output_folder, f"Qubit_{qubit_index}_temperatureSSdata_{timestamp}.h5")
    with h5py.File(h5_filename, 'w') as f:
        f.create_dataset("fidelity", data=fid)
        f.create_dataset("gain", data=gain)
        f.create_dataset("frequency", data=frequency)
        f.create_dataset("qubit_frequency", data=qubit_frequency)
        f.create_dataset("ground_iq_data", data=iq_list_g)
        f.create_dataset("excited_iq_data", data=iq_list_e)
        f.create_dataset("ig_new", data=ig_new)  # Save rotated ground-state I data
        f.create_dataset("ie_new", data=ie_new)  # Save rotated excited-state I data
        f.create_dataset("temperature_kelvin", data=temperature_k)
        f.create_dataset("temperature_millikelvin", data=temperature_mk)
        f.attrs['ground_state_population'] = ground_state_population
        f.attrs['excited_state_population_overlap'] = excited_state_population_overlap
        f.attrs['rotation_angle'] = rotation_angle

    print(f"Data and temperature saved for Qubit {qubit_index} to {h5_filename}")
    print(f"Temperature for Qubit {qubit_index}: {temperature_mk:.2f} mK")
    del experiment
