import numpy as np
import h5py
import os
import glob
import datetime

# Directory where your HDF5 files are stored
outerFolder1 = os.path.join("/home/nexusadmin/qick/NEXUS_sandbox/Data/Run30/2025-01-22")
outerFolder = outerFolder1 + "/readout_opt/Gain_Freq_Sweeps/"
print('For files inside: ', outerFolder)

def find_max_fidelity(file_path):
    with h5py.File(file_path, "r") as f:
        # Load the results dataset and metadata
        results = np.array(f["results"])
        gain_range = f.attrs["gain_range"]
        freq_range = f.attrs["freq_range"]
        reference_frequency = f.attrs["reference_frequency"]
        freq_steps = f.attrs["freq_steps"]
        gain_steps = f.attrs["gain_steps"]
        optimal_length = f.attrs["optimal_length"]

        # Calculate step sizes
        freq_step_size = (freq_range[1] - freq_range[0]) / freq_steps
        gain_step_size = (gain_range[1] - gain_range[0]) / gain_steps

        # Find the maximum fidelity value
        max_fidelity = np.max(results)

        # Find all indices where the maximum fidelity occurs
        max_indices = np.argwhere(results == max_fidelity)

        # Collect all configurations with the highest fidelity
        configurations = []
        for index in max_indices:
            freq_idx, gain_idx = index
            gain = gain_range[0] + gain_idx * gain_step_size
            freq_offset = (freq_range[0] - reference_frequency) + freq_idx * freq_step_size #this is doing freq - reference_frequency
            configurations.append((max_fidelity, gain, freq_offset, optimal_length))

        return max_fidelity, configurations

# Loop through each HDF5 file in the folder for each qubit
# for qubit_index in range(1, 5):
#     # Search for files matching the pattern with any timestamp
#     file_pattern = os.path.join(outerFolder, f"*_Qubit_{qubit_index}_*.h5")
#     # file_pattern = os.path.join(outerFolder, f"{str(datetime.date.today())}_Gain_Freq_Sweep_Qubit_{qubit_index}_*.h5")
#     # print("file pattern = " + file_pattern)
#     file_list = glob.glob(file_pattern)
#
#     if file_list:
#         # Initialize variables to track the maximum fidelity and configurations
#         overall_max_fidelity = -np.inf
#         all_configurations = []
#
#         # Iterate over each file and collect configurations tied for max fidelity
#         for file_path in file_list:
#             max_fidelity, configurations = find_max_fidelity(file_path)
#
#             # If this file's max fidelity is higher, reset the list of configurations
#             if max_fidelity > overall_max_fidelity:
#                 overall_max_fidelity = max_fidelity
#                 all_configurations = configurations
#             # If it matches the current max, add these configurations as well
#             elif max_fidelity == overall_max_fidelity:
#                 all_configurations.extend(configurations)
#
#         # Print all configurations with the highest fidelity
#         print(f"Qubit {qubit_index}:")
#         for max_fidelity, max_gain, max_freq_offset, optimal_length in all_configurations:
#             print(f"  Max Fidelity: {overall_max_fidelity:.4f}")
#             print(f"  Optimal Readout Length: {optimal_length:.4f} us")
#             print(f"  Optimal Gain: {max_gain:.4f} a.u.")
#             print(f"  Optimal Frequency Offset: {max_freq_offset:.4f} MHz\n")
#     else:
#         print(f"File for Qubit {qubit_index} not found.")

#-------------------------------When you want results below a specific res gain threshold for each resonator-------------------

def find_configurations_below_threshold(file_path, threshold):
    """
    Returns *all* (fidelity, gain, freq_offset, optimal_length) points for which
    the gain is strictly below the given threshold.
    """
    with h5py.File(file_path, "r") as f:
        # Read results as a Python object. If 'results' is truly a list-of-lists,
        # do something like:
        results = f["results"][()]  # This loads the dataset into a NumPy array
        results = results.tolist()  # Convert to list of lists if it's not already.

        # Retrieve metadata
        gain_range = f.attrs["gain_range"]
        freq_range = f.attrs["freq_range"]
        reference_frequency = f.attrs["reference_frequency"]
        freq_steps = f.attrs["freq_steps"]
        gain_steps = f.attrs["gain_steps"]
        optimal_length = f.attrs["optimal_length"]

    # Calculate step sizes
    freq_step_size = (freq_range[1] - freq_range[0]) / freq_steps
    gain_step_size = (gain_range[1] - gain_range[0]) / gain_steps

    valid_configurations = []

    # Outer loop: each row in 'results' corresponds to one freq_idx
    for freq_idx, row_of_fidelities in enumerate(results):
        # Inner loop: each element in the row corresponds to one gain_idx
        for gain_idx, fidelity in enumerate(row_of_fidelities):
            # Reconstruct the actual gain from 'gain_idx'
            gain = gain_range[0] + gain_idx * gain_step_size

            if gain < threshold:
                # Compute frequency offset from 'freq_idx'
                freq_offset = (freq_range[0] - reference_frequency) + freq_idx * freq_step_size

                valid_configurations.append((fidelity, gain, freq_offset, optimal_length))

    return valid_configurations

#Now getting results
punchout_thresholds = [0.25, 0.20, 0.25, 0.35]  # Set the threshold values for resonator gain
for qubit_index in range(1, 5):
    file_pattern = os.path.join(outerFolder, f"*_Qubit_{qubit_index}_*.h5")
    file_list = glob.glob(file_pattern)

    if not file_list:
        print(f"File(s) for Qubit {qubit_index} not found.")
        continue

    all_configs = []
    threshold = punchout_thresholds[qubit_index - 1]
    for file_path in file_list:
        valid_points = find_configurations_below_threshold(file_path, threshold)
        all_configs.extend(valid_points)

    if not all_configs:
        print(f"Qubit {qubit_index}: No valid points below threshold {threshold}.")
    else:
        # Find the maximum fidelity among the valid configs for this qubit
        best_fidelity = max(cfg[0] for cfg in all_configs)
        best_cfgs = [cfg for cfg in all_configs if cfg[0] == best_fidelity]

        print(f"Qubit {qubit_index} (Gain threshold: {threshold}):")
        for fidelity, gain, freq_offset, optimal_length in best_cfgs:
            print(f"  Fidelity: {fidelity:.4f}")
            print(f"  Gain: {gain:.4f}")
            print(f"  Freq Offset: {freq_offset:.4f}")
            print(f"  Opt Length: {optimal_length:.4f}\n")