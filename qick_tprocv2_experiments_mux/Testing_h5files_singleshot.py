import h5py
import numpy as np
import matplotlib.pyplot as plt

# Path to the HDF5 file for a specific qubit (replace with actual file path)
h5_filename = '/data/QICK_data/6transmon_run4a/2024-11-17/SingleShot_Test/qubit_4_data_20241117_205640.h5'  # Update with your file path


# Lists to store pulse lengths, average fidelities, and RMS fidelities
pulse_lengths = []
avg_fidelities = []
rms_fidelities = []

# Load data from the HDF5 file
with h5py.File(h5_filename, 'r') as h5_file:
    qubit_group = h5_file[f"Qubit_4"]

    # Iterate over each length to retrieve the avg_fidelity and rms_fidelity data
    for length_key in qubit_group.keys():
        length_group = qubit_group[length_key]

        # Store the pulse length
        pulse_length = float(length_key.split('_')[1])
        pulse_lengths.append(pulse_length)

        # Retrieve average fidelity and RMS fidelity
        avg_fid = length_group['avg_fidelity'][()]
        rms_fid = length_group['rms_fidelity'][()]

        avg_fidelities.append(avg_fid)
        rms_fidelities.append(rms_fid)

# Sort pulse lengths and corresponding fidelity data
sorted_indices = np.argsort(pulse_lengths)
pulse_lengths = np.array(pulse_lengths)[sorted_indices]
avg_fidelities = np.array(avg_fidelities)[sorted_indices]
rms_fidelities = np.array(rms_fidelities)[sorted_indices]

# Find the maximum average fidelity and corresponding length
max_fidelity = max(avg_fidelities[:10])
max_fid_index = avg_fidelities.tolist().index(max_fidelity)
max_length = pulse_lengths[max_fid_index]
print('first max length: ', max_length)

max_fidelity2 = max(avg_fidelities[:22])
max_fid_index2 = avg_fidelities.tolist().index(max_fidelity2)
max_length2 = pulse_lengths[max_fid_index2]
print('second max length: ', max_length2)

# Plot the average fidelity vs. pulse length with error bars
plt.figure()
plt.errorbar(pulse_lengths, avg_fidelities, yerr=rms_fidelities, fmt='-o', color='black', capsize=2)
#plt.axvline(x=max_length, linestyle="--", color="red")
#plt.axvline(x=max_length2, linestyle="--", color="green")
#plt.text(max_length + 0.1, max_fidelity-0.2, f'max length {max_length:.2f}', color='red')
#plt.text(max_length2 + 0.1, max_fidelity2-0.2, f'max length {max_length2:.2f}', color='green')
plt.xlabel('Readout and Pulse Length')
plt.ylabel('Fidelity')
plt.title('Fidelity vs. Pulse Length Q4')
plt.show()
plt.close()
