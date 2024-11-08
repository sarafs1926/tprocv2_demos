import numpy as np
import matplotlib.pyplot as plt
import os
from section_005_single_shot_ge import GainFrequencySweep
import datetime

date_str = str(datetime.date.today())
outerFolder = f"/data/QICK_data/6transmon_run4a/{date_str}/readout_opt/Gain_Freq_Sweeps/"

# Ensure the output folder exists
os.makedirs(outerFolder, exist_ok=True)

# Reference frequencies for each resonator in MHz
res_freq_ge = [6191.479, 6216.060, 6292.301, 6405.83, 6433.059, 6468.901]

# optimal_lengths = [5.59, 9.66, 7.62, 9.66, 5.59, 9.66]
optimal_lengths = [5.59, 5, 7.62, 6.67, 5.59, 9.66]

# Define sweeping parameters
gain_range = [0.5,1]  # Gain range in a.u.
freq_steps = 40
gain_steps = 10

# for QubitIndex in range(6):
QubitIndex= 3
print(f'Starting Qubit {QubitIndex + 1} measurements.')
# Select the reference frequency for the current qubit
reference_frequency = res_freq_ge[QubitIndex]
freq_range = [reference_frequency-1, reference_frequency + 1]  # Frequency range in MHz

sweep = GainFrequencySweep(QubitIndex,  optimal_lengths=optimal_lengths, output_folder=outerFolder)
results = sweep.run_sweep(freq_range, gain_range, freq_steps, gain_steps)
results = np.array(results)

''''
commented out for debugging
# Save results and metadata in an HDF5 file
h5_file = f"{outerFolder}Gain_Freq_Sweep_Qubit_{QubitIndex + 1}.h5"
with h5py.File(h5_file, "w") as f:
    # Store the data
    f.create_dataset("results", data=results)
    # Store metadata
    f.attrs["gain_range"] = gain_range
    f.attrs["freq_range"] = freq_range
    f.attrs["reference_frequency"] = reference_frequency
    f.attrs["freq_steps"] = freq_steps
    f.attrs["gain_steps"] = gain_steps

print(f"Saved data for Qubit {QubitIndex + 1} to {h5_file}")
'''

plt.imshow(results, aspect='auto',
           extent=[gain_range[0], gain_range[1], freq_range[0] - reference_frequency, freq_range[1] - reference_frequency],
           origin='lower')
plt.colorbar(label="Fidelity")
plt.xlabel("Readout pulse gain (a.u.)")  # Gain on x-axis
plt.ylabel("Readout frequency offset (MHz)")  # Frequency on y-axis
plt.title(f"Gain-Frequency Sweep for Qubit {QubitIndex + 1}")
plt.show()
# Save the plot
file= f"{outerFolder}Gain_Freq_Sweep_Qubit_{QubitIndex + 1}.png"
plt.savefig(file, dpi=600, bbox_inches='tight')
plt.close()  # Close the plot to free up memory

print(f"Saved plot for Qubit {QubitIndex + 1} to {file}")