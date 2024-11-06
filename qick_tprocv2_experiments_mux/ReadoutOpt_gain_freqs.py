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

optimal_lengths = [5.59, 9.66, 7.62, 9.66, 5.59, 9.66]

# Define sweeping parameters
gain_range = [0.5, 1.0]  # Gain range in a.u.
freq_steps = 45
gain_steps = 10

for QubitIndex in range(6):
    print(f'Starting Qubit {QubitIndex + 1} measurements.')
    # Select the reference frequency for the current qubit
    reference_frequency = res_freq_ge[QubitIndex]
    freq_range = [reference_frequency - 2, reference_frequency + 2]  # Frequency range in MHz

    sweep = GainFrequencySweep(QubitIndex,  optimal_lengths=optimal_lengths, output_folder=outerFolder)
    results = sweep.run_sweep(freq_range, gain_range, freq_steps, gain_steps)
    results = np.array(results)

    # Plotting with frequency offset
    plt.imshow(results, aspect='auto',
               extent=[freq_range[0] - reference_frequency, freq_range[1] - reference_frequency, gain_range[0], gain_range[1]],
               origin='lower')
    plt.colorbar(label="Fidelity")
    plt.xlabel("Readout frequency offset (MHz)")
    plt.ylabel("Readout pulse gain (a.u.)")
    plt.title(f"Gain-Frequency Sweep for Qubit {QubitIndex + 1}")
    #plt.show()
    # Save the plot
    file= f"{outerFolder}Gain_Freq_Sweep_Qubit_{QubitIndex + 1}.png"
    plt.savefig(file, dpi=600, bbox_inches='tight')
    plt.close()  # Close the plot to free up memory

    print(f"Saved plot for Qubit {QubitIndex + 1} to {file}")