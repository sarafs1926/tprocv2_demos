from section_005_single_shot_ge import SingleShot
from section_005_single_shot_ge import GainFrequencySweep
# from system_config import *
from system_config import QICK_experiment
import numpy as np
import matplotlib.pyplot as plt
import os
import h5py
import datetime
import time
from expt_config import *
from qick.asm_v2 import AveragerProgramV2
from tqdm import tqdm
from build_state import *
from expt_config import *


def create_folder_if_not_exists(folder_path):
    """Creates a folder at the given path if it doesn't already exist."""
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)


# Where to save data
prefix = str(datetime.date.today())
output_folder = "/data/QICK_data/6transmon_run4a/" + prefix + "/SingleShot_Test/"
create_folder_if_not_exists(output_folder)

n = 1  # Number of rounds
n_loops = 5  # Number of repetitions per length to average
tot_num_of_qubits = 6

# List of qubits and pulse lengths to measure
list_of_all_qubits = list(range(tot_num_of_qubits))
Qs = [0,1,2,3,4,5]

optimal_lengths = [None] * 6 # creates list where the script will be storing the optimal readout lengths for each qubit. We currently have 6 qubits in total.

#res_gain = [0.7, 0.9, 0.7, 0.7, 0.7, 0.9, 0.9]
res_gain = [1.0]*6

#lengs = np.linspace(0.5, 5, 19)  # increments of 0.25
lengs = np.linspace(0.5, 7, 27) # increments of 0.25

for QubitIndex in Qs:
    QubitIndex = int(QubitIndex)  # Ensure QubitIndex is an integer

    avg_fids = []
    rms_fids = []

    # Create a single HDF5 file for each qubit
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    # print("Output folder:", output_folder)
    # print("Type of output_folder:", type(output_folder))
    h5_filename = os.path.join(output_folder, f"qubit_{QubitIndex + 1}_data_{timestamp}.h5")
    with h5py.File(h5_filename, 'w') as h5_file:
        # Top-level group for the qubit
        qubit_group = h5_file.create_group(f"Qubit_{QubitIndex + 1}")

        # Iterate over rounds
        # for j in range(1, n + 1): #rounds
        # Create a group for each round
        # round_group = qubit_group.create_group(f"Round_{j}")
        fids = []  # Store fidelity values for each loop
        ground_iq_data = []  # Store ground state IQ data for each loop
        excited_iq_data = []  # Store excited state IQ data for each loop

        # Iterate over each readout pulse length
        for leng in lengs:
            # Subgroup for each readout length within the round
            length_group = qubit_group.create_group(f"Length_{leng}")

            for k in range(n_loops):  # loops for each read out length
                # ------------------------Single Shot-------------------------
                # Initialize experiment for each loop iteration
                #experiment = QICK_experiment(output_folder)
                experiment = QICK_experiment(output_folder, DAC_attenuator1=10, DAC_attenuator2=5, ADC_attenuator=10)
                # Set specific configuration values for each iteration
                experiment.readout_cfg['res_length'] = leng  # Set the current readout pulse length

                # Set gain for the current qubit
                gain = res_gain[QubitIndex]
                res_gains = experiment.mask_gain_res(QubitIndex, gain)  # Set gain for current qubit only

                experiment.readout_cfg['res_gain_ge'] = res_gains

                # ss = SingleShot(QubitIndex,list_of_all_qubits, output_folder, k, round(leng, 3)) #Old way
                ss = SingleShot(QubitIndex,list_of_all_qubits, output_folder, experiment, round_num=k, save_figs=False)  # New way
                fid, angle, iq_list_g, iq_list_e = ss.run(experiment.soccfg, experiment.soc)
                fids.append(fid)
                print(f'FID (round {k}) = {fid}')

                # Append IQ data for each loop
                ground_iq_data.append(iq_list_g)
                excited_iq_data.append(iq_list_e)

                # Save individual fidelity and IQ data for this loop
                loop_group = length_group.create_group(f"Loop_{k + 1}")
                loop_group.create_dataset("fidelity", data=fid)
                loop_group.create_dataset("ground_iq_data", data=iq_list_g)
                loop_group.create_dataset("excited_iq_data", data=iq_list_e)

                del experiment

            # Calculate average and RMS for fidelities across loops
            avg_fid = np.mean(fids)
            rms_fid = np.std(fids)
            avg_fids.append(avg_fid)
            rms_fids.append(rms_fid)

            # Calculate average IQ data across all loops
            avg_ground_iq = np.mean(ground_iq_data, axis=0)
            avg_excited_iq = np.mean(excited_iq_data, axis=0)

            # Save the averages and RMS to the HDF5 file for this length
            length_group.create_dataset("avg_fidelity", data=avg_fid)
            length_group.create_dataset("rms_fidelity", data=rms_fid)
            length_group.create_dataset("avg_ground_iq_data", data=avg_ground_iq)
            length_group.create_dataset("avg_excited_iq_data", data=avg_excited_iq)

            fids.clear()
            ground_iq_data.clear()
            excited_iq_data.clear()

    avg_max = max(avg_fids[:10])
    avg_max = max(avg_fids)
    avg_max_index = avg_fids.index(avg_max)
    max_len = lengs[avg_max_index]
    optimal_lengths[QubitIndex] = max_len
    
    # # Compute the first and second derivatives
    # first_derivative = np.gradient(avg_fids, lengs)
    # second_derivative = np.gradient(first_derivative, lengs)
    # 
    # # Find the index of the maximum second derivative (absolute value)
    # corner_index = np.argmax(np.abs(second_derivative))
    # max_len_corner = lengs[corner_index]
    # 
    # # Save the optimal length corresponding to the corner
    # optimal_lengths[QubitIndex] = max_len_corner
    

    # Plot the average fidelity vs. pulse length with error bars for each qubit
    plt.figure()
    plt.errorbar(lengs, avg_fids, yerr=rms_fids, fmt='-o', color='black')
    plt.axvline(x=max_len, linestyle="--", color="red")
    plt.text(max_len + 0.1, avg_fids[0], f'{max_len:.4f}', color='red')
    plt.xlabel('Readout and Pulse Length')
    plt.ylabel('Fidelity')
    plt.title(f'Avg Fidelity vs. Readout and Pulse Length for Qubit {QubitIndex + 1}, ({n_loops} repetitions)')
    plt.savefig(os.path.join(output_folder, f'fidelity_Q{QubitIndex + 1}_{timestamp}.png'), dpi=300)
    plt.close()

    del avg_fids, rms_fids, avg_ground_iq, avg_excited_iq, loop_group, length_group

########################################################################################################################
exit() #use this if you only want to run the first half of this script
optimal_lengths = [] #use when you are running the second half of this code separately
date_str = str(datetime.date.today())
outerFolder = f"/data/QICK_data/6transmon_run4a/{date_str}/readout_opt/Gain_Freq_Sweeps/"

# Ensure the output folder exists
os.makedirs(outerFolder, exist_ok=True)

# Reference frequencies for each resonator in MHz
res_freq_ge = [6191.419, 6216.1, 6292.361, 6405.77, 6432.759, 6468.481] #new
#res_freq_ge = [6191.439, 6216.0, 6292.261, 6405.79, 6432.899, 6468.501] old

# Define sweeping parameters
gain_range = [0.5,1]  # Gain range in a.u.
freq_steps = 30
gain_steps = 10

for QubitIndex in Qs:
    print(f'Starting Qubit {QubitIndex + 1} measurements.')
    # Select the reference frequency for the current qubit
    reference_frequency = res_freq_ge[QubitIndex]
    freq_range = [reference_frequency-1, reference_frequency + 1]  # Frequency range in MHz

    sweep = GainFrequencySweep(QubitIndex,list_of_all_qubits,  optimal_lengths=optimal_lengths, output_folder=outerFolder)
    results = sweep.run_sweep(freq_range, gain_range, freq_steps, gain_steps)
    results = np.array(results)

    # Save results and metadata in an HDF5 file
    timestamp = time.strftime("%H%M%S")
    h5_file = f"{outerFolder}Gain_Freq_Sweep_Qubit_{QubitIndex + 1}_{timestamp}.h5"
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

    plt.imshow(results, aspect='auto',
               extent=[gain_range[0], gain_range[1], freq_range[0] - reference_frequency, freq_range[1] - reference_frequency],
               origin='lower')
    plt.colorbar(label="Fidelity")
    plt.xlabel("Readout pulse gain (a.u.)")  # Gain on x-axis
    plt.ylabel("Readout frequency offset (MHz)")  # Frequency on y-axis
    plt.title(f"Gain-Frequency Sweep for Qubit {QubitIndex + 1}")
    #plt.show()
    # Save the plot
    file= f"{outerFolder}Gain_Freq_Sweep_Qubit_{QubitIndex + 1}_{timestamp}.png"
    plt.savefig(file, dpi=600, bbox_inches='tight')
    plt.close()  # Close the plot to free up memory
    del results, sweep
    print(f"Saved plot for Qubit {QubitIndex + 1} to {file}")
