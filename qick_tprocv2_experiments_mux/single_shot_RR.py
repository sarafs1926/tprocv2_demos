from section_005_single_shot_ge import SingleShot
from system_config import *
import numpy as np
import matplotlib.pyplot as plt
import os
import h5py
import datetime

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

# List of qubits and pulse lengths to measure
Qs = [0,1,2,3,4,5]  # already ran the first qubit, doing the other 5 now

lengs = np.linspace(0.5, 30, 30)  # Range of pulse lengths to measure

for QubitIndex in Qs:
    avg_fids = []
    rms_fids = []

    # Create a single HDF5 file for each qubit
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    h5_filename = os.path.join(output_folder, f"qubit_{QubitIndex + 1}_data_{timestamp}.h5")
    with h5py.File(h5_filename, 'w') as h5_file:
        # Top-level group for the qubit
        qubit_group = h5_file.create_group(f"Qubit_{QubitIndex + 1}")

        # Iterate over rounds
        #for j in range(1, n + 1): #rounds
        # Create a group for each round
        #round_group = qubit_group.create_group(f"Round_{j}")
        fids = []  # Store fidelity values for each loop
        ground_iq_data = []  # Store ground state IQ data for each loop
        excited_iq_data = []  # Store excited state IQ data for each loop

        # Iterate over each readout pulse length
        for leng in lengs:
            # Subgroup for each readout length within the round
            length_group = qubit_group.create_group(f"Length_{leng}")

            for k in range(n_loops): #loops for each read out length
                # ------------------------Single Shot-------------------------
                ss = SingleShot(QubitIndex, output_folder, k, round(leng, 3))
                fid, angle, iq_list_g, iq_list_e = ss.run(soccfg, soc)
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
    avg_max_index = avg_fids.index(avg_max)
    max_len = lengs[avg_max_index]

    # Plot the average fidelity vs. pulse length with error bars for each round
    plt.figure()
    plt.errorbar(lengs, avg_fids, yerr=rms_fids, fmt='-o', color='black')
    plt.axvline(x=max_len, linestyle="--", color="red")
    plt.text(max_len + 1, avg_fids[0], f'max length {max_len}')
    plt.xlabel('Readout and Pulse Length')
    plt.ylabel('Fidelity')
    plt.title(f'Avg Fidelity vs. Readout and Pulse Length for Qubit {QubitIndex + 1}, ({n_loops} repetitions)')
    plt.savefig(os.path.join(output_folder, f'fidelity_Q{QubitIndex + 1}_{timestamp}.png'), dpi=300)
    plt.close()

#print(f"Data saved for all qubits in {output_folder}")
