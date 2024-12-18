from section_005_single_shot_ge import SingleShot
#from system_config import *
from system_config import QICK_experiment
import numpy as np
import matplotlib.pyplot as plt
import os
import h5py
import datetime
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
output_folder = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{prefix}/SingleShot_Test" #"/data/QICK_data/6transmon_run4a/" + prefix + "/SingleShot_Test/"
create_folder_if_not_exists(output_folder)

n = 1  # Number of rounds
n_loops = 5  # Number of repetitions per length to average

# List of qubits and pulse lengths to measure
Qs = [0]

lengs = np.linspace(0.5, 5, 19) # increments of 0.25

for QubitIndex in Qs:
    QubitIndex = int(QubitIndex)  # Ensure QubitIndex is an integer

    avg_fids = []
    rms_fids = []

    # Create a single HDF5 file for each qubit
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    #print("Output folder:", output_folder)
    #print("Type of output_folder:", type(output_folder))
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
                # Initialize experiment for each loop iteration
                experiment = QICK_experiment(output_folder)
                # Set specific configuration values for each iteration
                experiment.readout_cfg['res_length'] = leng  # Set the current readout pulse length

                # Set gain for the current qubit
                gain = 1
                res_gains = experiment.mask_gain_res(QubitIndex, gain)  # Set gain for current qubit only

                experiment.readout_cfg['res_gain_ge'] = res_gains

                """
                exp_cfg = expt_cfg["Readout_Optimization"]
                q_config = all_qubit_state(experiment)
                Qubit = 'Q' + str(QubitIndex)
                
                # Combine specific qubit configuration with experiment-specific settings
                config = {**q_config[Qubit], **exp_cfg}"""
                #print(f"Single Shot configuration:", config)

                #ss = SingleShot(QubitIndex, output_folder, k, round(leng, 3)) #Old way
                ss = SingleShot(QubitIndex, output_folder, experiment, round_num=k, save_figs=True)#New way
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


    #avg_max = max(avg_fids[:10])
    avg_max = max(avg_fids)
    avg_max_index = avg_fids.index(avg_max)
    max_len = lengs[avg_max_index]

    # Plot the average fidelity vs. pulse length with error bars for each round
    plt.figure()
    plt.errorbar(lengs, avg_fids, yerr=rms_fids, fmt='-o', color='black')
    plt.axvline(x=max_len, linestyle="--", color="red")
    plt.text(max_len + 0.5, avg_fids[0], f'max length {max_len:.4f}', color='red')
    plt.xlabel('Readout and Pulse Length')
    plt.ylabel('Fidelity')
    plt.title(f'Avg Fidelity vs. Readout and Pulse Length for Qubit {QubitIndex + 1}, ({n_loops} repetitions)')
    plt.savefig(os.path.join(output_folder, f'fidelity_Q{QubitIndex + 1}_{timestamp}.png'), dpi=300)
    plt.close()

#print(f"Data saved for all qubits in {output_folder}")
