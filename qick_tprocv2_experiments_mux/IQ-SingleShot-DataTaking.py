import numpy as np
import h5py
import os
from section_005_single_shot_ge import SingleShot
from system_config import QICK_experiment
import datetime
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt


# Output folder configuration
date_str = str(datetime.date.today())
output_folder = f"/data/QICK_data/6transmon_run4a/{date_str}/Qubit_Temps/"
os.makedirs(output_folder, exist_ok=True)
tot_num_of_qubits = 6

list_of_all_qubits = list(range(tot_num_of_qubits))

# Optimized parameters for qubits 3 and 4

#For ADC 10 and DAC 2,10
# qubits = {
#     3: {"length": 2.5, "gain": 0.8, "freq_offset": -0.2667, "reference_freq": 6292.361, "qubit_freq": 4156.53},
#     4: {"length": 4.5, "gain": 0.95, "freq_offset": -0.1333, "reference_freq": 6405.77, "qubit_freq": 4459.20}
# }
#For ADC 20 and DAC 2,10
# qubits = {
#     3: {"length": 2.25, "gain": 0.8, "freq_offset": -0.2000, "reference_freq": 6292.361, "qubit_freq": 4156.53},
#     4: {"length": 4.0, "gain": 0.95, "freq_offset": -0.2667, "reference_freq": 6405.77, "qubit_freq": 4459.20}
# }

qubits = { #Note: gain and freq_offset were not optimized for this configuration, only readout length was optimized.
     3: {"length": 2.5, "gain": 1.0, "freq_offset": -0.0, "reference_freq": 6292.361, "qubit_freq": 4156.53},
     4: {"length": 4.0, "gain": 1.0, "freq_offset": -0.0, "reference_freq": 6405.77, "qubit_freq": 4459.20}
 }

# Loop through specified qubits and apply settings
for qubit_index, params in qubits.items():
    length = params["length"]
    gain = params["gain"]
    frequency = params["reference_freq"] + params["freq_offset"] #resonator offset frequency
    qubit_frequency = params["qubit_freq"]

    print(f"Processing Qubit {qubit_index} with Res_Length {length}, Res_Gain {gain}, Res Freq {frequency}, and Qubit Frequency {qubit_frequency}")
#For ADC 30 and DAC 2,10
    # Initialize experiment
    QubitIndex = qubit_index - 1
    #experiment = QICK_experiment(output_folder)
    experiment = QICK_experiment(output_folder, DAC_attenuator1=2, DAC_attenuator2=10, ADC_attenuator=30)
    experiment.readout_cfg['res_length'] = length
    experiment.readout_cfg['res_freq_ge'][QubitIndex] = frequency #resonator freq
    res_gains = experiment.mask_gain_res(QubitIndex, gain)
    experiment.readout_cfg['res_gain_ge'] = res_gains

    # Run the single-shot experiment
    ss = SingleShot(QubitIndex,list_of_all_qubits, output_folder, experiment, round_num=0, save_figs=False)
    fid, angle, iq_list_g, iq_list_e = ss.run(experiment.soccfg, experiment.soc)

    I_g = iq_list_g[QubitIndex][0].T[0]  # Ground-state I data
    Q_g = iq_list_g[QubitIndex][0].T[1]
    I_e = iq_list_e[QubitIndex][0].T[0]
    Q_e = iq_list_e[QubitIndex][0].T[1]

    # Call hist_ssf to get rotated I data
    fid, threshold, rotation_angle, ig_new, ie_new = ss.hist_ssf(
        data=[I_g, Q_g, I_e, Q_e], cfg=ss.config, plot=True)

    # Save results in HDF5
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    h5_filename = os.path.join(output_folder, f"Qubit_{qubit_index}_temperatureSSdata_{timestamp}.h5")
    with h5py.File(h5_filename, 'w') as f:
        f.create_dataset("fidelity", data=fid)
        f.create_dataset("gain", data=gain)
        f.create_dataset("offset_res_frequency", data=frequency)
        f.create_dataset("qubit_frequency", data=qubit_frequency)
        f.create_dataset("ground_iq_data", data=iq_list_g)
        f.create_dataset("excited_iq_data", data=iq_list_e)
        f.create_dataset("ig_new", data=ig_new)  # Save rotated ground-state I data
        f.create_dataset("ie_new", data=ie_new)  # Save rotated excited-state I data
        f.attrs['rotation_angle'] = rotation_angle

    print(f"Data saved for Qubit {qubit_index} to {h5_filename}")
    del experiment
