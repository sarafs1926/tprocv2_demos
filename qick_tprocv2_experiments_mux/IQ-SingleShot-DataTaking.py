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
    res_gains = experiment.set_gain_filter_ge(QubitIndex, gain)
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
