import numpy as np
from section_002_res_spec_ge_mux import ResonanceSpectroscopy
from section_004_qubit_spec_ge import QubitSpectroscopy
from section_006_amp_rabi_ge import AmplitudeRabiExperiment
from section_007_T1_ge import T1Measurement
from section_008_save_data_to_h5 import Data_H5
from section_005_single_shot_ge import SingleShot
from section_009_T2R_ge import T2RMeasurement
from section_010_T2E_ge import T2EMeasurement
#from expt_config import *
import glob
import re
import datetime
import ast
import os
import sys
sys.path.append(os.path.abspath("/home/quietuser/Documents/GitHub/tprocv2_demos/qick_tprocv2_experiments_mux/"))


date = '2024-12-10'
outerFolder = "/data/QICK_data/6transmon_run5/" + date + "/"
outerFolder_save_plots = "/data/QICK_data/6transmon_run5/" + date + "_plots/"

save_figs = True
fit_saved = False
signal = 'None'
figure_quality = 100 #ramp this up to like 500 for presentation plots

loader_config_instance = Data_H5(outerFolder)
sys_config = loader_config_instance.load_config('sys_config.h5') #_batch2
del loader_config_instance

loader_config_instance = Data_H5(outerFolder)
exp_config = loader_config_instance.load_config('expt_cfg.h5') #_batch2
del loader_config_instance

import numpy as np
import h5py
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
import math


def calculate_qubit_temperature(frequency_mhz, ground_state_population, excited_state_population):
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    h = 6.62607015e-34  # Planck's constant in J·s
    frequency_hz = frequency_mhz * 1e6
    T = (h * frequency_hz) / (k_B * np.log(ground_state_population / excited_state_population))
    return T


def fit_double_gaussian_with_full_coverage(iq_data):
    gmm = GaussianMixture(n_components=2)
    gmm.fit(iq_data.reshape(-1, 1))

    means = gmm.means_.flatten()
    covariances = np.sqrt(gmm.covariances_).flatten()
    weights = gmm.weights_

    ground_gaussian = np.argmin(means)
    excited_gaussian = 1 - ground_gaussian

    # Generate x values to approximate the crossing point
    x_vals = np.linspace(means[ground_gaussian] - 3 * covariances[ground_gaussian],
                         means[excited_gaussian] + 3 * covariances[excited_gaussian], 1000)

    # Calculate Gaussian fits for each x value
    ground_gaussian_fit = weights[ground_gaussian] * (1 / (np.sqrt(2 * np.pi) * covariances[ground_gaussian])) * np.exp(
        -0.5 * ((x_vals - means[ground_gaussian]) / covariances[ground_gaussian]) ** 2)
    excited_gaussian_fit = weights[excited_gaussian] * (
                1 / (np.sqrt(2 * np.pi) * covariances[excited_gaussian])) * np.exp(
        -0.5 * ((x_vals - means[excited_gaussian]) / covariances[excited_gaussian]) ** 2)

    # Find the x value where the two Gaussian functions are closest
    crossing_point = x_vals[np.argmin(np.abs(ground_gaussian_fit - excited_gaussian_fit))]

    labels = gmm.predict(iq_data.reshape(-1, 1))

    ground_data = iq_data[(labels == ground_gaussian) & (iq_data < crossing_point)]
    excited_data = iq_data[(labels == excited_gaussian) & (iq_data > crossing_point)]

    ground_state_population = len(ground_data) / len(iq_data)
    excited_state_population_overlap = len(excited_data) / len(iq_data)

    return ground_state_population, excited_state_population_overlap, gmm, means, covariances, weights, crossing_point, ground_gaussian, excited_gaussian

def process_string_of_nested_lists(data):
    # Remove extra whitespace and non-numeric characters.
    data = re.sub(r'\s*\[(\s*.*?\s*)\]\s*', r'[\1]', data)
    data = data.replace('[ ', '[')
    data = data.replace('[ ', '[')
    data = data.replace('[ ', '[')

    cleaned_data = ''.join(c for c in data if c.isdigit() or c in ['-', '.', ' ', 'e', '[', ']'])
    pattern = r'\[(.*?)\]'  # Regular expression to match data within brackets
    matches = re.findall(pattern, cleaned_data)
    result = []
    for match in matches:
        numbers = [float(x.strip('[').strip(']').replace("'", "").replace(" ", "").replace("  ", "")) for x in match.split()] # Convert strings to integers
        result.append(numbers)

    return result


def process_h5_data(data):
    # Check if the data is a byte string; decode if necessary.
    if isinstance(data, bytes):
        data_str = data.decode()
    elif isinstance(data, str):
        data_str = data
    else:
        raise ValueError("Unsupported data type. Data should be bytes or string.")

    # Remove extra whitespace and non-numeric characters.
    cleaned_data = ''.join(c for c in data_str if c.isdigit() or c in ['-', '.', ' ', 'e'])

    # Split into individual numbers, removing empty strings.
    numbers = [float(x) for x in cleaned_data.split() if x]
    return numbers

def string_to_float_list(input_string):
    try:
        # Remove 'np.float64()' parts
        cleaned_string = input_string.replace('np.float64(', '').replace(')', '')

        # Use ast.literal_eval for safe evaluation
        float_list = ast.literal_eval(cleaned_string)

        # Check if all elements are floats (or can be converted to floats)
        return [float(x) for x in float_list]
    except (ValueError, SyntaxError, TypeError):
        print("Error: Invalid input string format.  It should be a string representation of a list of numbers.")
        return None

# ----------------------------------------------Load/Plot/Save QSpec------------------------------------
outerFolder_expt = outerFolder + "/Data_h5/QSpec_ge/"
h5_files = glob.glob(os.path.join(outerFolder_expt, "*.h5"))

for h5_file in h5_files:
    save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]
    H5_class_instance = Data_H5(h5_file)
    load_data = H5_class_instance.load_from_h5(data_type=  'QSpec', save_r = int(save_round))

    populated_keys = []
    for q_key in load_data['QSpec']:
        # Access 'Dates' for the current q_key
        dates_list = load_data['QSpec'][q_key].get('Dates', [[]])

        # Check if any entry in 'Dates' is not NaN
        if any(
                not np.isnan(date)
                for date in dates_list[0]  # Iterate over the first batch of dates
        ):
            populated_keys.append(q_key)

    print(f"Populated keys: {populated_keys}")

    for q_key in populated_keys:
        for dataset in range(len(load_data['QSpec'][q_key].get('Dates', [])[0])):
            date = datetime.datetime.fromtimestamp(load_data['QSpec'][q_key].get('Dates', [])[0][dataset])
            I = process_h5_data(load_data['QSpec'][q_key].get('I', [])[0][dataset].decode())
            Q = process_h5_data(load_data['QSpec'][q_key].get('Q', [])[0][dataset].decode())
            #I_fit = load_data['QSpec'][q_key].get('I Fit', [])[0][dataset]
            #Q_fit = load_data['QSpec'][q_key].get('Q Fit', [])[0][dataset]
            freqs = process_h5_data(load_data['QSpec'][q_key].get('Frequencies', [])[0][dataset].decode())
            round_num = load_data['QSpec'][q_key].get('Round Num', [])[0][dataset]
            batch_num = load_data['QSpec'][q_key].get('Batch Num', [])[0][dataset]

            # if len(I)>0:
            #
            #     qspec_class_instance = QubitSpectroscopy(q_key, outerFolder_save_plots, round_num, signal, save_figs)
            #     q_spec_cfg = ast.literal_eval(exp_config['qubit_spec_ge'].decode())
            #     qspec_class_instance.plot_results(I, Q, freqs, q_spec_cfg, figure_quality)
            #     del qspec_class_instance

    del H5_class_instance

# ------------------------------------------------Load/Plot/Save SS---------------------------------------
outerFolder_expt = outerFolder + "/Data_h5/SS_ge/"
h5_files = glob.glob(os.path.join(outerFolder_expt, "*.h5"))

for h5_file in h5_files:

    save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]
    H5_class_instance = Data_H5(h5_file)
    load_data = H5_class_instance.load_from_h5(data_type=  'SS', save_r = int(save_round))

    populated_keys = []
    for q_key in load_data['SS']:
        # Access 'Dates' for the current q_key
        dates_list = load_data['SS'][q_key].get('Dates', [[]])

        # Check if any entry in 'Dates' is not NaN
        if any(
                not np.isnan(date)
                for date in dates_list[0]  # Iterate over the first batch of dates
        ):
            populated_keys.append(q_key)

    print(f"Populated keys: {populated_keys}")

    for q_key in populated_keys:
        for dataset in range(len(load_data['SS'][q_key].get('Dates', [])[0])):
            date= datetime.datetime.fromtimestamp(load_data['SS'][q_key].get('Dates', [])[0][dataset])
            angle = load_data['SS'][q_key].get('Angle', [])[0][dataset]
            fidelity = load_data['SS'][q_key].get('Fidelity', [])[0][dataset]
            I_g = process_h5_data(load_data['SS'][q_key].get('I_g', [])[0][dataset].decode())
            Q_g = process_h5_data(load_data['SS'][q_key].get('Q_g', [])[0][dataset].decode())
            I_e = process_h5_data(load_data['SS'][q_key].get('I_e', [])[0][dataset].decode())
            Q_e = process_h5_data(load_data['SS'][q_key].get('Q_e', [])[0][dataset].decode())
            round_num = load_data['SS'][q_key].get('Round Num', [])[0][dataset]
            batch_num = load_data['SS'][q_key].get('Batch Num', [])[0][dataset]

            I_g = np.array(I_g)
            Q_g = np.array(Q_g)
            I_e = np.array(I_e)
            Q_e = np.array(Q_e)

            if len(Q_g)>0:
                ss_class_instance = SingleShot(q_key, outerFolder_save_plots, round_num, save_figs)
                ss_cfg = ast.literal_eval(exp_config['Readout_Optimization'].decode())
                #ss_class_instance.hist_ssf(data=[I_g, Q_g, I_e, Q_e], cfg=ss_cfg, plot=True)
                fid, threshold, rotation_angle, ig_new, ie_new = ss_class_instance.hist_ssf(
                    data=[I_g, Q_g, I_e, Q_e], cfg=ss_cfg, plot=False)

                qubit_frequency =



                del ss_class_instance

    del H5_class_instance
