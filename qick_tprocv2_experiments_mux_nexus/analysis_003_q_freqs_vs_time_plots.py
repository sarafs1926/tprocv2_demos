import numpy as np

from section_002_res_spec_ge_mux import ResonanceSpectroscopy
from section_004_qubit_spec_ge_Franken import QubitSpectroscopy
from section_006_amp_rabi_ge import AmplitudeRabiExperiment
from section_007_T1_ge import T1Measurement
from section_008_save_data_to_h5 import Data_H5
from section_009_T2R_ge import T2RMeasurement
from section_010_T2E_ge import T2EMeasurement
#from expt_config import *
import glob
import re
import datetime
import ast
import os
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit

top_folder_dates = ['2024-11-21', '2024-11-21']
final_figure_quality = 50

#---------------------------------------get data--------------------------------
save_figs = False
fit_saved = False
signal = 'None'
figure_quality = 100 #ramp this up to like 500 for presentation plots

#---------definitions---------
def create_folder_if_not_exists(folder):
    """Creates a folder at the given path if it doesn't already exist."""
    if not os.path.exists(folder):
        os.makedirs(folder)

def exponential(x, a, b, c, d):
    return a * np.exp(-(x - b) / c) + d

def optimal_bins(data):
    n = len(data)
    if n == 0:
        return {}
    # Sturges' Rule
    sturges_bins = int(np.ceil(np.log2(n) + 1))
    return sturges_bins

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

# ----------Load/get data------------------------
qubit_frequencies = {i: [] for i in range(6)}
rounds = []
reps = []
file_names = []
date_times = {i: [] for i in range(6)}
mean_values = {}
show_legends = False

for folder_date in top_folder_dates:
    outerFolder = "/data/QICK_data/6transmon_run4a/" + folder_date + "/"
    outerFolder_save_plots = "/data/QICK_data/6transmon_run4a/" + folder_date + "_plots/"

    loader_config_instance = Data_H5(outerFolder)
    sys_config = loader_config_instance.load_config('sys_config_batch2.h5')
    del loader_config_instance

    loader_config_instance = Data_H5(outerFolder)
    exp_config = loader_config_instance.load_config('expt_cfg_batch2.h5')
    del loader_config_instance

    # ------------------------------------------Load/Plot/Save Q Spec------------------------------------
    outerFolder_expt = outerFolder + "/Data_h5/QSpec_ge/"
    h5_files = glob.glob(os.path.join(outerFolder_expt, "*.h5"))

    for h5_file in h5_files:
        save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]
        H5_class_instance = Data_H5(h5_file)
        load_data = H5_class_instance.load_from_h5(data_type='QSpec', save_r=int(save_round))

        for q_key in load_data['QSpec']:
            for dataset in range(len(load_data['QSpec'][q_key].get('Dates', [])[0])):
                date = datetime.datetime.fromtimestamp(load_data['QSpec'][q_key].get('Dates', [])[0][dataset])
                I = process_h5_data(load_data['QSpec'][q_key].get('I', [])[0][dataset].decode())
                Q = process_h5_data(load_data['QSpec'][q_key].get('Q', [])[0][dataset].decode())
                # I_fit = load_data['QSpec'][q_key].get('I Fit', [])[0][dataset]
                # Q_fit = load_data['QSpec'][q_key].get('Q Fit', [])[0][dataset]
                freqs = process_h5_data(load_data['QSpec'][q_key].get('Frequencies', [])[0][dataset].decode())
                round_num = load_data['QSpec'][q_key].get('Round Num', [])[0][dataset]
                batch_num = load_data['QSpec'][q_key].get('Batch Num', [])[0][dataset]

                if len(I) > 0:
                    qspec_class_instance = QubitSpectroscopy(q_key, outerFolder_save_plots, round_num, signal,
                                                             save_figs)
                    q_spec_cfg = ast.literal_eval(exp_config['qubit_spec_ge'].decode())
                    largest_amp_curve_mean, I_fit, Q_fit = qspec_class_instance.plot_results(I, Q, freqs, q_spec_cfg, figure_quality)

                    qubit_frequencies[q_key].extend([largest_amp_curve_mean])
                    date_times[q_key].extend([date.strftime("%Y-%m-%d %H:%M:%S")])

                    del qspec_class_instance

        del H5_class_instance



#---------------------------------plot-----------------------------------------------------
analysis_folder = "/data/QICK_data/6transmon_run4a/benchmark_analysis_plots/"
create_folder_if_not_exists(analysis_folder)
analysis_folder = "/data/QICK_data/6transmon_run4a/benchmark_analysis_plots/features_vs_time/"
create_folder_if_not_exists(analysis_folder)

font = 14
titles = [f"Qubit {i+1}" for i in range(6)]
colors = ['orange','blue','purple','green','brown','pink']
fig, axes = plt.subplots(2, 3, figsize=(12, 8))
plt.title('Qubit Frequencies vs Time',fontsize = font)
axes = axes.flatten()
titles = [f"Qubit {i + 1}" for i in range(6)]
from datetime import datetime
for i, ax in enumerate(axes):

    ax.set_title(titles[i], fontsize = font)

    x = date_times[i]
    y = qubit_frequencies[i]

    # Convert strings to datetime objects.
    datetime_objects = [datetime.strptime(date_string, "%Y-%m-%d %H:%M:%S") for date_string in x]

    # Combine datetime objects and y values into a list of tuples and sort by datetime.
    combined = list(zip(datetime_objects, y))
    combined.sort(reverse=True, key=lambda x: x[0])

    # Unpack them back into separate lists, in order from latest to most recent.
    sorted_x, sorted_y = zip(*combined)
    ax.scatter(sorted_x, sorted_y, color=colors[i])

    sorted_x = np.asarray(sorted(x))

    num_points = 5
    indices = np.linspace(0, len(sorted_x) - 1, num_points, dtype=int)

    # Set new x-ticks using the datetime objects at the selected indices
    ax.set_xticks(sorted_x[indices])
    ax.set_xticklabels([dt for dt in sorted_x[indices]], rotation=45)

    ax.scatter(x, y, color=colors[i])
    if show_legends:
        ax.legend(edgecolor='black')
    ax.set_xlabel('Time (Days)', fontsize=font-2)
    ax.set_ylabel('Qubit Frequency (MHz)', fontsize=font-2)
    ax.tick_params(axis='both', which='major', labelsize=8)

plt.tight_layout()
plt.savefig(analysis_folder + 'Q_Freqs.png', transparent=True, dpi=final_figure_quality)

#plt.show()
