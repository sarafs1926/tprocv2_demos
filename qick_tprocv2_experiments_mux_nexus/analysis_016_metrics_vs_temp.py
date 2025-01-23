import numpy as np
import os
import sys
sys.path.append(os.path.abspath("/home/nexusadmin/Documents/GitHub/tprocv2_demos/qick_tprocv2_experiments_mux_nexus/"))

from section_002_res_spec_ge_mux import ResonanceSpectroscopy
from section_004_qubit_spec_ge import QubitSpectroscopy
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
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit
import pandas as pd

class GetThermData:
    def __init__(self, folder):
        self.folder = folder

    def run(self):
        from datetime import datetime
        unix_dates = []
        mcp2_temps = []
        magcan_temps = []

        for file_name in os.listdir(self.folder):
            if file_name.endswith('.csv'):
                file_path = os.path.join(self.folder, file_name)
                data = pd.read_csv(file_path)

                if {'Unix Date', 'MCP2 Temp (mK)', 'Mag Can Temp (mK)'}.issubset(data.columns):
                    unix_dates.extend(data['Unix Date'].tolist())
                    mcp2_temps.extend(data['MCP2 Temp (mK)'].tolist())
                    magcan_temps.extend(data['Mag Can Temp (mK)'].tolist())
        datetime_dates = []
        for i in unix_dates:
            datetime_dates.append(datetime.utcfromtimestamp(i))
        return datetime_dates, mcp2_temps, magcan_temps

class ResonatorFreqVsTemp:
    def __init__(self, figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs, fit_saved,
                 signal, run_name, exp_config):
        self.figure_quality = figure_quality
        self.number_of_qubits = number_of_qubits
        self.save_figs = save_figs
        self.fit_saved = fit_saved
        self.signal = signal
        self.run_name = run_name
        self.top_folder_dates = top_folder_dates
        self.final_figure_quality = final_figure_quality
        self.exp_config = exp_config

    def datetime_to_unix(self, dt):
        # Convert to Unix timestamp
        unix_timestamp = int(dt.timestamp())
        return unix_timestamp

    def create_folder_if_not_exists(self, folder):
        """Creates a folder at the given path if it doesn't already exist."""
        if not os.path.exists(folder):
            os.makedirs(folder)

    def exponential(self, x, a, b, c, d):
        return a * np.exp(-(x - b) / c) + d

    def optimal_bins(self, data):
        n = len(data)
        if n == 0:
            return {}
        # Sturges' Rule
        sturges_bins = int(np.ceil(np.log2(n) + 1))
        return sturges_bins

    def process_string_of_nested_lists(self, data):
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


    def process_h5_data(self, data):
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

    def string_to_float_list(self, input_string):
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

    def run(self):
        import datetime
        # ----------Load/get data------------------------
        resonator_centers = {i: [] for i in range(self.number_of_qubits)}
        rounds = []
        reps = []
        file_names = []
        date_times = {i: [] for i in range(self.number_of_qubits)}
        mean_values = {}

        for folder_date in self.top_folder_dates:
            outerFolder = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/" + folder_date + "/"
            outerFolder_save_plots = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/" + folder_date + "_plots/"

            # ------------------------------------------Load/Plot/Save Res Spec------------------------------------
            outerFolder_expt = outerFolder + "/Data_h5/Res_ge/"
            h5_files = glob.glob(os.path.join(outerFolder_expt, "*.h5"))

            for h5_file in h5_files:
                save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]
                H5_class_instance = Data_H5(h5_file)
                # H5_class_instance.print_h5_contents(h5_file)
                load_data = H5_class_instance.load_from_h5(data_type='Res', save_r=int(save_round))

                # just look at this resonator data, should have batch_num of arrays in each one
                # right now the data writes the same thing batch_num of times, so it will do the same 5 datasets 5 times, until you fix this just grab the first one (All 5)
                for q_key in load_data['Res']:
                    # print("all batch_num datasets------------------------", load_data['Res'][q_key].get('Amps', [])[0])
                    # print("one dataset------------------------",load_data['Res'][q_key].get('Amps', [])[0][0].decode())
                    # go through each dataset in the batch and plot
                    for dataset in range(len(load_data['Res'][q_key].get('Dates', [])[0])):
                        if 'nan' in str(load_data['Res'][q_key].get('Dates', [])[0][dataset]):
                            continue

                        date = datetime.datetime.fromtimestamp(
                            load_data['Res'][q_key].get('Dates', [])[0][dataset])  # single date per dataset

                        freq_pts = self.process_h5_data(load_data['Res'][q_key].get('freq_pts', [])[0][
                                                       dataset].decode())  # comes in as an array but put into a byte string, need to convert to list
                        freq_center = self.process_h5_data(load_data['Res'][q_key].get('freq_center', [])[0][
                                                          dataset].decode())  # comes in as an array but put into a string, need to convert to list
                        freqs_found = self.string_to_float_list(load_data['Res'][q_key].get('Found Freqs', [])[0][
                                                               dataset].decode())  # comes in as a list of floats in string format, need to convert
                        amps = self.process_string_of_nested_lists(
                            load_data['Res'][q_key].get('Amps', [])[0][dataset].decode())  # list of lists
                        round_num = load_data['Res'][q_key].get('Round Num', [])[0][dataset]  # already a float
                        batch_num = load_data['Res'][q_key].get('Batch Num', [])[0][dataset]

                        if len(freq_pts) > 0:
                            res_class_instance = ResonanceSpectroscopy(q_key, outerFolder_save_plots, round_num, self.save_figs)
                            res_spec_cfg = ast.literal_eval(self.exp_config['res_spec'].decode())
                            res_freqs = res_class_instance.get_results(freq_pts, freq_center, amps)

                            resonator_centers[q_key].extend([res_freqs[q_key]])
                            date_times[q_key].extend([date.strftime("%Y-%m-%d %H:%M:%S")])

                            del res_class_instance

                del H5_class_instance
        return date_times, resonator_centers

    def plot(self, date_times, resonator_centers, mcp2_dates, mcp2_temps, show_legends):
        #---------------------------------plot-----------------------------------------------------
        analysis_folder = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/benchmark_analysis_plots/"
        self.create_folder_if_not_exists(analysis_folder)
        analysis_folder = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/benchmark_analysis_plots/features_vs_time/"
        self.create_folder_if_not_exists(analysis_folder)

        font = 14
        colors = ['orange','blue','purple','green','brown','pink']
        fig, axes = plt.subplots(2, 3, figsize=(12, 8))
        plt.title('Resonator centers vs Time',fontsize = font)
        axes = axes.flatten()
        titles = [f"Res {i + 1}" for i in range(self.number_of_qubits)]
        from datetime import datetime
        for i, ax in enumerate(axes):

            ax.set_title(titles[i], fontsize=font)

            x = date_times[i]
            y = resonator_centers[i]

            datetime_objects = [datetime.strptime(date_string, "%Y-%m-%d %H:%M:%S") for date_string in x]

            combined = list(zip(datetime_objects, y))
            combined.sort(reverse=True, key=lambda item: item[0])
            sorted_x, sorted_y = zip(*combined)

            sorted_x = np.array(sorted_x)
            sorted_y = np.array(sorted_y)

            combined_mcp2 = list(zip(mcp2_dates, mcp2_temps))
            combined_mcp2.sort(reverse=True, key=lambda item: item[0])
            sorted_x_mcp2, sorted_y_mcp2 = zip(*combined_mcp2)

            sorted_x_mcp2 = np.array(sorted_x_mcp2)
            sorted_y_mcp2 = np.array(sorted_y_mcp2)

            ax.scatter(sorted_x, sorted_y, color=colors[i])

            num_points = 5
            indices = np.linspace(0, len(sorted_x) - 1, num_points, dtype=int)
            ax.set_xticks(sorted_x[indices])
            ax.set_xticklabels([dt for dt in sorted_x[indices]], rotation=45)

            min_datetime = sorted_x.min()
            max_datetime = sorted_x.max()

            valid_indices = (sorted_x_mcp2 >= min_datetime) & (sorted_x_mcp2 <= max_datetime)
            filtered_x_mcp2 = sorted_x_mcp2[valid_indices]
            filtered_y_mcp2 = sorted_y_mcp2[valid_indices]

            ax.set_xlabel('Time (Days)', fontsize=font - 2)
            ax.set_ylabel('Resonator Center (MHz)', fontsize=font - 2)
            ax.tick_params(axis='both', which='major', labelsize=8)

            ax2 = ax.twinx()
            ax2.scatter(filtered_x_mcp2, filtered_y_mcp2, color='darkblue')
            ax2.set_ylabel('Temperature (mK)', fontsize=font - 2, color='darkblue')
            ax2.tick_params(axis='y', labelcolor='darkblue', labelsize=10)
            
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(analysis_folder + 'Res_Centers_vs_time_thermometry.pdf', transparent=True, dpi=self.final_figure_quality)

        #plt.show()
class QubitFreqsVsTemp:
    def __init__(self, figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs, fit_saved,
                 signal, run_name, exp_config):
        self.save_figs = save_figs
        self.fit_saved = fit_saved
        self.signal = signal
        self.figure_quality = figure_quality
        self.run_name = run_name
        self.number_of_qubits = number_of_qubits
        self.final_figure_quality = final_figure_quality
        self.top_folder_dates = top_folder_dates
        self.exp_config = exp_config

    def datetime_to_unix(self, dt):
        # Convert to Unix timestamp
        unix_timestamp = int(dt.timestamp())
        return unix_timestamp

    def unix_to_datetime(self, unix_timestamp):
        # Convert the Unix timestamp to a datetime object
        dt = datetime.fromtimestamp(unix_timestamp)
        return dt

    def create_folder_if_not_exists(self, folder):
        """Creates a folder at the given path if it doesn't already exist."""
        if not os.path.exists(folder):
            os.makedirs(folder)

    def exponential(self, x, a, b, c, d):
        return a * np.exp(-(x - b) / c) + d

    def optimal_bins(self, data):
        n = len(data)
        if n == 0:
            return {}
        # Sturges' Rule
        sturges_bins = int(np.ceil(np.log2(n) + 1))
        return sturges_bins

    def process_string_of_nested_lists(self, data):
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

    def process_h5_data(self, data):
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

    def string_to_float_list(self, input_string):
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

    def run(self):
        import datetime

        qubit_frequencies = {i: [] for i in range(self.number_of_qubits)}
        rounds = []
        reps = []
        file_names = []
        date_times = {i: [] for i in range(self.number_of_qubits)}
        mean_values = {}
        for folder_date in self.top_folder_dates:
            outerFolder = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/" + folder_date + "/"
            outerFolder_save_plots = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/" + folder_date + "_plots/"

            # ------------------------------------------Load/Plot/Save Q Spec------------------------------------
            outerFolder_expt = outerFolder + "/Data_h5/QSpec_ge/"
            h5_files = glob.glob(os.path.join(outerFolder_expt, "*.h5"))

            for h5_file in h5_files:
                save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]

                H5_class_instance = Data_H5(h5_file)
                load_data = H5_class_instance.load_from_h5(data_type='QSpec', save_r=int(save_round))

                for q_key in load_data['QSpec']:
                    for dataset in range(len(load_data['QSpec'][q_key].get('Dates', [])[0])):
                        if 'nan' in str(load_data['QSpec'][q_key].get('Dates', [])[0][dataset]):
                            continue
                        date = datetime.datetime.fromtimestamp(load_data['QSpec'][q_key].get('Dates', [])[0][dataset])
                        I = self.process_h5_data(load_data['QSpec'][q_key].get('I', [])[0][dataset].decode())
                        Q = self.process_h5_data(load_data['QSpec'][q_key].get('Q', [])[0][dataset].decode())
                        # I_fit = load_data['QSpec'][q_key].get('I Fit', [])[0][dataset]
                        # Q_fit = load_data['QSpec'][q_key].get('Q Fit', [])[0][dataset]
                        freqs = self.process_h5_data(load_data['QSpec'][q_key].get('Frequencies', [])[0][dataset].decode())
                        round_num = load_data['QSpec'][q_key].get('Round Num', [])[0][dataset]
                        batch_num = load_data['QSpec'][q_key].get('Batch Num', [])[0][dataset]

                        if len(I) > 0:
                            qspec_class_instance = QubitSpectroscopy(q_key, outerFolder_save_plots, round_num, self.signal,
                                                                     self.save_figs)
                            q_spec_cfg = ast.literal_eval(self.exp_config['qubit_spec_ge'].decode())
                            largest_amp_curve_mean, I_fit, Q_fit = qspec_class_instance.get_results(I, Q, freqs)

                            qubit_frequencies[q_key].extend([largest_amp_curve_mean])
                            date_times[q_key].extend([date.strftime("%Y-%m-%d %H:%M:%S")])

                            del qspec_class_instance

                del H5_class_instance
        return date_times, qubit_frequencies

    def plot(self,date_times, qubit_frequencies, mcp2_dates, mcp2_temps,show_legends):
        #---------------------------------plot-----------------------------------------------------
        analysis_folder = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/benchmark_analysis_plots/"
        self.create_folder_if_not_exists(analysis_folder)
        analysis_folder = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/benchmark_analysis_plots/features_vs_time/"
        self.create_folder_if_not_exists(analysis_folder)

        font = 14
        titles = [f"Qubit {i+1}" for i in range(self.number_of_qubits)]
        colors = ['orange','blue','purple','green','brown','pink']
        fig, axes = plt.subplots(2, 3, figsize=(12, 8))
        plt.title('Qubit Frequencies vs Time',fontsize = font)
        axes = axes.flatten()
        titles = [f"Qubit {i + 1}" for i in range(self.number_of_qubits)]
        from datetime import datetime
        for i, ax in enumerate(axes):

            ax.set_title(titles[i], fontsize = font)

            x = date_times[i]
            y = qubit_frequencies[i]

            datetime_objects = [datetime.strptime(date_string, "%Y-%m-%d %H:%M:%S") for date_string in x]

            combined = list(zip(datetime_objects, y))
            combined.sort(reverse=True, key=lambda item: item[0])
            sorted_x, sorted_y = zip(*combined)

            sorted_x = np.array(sorted_x)
            sorted_y = np.array(sorted_y)

            combined_mcp2 = list(zip(mcp2_dates, mcp2_temps))
            combined_mcp2.sort(reverse=True, key=lambda item: item[0])
            sorted_x_mcp2, sorted_y_mcp2 = zip(*combined_mcp2)

            sorted_x_mcp2 = np.array(sorted_x_mcp2)
            sorted_y_mcp2 = np.array(sorted_y_mcp2)

            ax.scatter(sorted_x, sorted_y, color=colors[i])

            num_points = 5
            indices = np.linspace(0, len(sorted_x) - 1, num_points, dtype=int)
            ax.set_xticks(sorted_x[indices])
            ax.set_xticklabels([dt for dt in sorted_x[indices]], rotation=45)

            min_datetime = sorted_x.min()
            max_datetime = sorted_x.max()

            valid_indices = (sorted_x_mcp2 >= min_datetime) & (sorted_x_mcp2 <= max_datetime)
            filtered_x_mcp2 = sorted_x_mcp2[valid_indices]
            filtered_y_mcp2 = sorted_y_mcp2[valid_indices]


            ax.set_xlabel('Time (Days)', fontsize=font-2)
            ax.set_ylabel('Qubit Frequency (MHz)', fontsize=font-2)
            ax.tick_params(axis='both', which='major', labelsize=8)

            ax2 = ax.twinx()
            ax2.scatter(filtered_x_mcp2, filtered_y_mcp2, color='darkblue')
            ax2.set_ylabel('Temperature (mK)', fontsize=font - 2, color='darkblue')
            ax2.tick_params(axis='y', labelcolor='darkblue', labelsize=10)

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])

        plt.savefig(analysis_folder + 'Q_Freqs_vs_time_thermometry.pdf', transparent=True, dpi=self.final_figure_quality)

        #plt.show()

class PiAmpsVsTemp:
    def __init__(self, figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs, fit_saved,
                 signal, run_name, exp_config):
        self.save_figs = save_figs
        self.fit_saved = fit_saved
        self.signal = signal
        self.figure_quality = figure_quality
        self.run_name = run_name
        self.number_of_qubits = number_of_qubits
        self.final_figure_quality = final_figure_quality
        self.top_folder_dates = top_folder_dates
        self.exp_config = exp_config

    def datetime_to_unix(self, dt):
        # Convert to Unix timestamp
        unix_timestamp = int(dt.timestamp())
        return unix_timestamp

    def unix_to_datetime(self, unix_timestamp):
        # Convert the Unix timestamp to a datetime object
        dt = datetime.fromtimestamp(self, unix_timestamp)
        return dt

    def create_folder_if_not_exists(self, folder):
        """Creates a folder at the given path if it doesn't already exist."""
        if not os.path.exists(folder):
            os.makedirs(folder)

    def exponential(self, x, a, b, c, d):
        return a * np.exp(-(x - b) / c) + d

    def optimal_bins(self, data):
        n = len(data)
        if n == 0:
            return {}
        # Sturges' Rule
        sturges_bins = int(np.ceil(np.log2(n) + 1))
        return sturges_bins

    def process_string_of_nested_lists(self, data):
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

    def process_h5_data(self, data):
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

    def string_to_float_list(self, input_string):
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

    def run(self, plot_depths=False):
        import datetime
        # ----------Load/get data------------------------
        pi_amps = {i: [] for i in range(self.number_of_qubits)}
        depths = {i: [] for i in range(self.number_of_qubits)}
        rounds = []
        reps = []
        file_names = []
        date_times = {i: [] for i in range(self.number_of_qubits)}
        mean_values = {}
        for folder_date in self.top_folder_dates:
            outerFolder = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/" + folder_date + "/"
            outerFolder_save_plots = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/" + folder_date + "_plots/"

            # ------------------------------------------------Load/Plot/Save Rabi---------------------------------------
            outerFolder_expt = outerFolder + "/Data_h5/Rabi_ge/"
            h5_files = glob.glob(os.path.join(outerFolder_expt, "*.h5"))
            for h5_file in h5_files:
                save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]
                H5_class_instance = Data_H5(h5_file)
                load_data = H5_class_instance.load_from_h5(data_type='Rabi', save_r=int(save_round))

                for q_key in load_data['Rabi']:
                    for dataset in range(len(load_data['Rabi'][q_key].get('Dates', [])[0])):
                        if 'nan' in str(load_data['Rabi'][q_key].get('Dates', [])[0][dataset]):
                            continue
                        date = datetime.datetime.fromtimestamp(load_data['Rabi'][q_key].get('Dates', [])[0][dataset])
                        I = self.process_h5_data(load_data['Rabi'][q_key].get('I', [])[0][dataset].decode())
                        Q = self.process_h5_data(load_data['Rabi'][q_key].get('Q', [])[0][dataset].decode())
                        gains = self.process_h5_data(load_data['Rabi'][q_key].get('Gains', [])[0][dataset].decode())
                        # fit = load_data['Rabi'][q_key].get('Fit', [])[0][dataset]
                        round_num = load_data['Rabi'][q_key].get('Round Num', [])[0][dataset]
                        batch_num = load_data['Rabi'][q_key].get('Batch Num', [])[0][dataset]

                        if len(I) > 0:
                            rabi_class_instance = AmplitudeRabiExperiment(q_key, outerFolder_save_plots, round_num,
                                                                          self.signal, self.save_figs)
                            rabi_cfg = ast.literal_eval(self.exp_config['power_rabi_ge'].decode())
                            I = np.asarray(I)
                            Q = np.asarray(Q)
                            gains = np.asarray(gains)
                            if plot_depths:
                                best_signal_fit, pi_amp, depth = rabi_class_instance.get_results(I, Q, gains,
                                                                                                 grab_depths=True)
                                depths[q_key].extend([depth])
                            else:
                                best_signal_fit, pi_amp = rabi_class_instance.get_results(I, Q, gains)

                            pi_amps[q_key].extend([pi_amp])
                            date_times[q_key].extend([date.strftime("%Y-%m-%d %H:%M:%S")])

                            del rabi_class_instance
                del H5_class_instance
        if plot_depths:
            return date_times, pi_amps, depths
        else:
            return date_times, pi_amps

    def plot(self, date_times, pi_amps, mcp2_dates, mcp2_temps, show_legends):
        #---------------------------------plot-----------------------------------------------------
        analysis_folder = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/benchmark_analysis_plots/"
        self.create_folder_if_not_exists(analysis_folder)
        analysis_folder = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/benchmark_analysis_plots/features_vs_time/"
        self.create_folder_if_not_exists(analysis_folder)

        font = 14
        titles = [f"Qubit {i+1}" for i in range(self.number_of_qubits)]
        colors = ['orange','blue','purple','green','brown','pink']
        fig, axes = plt.subplots(2, 3, figsize=(12, 8))
        plt.title('Pi Amplitudes vs Time',fontsize = font)
        axes = axes.flatten()
        titles = [f"Qubit {i + 1}" for i in range(self.number_of_qubits)]
        from datetime import datetime
        for i, ax in enumerate(axes):

            ax.set_title(titles[i], fontsize = font)

            x = date_times[i]
            y = pi_amps[i]

            datetime_objects = [datetime.strptime(date_string, "%Y-%m-%d %H:%M:%S") for date_string in x]

            combined = list(zip(datetime_objects, y))
            combined.sort(reverse=True, key=lambda item: item[0])
            sorted_x, sorted_y = zip(*combined)

            sorted_x = np.array(sorted_x)
            sorted_y = np.array(sorted_y)

            combined_mcp2 = list(zip(mcp2_dates, mcp2_temps))
            combined_mcp2.sort(reverse=True, key=lambda item: item[0])
            sorted_x_mcp2, sorted_y_mcp2 = zip(*combined_mcp2)

            sorted_x_mcp2 = np.array(sorted_x_mcp2)
            sorted_y_mcp2 = np.array(sorted_y_mcp2)

            ax.scatter(sorted_x, sorted_y, color=colors[i])

            num_points = 5
            indices = np.linspace(0, len(sorted_x) - 1, num_points, dtype=int)
            ax.set_xticks(sorted_x[indices])
            ax.set_xticklabels([dt for dt in sorted_x[indices]], rotation=45)

            min_datetime = sorted_x.min()
            max_datetime = sorted_x.max()

            valid_indices = (sorted_x_mcp2 >= min_datetime) & (sorted_x_mcp2 <= max_datetime)
            filtered_x_mcp2 = sorted_x_mcp2[valid_indices]
            filtered_y_mcp2 = sorted_y_mcp2[valid_indices]

            ax.set_xlabel('Time (Days)', fontsize=font-2)
            ax.set_ylabel('Pi Amp (a.u.)', fontsize=font-2)
            ax.tick_params(axis='both', which='major', labelsize=8)

            ax2 = ax.twinx()
            ax2.scatter(filtered_x_mcp2, filtered_y_mcp2, color='darkblue')
            ax2.set_ylabel('Temperature (mK)', fontsize=font - 2, color='darkblue')
            ax2.tick_params(axis='y', labelcolor='darkblue', labelsize=10)

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(analysis_folder + 'Pi_Amps_vs_time_thermometry.pdf', transparent=True, dpi=self.final_figure_quality)

        #plt.show()

class T1VsTemp:
    def __init__(self, figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs, fit_saved,
                 signal, run_name, exp_config):
        self.save_figs = save_figs
        self.fit_saved = fit_saved
        self.signal = signal
        self.figure_quality = figure_quality
        self.run_name = run_name
        self.number_of_qubits = number_of_qubits
        self.final_figure_quality = final_figure_quality
        self.top_folder_dates = top_folder_dates
        self.exp_config = exp_config

    def datetime_to_unix(self, dt):
        # Convert to Unix timestamp
        unix_timestamp = int(dt.timestamp())
        return unix_timestamp

    def unix_to_datetime(self, unix_timestamp):
        # Convert the Unix timestamp to a datetime object
        dt = datetime.fromtimestamp(unix_timestamp)
        return dt

    def create_folder_if_not_exists(self, folder):
        """Creates a folder at the given path if it doesn't already exist."""
        if not os.path.exists(folder):
            os.makedirs(folder)

    def exponential(self, x, a, b, c, d):
        return a * np.exp(-(x - b) / c) + d

    def optimal_bins(self, data):
        n = len(data)
        if n == 0:
            return {}
        # Sturges' Rule
        sturges_bins = int(np.ceil(np.log2(n) + 1))
        return sturges_bins

    def process_string_of_nested_lists(self, data):
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


    def process_h5_data(self, data):
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

    def string_to_float_list(self, input_string):
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

    def run(self):
        import datetime

        # ----------Load/get data------------------------
        t1_vals = {i: [] for i in range(self.number_of_qubits)}
        rounds = []
        reps = []
        file_names = []
        date_times = {i: [] for i in range(self.number_of_qubits)}
        mean_values = {}

        for folder_date in self.top_folder_dates:
            outerFolder = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/" + folder_date + "/"
            outerFolder_save_plots = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/" + folder_date + "_plots/"

            # ------------------------------------------------Load/Plot/Save T1----------------------------------------------
            outerFolder_expt = outerFolder + "/Data_h5/T1_ge/"
            h5_files = glob.glob(os.path.join(outerFolder_expt, "*.h5"))

            for h5_file in h5_files:

                save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]
                H5_class_instance = Data_H5(h5_file)
                load_data = H5_class_instance.load_from_h5(data_type='T1', save_r=int(save_round))

                for q_key in load_data['T1']:
                    for dataset in range(len(load_data['T1'][q_key].get('Dates', [])[0])):
                        if 'nan' in str(load_data['T1'][q_key].get('Dates', [])[0][dataset]):
                            continue
                        # T1 = load_data['T1'][q_key].get('T1', [])[0][dataset]
                        # errors = load_data['T1'][q_key].get('Errors', [])[0][dataset]
                        date = datetime.datetime.fromtimestamp(load_data['T1'][q_key].get('Dates', [])[0][dataset])
                        I = self.process_h5_data(load_data['T1'][q_key].get('I', [])[0][dataset].decode())
                        Q = self.process_h5_data(load_data['T1'][q_key].get('Q', [])[0][dataset].decode())
                        delay_times = self.process_h5_data(load_data['T1'][q_key].get('Delay Times', [])[0][dataset].decode())
                        # fit = load_data['T1'][q_key].get('Fit', [])[0][dataset]
                        round_num = load_data['T1'][q_key].get('Round Num', [])[0][dataset]
                        batch_num = load_data['T1'][q_key].get('Batch Num', [])[0][dataset]

                        if len(I) > 0:
                            T1_class_instance = T1Measurement(q_key, outerFolder_save_plots, round_num, self.signal, self.save_figs,
                                                              fit_data=True)
                            T1_spec_cfg = ast.literal_eval(self.exp_config['T1_ge'].decode())
                            q1_fit_exponential, T1_err, T1_est, plot_sig = T1_class_instance.t1_fit(I, Q, delay_times)
                            if T1_est < 0:
                                print("The value is negative, continuing...")
                                continue
                            if T1_est > 1000:
                                print("The value is above 1000 us, this is a bad fit, continuing...")
                                continue
                            t1_vals[q_key].extend([T1_est])
                            date_times[q_key].extend([date.strftime("%Y-%m-%d %H:%M:%S")])

                            del T1_class_instance

                del H5_class_instance
        return date_times, t1_vals

    def plot(self, date_times, t1_vals, mcp2_dates, mcp2_temps, show_legends):
        #---------------------------------plot-----------------------------------------------------
        analysis_folder = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/benchmark_analysis_plots/"
        self.create_folder_if_not_exists(analysis_folder)
        analysis_folder = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/benchmark_analysis_plots/features_vs_time/"
        self.create_folder_if_not_exists(analysis_folder)

        font = 14
        titles = [f"Qubit {i+1}" for i in range(self.number_of_qubits)]
        colors = ['orange','blue','purple','green','brown','pink']
        fig, axes = plt.subplots(2, 3, figsize=(12, 8))
        plt.title('T1 Values vs Time',fontsize = font)
        axes = axes.flatten()
        titles = [f"Qubit {i + 1}" for i in range(self.number_of_qubits)]
        from datetime import datetime
        for i, ax in enumerate(axes):

            ax.set_title(titles[i], fontsize = font)

            x = date_times[i]
            y = t1_vals[i]

            datetime_objects = [datetime.strptime(date_string, "%Y-%m-%d %H:%M:%S") for date_string in x]

            combined = list(zip(datetime_objects, y))
            combined.sort(reverse=True, key=lambda item: item[0])
            sorted_x, sorted_y = zip(*combined)

            sorted_x = np.array(sorted_x)
            sorted_y = np.array(sorted_y)

            combined_mcp2 = list(zip(mcp2_dates, mcp2_temps))
            combined_mcp2.sort(reverse=True, key=lambda item: item[0])
            sorted_x_mcp2, sorted_y_mcp2 = zip(*combined_mcp2)

            sorted_x_mcp2 = np.array(sorted_x_mcp2)
            sorted_y_mcp2 = np.array(sorted_y_mcp2)

            ax.scatter(sorted_x, sorted_y, color=colors[i])

            num_points = 5
            indices = np.linspace(0, len(sorted_x) - 1, num_points, dtype=int)
            ax.set_xticks(sorted_x[indices])
            ax.set_xticklabels([dt for dt in sorted_x[indices]], rotation=45)

            min_datetime = sorted_x.min()
            max_datetime = sorted_x.max()

            valid_indices = (sorted_x_mcp2 >= min_datetime) & (sorted_x_mcp2 <= max_datetime)
            filtered_x_mcp2 = sorted_x_mcp2[valid_indices]
            filtered_y_mcp2 = sorted_y_mcp2[valid_indices]

            ax.set_xlabel('Time (Days)', fontsize=font-2)
            ax.set_ylabel('T1 (us)', fontsize=font-2)
            ax.tick_params(axis='both', which='major', labelsize=8)

            ax2 = ax.twinx()
            ax2.scatter(filtered_x_mcp2, filtered_y_mcp2, color='darkblue')
            ax2.set_ylabel('Temperature (mK)', fontsize=font - 2, color='darkblue')
            ax2.tick_params(axis='y', labelcolor='darkblue', labelsize=10)

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(analysis_folder + 'T1_vals_vs_time_thermometry.pdf', transparent=True, dpi=self.final_figure_quality)

        #plt.show()

class T2rVsTemp:
    def __init__(self, figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs, fit_saved,
                 signal, run_name, exp_config):
        self.save_figs = save_figs
        self.fit_saved = fit_saved
        self.signal = signal
        self.figure_quality = figure_quality
        self.run_name = run_name
        self.number_of_qubits = number_of_qubits
        self.final_figure_quality = final_figure_quality
        self.top_folder_dates = top_folder_dates
        self.exp_config = exp_config

    def datetime_to_unix(self, dt):
        # Convert to Unix timestamp
        unix_timestamp = int(dt.timestamp())
        return unix_timestamp

    def unix_to_datetime(self, unix_timestamp):
        # Convert the Unix timestamp to a datetime object
        dt = datetime.fromtimestamp(unix_timestamp)
        return dt

    def create_folder_if_not_exists(self, folder):
        """Creates a folder at the given path if it doesn't already exist."""
        if not os.path.exists(folder):
            os.makedirs(folder)

    def exponential(self, x, a, b, c, d):
        return a * np.exp(-(x - b) / c) + d

    def optimal_bins(self, data):
        n = len(data)
        if n == 0:
            return {}
        # Sturges' Rule
        sturges_bins = int(np.ceil(np.log2(n) + 1))
        return sturges_bins

    def process_string_of_nested_lists(self, data):
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

    def process_h5_data(self, data):
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

    def string_to_float_list(self, input_string):
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

    def run(self):
        import datetime
        # ----------Load/get data------------------------
        t2_vals = {i: [] for i in range(self.number_of_qubits)}
        rounds = []
        reps = []
        file_names = []
        date_times = {i: [] for i in range(self.number_of_qubits)}
        mean_values = {}

        for folder_date in self.top_folder_dates:
            outerFolder = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/" + folder_date + "/"
            outerFolder_save_plots = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/" + folder_date + "_plots/"

            # -------------------------------------------------------Load/Plot/Save T2------------------------------------------
            outerFolder_expt = outerFolder + "/Data_h5/T2_ge/"
            h5_files = glob.glob(os.path.join(outerFolder_expt, "*.h5"))

            for h5_file in h5_files:
                save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]
                H5_class_instance = Data_H5(h5_file)
                load_data = H5_class_instance.load_from_h5(data_type='T2', save_r=int(save_round))

                for q_key in load_data['T2']:
                    for dataset in range(len(load_data['T2'][q_key].get('Dates', [])[0])):
                        if 'nan' in str(load_data['T2'][q_key].get('Dates', [])[0][dataset]):
                            continue
                        # T2 = load_data['T2'][q_key].get('T2', [])[0][dataset]
                        # errors = load_data['T2'][q_key].get('Errors', [])[0][dataset]
                        date = datetime.datetime.fromtimestamp(load_data['T2'][q_key].get('Dates', [])[0][dataset])
                        I = self.process_h5_data(load_data['T2'][q_key].get('I', [])[0][dataset].decode())
                        Q = self.process_h5_data(load_data['T2'][q_key].get('Q', [])[0][dataset].decode())
                        delay_times = self.process_h5_data(load_data['T2'][q_key].get('Delay Times', [])[0][dataset].decode())
                        # fit = load_data['T2'][q_key].get('Fit', [])[0][dataset]
                        round_num = load_data['T2'][q_key].get('Round Num', [])[0][dataset]
                        batch_num = load_data['T2'][q_key].get('Batch Num', [])[0][dataset]

                        if len(I) > 0:
                            T2_class_instance = T2RMeasurement(q_key, outerFolder_save_plots, round_num, self.signal,
                                                               self.save_figs, fit_data=True)
                            try:
                                fitted, t2r_est, t2r_err, plot_sig = T2_class_instance.t2_fit(delay_times, I, Q)
                            except:
                                continue
                            T2_cfg = ast.literal_eval(self.exp_config['Ramsey_ge'].decode())
                            if t2r_est < 0:
                                print("The value is negative, continuing...")
                                continue
                            if t2r_est > 1000:
                                print("The value is above 1000 us, this is a bad fit, continuing...")
                                continue
                            t2_vals[q_key].extend([t2r_est])
                            date_times[q_key].extend([date.strftime("%Y-%m-%d %H:%M:%S")])

                            del T2_class_instance

                del H5_class_instance
        return date_times, t2_vals

    def plot(self, date_times, t2_vals, mcp2_dates, mcp2_temps, show_legends):
        #---------------------------------plot-----------------------------------------------------
        analysis_folder = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/benchmark_analysis_plots/"
        self.create_folder_if_not_exists(analysis_folder)
        analysis_folder = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/benchmark_analysis_plots/features_vs_time/"
        self.create_folder_if_not_exists(analysis_folder)

        font = 14
        titles = [f"Qubit {i+1}" for i in range(self.number_of_qubits)]
        colors = ['orange','blue','purple','green','brown','pink']
        fig, axes = plt.subplots(2, 3, figsize=(12, 8))
        plt.title('T2 Values vs Time',fontsize = font)
        axes = axes.flatten()
        titles = [f"Qubit {i + 1}" for i in range(self.number_of_qubits)]
        from datetime import datetime
        for i, ax in enumerate(axes):

            ax.set_title(titles[i], fontsize = font)

            x = date_times[i]
            y = t2_vals[i]

            datetime_objects = [datetime.strptime(date_string, "%Y-%m-%d %H:%M:%S") for date_string in x]

            combined = list(zip(datetime_objects, y))
            combined.sort(reverse=True, key=lambda item: item[0])
            sorted_x, sorted_y = zip(*combined)

            sorted_x = np.array(sorted_x)
            sorted_y = np.array(sorted_y)

            combined_mcp2 = list(zip(mcp2_dates, mcp2_temps))
            combined_mcp2.sort(reverse=True, key=lambda item: item[0])
            sorted_x_mcp2, sorted_y_mcp2 = zip(*combined_mcp2)

            sorted_x_mcp2 = np.array(sorted_x_mcp2)
            sorted_y_mcp2 = np.array(sorted_y_mcp2)

            ax.scatter(sorted_x, sorted_y, color=colors[i])

            num_points = 5
            indices = np.linspace(0, len(sorted_x) - 1, num_points, dtype=int)
            ax.set_xticks(sorted_x[indices])
            ax.set_xticklabels([dt for dt in sorted_x[indices]], rotation=45)

            min_datetime = sorted_x.min()
            max_datetime = sorted_x.max()

            valid_indices = (sorted_x_mcp2 >= min_datetime) & (sorted_x_mcp2 <= max_datetime)
            filtered_x_mcp2 = sorted_x_mcp2[valid_indices]
            filtered_y_mcp2 = sorted_y_mcp2[valid_indices]

            ax.set_xlabel('Time (Days)', fontsize=font-2)
            ax.set_ylabel('T2 (us)', fontsize=font-2)
            ax.tick_params(axis='both', which='major', labelsize=8)

            ax2 = ax.twinx()
            ax2.scatter(filtered_x_mcp2, filtered_y_mcp2, color='darkblue')
            ax2.set_ylabel('Temperature (mK)', fontsize=font - 2, color='darkblue')
            ax2.tick_params(axis='y', labelcolor='darkblue', labelsize=10)

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(analysis_folder + 'T2_vals_vs_time_thermometry.pdf', transparent=True, dpi=self.final_figure_quality)

        #plt.show()

class T2eVsTemp:
    def __init__(self, figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs, fit_saved,
                 signal, run_name, exp_config):
        self.save_figs = save_figs
        self.fit_saved = fit_saved
        self.signal = signal
        self.figure_quality = figure_quality
        self.run_name = run_name
        self.number_of_qubits = number_of_qubits
        self.final_figure_quality = final_figure_quality
        self.top_folder_dates = top_folder_dates
        self.exp_config = exp_config

    def datetime_to_unix(self, dt):
        # Convert to Unix timestamp
        unix_timestamp = int(dt.timestamp())
        return unix_timestamp

    def unix_to_datetime(self, unix_timestamp):
        # Convert the Unix timestamp to a datetime object
        dt = datetime.fromtimestamp(unix_timestamp)
        return dt

    def create_folder_if_not_exists(self, folder):
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

    def process_string_of_nested_lists(self, data):
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

    def process_h5_data(self, data):
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

    def string_to_float_list(self, input_string):
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

    def run(self):
        import datetime

        # ----------Load/get data------------------------
        t2e_vals = {i: [] for i in range(self.number_of_qubits)}
        rounds = []
        reps = []
        file_names = []
        date_times = {i: [] for i in range(self.number_of_qubits)}
        mean_values = {}

        for folder_date in self.top_folder_dates:
            outerFolder = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/" + folder_date + "/"
            outerFolder_save_plots = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/" + folder_date + "_plots/"

            # -------------------------------------------------------Load/Plot/Save T2E------------------------------------------
            outerFolder_expt = outerFolder + "/Data_h5/T2E_ge/"
            h5_files = glob.glob(os.path.join(outerFolder_expt, "*.h5"))

            for h5_file in h5_files:
                save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]
                H5_class_instance = Data_H5(h5_file)
                load_data = H5_class_instance.load_from_h5(data_type='T2E', save_r=int(save_round))

                for q_key in load_data['T2E']:
                    for dataset in range(len(load_data['T2E'][q_key].get('Dates', [])[0])):
                        if 'nan' in str(load_data['T2E'][q_key].get('Dates', [])[0][dataset]):
                            continue
                        # T2 = load_data['T2E'][q_key].get('T2', [])[0][dataset]
                        # errors = load_data['T2E'][q_key].get('Errors', [])[0][dataset]
                        date = datetime.datetime.fromtimestamp(load_data['T2E'][q_key].get('Dates', [])[0][dataset])
                        I = self.process_h5_data(load_data['T2E'][q_key].get('I', [])[0][dataset].decode())
                        Q = self.process_h5_data(load_data['T2E'][q_key].get('Q', [])[0][dataset].decode())
                        delay_times = self.process_h5_data(load_data['T2E'][q_key].get('Delay Times', [])[0][dataset].decode())
                        # fit = load_data['T2E'][q_key].get('Fit', [])[0][dataset]
                        round_num = load_data['T2E'][q_key].get('Round Num', [])[0][dataset]
                        batch_num = load_data['T2E'][q_key].get('Batch Num', [])[0][dataset]

                        if len(I) > 0:
                            T2E_class_instance = T2EMeasurement(q_key, outerFolder_save_plots, round_num, self.signal, self.save_figs,
                                                               fit_data=True)
                            try:
                                fitted, t2e_est, t2e_err, plot_sig = T2E_class_instance.t2_fit(delay_times, I, Q)
                            except:
                                continue
                            T2E_cfg = ast.literal_eval(self.exp_config['SpinEcho_ge'].decode())
                            if t2e_est < 0:
                                print("The value is negative, continuing...")
                                continue
                            if t2e_est > 1000:
                                print("The value is above 1000 us, this is a bad fit, continuing...")
                                continue
                            t2e_vals[q_key].extend([t2e_est])
                            date_times[q_key].extend([date.strftime("%Y-%m-%d %H:%M:%S")])

                            del T2E_class_instance
                del H5_class_instance
        return date_times, t2e_vals

    def plot(self, date_times, t2e_vals, mcp2_dates, mcp2_temps, show_legends):
        #---------------------------------plot-----------------------------------------------------
        analysis_folder = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/benchmark_analysis_plots/"
        self.create_folder_if_not_exists(analysis_folder)
        analysis_folder = f"/home/nexusadmin/qick/NEXUS_sandbox/Data/{self.run_name}/benchmark_analysis_plots/features_vs_time/"
        self.create_folder_if_not_exists(analysis_folder)

        font = 14
        titles = [f"Qubit {i+1}" for i in range(self.number_of_qubits)]
        colors = ['orange','blue','purple','green','brown','pink']
        fig, axes = plt.subplots(2, 3, figsize=(12, 8))
        plt.title('T2E Values vs Time',fontsize = font)
        axes = axes.flatten()
        titles = [f"Qubit {i + 1}" for i in range(self.number_of_qubits)]
        from datetime import datetime
        for i, ax in enumerate(axes):

            ax.set_title(titles[i], fontsize = font)

            x = date_times[i]
            y = t2e_vals[i]

            # Convert strings to datetime objects.
            datetime_objects = [datetime.strptime(date_string, "%Y-%m-%d %H:%M:%S") for date_string in x]

            datetime_objects = [datetime.strptime(date_string, "%Y-%m-%d %H:%M:%S") for date_string in x]

            combined = list(zip(datetime_objects, y))
            combined.sort(reverse=True, key=lambda item: item[0])
            sorted_x, sorted_y = zip(*combined)

            sorted_x = np.array(sorted_x)
            sorted_y = np.array(sorted_y)

            combined_mcp2 = list(zip(mcp2_dates, mcp2_temps))
            combined_mcp2.sort(reverse=True, key=lambda item: item[0])
            sorted_x_mcp2, sorted_y_mcp2 = zip(*combined_mcp2)

            sorted_x_mcp2 = np.array(sorted_x_mcp2)
            sorted_y_mcp2 = np.array(sorted_y_mcp2)

            ax.scatter(sorted_x, sorted_y, color=colors[i])

            num_points = 5
            indices = np.linspace(0, len(sorted_x) - 1, num_points, dtype=int)
            ax.set_xticks(sorted_x[indices])
            ax.set_xticklabels([dt for dt in sorted_x[indices]], rotation=45)

            min_datetime = sorted_x.min()
            max_datetime = sorted_x.max()

            valid_indices = (sorted_x_mcp2 >= min_datetime) & (sorted_x_mcp2 <= max_datetime)
            filtered_x_mcp2 = sorted_x_mcp2[valid_indices]
            filtered_y_mcp2 = sorted_y_mcp2[valid_indices]

            ax.set_xlabel('Time (Days)', fontsize=font-2)
            ax.set_ylabel('T2E (us)', fontsize=font-2)
            ax.tick_params(axis='both', which='major', labelsize=8)

            ax2 = ax.twinx()
            ax2.scatter(filtered_x_mcp2, filtered_y_mcp2, color='darkblue')
            ax2.set_ylabel('Temperature (mK)', fontsize=font - 2, color='darkblue')
            ax2.tick_params(axis='y', labelcolor='darkblue', labelsize=10)

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(analysis_folder + 'T2E_vals_vs_time_thermometry.pdf', transparent=True, dpi=self.final_figure_quality)

        #plt.show()
