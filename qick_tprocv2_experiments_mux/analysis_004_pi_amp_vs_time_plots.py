import numpy as np
import os
import sys
sys.path.append(os.path.abspath("/home/quietuser/Documents/GitHub/tprocv2_demos/qick_tprocv2_experiments_mux/"))

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
import os
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit

class PiAmpsVsTime:
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

    def run(self, plot_depths = False):
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
            outerFolder = f"/data/QICK_data/{self.run_name}/" + folder_date + "/"
            outerFolder_save_plots = f"/data/QICK_data/{self.run_name}/" + folder_date + "_plots/"

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

    def plot(self, date_times, pi_amps, show_legends):
        #---------------------------------plot-----------------------------------------------------
        analysis_folder = f"/data/QICK_data/{self.run_name}/benchmark_analysis_plots/"
        self.create_folder_if_not_exists(analysis_folder)
        analysis_folder = f"/data/QICK_data/{self.run_name}/benchmark_analysis_plots/features_vs_time/"
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
            ax.set_ylabel('Pi Amp (a.u.)', fontsize=font-2)
            ax.tick_params(axis='both', which='major', labelsize=8)

        plt.tight_layout()
        plt.savefig(analysis_folder + 'Pi_Amps.pdf', transparent=True, dpi=self.final_figure_quality)

        #plt.show()

    def plot_vs_signal_depth(self, date_times, pi_amps, depths, show_legends):
        #---------------------------------plot-----------------------------------------------------
        analysis_folder = f"/data/QICK_data/{self.run_name}/benchmark_analysis_plots/"
        self.create_folder_if_not_exists(analysis_folder)
        analysis_folder = f"/data/QICK_data/{self.run_name}/benchmark_analysis_plots/other/"
        self.create_folder_if_not_exists(analysis_folder)

        font = 14
        titles = [f"Qubit {i+1}" for i in range(self.number_of_qubits)]
        colors = ['orange','blue','purple','green','brown','pink']
        fig, axes = plt.subplots(2, 3, figsize=(12, 8))
        plt.title('Pi Amplitudes vs Signal Depths',fontsize = font)
        axes = axes.flatten()
        titles = [f"Qubit {i + 1}" for i in range(self.number_of_qubits)]
        from datetime import datetime

        for i, ax in enumerate(axes):

            ax.set_title(titles[i], fontsize = font)

            x = depths[i]
            y = pi_amps[i]

            ax.scatter(x, y, color=colors[i])
            if show_legends:
                ax.legend(edgecolor='black')
            ax.set_xlabel('Signal Depth', fontsize=font-2)
            ax.set_ylabel('Pi Amp (a.u.)', fontsize=font-2)
            ax.tick_params(axis='both', which='major', labelsize=8)

        plt.tight_layout()
        plt.savefig(analysis_folder + 'Pi_Amps_vs_depth.pdf', transparent=True, dpi=self.final_figure_quality)

        #plt.show()

    def plot_signal_depth_vs_time(self, date_times, pi_amps, depths, show_legends):
        #---------------------------------plot-----------------------------------------------------
        analysis_folder = f"/data/QICK_data/{self.run_name}/benchmark_analysis_plots/"
        self.create_folder_if_not_exists(analysis_folder)
        analysis_folder = f"/data/QICK_data/{self.run_name}/benchmark_analysis_plots/other/"
        self.create_folder_if_not_exists(analysis_folder)

        font = 14
        titles = [f"Qubit {i+1}" for i in range(self.number_of_qubits)]
        colors = ['orange','blue','purple','green','brown','pink']
        fig, axes = plt.subplots(2, 3, figsize=(12, 8))
        plt.title('Signal Depths vs Time',fontsize = font)
        axes = axes.flatten()
        titles = [f"Qubit {i + 1}" for i in range(self.number_of_qubits)]
        from datetime import datetime

        for i, ax in enumerate(axes):

            ax.set_title(titles[i], fontsize = font)

            x = date_times[i]
            y = depths[i]

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
            ax.set_ylabel('Signal Depth (a.u.)', fontsize=font-2)
            ax.tick_params(axis='both', which='major', labelsize=8)

        plt.tight_layout()
        plt.savefig(analysis_folder + 'depths_vs_time.pdf', transparent=True, dpi=self.final_figure_quality)

        #plt.show()

    def plot_vs_temps(self, date_times, pi_amps, temps, show_legends):
        #---------------------------------plot-----------------------------------------------------
        analysis_folder = f"/data/QICK_data/{self.run_name}/benchmark_analysis_plots/"
        self.create_folder_if_not_exists(analysis_folder)
        analysis_folder = f"/data/QICK_data/{self.run_name}/benchmark_analysis_plots/other/"
        self.create_folder_if_not_exists(analysis_folder)

        font = 14
        titles = [f"Qubit {i+1}" for i in range(self.number_of_qubits)]
        colors = ['orange','blue','purple','green','brown','pink']
        fig, axes = plt.subplots(2, 3, figsize=(12, 8))
        plt.title('Pi Amplitudes vs Qubit temps',fontsize = font)
        axes = axes.flatten()
        titles = [f"Qubit {i + 1}" for i in range(self.number_of_qubits)]
        from datetime import datetime

        for i, ax in enumerate(axes):

            ax.set_title(titles[i], fontsize = font)
            x = temps[i]
            y = pi_amps[i]

            ax.scatter(x, y, color=colors[i])
            if show_legends:
                ax.legend(edgecolor='black')
            ax.set_xlabel('Qubit temp (mK)', fontsize=font-2)
            ax.set_ylabel('Pi Amp (a.u)', fontsize=font-2)
            ax.tick_params(axis='both', which='major', labelsize=8)

        plt.tight_layout()
        plt.savefig(analysis_folder + 'Pi_Amps_vs_qtemp.pdf', transparent=True, dpi=self.final_figure_quality)

        #plt.show()

    def plot_vs_ssf(self, date_times, pi_amps, ssf, show_legends):
        #---------------------------------plot-----------------------------------------------------
        analysis_folder = f"/data/QICK_data/{self.run_name}/benchmark_analysis_plots/"
        self.create_folder_if_not_exists(analysis_folder)
        analysis_folder = f"/data/QICK_data/{self.run_name}/benchmark_analysis_plots/other/"
        self.create_folder_if_not_exists(analysis_folder)

        font = 14
        titles = [f"Qubit {i+1}" for i in range(self.number_of_qubits)]
        colors = ['orange','blue','purple','green','brown','pink']
        fig, axes = plt.subplots(2, 3, figsize=(12, 8))
        plt.title('Pi Amplitudes vs SSF',fontsize = font)
        axes = axes.flatten()
        titles = [f"Qubit {i + 1}" for i in range(self.number_of_qubits)]
        from datetime import datetime

        for i, ax in enumerate(axes):

            ax.set_title(titles[i], fontsize = font)
            x = ssf[i]
            y = pi_amps[i]

            ax.scatter(x, y, color=colors[i])
            if show_legends:
                ax.legend(edgecolor='black')
            ax.set_xlabel('SSF', fontsize=font-2)
            ax.set_ylabel('Pi Amp (a.u)', fontsize=font-2)
            ax.tick_params(axis='both', which='major', labelsize=8)

        plt.tight_layout()
        plt.savefig(analysis_folder + 'Pi_Amps_vs_ssf.pdf', transparent=True, dpi=self.final_figure_quality)

        #plt.show()