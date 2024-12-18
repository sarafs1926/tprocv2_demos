import numpy as np
import os
import sys
from section_002_res_spec_ge_mux import ResonanceSpectroscopy
from section_004_qubit_spec_ge import QubitSpectroscopy
from section_006_amp_rabi_ge import AmplitudeRabiExperiment
from section_007_T1_ge import T1Measurement
from section_008_save_data_to_h5 import Data_H5
from section_009_T2R_ge import T2RMeasurement
from section_010_T2E_ge import T2EMeasurement
import glob
import re
import datetime
import ast
import os
import matplotlib.pyplot as plt
from scipy.stats import norm
import json
import h5py
from scipy.optimize import curve_fit

class T2eHistCumulErrPlots:
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
        # ----------Load/get data from T2E------------------------
        t2e_vals = {i: [] for i in range(self.number_of_qubits)}
        t2e_errs = {i: [] for i in range(self.number_of_qubits)}
        qubit_for_this_index = []
        rounds = []
        reps = []
        file_names = []
        dates = {i: [] for i in range(self.number_of_qubits)}

        for folder_date in self.top_folder_dates:
            outerFolder = f"/data/QICK_data/{self.run_name}/" + folder_date + "/"
            outerFolder_save_plots = f"/data/QICK_data/{self.run_name}/" + folder_date + "_plots/"

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
                        # T2E = load_data['T2E'][q_key].get('T2E', [])[0][dataset]
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
                                fitted, T2E, T2E_err, plot_sig = T2E_class_instance.t2_fit(delay_times, I, Q)
                            except Exception as e:
                                print("Error fitting:", e)
                                continue

                            T2E_cfg = ast.literal_eval(self.exp_config['SpinEcho_ge'].decode())
                            if T2E < 0:
                                print("The value is negative, continuing...")
                                continue
                            if T2E > 1000:
                                print("The value is above 1000 us, this is a bad fit, continuing...")
                                continue
                            t2e_vals[q_key].extend([T2E])  # Store T1 values
                            t2e_errs[q_key].extend([T2E_err])  # Store T1 error values
                            dates[q_key].extend([date.strftime("%Y-%m-%d %H:%M:%S")])  # Decode bytes to string

                            del T2E_class_instance

                del H5_class_instance
        return dates, t2e_vals, t2e_errs

    def plot(self, dates, t2e_vals, t2e_errs, show_legends):
        #---------------------------------plot-----------------------------------------------------
        analysis_folder = f"/data/QICK_data/{self.run_name}/benchmark_analysis_plots/"
        self.create_folder_if_not_exists(analysis_folder)
        analysis_folder = f"/data/QICK_data/{self.run_name}/benchmark_analysis_plots/T2E/"
        self.create_folder_if_not_exists(analysis_folder)

        fig, axes = plt.subplots(2, 3, figsize=(12, 8))
        axes = axes.flatten()
        font = 14
        titles = [f"Qubit {i+1}" for i in range(self.number_of_qubits)]
        gaussian_xvals =  {i: [] for i in range(0, self.number_of_qubits)}
        gaussian_yvals =  {i: [] for i in range(0, self.number_of_qubits)}
        gaussian_colors = {i: [] for i in range(0, self.number_of_qubits)}
        gaussian_dates = {i: [] for i in range(0, self.number_of_qubits)}
        mean_values = {}
        std_values = {}
        colors = ['orange','blue','purple','green','brown','pink']
        for i, ax in enumerate(axes):
            if len(dates[i])>1:
                date_label = dates[i][0]
            else:
                date_label = ''


            if len(t2e_vals[i]) >1:
                optimal_bin_num = self.optimal_bins(t2e_vals[i])

                # Fit a Gaussian to the raw data instead of the histogram
                # get the mean and standard deviation of the data
                mu_1, std_1 = norm.fit(t2e_vals[i])
                mean_values[f"Qubit {i + 1}"] = mu_1  # Store the mean value for each qubit
                std_values[f"Qubit {i + 1}"] = std_1

                # Generate x values for plotting a gaussian based on this mean and standard deviation
                x_1 = np.linspace(min(t2e_vals[i]), max(t2e_vals[i]), optimal_bin_num)
                p_1 = norm.pdf(x_1, mu_1, std_1)

                # Calculate histogram data for t1_vals[i]
                hist_data_1, bins_1 = np.histogram(t2e_vals[i], bins=optimal_bin_num)
                bin_centers_1 = (bins_1[:-1] + bins_1[1:]) / 2

                # Scale the Gaussian curve to match the histogram
                # the gaussian height natrually doesnt match the bin heights in the histograms
                # np.diff(bins_1)  calculates the width of each bin by taking the difference between bin edges
                # the total counts are in hist_data_1.sum()
                # to scale, multiply data gaussian by bin width to convert the probability density to probability within each bin
                # then multiply by the total count to scale the probability to match the overall number of datapoints
                # https://mathematica.stackexchange.com/questions/262314/fit-function-to-histogram
                # https://stackoverflow.com/questions/23447262/fitting-a-gaussian-to-a-histogram-with-matplotlib-and-numpy-wrong-y-scaling
                ax.plot(x_1, p_1 * (np.diff(bins_1) * hist_data_1.sum()), 'b--', linewidth=2, color=colors[i])

                # Plot histogram and Gaussian fit for t1_vals[i]
                ax.hist(t2e_vals[i], bins=optimal_bin_num, alpha=0.7, color=colors[i], edgecolor='black', label=date_label)

                #make a fuller gaussian to make smoother lotting for cumulative plot
                x_1_full = np.linspace(min(t2e_vals[i]), max(t2e_vals[i]), 2000)
                p_1_full = norm.pdf(x_1_full, mu_1, std_1)

                gaussian_xvals[i].append(x_1_full)
                gaussian_yvals[i].append(p_1_full )
                gaussian_colors[i].append(colors[i])
                gaussian_dates[i].append(date_label)

                #rough start at errors:
                #counts, bin_edges, _ = ax.hist(t1_vals[i], bins=20, alpha=0.7, color='blue', edgecolor='black')
                #bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
                #bin_errors = [np.sqrt(np.sum(t1_errs[i])) for _ in range(len(bin_centers))]  #using error propogation through the sum of valuse in each bin
                #ax.errorbar(bin_centers, counts, yerr=bin_errors, fmt='o', color='red', ecolor='black', capsize=3, linestyle='None')
                if show_legends:
                    ax.legend()
                ax.set_title(titles[i] + f" $\mu$: {mu_1:.2f} $\sigma$:{std_1:.2f}",fontsize = font)
                ax.set_xlabel('T2E (µs)',fontsize = font)
                ax.set_ylabel('Frequency',fontsize = font)
                ax.tick_params(axis='both', which='major', labelsize=font)

        plt.tight_layout()
        plt.savefig( analysis_folder + 'hists.pdf', transparent=True, dpi=self.final_figure_quality)

        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        plt.title('Cumulative Distribution',fontsize = font)
        for i in range(0, len(t2e_vals)):
            if len(dates[i])>1:
                date_label = dates[i][0]
            else:
                date_label = ''

            if len(t2e_vals[i]) > 1:
                t1_vals_sorted = np.sort(t2e_vals[i])
                len_samples = len(t1_vals_sorted)
                var = np.linspace(1,len_samples,len_samples)/ len_samples

                cumulative_gaussian = np.cumsum(gaussian_yvals[i][0]) / np.sum(gaussian_yvals[i][0])
                ax.scatter(t1_vals_sorted,var,color = colors[i], label = f'Q{i+1}', s = 5)
                ax.plot(gaussian_xvals[i][0], cumulative_gaussian, color=colors[i], label='Gauss Fit ' + f'Q{i + 1}',
                        linestyle='--')
                ax.tick_params(axis='both', which='major', labelsize=font)
        #ax.set_title('')
        ax.set_xlabel('T2E (us)',fontsize = font)
        ax.set_ylabel('Cumulative Distribution',fontsize = font)
        ax.loglog()
        ax.legend(edgecolor='black')
        #ax.set_xlim(10**0, 10**3)
        #ax.set_ylim(10 ** -7, 10 ** 0) #to compare to johns plot, need to adjust a little
        plt.tight_layout()
        plt.savefig(analysis_folder + 'cumulative.pdf', transparent=True, dpi=self.final_figure_quality)



        fig, axes = plt.subplots(2, 3, figsize=(12, 8))
        plt.title('Fit Error vs T2E Time',fontsize = font)
        axes = axes.flatten()
        titles = [f"Qubit {i + 1}" for i in range(self.number_of_qubits)]
        for i, ax in enumerate(axes):

            if len(dates[i])>1:
                date_label = dates[i][0]
            else:
                date_label = ''
            ax.set_title(titles[i], fontsize = font)
            ax.scatter(t2e_vals[i], t2e_errs[i], label = date_label, color = colors[i])
            if show_legends:
                ax.legend(edgecolor='black')
            ax.set_xlabel('T2E (us)', fontsize = font)
            ax.set_ylabel('Fit error (us)', fontsize = font)
            ax.tick_params(axis='both', which='major', labelsize=font)
        plt.tight_layout()
        plt.savefig(analysis_folder + 'errs.pdf', transparent=True, dpi=self.final_figure_quality)
        #plt.show()
        return std_values, mean_values