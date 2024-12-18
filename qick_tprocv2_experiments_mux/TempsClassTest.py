from section_004_qubit_spec_ge import QubitSpectroscopy
from section_008_save_data_to_h5 import Data_H5
from section_005_single_shot_ge import SingleShot
import glob
import re
import datetime
import ast
import os
import sys
from matplotlib.dates import DateFormatter
import numpy as np
import h5py
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
import math

sys.path.append(os.path.abspath("/home/quietuser/Documents/GitHub/tprocv2_demos/qick_tprocv2_experiments_mux/"))

class Temperature_Results:
    def __init__(self, date, figure_quality, save_figs, fit_saved, signal, run_name, number_of_qubits, outerFolder,
                 outerFolder_save_plots):
        #self.exp_cfg = None
        self.date = date
        self.figure_quality = figure_quality
        self.save_figs = save_figs
        self.fit_saved = fit_saved
        self.signal = signal
        self.run_name = run_name
        self.number_of_qubits = number_of_qubits
        self.outerFolder = outerFolder
        self.outerFolder_save_plots = outerFolder_save_plots

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
            numbers = [float(x.strip('[').strip(']').replace("'", "").replace(" ", "").replace("  ", "")) for x in
                       match.split()]  # Convert strings to integers
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
        loader_config_instance = Data_H5(self.outerFolder)
        sys_config = loader_config_instance.load_config('sys_config.h5')  # _batch2
        del loader_config_instance

        loader_config_instance = Data_H5(self.outerFolder)
        exp_config = loader_config_instance.load_config('expt_cfg.h5')  # _batch2
        self.exp_cfg = exp_config

        #self.load_plot_save_ss(exp_config)
        del loader_config_instance


    def calculate_qubit_temperature(self, frequency_mhz, ground_state_population, excited_state_population):
        k_B = 1.380649e-23  # Boltzmann constant in J/K
        h = 6.62607015e-34  # Planck's constant in J·s
        frequency_hz = frequency_mhz * 1e6
        T = (h * frequency_hz) / (k_B * np.log(ground_state_population / excited_state_population))
        return T

    def fit_double_gaussian_with_full_coverage(self, iq_data):
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
        ground_gaussian_fit = weights[ground_gaussian] * (
                    1 / (np.sqrt(2 * np.pi) * covariances[ground_gaussian])) * np.exp(
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

        return ground_state_population, excited_state_population_overlap, gmm, means, covariances, weights, crossing_point, ground_gaussian, excited_gaussian, ground_data, excited_data, iq_data

    # ----------------------------------------------Load/Plot/Save QSpec------------------------------------
    def load_plot_save_q_spec(self):
        self.outerFolder_expt = self.outerFolder + "/Data_h5/QSpec_ge/"
        h5_files = glob.glob(os.path.join(self.outerFolder_expt, "*.h5"))
        qubit_frequencies = []
        for h5_file in h5_files:
            save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]
            H5_class_instance = Data_H5(h5_file)
            load_data = H5_class_instance.load_from_h5(data_type='QSpec', save_r=int(save_round))

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
                    I = self.process_h5_data(load_data['QSpec'][q_key].get('I', [])[0][dataset].decode())
                    Q = self.process_h5_data(load_data['QSpec'][q_key].get('Q', [])[0][dataset].decode())
                    # I_fit = load_data['QSpec'][q_key].get('I Fit', [])[0][dataset]
                    # Q_fit = load_data['QSpec'][q_key].get('Q Fit', [])[0][dataset]
                    freqs = self.process_h5_data(load_data['QSpec'][q_key].get('Frequencies', [])[0][dataset].decode())
                    round_num = load_data['QSpec'][q_key].get('Round Num', [])[0][dataset]
                    batch_num = load_data['QSpec'][q_key].get('Batch Num', [])[0][dataset]

                    if len(I) > 0:
                        qspec_class_instance = QubitSpectroscopy(q_key, self.outerFolder_save_plots, round_num,
                                                                 self.signal, self.save_figs)

                        # Ensure exp_cfg is loaded
                        if not hasattr(self, 'exp_cfg') or self.exp_cfg is None:
                            loader_config_instance = Data_H5(self.outerFolder)
                            self.exp_cfg = loader_config_instance.load_config('expt_cfg.h5')
                            del loader_config_instance

                        q_spec_cfg = ast.literal_eval(self.exp_cfg['qubit_spec_ge'].decode())
                        # print('q_spec_cfg: ', q_spec_cfg)
                        #qspec_class_instance.plot_results(I, Q, freqs, q_spec_cfg, self.figure_quality)
                        largest_amp_curve_mean, I_fit, Q_fit = qspec_class_instance.get_results(I, Q, freqs)

                        qubit_frequencies.append({
                            'h5_file': h5_file,
                            'q_key': q_key,
                            'date': date,
                            'largest_amp_curve_mean': largest_amp_curve_mean
                        })

                        del qspec_class_instance

            del H5_class_instance

        print('qubit_frequencies: ', qubit_frequencies)
        return qubit_frequencies
    # ------------------------------------------------Load/Plot/Save SS---------------------------------------
    def load_plot_save_ss(self):
        qubit_frequencies = self.load_plot_save_q_spec()
        temperature_folder = os.path.join(self.outerFolder, "Temperatures")
        if not os.path.exists(temperature_folder):
            os.makedirs(temperature_folder)

        outerFolder_expt = self.outerFolder + "/Data_h5/SS_ge/"
        h5_files = glob.glob(os.path.join(self.outerFolder_expt, "*.h5"))
        # Initialize a dictionary to store temperatures for each qubit
        qubit_temperatures = {i: [] for i in range(self.number_of_qubits)}  # Assuming there are 6 qubits


        for h5_file in h5_files:

            save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]
            H5_class_instance = Data_H5(h5_file)
            load_data = H5_class_instance.load_from_h5(data_type='SS', save_r=int(save_round))

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
                qubit_folder = os.path.join(temperature_folder, f"Qubit_{q_key + 1}")
                if not os.path.exists(qubit_folder):
                    os.makedirs(qubit_folder)
                for dataset in range(len(load_data['SS'][q_key].get('Dates', [])[0])):
                    timestamp = load_data['SS'][q_key].get('Dates', [])[0][dataset]
                    date = datetime.datetime.fromtimestamp(load_data['SS'][q_key].get('Dates', [])[0][dataset])
                    angle = load_data['SS'][q_key].get('Angle', [])[0][dataset]
                    fidelity = load_data['SS'][q_key].get('Fidelity', [])[0][dataset]
                    I_g = self.process_h5_data(load_data['SS'][q_key].get('I_g', [])[0][dataset].decode())
                    Q_g = self.process_h5_data(load_data['SS'][q_key].get('Q_g', [])[0][dataset].decode())
                    I_e = self.process_h5_data(load_data['SS'][q_key].get('I_e', [])[0][dataset].decode())
                    Q_e = self.process_h5_data(load_data['SS'][q_key].get('Q_e', [])[0][dataset].decode())
                    round_num = load_data['SS'][q_key].get('Round Num', [])[0][dataset]
                    batch_num = load_data['SS'][q_key].get('Batch Num', [])[0][dataset]

                    I_g = np.array(I_g)
                    Q_g = np.array(Q_g)
                    I_e = np.array(I_e)
                    Q_e = np.array(Q_e)

                    if len(Q_g) > 0:
                        ss_class_instance = SingleShot(q_key, self.outerFolder_save_plots, round_num, self.save_figs)
                        ss_cfg = ast.literal_eval(self.exp_cfg['Readout_Optimization'].decode())
                        #ss_class_instance.hist_ssf(data=[I_g, Q_g, I_e, Q_e], cfg=ss_cfg, plot=True)
                        fid, threshold, rotation_angle, ig_new, ie_new = ss_class_instance.hist_ssf(
                            data=[I_g, Q_g, I_e, Q_e], cfg=ss_cfg, plot=False)

                        # Extract the relevant portion of the file name for matching
                        base_h5_file = "_".join(
                            h5_file.split('/')[-1].split('_')[:2])  # Extract up to 2024-12-11_11-45-27
                        # print(f"Base H5 File for Matching: {base_h5_file}")

                        qubit_frequency = [
                            entry['largest_amp_curve_mean']
                            for entry in qubit_frequencies
                            if
                            "_".join(entry['h5_file'].split('/')[-1].split('_')[:2]) == base_h5_file and entry[
                                'q_key'] == q_key
                        ]
                        if len(qubit_frequency) == 0:
                            print(f"No match found for h5_file: {base_h5_file}, q_key: {q_key}. Skipping.")
                            continue

                        qubit_frequency = qubit_frequency[0]  # there should only be one value inside this list
                        ground_state_population, excited_state_population_overlap, gmm, means, covariances, weights, crossing_point, ground_gaussian, excited_gaussian, ground_data, excited_data, iq_data = self.fit_double_gaussian_with_full_coverage(
                            ig_new)
                        temperature_k = self.calculate_qubit_temperature(qubit_frequency, ground_state_population,
                                                                    excited_state_population_overlap)
                        temperature_mk = temperature_k * 1e3

                        print(f"Ground state population: {ground_state_population}")
                        print(f"Excited state (leakage) population: {excited_state_population_overlap}")
                        print(f"Qubit {q_key + 1} Temperature: {temperature_mk:.2f} mK", "\n")

                        qubit_temperatures[q_key].append((temperature_mk, timestamp))  # save temps for each qubit

                        # Plotting double gaussian distributions and fitting
                        xlims = [np.min(ig_new), np.max(ig_new)]
                        plt.figure(figsize=(10, 6))

                        # Plot histogram for `ig_new`
                        steps = 3000
                        numbins = round(math.sqrt(steps))
                        n, bins, _ = plt.hist(ig_new, bins=numbins, range=xlims, density=False, alpha=0.5,
                                              label='Histogram of $I_g$',
                                              color='gray')
                        # print(numbins)
                        # Use the midpoints of bins to create boolean masks
                        bin_centers = (bins[:-1] + bins[1:]) / 2
                        ground_region = (bin_centers < crossing_point)
                        excited_region = (bin_centers >= crossing_point)

                        # Calculate scaling factors for each region
                        scaling_factor_ground = max(n[ground_region]) / max(
                            (weights[ground_gaussian] / (np.sqrt(2 * np.pi) * covariances[ground_gaussian])) * np.exp(
                                -0.5 * ((bin_centers[ground_region] - means[ground_gaussian]) / covariances[
                                    ground_gaussian]) ** 2))

                        scaling_factor_excited = max(n[excited_region]) / max(
                            (weights[excited_gaussian] / (np.sqrt(2 * np.pi) * covariances[excited_gaussian])) * np.exp(
                                -0.5 * ((bin_centers[excited_region] - means[excited_gaussian]) / covariances[
                                    excited_gaussian]) ** 2))

                        # Generate x values for plotting Gaussian components
                        x = np.linspace(xlims[0], xlims[1], 1000)
                        ground_gaussian_fit = scaling_factor_ground * (
                                weights[ground_gaussian] / (
                                    np.sqrt(2 * np.pi) * covariances[ground_gaussian])) * np.exp(
                            -0.5 * ((x - means[ground_gaussian]) / covariances[ground_gaussian]) ** 2)
                        excited_gaussian_fit = scaling_factor_excited * (
                                weights[excited_gaussian] / (
                                    np.sqrt(2 * np.pi) * covariances[excited_gaussian])) * np.exp(
                            -0.5 * ((x - means[excited_gaussian]) / covariances[excited_gaussian]) ** 2)

                        plt.plot(x, ground_gaussian_fit, label='Ground Gaussian Fit', color='blue', linewidth=2)
                        plt.plot(x, excited_gaussian_fit, label='Excited (leakage) Gaussian Fit', color='red',
                                 linewidth=2)
                        plt.axvline(crossing_point, color='black', linestyle='--', linewidth=1,
                                    label=f'Crossing Point ({crossing_point:.2f})')

                        # Add shading for ground and excited state regions
                        x_vals = np.linspace(np.min(ig_new), np.max(ig_new), 1000)

                        # Add shading for ground_data points
                        plt.hist(
                            ground_data, bins=numbins, range=[np.min(ig_new), np.max(ig_new)], density=False,
                            alpha=0.5, color="blue", label="Ground Data Region", zorder=2
                        )

                        # Add shading for excited_data points
                        plt.hist(
                            excited_data, bins=numbins, range=[np.min(ig_new), np.max(ig_new)], density=False,
                            alpha=0.5, color="red", label="Excited Data Region", zorder=3
                        )

                        # plt.hist(
                        #     iq_data, bins=numbins, range=[np.min(ig_new), np.max(ig_new)], density=False,
                        #     alpha=0.2, color="green", label="All IQ Data Region", zorder=1
                        # )

                        plt.title(
                            f"Fidelity Histogram and Double Gaussian Fit ; Qubit {q_key + 1}; Fidelity = {fidelity * 100:.2f}%")
                        plt.xlabel('$I_g$', fontsize=14)
                        # plt.ylabel('Probability Density', fontsize=14)
                        plt.legend()
                        # plt.show()

                        # Save the plot to the Temperatures folder
                        plot_filename = os.path.join(qubit_folder,
                                                     f"Qubit{q_key + 1}_Fidelityhist_gaussianfit_Dataset{dataset}_{datetime.datetime.now().strftime('%Y%m%d%H%M%S')}.png")
                        plt.savefig(plot_filename)
                        # print(f"Plot saved to {plot_filename}")
                        plt.close()

                        # Cleanup
                        del ig_new, ground_gaussian_fit, excited_gaussian_fit, n, bins, x, gmm, \
                            bin_centers, ground_region, excited_region, scaling_factor_ground, \
                            scaling_factor_excited, weights, means, covariances, ground_gaussian, \
                            excited_gaussian, crossing_point

                        del ss_class_instance

            del H5_class_instance
            return qubit_temperatures

    def temperature_plots(self, qubit_temperatures):
        # ------------------------------------------Histograms--------------------------------------------------------

        temperature_folder = os.path.join(self.outerFolder, "Temperatures")
        if not os.path.exists(temperature_folder):
            os.makedirs(temperature_folder)

        temp_histograms_folder = os.path.join(temperature_folder, "Temps_Histograms")
        if not os.path.exists(temp_histograms_folder):
            os.makedirs(temp_histograms_folder)

        # After the loop is complete, create the subplots for temperature histograms
        plt.figure(figsize=(15, 10))
        colors = ['orange', 'blue', 'purple', 'green', 'brown', 'pink']
        for qubit_id, temperatures in qubit_temperatures.items():
            temperatures = [temp[0] for temp in qubit_temperatures[qubit_id]]
            plt.subplot(2, 3, qubit_id + 1)
            plt.hist(temperatures, bins=20, color=colors[qubit_id], alpha=0.7, edgecolor='black')
            plt.title(f"Qubit {qubit_id + 1} Temperature Distribution")
            plt.xlabel("Temperature (mK)")
            plt.ylabel("Frequency")
            plt.grid(alpha=0.3)

        # Adjust layout and save the figure
        plt.tight_layout()
        hist_plot_filename = os.path.join(temp_histograms_folder,
                                          f"Temperature_Distributions_{datetime.datetime.now().strftime('%Y%m%d%H%M%S')}.png")
        plt.savefig(hist_plot_filename)
        print(f"Temperature histograms saved to {hist_plot_filename}")

        # Cleanup
        plt.close()

        # ------------------------------------------Scatter Plots---------------------------------------------------
        # Create a new folder to store temperature scatter plots
        temp_scatter_folder = os.path.join(temperature_folder, "Temps_Scatter")
        if not os.path.exists(temp_scatter_folder):
            os.makedirs(temp_scatter_folder)

        # Create subplots for temperatures over time
        plt.figure(figsize=(15, 10))
        colors = ['orange', 'blue', 'purple', 'green', 'brown', 'pink']
        for qubit_id, data in qubit_temperatures.items():
            if data:  # Check if there's data for the qubit
                temperatures, timestamps = zip(*data)  # Unpack temperatures and timestamps
                timestamps = [datetime.datetime.fromtimestamp(ts) for ts in
                              timestamps]  # Convert timestamps to datetime

                ax = plt.subplot(2, 3, qubit_id + 1)  # Capture the subplot into 'ax'
                ax.scatter(timestamps, temperatures, color=colors[qubit_id], alpha=0.7, edgecolor='black')
                ax.set_title(f"Qubit {qubit_id + 1} Temperature Over Time")
                ax.set_xlabel("Time")
                ax.set_ylabel("Temperature (mK)")
                ax.grid(alpha=0.3)

                # Format x-axis to show full timestamps
                date_formatter = DateFormatter('%m-%d %H:%M')  # Customize date and time format
                ax.xaxis.set_major_formatter(date_formatter)
                plt.xticks(rotation=45)

        # Adjust layout and save the figure
        plt.tight_layout()
        scatter_plot_filename = os.path.join(temp_scatter_folder,
                                             f"Temperature_Over_Time_{datetime.datetime.now().strftime('%Y%m%d%H%M%S')}.png")
        plt.savefig(scatter_plot_filename)
        print(f"Temperature scatter plots saved to {scatter_plot_filename}")

        # Cleanup
        plt.close()


