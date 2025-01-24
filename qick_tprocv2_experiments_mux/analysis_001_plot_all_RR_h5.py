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

class PlotAllRR:
    def __init__(self, list_of_all_qubits, date, figure_quality, save_figs, fit_saved, signal, run_name, number_of_qubits, outerFolder,
                 outerFolder_save_plots, exp_config):
        self.date = date
        self.list_of_all_qubits = list_of_all_qubits
        self.figure_quality = figure_quality
        self.save_figs = save_figs
        self.fit_saved = fit_saved
        self.signal = signal
        self.run_name = run_name
        self.number_of_qubits = number_of_qubits
        self.outerFolder = outerFolder
        self.outerFolder_save_plots = outerFolder_save_plots
        self.exp_config = exp_config

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
    
    def run(self, plot_res_spec = True, plot_q_spec = True, plot_rabi = True, plot_ss = True, plot_t1 = True,
            plot_t2r = True, plot_t2e = True):

        if plot_res_spec:
            self.load_plot_save_res_spec(self.exp_config)
        if plot_q_spec:
            self.load_plot_save_q_spec(self.exp_config)
        if plot_rabi:
            self.load_plot_save_rabi(self.exp_config)
        if plot_ss:
            self.load_plot_save_ss(self.exp_config)
        if plot_t1:
            self.load_plot_save_t1(self.exp_config)
        if plot_t2r:
            self.load_plot_save_t2r(self.exp_config)
        if plot_t2e:
            self.load_plot_save_t2e(self.exp_config)
        

    def load_plot_save_res_spec(self, exp_config):
        # ------------------------------------------Load/Plot/Save Res Spec------------------------------------
        outerFolder_expt = self.outerFolder + "/Data_h5/Res_ge/"
        h5_files = glob.glob(os.path.join(outerFolder_expt, "*.h5"))
        
        for h5_file in h5_files:
            save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]
            H5_class_instance = Data_H5(h5_file)
            #H5_class_instance.print_h5_contents(h5_file)
            load_data = H5_class_instance.load_from_h5(data_type=  'Res', save_r = int(save_round))
        
            #just look at this resonator data, should have batch_num of arrays in each one
            #right now the data writes the same thing batch_num of times, so it will do the same 5 datasets 5 times, until you fix this just grab the first one (All 5)
        
            populated_keys = []
            for q_key in load_data['Res']:
                # Access 'Dates' for the current q_key
                dates_list = load_data['Res'][q_key].get('Dates', [[]])
        
                # Check if any entry in 'Dates' is not NaN
                if any(
                        not np.isnan(date)
                        for date in dates_list[0]  # Iterate over the first batch of dates
                ):
                    populated_keys.append(q_key)
        
            for q_key in populated_keys:
                #go through each dataset in the batch and plot
                for dataset in range(len(load_data['Res'][q_key].get('Dates', [])[0])):
                    date = datetime.datetime.fromtimestamp(load_data['Res'][q_key].get('Dates', [])[0][dataset])   #single date per dataset
                    freq_pts = self.process_h5_data(load_data['Res'][q_key].get('freq_pts', [])[0][dataset].decode())   # comes in as an array but put into a byte string, need to convert to list
                    freq_center = self.process_h5_data(load_data['Res'][q_key].get('freq_center', [])[0][dataset].decode()) # comes in as an array but put into a string, need to convert to list
                    freqs_found = self.string_to_float_list(load_data['Res'][q_key].get('Found Freqs', [])[0][dataset].decode()) #comes in as a list of floats in string format, need to convert
                    amps =  self.process_string_of_nested_lists(load_data['Res'][q_key].get('Amps', [])[0][dataset].decode())  #list of lists
                    round_num = load_data['Res'][q_key].get('Round Num', [])[0][dataset] #already a float
                    batch_num = load_data['Res'][q_key].get('Batch Num', [])[0][dataset]
                    freq_pts_data = load_data['Res'][q_key].get('freq_pts', [])[0][dataset].decode()
        
                    # Replace whitespace between numbers with commas to make it a valid list
                    formatted_str = freq_pts_data.replace('  ', ',').replace('\n', '')
                    formatted_str = formatted_str.replace(' ', ',').replace('\n', '')
                    formatted_str = formatted_str.replace(',]', ']').replace('\n', '')
                    formatted_str = formatted_str.replace('],[', '],[')
                    formatted_str = re.sub(r",,", ",", formatted_str)
                    formatted_str = re.sub(r",\s*([\]])", r"\1", formatted_str)
                    formatted_str = re.sub(r"(\d+)\.,", r"\1.0,",
                                           formatted_str)  # Fix malformed floating-point numbers (e.g., '5829.,' -> '5829.0')
                    # Convert to NumPy array
                    freq_points = np.array(eval(formatted_str))
        
                    if len(freq_pts) > 0:
                        res_class_instance = ResonanceSpectroscopy(q_key, self.outerFolder_save_plots, round_num, self.save_figs)
                        res_spec_cfg = ast.literal_eval(self.exp_config['res_spec'].decode())
                        res_class_instance.plot_results(freq_points, freq_center, amps, res_spec_cfg, self.figure_quality)
                        del res_class_instance
        
            del H5_class_instance

    def load_plot_save_q_spec(self, exp_config):
        # ----------------------------------------------Load/Plot/Save QSpec------------------------------------
        outerFolder_expt = self.outerFolder + "/Data_h5/QSpec_ge/"
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
        
            for q_key in populated_keys:
                for dataset in range(len(load_data['QSpec'][q_key].get('Dates', [])[0])):
                    date = datetime.datetime.fromtimestamp(load_data['QSpec'][q_key].get('Dates', [])[0][dataset])
                    I = self.process_h5_data(load_data['QSpec'][q_key].get('I', [])[0][dataset].decode())
                    Q = self.process_h5_data(load_data['QSpec'][q_key].get('Q', [])[0][dataset].decode())
                    #I_fit = load_data['QSpec'][q_key].get('I Fit', [])[0][dataset]
                    #Q_fit = load_data['QSpec'][q_key].get('Q Fit', [])[0][dataset]
                    freqs = self.process_h5_data(load_data['QSpec'][q_key].get('Frequencies', [])[0][dataset].decode())
                    round_num = load_data['QSpec'][q_key].get('Round Num', [])[0][dataset]
                    batch_num = load_data['QSpec'][q_key].get('Batch Num', [])[0][dataset]
        
                    if len(I)>0:
        
                        qspec_class_instance = QubitSpectroscopy(q_key, self.outerFolder_save_plots, round_num, self.signal, self.save_figs)
                        q_spec_cfg = ast.literal_eval(self.exp_config['qubit_spec_ge'].decode())
                        #print('q_spec_cfg: ', q_spec_cfg)
                        qspec_class_instance.plot_results(I, Q, freqs, q_spec_cfg, self.figure_quality)
                        del qspec_class_instance
        
            del H5_class_instance

    def load_plot_save_rabi(self, exp_config):
        # ------------------------------------------------Load/Plot/Save Rabi---------------------------------------
        outerFolder_expt = self.outerFolder + "/Data_h5/Rabi_ge/"
        h5_files = glob.glob(os.path.join(outerFolder_expt, "*.h5"))
        
        for h5_file in h5_files:
        
            save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]
            H5_class_instance = Data_H5(h5_file)
            load_data = H5_class_instance.load_from_h5(data_type=  'Rabi', save_r = int(save_round))
        
            populated_keys = []
            for q_key in load_data['Rabi']:
                # Access 'Dates' for the current q_key
                dates_list = load_data['Rabi'][q_key].get('Dates', [[]])
        
                # Check if any entry in 'Dates' is not NaN
                if any(
                        not np.isnan(date)
                        for date in dates_list[0]  # Iterate over the first batch of dates
                ):
                    populated_keys.append(q_key)
        
            for q_key in populated_keys:
                for dataset in range(len(load_data['Rabi'][q_key].get('Dates', [])[0])):
                    date= datetime.datetime.fromtimestamp(load_data['Rabi'][q_key].get('Dates', [])[0][dataset])
                    I = self.process_h5_data(load_data['Rabi'][q_key].get('I', [])[0][dataset].decode())
                    Q = self.process_h5_data(load_data['Rabi'][q_key].get('Q', [])[0][dataset].decode())
                    gains = self.process_h5_data(load_data['Rabi'][q_key].get('Gains', [])[0][dataset].decode())
                    #fit = load_data['Rabi'][q_key].get('Fit', [])[0][dataset]
                    round_num = load_data['Rabi'][q_key].get('Round Num', [])[0][dataset]
                    batch_num = load_data['Rabi'][q_key].get('Batch Num', [])[0][dataset]
        
                    if len(I)>0:
        
                        rabi_class_instance = AmplitudeRabiExperiment(q_key, self.outerFolder_save_plots, round_num, self.signal, self.save_figs)
                        rabi_cfg = ast.literal_eval(self.exp_config['power_rabi_ge'].decode())
                        I = np.asarray(I)
                        Q = np.asarray(Q)
                        gains = np.asarray(gains)
                        rabi_class_instance.plot_results(I, Q, gains, rabi_cfg, self.figure_quality)
                        del rabi_class_instance
        
            del H5_class_instance

    def load_plot_save_ss(self, exp_config):
        
        # ------------------------------------------------Load/Plot/Save SS---------------------------------------
        outerFolder_expt = self.outerFolder + "/Data_h5/SS_ge/"
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
        
            for q_key in populated_keys:
                for dataset in range(len(load_data['SS'][q_key].get('Dates', [])[0])):
                    date= datetime.datetime.fromtimestamp(load_data['SS'][q_key].get('Dates', [])[0][dataset])
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
        
                    if len(Q_g)>0:
                        ss_class_instance = SingleShot(q_key, self.outerFolder_save_plots, round_num, self.save_figs)
                        ss_cfg = ast.literal_eval(self.exp_config['Readout_Optimization'].decode())
                        ss_class_instance.hist_ssf(data=[I_g, Q_g, I_e, Q_e], cfg=ss_cfg, plot=True)
                        del ss_class_instance
        
            del H5_class_instance

    def load_plot_save_t1(self, exp_config):
        # ------------------------------------------------Load/Plot/Save T1----------------------------------------------
        outerFolder_expt = self.outerFolder + "/Data_h5/T1_ge/"
        h5_files = glob.glob(os.path.join(outerFolder_expt, "*.h5"))
        
        for h5_file in h5_files:
        
            save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]
            H5_class_instance = Data_H5(h5_file)
            load_data = H5_class_instance.load_from_h5(data_type=  'T1', save_r = int(save_round))
        
            populated_keys = []
            for q_key in load_data['T1']:
                # Access 'Dates' for the current q_key
                dates_list = load_data['T1'][q_key].get('Dates', [[]])
        
                # Check if any entry in 'Dates' is not NaN
                if any(
                        not np.isnan(date)
                        for date in dates_list[0]  # Iterate over the first batch of dates
                ):
                    populated_keys.append(q_key)
        
            for q_key in populated_keys:
                for dataset in range(len(load_data['T1'][q_key].get('Dates', [])[0])):
                    #T1 = load_data['T1'][q_key].get('T1', [])[0][dataset]
                    #errors = load_data['T1'][q_key].get('Errors', [])[0][dataset]
                    date= datetime.datetime.fromtimestamp(load_data['T1'][q_key].get('Dates', [])[0][dataset])
                    I = self.process_h5_data(load_data['T1'][q_key].get('I', [])[0][dataset].decode())
                    Q = self.process_h5_data(load_data['T1'][q_key].get('Q', [])[0][dataset].decode())
                    delay_times = self.process_h5_data(load_data['T1'][q_key].get('Delay Times', [])[0][dataset].decode())
                    #fit = load_data['T1'][q_key].get('Fit', [])[0][dataset]
                    round_num = load_data['T1'][q_key].get('Round Num', [])[0][dataset]
                    batch_num = load_data['T1'][q_key].get('Batch Num', [])[0][dataset]
        
                    if len(I)>0:
                        T1_class_instance = T1Measurement(q_key, self.outerFolder_save_plots, round_num, self.signal, self.save_figs, fit_data = True)
                        T1_spec_cfg = ast.literal_eval(self.exp_config['T1_ge'].decode())
                        T1_class_instance.plot_results(I, Q, delay_times, date, T1_spec_cfg, self.figure_quality)
                        del T1_class_instance
        
            del H5_class_instance

    def load_plot_save_t2r(self, exp_config):
        # -------------------------------------------------------Load/Plot/Save T2R------------------------------------------
        outerFolder_expt = self.outerFolder + "/Data_h5/T2_ge/"
        h5_files = glob.glob(os.path.join(outerFolder_expt, "*.h5"))
        
        for h5_file in h5_files:
            save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]
            H5_class_instance = Data_H5(h5_file)
            load_data = H5_class_instance.load_from_h5(data_type=  'T2', save_r = int(save_round))
        
            populated_keys = []
            for q_key in load_data['T2']:
                # Access 'Dates' for the current q_key
                dates_list = load_data['T2'][q_key].get('Dates', [[]])
        
                # Check if any entry in 'Dates' is not NaN
                if any(
                        not np.isnan(date)
                        for date in dates_list[0]  # Iterate over the first batch of dates
                ):
                    populated_keys.append(q_key)
        
            for q_key in populated_keys:
                for dataset in range(len(load_data['T2'][q_key].get('Dates', [])[0])):
                    #T2 = load_data['T2'][q_key].get('T2', [])[0][dataset]
                    #errors = load_data['T2'][q_key].get('Errors', [])[0][dataset]
                    date = datetime.datetime.fromtimestamp(load_data['T2'][q_key].get('Dates', [])[0][dataset])
                    I = self.process_h5_data(load_data['T2'][q_key].get('I', [])[0][dataset].decode())
                    Q = self.process_h5_data(load_data['T2'][q_key].get('Q', [])[0][dataset].decode())
                    delay_times = self.process_h5_data(load_data['T2'][q_key].get('Delay Times', [])[0][dataset].decode())
                    #fit = load_data['T2'][q_key].get('Fit', [])[0][dataset]
                    round_num = load_data['T2'][q_key].get('Round Num', [])[0][dataset]
                    batch_num = load_data['T2'][q_key].get('Batch Num', [])[0][dataset]
        
                    if len(I) > 0:
                        T2_class_instance = T2RMeasurement(q_key, self.outerFolder_save_plots, round_num, self.signal, self.save_figs, fit_data = True)
                        try:
                            fitted, t2r_est, t2r_err, plot_sig = T2_class_instance.t2_fit(delay_times, I, Q)
                        except Exception as e:
                            print('Fit didnt work due to error: ', e)
                            continue
                        T2_cfg = ast.literal_eval(self.exp_config['Ramsey_ge'].decode())
                        T2_class_instance.plot_results(I, Q, delay_times, date, fitted, t2r_est, t2r_err, plot_sig, config = T2_cfg, fig_quality=self.figure_quality)
                        del T2_class_instance
        
            del H5_class_instance

    def load_plot_save_t2e(self, exp_config):
        # -------------------------------------------------------Load/Plot/Save T2E------------------------------------------
        outerFolder_expt = self.outerFolder + "/Data_h5/T2E_ge/"
        h5_files = glob.glob(os.path.join(outerFolder_expt, "*.h5"))
        
        for h5_file in h5_files:
            save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]
            H5_class_instance = Data_H5(h5_file)
            load_data = H5_class_instance.load_from_h5(data_type=  'T2E', save_r = int(save_round))
            populated_keys = []
            for q_key in load_data['T2E']:
                # Access 'Dates' for the current q_key
                dates_list = load_data['T2E'][q_key].get('Dates', [[]])
        
                # Check if any entry in 'Dates' is not NaN
                if any(
                        not np.isnan(date)
                        for date in dates_list[0]  # Iterate over the first batch of dates
                ):
                    populated_keys.append(q_key)
        
            for q_key in populated_keys:
                for dataset in range(len(load_data['T2E'][q_key].get('Dates', [])[0])):
                    #T2 = load_data['T2E'][q_key].get('T2', [])[0][dataset]
                    #errors = load_data['T2E'][q_key].get('Errors', [])[0][dataset]
                    date = datetime.datetime.fromtimestamp(load_data['T2E'][q_key].get('Dates', [])[0][dataset])
                    I = self.process_h5_data(load_data['T2E'][q_key].get('I', [])[0][dataset].decode())
                    Q = self.process_h5_data(load_data['T2E'][q_key].get('Q', [])[0][dataset].decode())
                    delay_times = self.process_h5_data(load_data['T2E'][q_key].get('Delay Times', [])[0][dataset].decode())
                    #fit = load_data['T2E'][q_key].get('Fit', [])[0][dataset]
                    round_num = load_data['T2E'][q_key].get('Round Num', [])[0][dataset]
                    batch_num = load_data['T2E'][q_key].get('Batch Num', [])[0][dataset]
        
                    if len(I) > 0:
                        T2E_class_instance = T2EMeasurement(q_key, self.outerFolder_save_plots, round_num, self.signal, self.save_figs, fit_data = True)
                        try:
                            fitted, t2e_est, t2e_err, plot_sig = T2E_class_instance.t2_fit(delay_times, I, Q)
                        except Exception as e:
                            print('Fit didnt work due to error: ', e)
                            continue
                        T2E_cfg = ast.literal_eval(self.exp_config['SpinEcho_ge'].decode())
                        T2E_class_instance.plot_results(I, Q, delay_times, date, fitted, t2e_est, t2e_err, plot_sig, config = T2E_cfg, fig_quality=self.figure_quality)
                        del T2E_class_instance
        
            del H5_class_instance