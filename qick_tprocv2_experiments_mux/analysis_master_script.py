from syspurpose.files import three_way_merge

from analysis_001_plot_all_RR_h5 import PlotAllRR
from analysis_002_res_centers_vs_time_plots import ResonatorFreqVsTime
from analysis_003_q_freqs_vs_time_plots import QubitFreqsVsTime
from analysis_004_pi_amp_vs_time_plots import PiAmpsVsTime
from analysis_006_T1_vs_time_plots import T1VsTime
from analysis_005_Qtemp_vs_time_plots import QTempsVsTime
from analysis_007_T2R_vs_time_plots import T2rVsTime
from analysis_008_T2E_vs_time_plots import T2eVsTime
from analysis_009_T1_hist_cumul_err_plots import T1HistCumulErrPlots
from analysis_010_T2R_hist_cumul_err_plots import T2rHistCumulErrPlots
from analysis_011_T2E_hist_cumul_err_plots import T2eHistCumulErrPlots
from analysis_012_save_run_data import SaveRunData
from analysis_013_update_saved_run_data_notes import UpdateNote
from analysis_014_temperature_calcsandplots import TempCalcAndPlots
from analysis_015_plot_all_run_stats import CompareRuns
from analysis_016_metrics_vs_temp import (ResonatorFreqVsTemp, QubitFreqsVsTemp,
                                          PiAmpsVsTemp, T1VsTemp, T2rVsTemp, T2eVsTemp)
from section_008_save_data_to_h5 import Data_H5
from analysis_000_load_configs import LoadConfigs

###################################################### Set These #######################################################
save_figs = True
fit_saved = False
show_legends = False
signal = 'None'
number_of_qubits = 6
run_number = 2 #starting from first run with qubits
figure_quality = 100 #ramp this up to like 500 for presentation plots
final_figure_quality = 500
run_name = '6transmon_run5'
run_notes = ('Added more eccosorb filters and a lpf on mxc before and after the device. Added thermometry '
             'next to the device') #please make it brief for the plot
top_folder_dates = ['2024-12-09', '2024-12-10', '2024-12-11', '2024-12-12', '2024-12-13', '2024-12-14', '2024-12-15',
                    '2024-12-16', '2024-12-17', '2024-12-18', '2024-12-19', '2024-12-20']
#top_folder_dates = ['2024-12-18', '2024-12-19']
###################################### 00: Load Configs for Plotting Titles ############################################
date = '2024-12-16'  #only plot all of the data for one date at a time because there is a lot
outerFolder = f"/data/QICK_data/{run_name}/" + date + "/"
config_loader = LoadConfigs(outerFolder)
sys_config, exp_config = config_loader.run()

# ################################################### 01: Plot Everything ################################################
# outerFolder_save_plots = f"/data/QICK_data/{run_name}/" + date + "_plots/"
# plotter = PlotAllRR(date, figure_quality, save_figs, fit_saved, signal, run_name, number_of_qubits, outerFolder,
#                  outerFolder_save_plots, exp_config)
# plotter.run(plot_res_spec = False, plot_q_spec = False, plot_rabi = False, plot_ss = True, plot_t1 = False,
#             plot_t2r = False, plot_t2e = False)
#
# ########################################## 02: Resonator Freqs vs Time Plots ###########################################
# res_spec_vs_time = ResonatorFreqVsTime(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates,
#                                        save_figs, fit_saved, signal, run_name, exp_config)
# date_times, res_freqs = res_spec_vs_time.run()
# res_spec_vs_time.plot(date_times, res_freqs, show_legends)
#
# ############################################ 03: Qubit Freqs vs Time Plots #############################################
# q_spec_vs_time = QubitFreqsVsTime(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates,
#                                   save_figs, fit_saved, signal, run_name, exp_config)
# date_times, q_freqs = q_spec_vs_time.run()
# q_spec_vs_time.plot(date_times, q_freqs, show_legends)
#
# ############################################## 04: Pi Amp vs Time Plots ###############################################
# pi_amps_vs_time = PiAmpsVsTime(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs,
#                               fit_saved,signal, run_name, exp_config)
# #can only have the 'depths' argument returned here if plot_depths=True, otherwise delete it
# #date_times, pi_amps, depths = pi_amps_vs_time.run(plot_depths=True)
# date_times, pi_amps = pi_amps_vs_time.run(plot_depths=False)
# pi_amps_vs_time.plot(date_times, pi_amps, show_legends)
# # pi_amps_vs_time.plot_vs_signal_depth(date_times, pi_amps, depths, show_legends)
# # pi_amps_vs_time.plot_signal_depth_vs_time(date_times, pi_amps, depths, show_legends)
# #
# # temps_class_obj = TempCalcAndPlots(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs, fit_saved,
# #                  signal, run_name, exp_config, outerFolder)
# #
# # temps, qubit_temp_dates = temps_class_obj.get_temps()
# # filtered_pi_amps = temps_class_obj.get_filtered_pi_amps(qubit_temp_dates, date_times, pi_amps)
# # pi_amps_vs_time.plot_vs_temps(date_times, filtered_pi_amps, temps, show_legends)
# #
# # ssf, qubit_ssf_dates = temps_class_obj.get_ssf()
# # filtered_pi_amps = temps_class_obj.get_filtered_pi_amps(qubit_ssf_dates, date_times, pi_amps)
# # pi_amps_vs_time.plot_vs_ssf(date_times, filtered_pi_amps, ssf, show_legends)
#
# ############# 05: Qubit Temp vs time (not working currently, use arianna code at bottom) #############################
# qtemp_vs_time = QTempsVsTime(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs,
#                                fit_saved,signal, run_name, exp_config)
#
# qubit_temp_dates, qubit_temperatures = qtemp_vs_time.run()
# qtemp_vs_time.plot(qubit_temp_dates, qubit_temperatures, show_legends)
#
# ################################################# 06: T1 vs Time Plots #################################################
# t1_vs_time = T1VsTime(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs, fit_saved,
#                  signal, run_name, exp_config)
# date_times, t1_vals = t1_vs_time.run()
# t1_vs_time.plot(date_times, t1_vals, show_legends)
#
# ################################################# 07: T2R vs Time Plots ################################################
# t2r_vs_time = T2rVsTime(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs, fit_saved,
#                  signal, run_name, exp_config)
# date_times, t2r_vals = t2r_vs_time.run()
# t2r_vs_time.plot(date_times, t2r_vals, show_legends)
#
# ################################################# 08: T2E vs Time Plots ################################################
# t2e_vs_time = T2eVsTime(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs, fit_saved,
#                  signal, run_name, exp_config)
# date_times, t2e_vals = t2e_vs_time.run()
# t2e_vs_time.plot(date_times, t2e_vals, show_legends)
#
# ############################################## 09: T1 hist/cumul/err Plots #############################################
# t1_distribution_plots = T1HistCumulErrPlots(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates,
#                                             save_figs, fit_saved, signal, run_name, run_notes, run_number, exp_config)
# dates, t1_vals, t1_errs = t1_distribution_plots.run()
# gaussian_dates, t1_std_values, t1_mean_values = t1_distribution_plots.plot(dates, t1_vals, t1_errs, show_legends)
#
# ############################################## 10: T2R hist/cumul/err Plots ############################################
# t2r_distribution_plots = T2rHistCumulErrPlots(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates,
#                                             save_figs, fit_saved, signal, run_name, exp_config)
# dates, t2r_vals, t2r_errs = t2r_distribution_plots.run(t1_vals)
# t2r_std_values, t2r_mean_values = t2r_distribution_plots.plot(dates, t2r_vals, t2r_errs, show_legends)
#
# ############################################## 11: T2E hist/cumul/err Plots ############################################
# t2e_distribution_plots = T2eHistCumulErrPlots(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates,
#                                             save_figs, fit_saved, signal, run_name, exp_config)
# dates, t2e_vals, t2e_errs = t2e_distribution_plots.run(t1_vals)
# t2e_std_values, t2e_mean_values = t2e_distribution_plots.plot(dates, t2e_vals, t2e_errs, show_legends)
#
# ############################ 12: Save the Key Statistics for This Run to Compare Later #################################
# #need to run 00, and 08-10 before this to get all of the variables
# saver = SaveRunData(run_number, run_notes)
# saver.run(gaussian_dates, t1_vals, t1_errs, t1_std_values, t1_mean_values, t2r_vals, t2r_errs, t2r_std_values,
#             t2r_mean_values, t2e_vals, t2e_errs, t2e_std_values, t2e_mean_values)
#
# # ################################## 13: Update Saved Run Notes For Comparison Plot ######################################
# # run_number_to_update = 2
# # new_run_notes = ("Added more eccosorb filters and a lpf on mxc before and after the device. Added thermometry "
# #                  "next to the device")
# # updater = UpdateNote(run_number_to_update, new_run_notes)
# # updater.run()
#
# ################################################ 14: Run Comparison Plots ##############################################
# run_number_list = [1,2]
# comparing_runs = CompareRuns(run_number_list)
# comparing_runs.run()

############################################### 15: Qubit Temperature Plots ############################################
temps_class_obj = TempCalcAndPlots(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs, fit_saved,
                 signal, run_name, exp_config, outerFolder)
temps_class_obj.run()

# ########################################## 16: Metrics Vs Temperature Plots ############################################
# res_spec_vs_temp = ResonatorFreqVsTemp(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates,
#                                        save_figs, fit_saved, signal, run_name, exp_config)
# date_times, res_freqs = res_spec_vs_temp.run()
# res_spec_vs_temp.plot(date_times, res_freqs, show_legends)
#
# q_spec_vs_temp = QubitFreqsVsTemp(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates,
#                                   save_figs, fit_saved, signal, run_name, exp_config)
# date_times, q_freqs = q_spec_vs_temp.run()
#
# pi_amps_vs_temp = PiAmpsVsTemp(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs,
#                                fit_saved,signal, run_name, exp_config)
# date_times, pi_amps = pi_amps_vs_temp.run()
# pi_amps_vs_temp.plot(date_times, pi_amps, show_legends)
#
# t1_vs_temp = T1VsTemp(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs, fit_saved,
#                  signal, run_name, exp_config)
# date_times, t1_vals = t1_vs_temp.run()
# t1_vs_temp.plot(date_times, t1_vals, show_legends)
#
# t2r_vs_temp = T2rVsTemp(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs, fit_saved,
#                  signal, run_name, exp_config)
# date_times, t2r_vals = t2r_vs_temp.run()
# t2r_vs_temp.plot(date_times, t2r_vals, show_legends)
#
# t2e_vs_temp = T2eVsTemp(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs, fit_saved,
#                  signal, run_name, exp_config)
# date_times, t2e_vals = t2e_vs_temp.run()
# t2e_vs_temp.plot(date_times, t2e_vals, show_legends)