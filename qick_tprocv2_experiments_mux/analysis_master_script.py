from analysis_001_plot_all_RR_h5 import PlotAllRR
from analysis_002_res_centers_vs_time_plots import ResonatorFreqVsTime
from analysis_003_q_freqs_vs_time_plots import QubitFreqsVsTime
from analysis_004_pi_amp_vs_time_plots import PiAmpsVsTime
from analysis_005_T1_vs_time_plots import T1VsTime
from analysis_006_T2R_vs_time_plots import T2rVsTime
from analysis_007_T2E_vs_time_plots import T2eVsTime
from analysis_008_T1_hist_cumul_err_plots import T1HistCumulErrPlots
from analysis_009_T2R_hist_cumul_err_plots import T2rHistCumulErrPlots
from analysis_010_T2E_hist_cumul_err_plots import T2eHistCumulErrPlots
from analysis_011_save_run_data import SaveRunData
from analysis_014_plot_all_run_stats import CompareRuns
from analysis_012_update_saved_run_data_notes import UpdateNote

####################################################### Configs ########################################################
save_figs = True
fit_saved = False
show_legends = False
signal = 'None'
number_of_qubits = 6
run_number = 2 #starting from first run with qubits
figure_quality = 100 #ramp this up to like 500 for presentation plots
final_figure_quality = 500
run_name = '6transmon_run5'
run_notes = 'Added eccosorb filters' #please make it brief for the plot
top_folder_dates = ['2024-12-10', '2024-12-11', '2024-12-12', '2024-12-13']

##################################################### Plot Everything ##################################################
date = '2024-12-10'
outerFolder = f"/data/QICK_data/{run_name}/" + date + "/"
outerFolder_save_plots = f"/data/QICK_data/{run_name}/" + date + "_plots/"
plotter = PlotAllRR(date, figure_quality, save_figs, fit_saved, signal, run_name, number_of_qubits, outerFolder,
                 outerFolder_save_plots)
plotter.run()

############################################ Resonator Freqs vs Time Plots #############################################
res_spec_vs_time = ResonatorFreqVsTime(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates,
                                       save_figs, fit_saved, signal, run_name)
date_times, res_freqs = res_spec_vs_time.run()
res_spec_vs_time.plot(date_times, res_freqs, show_legends)

############################################## Qubit Freqs vs Time Plots ###############################################
q_spec_vs_time = QubitFreqsVsTime(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates,
                                  save_figs, fit_saved, signal, run_name)
date_times, q_freqs = q_spec_vs_time.run()
q_spec_vs_time.plot(date_times, q_freqs, show_legends)

################################################# Pi Amp vs Time Plots #################################################
pi_amps_vs_time = PiAmpsVsTime(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs,
                               fit_saved,signal, run_name)
date_times, pi_amps = pi_amps_vs_time.run()
pi_amps_vs_time.plot(date_times, pi_amps, show_legends)

################################################### T1 vs Time Plots ###################################################
t1_vs_time = T1VsTime(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs, fit_saved,
                 signal, run_name)
date_times, t1_vals = t1_vs_time.run()
t1_vs_time.plot(date_times, t1_vals, show_legends)

################################################### T2R vs Time Plots ##################################################
t2r_vs_time = T2rVsTime(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs, fit_saved,
                 signal, run_name)
date_times, t2r_vals = t2r_vs_time.run()
t2r_vs_time.plot(date_times, t2r_vals, show_legends)

################################################### T2E vs Time Plots ##################################################
t2e_vs_time = T2eVsTime(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates, save_figs, fit_saved,
                 signal, run_name)
date_times, t2e_vals = t2e_vs_time.run()
t2e_vs_time.plot(date_times, t2e_vals, show_legends)

################################################ T1 hist/cumul/err Plots ###############################################
t1_distribution_plots = T1HistCumulErrPlots(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates,
                                            save_figs, fit_saved, signal, run_name, run_notes, run_number)
dates, t1_vals, t1_errs = t1_distribution_plots.run()
gaussian_dates, t1_std_values, t1_mean_values = t1_distribution_plots.plot(dates, t1_vals, t1_errs, show_legends)

################################################ T2R hist/cumul/err Plots ##############################################
t2r_distribution_plots = T2rHistCumulErrPlots(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates,
                                            save_figs, fit_saved, signal, run_name)
dates, t2r_vals, t2r_errs = t2r_distribution_plots.run()
t2r_std_values, t2r_mean_values = t2r_distribution_plots.plot(dates, t2r_vals, t2r_errs, show_legends)

################################################ T2E hist/cumul/err Plots ##############################################
t2e_distribution_plots = T2eHistCumulErrPlots(figure_quality, final_figure_quality, number_of_qubits, top_folder_dates,
                                            save_figs, fit_saved, signal, run_name)
dates, t2e_vals, t2e_errs = t2e_distribution_plots.run()
t2e_std_values, t2e_mean_values = t2e_distribution_plots.plot(dates, t2e_vals, t2e_errs, show_legends)

############################## Save the Key Statistics for This Run to Compare Later ###################################
SaveRunData(gaussian_dates, t1_vals, t1_errs, t1_std_values, t1_mean_values, t2r_vals, t2r_errs, t2r_std_values,
            t2r_mean_values, t2e_vals, t2e_errs, t2e_std_values, t2e_mean_values)

#################################### Update Saved Run Notes For Comparison Plot ########################################
run_number_to_update = 2
new_run_notes = ("Added more eccosorb filters and a lpf on mxc before and after the device. Added thermometry "
                 "next to the device")
updater = UpdateNote(run_number_to_update, new_run_notes)
updater.run()

################################################## Run Comparison Plots ################################################
run_number_list = [1,2]
comparing_runs = CompareRuns(run_number_list)
comparing_runs.run()