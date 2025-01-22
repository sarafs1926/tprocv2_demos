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
from analysis_016_metrics_vs_temp import (ResonatorFreqVsTemp, GetThermData, QubitFreqsVsTemp,
                                          PiAmpsVsTemp, T1VsTemp, T2rVsTemp, T2eVsTemp)
from section_008_save_data_to_h5 import Data_H5
from analysis_000_load_configs import LoadConfigs
from analysis_017_plot_metric_dependencies import PlotMetricDependencies

###################################################### Set These #######################################################
save_figs = True
fit_saved = False
show_legends = False
signal = 'None'
number_of_qubits = 6
run_number = 2 #starting from first run with qubits
figure_quality = 100 #ramp this up to like 500 for presentation plots
final_figure_quality = 200
run_name = '6transmon_run4'
# run_notes = ('Added more eccosorb filters and a lpf on mxc before and after the device. Added thermometry '
#              'next to the device') #please make it brief for the plot
top_folder_dates = ['2024-11-21', '2024-11-23', '2024-11-24']

date = '2024-11-23'  #only plot all of the data for one date at a time because there is a lot
outerFolder = f"/data/QICK_data/{run_name}/" + date + "/"
config_loader = LoadConfigs(outerFolder)
sys_config, exp_config = config_loader.run()
########################################################################################################################

#--------------------------------------- Load ALl data for this run ----------------------------------------------------
r = 1 #Will be analyzing data for this run. Note: run 1 = run 4a at QUIET and run 2 = run 5a at QUIET
run_number_list = [1]
run_stats_folder = f"run_stats/run{r}/"
filename = run_stats_folder + 'experiment_data.h5'
compare_runs = CompareRuns(run_number_list) #class instance
data = compare_runs.load_from_h5(filename)
#print(data.keys()) #if you want to see the type of data that is contained inside the file

#-----------------------------------Extracting relevant Q3 data for our plots-------------------------------------------

# Extract date/times and T1 values for qubit 3
date_times_t1_q3 = data['date_times_t1']['2']
t1_vals_q3       = data['t1_vals']['2']

date_times_pi_amps_q3 = data['date_times_pi_amps']['2']
pi_amps_q3 = data['pi_amps']['2']

q_freqs_q3 = data['q_freqs']['2']
date_times_q_spec_q3 = data['date_times_q_spec']['2']

#----------------------------------Q3 metrics vs time in one plot-------------------------------
plotter = PlotMetricDependencies(run_name, number_of_qubits, final_figure_quality)

#function originally was set up for qubit 1, ignore that it says q1, just provide the right data.
plotter.plot_q1_temp_and_t1(q1_temp_times= None, q1_temps= None, q1_t1_times=date_times_t1_q3, q1_t1_vals=t1_vals_q3, temp_label="Q3 Qubit Temp (mK)", t1_label="T1 (µs)",
                            magcan_dates = None, magcan_temps = None, magcan_label = "Mag Can Temp (mK)",
                            mcp2_dates = None, mcp2_temps = None, mcp2_label = "MCP2 Temp (mK)",
                            Q1_freqs = q_freqs_q3, Q1_dates_spec = date_times_q_spec_q3, qspec_label = "Frequency (MHz)",
                            date_times_pi_amps_Q1 = date_times_pi_amps_q3, pi_amps_Q1 = pi_amps_q3, pi_amps_label = "Pi Amp (a.u.)")
