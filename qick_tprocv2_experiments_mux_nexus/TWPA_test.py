import sys
import os
import numpy as np
sys.path.append(os.path.abspath("/home/nexusadmin/Documents/GitHub/tprocv2_demos/qick_tprocv2_experiments_mux_nexus/"))
from system_config import QICK_experiment
from section_003_punch_out_ge_mux import TWPASweep, TWPAConsistency
import datetime

outerFolder = os.path.join("/home/nexusadmin/qick/NEXUS_sandbox/Data/", str(datetime.date.today()))

experiment = QICK_experiment(outerFolder)

## Make Sure only 1 of the following 3 blocks in uncommented and check input values!

# #### If Running TWPA Power Sweep #####
# TWPA_sweep   = TWPASweep(outerFolder, experiment)
# pump_freq = 7.8238e9  #7.9159e9    #7.8613e9    #7.47e9
# start_power, stop_power, power_num_points = -13.0, -11.9, 11 #12, 12.5
# TWPA_sweep.run(experiment.soccfg, experiment.soc, start_power, stop_power, pump_freq, power_num_points,
#                power_sweep = True, plot_res_sweeps = True, plot_gains = True)
# del TWPA_sweep
# del experiment


#### If Running TWPA Freq Sweep #####
# TWPA_sweep   = TWPASweep(outerFolder, experiment)
# pump_power = -12.5
# start_freq, stop_freq, freq_num_points = 7.81e9, 7.83e9, 16 ##Don't go above 7.9159 GHz #7.811, 7.844
# #7.82 and 7.84 for -12.5
# # 7.82ish, 7.815, 7.83, for -12.0
# TWPA_sweep.run(experiment.soccfg, experiment.soc, start_freq, stop_freq, pump_power, freq_num_points,
#                power_sweep = False, plot_res_sweeps = True, plot_gains = True)
# del TWPA_sweep
# del experiment


# #### If Running TWPA Consistency Check #####
TWPA_cons = TWPAConsistency(outerFolder, experiment)
pump_power = -12.5
pump_freq = 7.829e9
num_reps = 7
TWPA_cons.run(experiment.soccfg, experiment.soc, pump_power, pump_freq, num_reps, plot_res_sweeps=True, plot_gains=True)
del TWPA_cons
del experiment