import sys
import os
import numpy as np
sys.path.append(os.path.abspath("/home/nexusadmin/Documents/GitHub/tprocv2_demos/qick_tprocv2_experiments_mux_nexus/"))
from system_config import QICK_experiment
from section_003_punch_out_ge_mux import TWPA_Sweep
import datetime

att_1=999
att_2=999

outerFolder = os.path.join("/home/nexusadmin/qick/NEXUS_sandbox/Data/", str(datetime.date.today()))

experiment = QICK_experiment(outerFolder)
TWPA_sweep   = TWPA_Sweep(outerFolder, experiment)

TWPA_freq = 7.47e9    #7.9159e9    #7.8613e9    #7.47e9
start_power, stop_power, num_points = -20, -20, 1
TWPA_sweep.run(experiment.soccfg, experiment.soc, start_power, stop_power, TWPA_freq, num_points, att_1, att_2, plot_Center_shift = True, plot_res_sweeps = True, plot_gains = True)

del TWPA_sweep
del experiment
