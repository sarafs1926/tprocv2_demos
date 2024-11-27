import sys
import os
import numpy as np
sys.path.append(os.path.abspath("/home/nexusadmin/Documents/GitHub/tprocv2_demos/qick_tprocv2_experiments_mux_nexus/"))
from system_config import QICK_experiment
from section_003_punch_out_ge_mux import PunchOut
import datetime

att_1=999
att_2=999

outerFolder = os.path.join("/home/nexusadmin/qick/NEXUS_sandbox/Data/", str(datetime.date.today()))

experiment = QICK_experiment(outerFolder)
punch_out   = PunchOut(outerFolder, experiment)

start_gain, stop_gain, num_points = 0.1, 1, 10
punch_out.run(experiment.soccfg, experiment.soc, start_gain, stop_gain, num_points, att_1, att_2, plot_Center_shift = True, plot_res_sweeps = True)

del punch_out
del experiment
