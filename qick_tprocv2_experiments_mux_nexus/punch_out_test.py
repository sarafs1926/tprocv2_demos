import sys
import os
import numpy as np
sys.path.append(os.path.abspath("/home/quietuser/Documents/GitHub/tprocv2_demos/qick_tprocv2_experiments_mux/"))
from system_config import QICK_experiment
from section_003_punch_out_ge_mux import PunchOut
import datetime

sweep_DAC_attenuator1 =[5] #np.linspace(5,20, 4)
sweep_DAC_attenuator2 =[10]#[15,20,25,30] #np.linspace(5,20,4)

outerFolder = "/data/QICK_data/6transmon_run4a/" + str(datetime.date.today()) + "/"
for att_1 in sweep_DAC_attenuator1:
    for att_2 in sweep_DAC_attenuator2:
        att_1 = round(att_1, 3)
        att_2 = round(att_2, 3)
        experiment = QICK_experiment(outerFolder, DAC_attenuator1 = att_1, DAC_attenuator2 = att_2)
        punch_out   = PunchOut(outerFolder, experiment)

        start_gain, stop_gain, num_points = 0.1, 1, 10
        punch_out.run(experiment.soccfg, experiment.soc, start_gain, stop_gain, num_points, att_1, att_2, plot_Center_shift = True, plot_res_sweeps = True)

        del punch_out
        del experiment
