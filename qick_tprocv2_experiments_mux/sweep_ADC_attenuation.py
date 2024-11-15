import sys
import os
import numpy as np
sys.path.append(os.path.abspath("/home/quietuser/Documents/GitHub/tprocv2_demos/qick_tprocv2_experiments_mux/"))
from system_config import QICK_experiment
from section_003_punch_out_ge_mux import PunchOut
import datetime

sweep_ADC_attenuator1 =[10,15,20,25,30] #np.linspace(5,20, 4)

outerFolder = "/data/QICK_data/6transmon_run4a/" + str(datetime.date.today()) + "/"
for att_1 in sweep_ADC_attenuator1:
    experiment = QICK_experiment(outerFolder, ADC_attenuator = att_1)
    #put Arianna style code here
    del experiment