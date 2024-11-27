import sys
import os
import numpy as np
sys.path.append(os.path.abspath("/home/quietuser/Documents/GitHub/tprocv2_demos/qick_tprocv2_experiments_mux/"))
from section_001_time_of_flight import TOFExperiment
from section_002_res_spec_ge_mux import ResonanceSpectroscopy
from section_005_single_shot_ge import SingleShot
from system_config import *

# N benchmark
n = 10
Qs = [3] #[0,1,2,3,4,5]
j = 0
angles = {i: [] for i in range(0, 6)}
angles_avg = {i: [] for i in range(0, 6)}
while j < n:
    j += 1
    for QubitIndex in Qs:

        # ---------------------TOF------------------------
        #tof = TOFExperiment(QubitIndex, outerFolder, j)
        #tof.run(soccfg, soc)

        #---------------------Res spec---------------------
        res_spec = ResonanceSpectroscopy(QubitIndex, outerFolder, j)
        res_freqs = res_spec.run(soccfg, soc)

        #-----------------Roll Signal into I---------------
        #get the average theta value, then use that to rotate the signal. Plug that value into system_config res_phase
        ss = SingleShot(QubitIndex, '/data', j, round(4, 3))
        fid, angle = ss.run(soccfg, soc)
        angles[QubitIndex].append(angle)

for i in range(0, 6):
    print(f'Q{i} angle samples: ', angles[i])
    print(f'Q{i} avg theta: ', np.average(angles[i]))
    angles_avg[i].append(np.average(angles[i]))