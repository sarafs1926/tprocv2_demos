import sys
import os
sys.path.append(os.path.abspath("/home/quietuser/Documents/GitHub/tprocv2_demos/qick_tprocv2_experiments_mux/"))
from section_001_time_of_flight import TOFExperiment
from section_002_res_spec_ge_mux import ResonanceSpectroscopy
from section_003_qubit_spec_ge import QubitSpectroscopy
from section_004_amp_rabi_ge import AmplitudeRabiExperiment
from section_006_T1_ge import T1Measurement
from system_config import *



# N benchmark
n = 100000
j = 0
Qs=[0,1,2,3,4,5]

while j < n:
    j += 1
    for QubitIndex in Qs:

        # ---------------------TOF---------------------
        #tof = TOFExperiment(QubitIndex, outerFolder, j)
        #tof.run(soccfg, soc)

        #---------------------Res spec---------------------
        res_spec = ResonanceSpectroscopy(QubitIndex, outerFolder, j)
        res_freqs = res_spec.run(soccfg, soc)

        #--------------------Qubit spec--------------------
        # Right now this does not return qubit frequency or update the config with the found values, do we want to change that?
        #only need to update the res_freqs here, it will update the imported config and will change all following classes
        q_spec = QubitSpectroscopy(QubitIndex, outerFolder, res_freqs, j)
        q_spec.run(soccfg, soc)

        #-----------------------Rabi-----------------------
        # Right now this does not have updated fitting, we need to make sure this fit works every time
        rabi = AmplitudeRabiExperiment(QubitIndex, outerFolder, j)
        rabi.run(soccfg, soc)

        #------------------------T1-------------------------
        # Also need to update the fit here. Maybe do custom fits for all three of these classes (QSpec/Rabi/T1)
        t1 = T1Measurement(QubitIndex, outerFolder, j)
        t1.run(soccfg, soc)