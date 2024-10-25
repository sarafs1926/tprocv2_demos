from section_001_time_of_flight import TOFExperiment
from section_002_res_spec_ge_mux import ResonanceSpectroscopy
from section_003_qubit_spec_ge import QubitSpectroscopy
from section_004_amp_rabi_ge import AmplitudeRabiExperiment
from section_005_T1_ge import T1Measurement
from system_config import *

# N benchmark
n = 1

expt_name = "T1_ge"
j = 0
while j < n:
    j += 1
    for QubitIndex in range(0,2):

        QubitIndex = 1 #only change when you figure out how to not set everything else to0 except q2 in mask for resonator


        tof = TOFExperiment(QubitIndex, outerFolder)
        #---------------------res spec---------------------
        res_spec = ResonanceSpectroscopy(QubitIndex, outerFolder)
        res_spec.run(soccfg, soc)
        #--------------------qubit spec--------------------
        q_spec = QubitSpectroscopy(QubitIndex,outerFolder)
        q_spec.run(soccfg, soc)
        #-----------------------rabi-----------------------
        rabi = AmplitudeRabiExperiment(QubitIndex, outerFolder)
        #------------------------t1-------------------------
        t1 = T1Measurement(QubitIndex, outerFolder, expt_name)
        t1.run(soccfg, soc)