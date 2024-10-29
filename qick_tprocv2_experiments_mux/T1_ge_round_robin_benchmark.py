from section_001_time_of_flight import TOFExperiment
from section_002_res_spec_ge_mux import ResonanceSpectroscopy
from section_003_qubit_spec_ge import QubitSpectroscopy
from section_004_amp_rabi_ge import AmplitudeRabiExperiment
from section_005_T1_ge import T1Measurement
from system_config import *

# N benchmark
n = 1
j = 0
while j < n:
    j += 1
    for QubitIndex in range(0,2):

        QubitIndex = 1 #only change when you figure out how to not set everything else to0 except q2 in mask for resonator

        # ---------------------TOF---------------------
        tof = TOFExperiment(QubitIndex, outerFolder)
        tof.run(soccfg, soc)

        #---------------------Res spec---------------------
        res_spec = ResonanceSpectroscopy(QubitIndex, outerFolder)
        config = res_spec.run(soccfg, soc)

        #--------------------Qubit spec--------------------
        # right now this does not return qubit frequency or update the config with the found values, do we want to change that?
        q_spec = QubitSpectroscopy(QubitIndex, outerFolder, config)
        config = q_spec.run(soccfg, soc)

        #-----------------------Rabi-----------------------
        #right now this does not have updated fitting, we need to make sure this fit works every time
        rabi = AmplitudeRabiExperiment(QubitIndex, outerFolder, config)
        config = rabi.run(soccfg, soc)

        #------------------------T1-------------------------
        t1 = T1Measurement(QubitIndex, outerFolder, config)
        config = t1.run(soccfg, soc)