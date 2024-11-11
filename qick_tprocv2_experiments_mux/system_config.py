from qick import *
import sys
import os
sys.path.append(os.path.abspath("/home/quietuser/Documents/GitHub/tprocv2_demos"))
from tprocv2_demos.qick_tprocv2_experiments_mux.socProxy import makeProxy
import os
import datetime
import numpy as np

class QICK_experiment:
    def __init__(self, folder):
        # Where do you want to save data
        self.outerFolder = folder
        self.create_folder_if_not_exists(self.outerFolder)

        # Make proxy to the QICK
        self.soc, self.soccfg = makeProxy()
        #print(self.soccfg)

        self.FSGEN_CH = 0
        self.MIXMUXGEN_CH = 4
        self.MUXRO_CH = [2, 3, 4, 5, 6, 7]

        # Qubit you want to work with
        self.QUBIT_INDEX = 5

        # Hardware Configuration
        self.hw_cfg = {
            # DAC
            "qubit_ch": [self.FSGEN_CH] * 6,  # Qubit Channel Port, Full-speed DAC
            "res_ch": [self.MIXMUXGEN_CH] * 6,  # Single Tone Readout Port, MUX DAC
            # "qubit_ch_ef": [GEN_CH5]*6, # Qubit ef Channel, Full-speed DAC
            "nqz_qubit": 1,
            "nqz_res": 2,
            # ADC
            "ro_ch": [self.MUXRO_CH] * 6  # MUX readout channel
        }

        # Readout Configuration
        self.readout_cfg = {
            "trig_time": 0.75,  # [Clock ticks] - get this value from TOF experiment

            # Changes related to the resonator output channel
            "mixer_freq": 6000,  # [MHz]
            "res_freq_ge": [6191.459, 6216.02, 6292.281, 6405.81, 6432.799, 6468.281],  # MHz, latest
            #"res_freq_ge": [6191.479, 6216.040, 6292.301, 6405.83, 6433.059, 6468.901],  # MHz, old
            # "res_gain_ge": [1] + [0]*5,
            "res_gain_ge": [1, 1, 1, 0.7, 1, 1],
            # set_res_gain_ge(QUBIT_INDEX), #utomatically sets all gains to zero except for the qubit you are observing
            # "res_gain_ge": [1,1,0.7,0.7,0.7,1], #[0.4287450656184295, 0.4903077560386716, 0.4903077560386716, 0.3941941738241592, 0.3941941738241592, 0.4903077560386716],  # DAC units
            # "res_freq_ef": [7149.44, 0, 0, 0, 0, 0], # [MHz]
            # "res_gain_ef": [0.6, 0, 0, 0, 0, 0], # [DAC units]
            "res_length": 9.0,  # [us] (1.0 for res spec)
            "res_phase": [0] * 3 + [-1.975688 * 180 / np.pi] + [-0.480887 * 180 / np.pi] + [0],
            # [-0.1006 *360/np.pi, -2.412527*360/np.pi, -1.821284*360/np.pi, -1.90962*360/np.pi, -0.566479*360/np.pi, -0.5941687*360/np.pi], # Rotation Angle From QICK Function, is the ang of 10 ss angles per qubit
            # "res_phase": [0]*6,#[-0.1006 *360/np.pi, -2.412527*360/np.pi, -1.821284*360/np.pi, -1.90962*360/np.pi, -0.566479*360/np.pi, -0.5941687*360/np.pi], # Rotation Angle From QICK Function, is the ang of 10 ss angles per qubit
            "ro_phase": [0, 0, 0, 0, 0, 0],  # Rotation Angle From QICK Function
            # "threshold": [0, 0, 0, 0, 0, 0], # Threshold for Distinguish g/e, from QICK Function
        }

        # Qubit Configuration
        self.qubit_cfg = {
            "qubit_freq_ge": [4184.13, 3821.142, 4156.58, 4459.19, 4471.10, 4997.87], #new
            #"qubit_freq_ge": [4184.15, 3821.156, 4156.88, 4459.12, 4471.18, 4998.04],  # Freqs of Qubit g/e Transition, old
            # 4184.131522750292, 3821.142737167706, 4156.57897361413, 4459.1949814835025, 4471.099595815802, 4997.870725348132 new (full values)
            "qubit_gain_ge": [0.05] * 6,
            # [0.4287450656184295, 0.4287450656184295, 0.4903077560386716, 0.6, 0.4903077560386716, 0.4287450656184295], # For spec pulse
            "qubit_length_ge": 15,  # [us] for spec Pulse
            "qubit_phase": 0,  # [deg]
            "sigma": [0.08, 0.15, 0.1, 0.08, 0.12, 0.13],  # [us] for Gaussian Pulse
            "pi_amp": [0.8, 0.78, 0.8, 0.8, 0.8, 0.8],
        }

    def create_folder_if_not_exists(self, folder):
        """Creates a folder at the given path if it doesn't already exist."""
        if not os.path.exists(folder):
            os.makedirs(folder)

    def set_gain_filter_ge(self, QUBIT_INDEX, IndexGain = 1, num_qubits=6):
        """Sets the gain for the selected qubit to 1, others to 0."""
        filtered_gain_ge = [0] * num_qubits  # Initialize all gains to 0
        if 0 <= QUBIT_INDEX < num_qubits: #makes sure you are within the range of options
            filtered_gain_ge[QUBIT_INDEX] = IndexGain  # Set the gain for the selected qubit
        return filtered_gain_ge




