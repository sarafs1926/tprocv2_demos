from qick import *
import sys
import os
sys.path.append(os.path.abspath("/home/quietuser/Documents/GitHub/tprocv2_demos"))
from tprocv2_demos.qick_tprocv2_experiments_mux.socProxy import makeProxy
import os
import datetime
import numpy as np

class QICK_experiment:
    def __init__(self, folder, DAC_attenuator1 = 5, DAC_attenuator2 = 15, ADC_attenuator = 15 ):
        # Where do you want to save data
        self.outerFolder = folder
        self.create_folder_if_not_exists(self.outerFolder)

        # attenuation settings
        self.DAC_attenuator1 = DAC_attenuator1
        self.DAC_attenuator2 = DAC_attenuator2
        self.ADC_attenuator = ADC_attenuator

        # Make proxy to the QICK
        self.soc, self.soccfg = makeProxy()
        print(self.soccfg)

        self.FSGEN_CH = 6 # 0 for "old QICK", 6 for RF board
        self.MIXMUXGEN_CH = 4 # Readout resonator DAC channel
        self.MUXRO_CH = [2, 3, 4, 5, 6, 7]
        self.MUXRO_CH_RF = 5  # New variable that we need for QICK box

        ### NEW for the RF board
        self.qubit_center_freq = 4400  # To be in the middle of the qubit freqs.
        self.res_center_freq = 6330  # To be in the middle of the res freqs. 3000-5000 see nothing,6000 and 7000 see something, 8000+ see nothing
        self.soc.rfb_set_gen_filter(self.MIXMUXGEN_CH, fc=self.res_center_freq / 1000, ftype='bandpass', bw=1.0)
        self.soc.rfb_set_gen_filter(self.FSGEN_CH, fc=self.qubit_center_freq / 1000, ftype='bandpass', bw=1.0)
        self.soc.rfb_set_ro_filter(self.MUXRO_CH_RF, fc=self.res_center_freq / 1000, ftype='bandpass', bw=1.0)
        # Set attenuator on DAC.
        self.soc.rfb_set_gen_rf(self.MIXMUXGEN_CH, self.DAC_attenuator1, self.DAC_attenuator2)  # Verified 30->25 see increased gain in loopback
        self.soc.rfb_set_gen_rf(self.FSGEN_CH, 5, 10)  # Verified 30->25 see increased gain in loopback
        # Set attenuator on ADC.
        ### IMPORTANT: set this to 30 and you get 60 dB of warm gain. Set to 0 and you get 90 dB of warm gain
        self.soc.rfb_set_ro_rf(self.MUXRO_CH_RF, self.ADC_attenuator)  # Verified 30->25 see increased gain in loopback


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
            "res_freq_ge": [6191.419, 6216.1, 6292.361, 6405.77, 6432.759, 6468.481],  # MHz, new
            #"res_freq_ge": [6191.439, 6216.0, 6292.261, 6405.79, 6432.899, 6468.501],  # MHz, old
            # "res_freq_ge": [6191.459, 6216.02, 6292.281, 6405.81, 6432.799, 6468.281],  # MHz, old
            #"res_freq_ge": [6191.479, 6216.040, 6292.301, 6405.83, 6433.059, 6468.901],  # MHz, old
            # "res_gain_ge": [1] + [0]*5,
            "res_gain_ge": [1, 1, 1, 1, 1, 1],
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
            "qubit_freq_ge": [4184.14, 3821.149, 4156.53, 4459.20, 4471.12, 4997.86],  # new
            #"qubit_freq_ge": [4184.14, 3821.144, 4156.57, 4459.19, 4471.12, 4997.86], #old
            #"qubit_freq_ge": [4184.13, 3821.142, 4156.58, 4459.19, 4471.10, 4997.87], #old
            #"qubit_freq_ge": [4184.15, 3821.156, 4156.88, 4459.12, 4471.18, 4998.04],  # Freqs of Qubit g/e Transition, old
            "qubit_gain_ge": [0.05] * 6,
            # [0.4287450656184295, 0.4287450656184295, 0.4903077560386716, 0.6, 0.4903077560386716, 0.4287450656184295], # For spec pulse
            "qubit_length_ge": 15,  # [us] for spec Pulse
            "qubit_phase": 0,  # [deg]
            "sigma": [0.08, 0.18, 0.14, 0.13, 0.18, 0.6],  # [us] for Gaussian Pulse
            "pi_amp": [0.92, 0.87, 0.75, 0.73, 0.77, 0.78],
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




