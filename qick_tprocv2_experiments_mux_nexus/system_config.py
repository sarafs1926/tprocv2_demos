from qick import *
import sys
import os
sys.path.append(os.path.abspath("/home/quietuser/Documents/GitHub/tprocv2_demos"))
from tprocv2_demos.qick_tprocv2_experiments_mux_nexus.socProxy import makeProxy
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

        self.FSGEN_CH = 8 # 0 for "old QICK", 6 for RF board
        self.MIXMUXGEN_CH = 4 # Readout resonator DAC channel
        self.MUXRO_CH = [2, 3, 4, 5]
        # self.MUXRO_CH_RF = 5  # New variable that we need for QICK box

        # self.TESTCH_DAC = 5 # loopback channel for RF board
        # self.TESTCH_ADC = 0  # loopback channel for RF board
        # self.TESTCH_ADC_RF = 4  # New variable that we need for QICK box

        # From mux_simultaneous
        # GEN_CH8 = 8
        # GEN_CH10 = 10
        # GEN_CH12 = 12
        # GEN_CH14 = 14
        # MIXMUXGEN_CH = 4
        # MUXRO_CH = [2, 3, 4, 5]
        # # Qubit you want to work with
        # QUBIT_INDEX = 0

        ### NEW for the RF board
        # self.qubit_center_freq = 4400  # To be in the middle of the qubit freqs.
        # self.res_center_freq = 6330  # To be in the middle of the res freqs. 3000-5000 see nothing,6000 and 7000 see something, 8000+ see nothing
        # self.soc.rfb_set_gen_filter(self.MIXMUXGEN_CH, fc=self.res_center_freq / 1000, ftype='bandpass', bw=1.0)
        # self.soc.rfb_set_gen_filter(self.FSGEN_CH, fc=self.qubit_center_freq / 1000, ftype='bandpass', bw=1.0)
        # self.soc.rfb_set_ro_filter(self.MUXRO_CH_RF, fc=self.res_center_freq / 1000, ftype='bandpass', bw=1.0)
        # # Set attenuator on DAC.
        # self.soc.rfb_set_gen_rf(self.MIXMUXGEN_CH, self.DAC_attenuator1, self.DAC_attenuator2)  # Verified 30->25 see increased gain in loopback
        # self.soc.rfb_set_gen_rf(self.FSGEN_CH, 5, 4)  # Verified 30->25 see increased gain in loopback
        # # Set attenuator on ADC.
        # ### IMPORTANT: set this to 30 and you get 60 dB of warm gain. Set to 0 and you get 90 dB of warm gain
        # self.soc.rfb_set_ro_rf(self.MUXRO_CH_RF, self.ADC_attenuator)  # Verified 30->25 see increased gain in loopback


        # Qubit you want to work with
        self.QUBIT_INDEX = 0

        # Hardware Configuration
        self.hw_cfg = {
            # DAC
            "qubit_ch": [self.FSGEN_CH] * 4,  # Qubit Channel Port, Full-speed DAC
            "res_ch": [self.MIXMUXGEN_CH] * 4,  # Single Tone Readout Port, MUX DAC
            # "qubit_ch_ef": [GEN_CH5]*6, # Qubit ef Channel, Full-speed DAC
            "nqz_qubit": 2,
            "nqz_res": 2,
            # ADC
            "ro_ch": [self.MUXRO_CH] * 4  # MUX readout channel
        }

        # Readout Configuration
        self.readout_cfg = {
            "trig_time": 0.75,  # [Clock ticks] - get this value from TOF experiment
            # Changes related to the resonator output channel
            "mixer_freq": 5500, # [MHz]
            "res_freq_ge": [6187.191, 5827.678, 6074.095, 5958.673], #MHz
            "res_gain_ge": [0.15]*4, #[1, 1, 1, 1],
            "res_length": 10,  # [us] (1.0 for res spec)
            "res_phase": [0] * 4,
            "ro_phase": [0] * 4,  # Rotation Angle From QICK Function
        }

        # Qubit Configuration
        self.qubit_cfg = {
            "qubit_mixer_freq": 4300,  # [MHz]
            "qubit_freq_ge": [4909, 4749.4, 4569, 4759],  # Freqs of Qubit g/e Transition
            "qubit_gain_ge": [0.2] * 4, # [0.05] * 4
            "qubit_length_ge": 20,  # [us] for spec Pulse
            "qubit_phase": 0,  # [deg]
             "sigma": [0.04, 0.04, 0.03, 0.04], #[0.08, 0.15, 0.11, 0.09], # TO DO CHANGE THIS (11/26)
            "pi_amp": [1.0, 0.93, 0.77, 0.8], # TO DO CHANGE THIS (11/26)

            # "qubit_freqs_ge": [4909, 4749.4, 4569, 4759],  # Freqs of Qubit g/e Transition
            # "qubit_gains_ge": [1] * 4,  # [0.05] * 4
            # "qubit_phases": [0] * 4,  # [deg]
        }

    def create_folder_if_not_exists(self, folder):
        """Creates a folder at the given path if it doesn't already exist."""
        if not os.path.exists(folder):
            os.makedirs(folder)

    def mask_gain_res(self, QUBIT_INDEX, IndexGain = 1, num_qubits=4):
        """Sets the gain for the selected qubit to 1, others to 0."""
        filtered_gain_ge = [0] * num_qubits  # Initialize all gains to 0
        if 0 <= QUBIT_INDEX < num_qubits: #makes sure you are within the range of options
            filtered_gain_ge[QUBIT_INDEX] = IndexGain  # Set the gain for the selected qubit
        return filtered_gain_ge




