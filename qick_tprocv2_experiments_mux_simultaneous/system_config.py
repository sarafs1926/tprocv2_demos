from qick import *
from tprocv2_demos.qick_tprocv2_experiments_mux_simultaneous.socProxy import makeProxy
import os
import datetime

def create_folder_if_not_exists(folder_path):
    """Creates a folder at the given path if it doesn't already exist."""
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

# Where do you want to save data
prefix = str(datetime.date.today())
outerFolder = "/home/nexusadmin/qick/NEXUS_sandbox/Data/" + prefix + "/"
create_folder_if_not_exists(outerFolder)

# Make proxy to the QICK
soc, soccfg = makeProxy()
print(soccfg)

GEN_CH8 = 8
GEN_CH10 = 10
GEN_CH12 = 12
GEN_CH14 = 14
MIXMUXGEN_CH = 4
MUXRO_CH = [2,3,4,5]
# Qubit you want to work with
QUBIT_INDEX = 0

# Hardware Configuration
hw_cfg={
    # DAC
    "qubit_ch": [GEN_CH8]*4, # Qubit Channel Port, Full-speed DAC
    "res_ch": [MIXMUXGEN_CH]*4, # Single Tone Readout Port, MUX DAC
    # "qubit_ch_ef": [GEN_CH5]*6, # Qubit ef Channel, Full-speed DAC
    "nqz_qubit": 2,
    "nqz_res": 2,
    # ADC
    "ro_ch": [MUXRO_CH] * 4 # MUX readout channel
       }

# Readout Configuration
readout_cfg={
    "trig_time": 0.75, # [Clock ticks] - get this value from TOF experiment
    "ro_length": 4.0, # [us]

    # Changes related to the resonator output channel
    "mixer_freq": 5550,  # [MHz]
    "res_freq_ge": [6187.191, 5827.678, 6074.095, 5958.673],
    # "res_gain_ge": [1] + [0]*5,
    "res_gain_ge": [0,1,0,0], #[0] + [1] + [0] * 2,
    # "res_gain_ge": [1,1,0.7,0.7,0.7,1], #[0.4287450656184295, 0.4903077560386716, 0.4903077560386716, 0.3941941738241592, 0.3941941738241592, 0.4903077560386716],  # DAC units
    # "res_freq_ef": [7149.44, 0, 0, 0, 0, 0], # [MHz]
    # "res_gain_ef": [0.6, 0, 0, 0, 0, 0], # [DAC units]
    "res_length": 4.0, # [us] (1.0 for res spec)
    "res_phase": [0, 0, 0, 0], # Rotation Angle From QICK Function
    "ro_phase": [0, 0, 0, 0],  # Rotation Angle From QICK Function
    # "threshold": [0, 0, 0, 0, 0, 0], # Threshold for Distinguish g/e, from QICK Function
    }

# Qubit Configuration
qubit_cfg={
    "qubit_freq_ge": [4909, 4749.4, 4569, 4759], # Freqs of Qubit g/e Transition
    "qubit_gain_ge": [1]*4, #[0.4287450656184295, 0.4287450656184295, 0.4903077560386716, 0.6, 0.4903077560386716, 0.4287450656184295], # For spec pulse
    "qubit_length_ge": 15, # [us] for spec Pulse
    "qubit_phase": 0, # [deg]
    "sigma": [0.08, 0.15, 0.1, 0.08], # [us] for Gaussian Pulse
    "pi_amp": [0.8, 0.74, 0.8, 0.8],
    "mixer_freq_q": 4300,  # [MHz]

}
