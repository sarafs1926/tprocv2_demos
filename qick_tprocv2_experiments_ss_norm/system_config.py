from qick import *
from qick.pyro import make_proxy

soc, soccfg = make_proxy(ns_host="192.168.1.147", ns_port=8000, proxy_name="rfsoc")
print(soccfg)

# Where do you want to save data
DATA_PATH = "M:/malab\People\Santi\data\candle_qubit"

# DAC Signal Generating Channels
GEN_CH0 = 0
GEN_CH1 = 1
GEN_CH2 = 2
GEN_CH3 = 3
GEN_CH4 = 4
GEN_CH5 = 5
GEN_CH6 = 6
GEN_CH7 = 7
GEN_CH8 = 8
GEN_CH9 = 9
GEN_CH10 = 10
GEN_CH11 = 11
GEN_CH12 = 12
GEN_CH13 = 13
GEN_CH14 = 14
GEN_CH15 = 15

# ADC Readout Channels
RO_CH0 = 0
RO_CH1 = 1
RO_CH2 = 2
RO_CH3 = 3
RO_CH4 = 4
RO_CH5 = 5
RO_CH6 = 6
RO_CH7 = 7
RO_CH8 = 8
RO_CH9 = 9

# Qubit you want to work with
QUBIT_INDEX = 0
# Perform single shot for ge calibration?
SS = True

# Hardware Configuration
hw_cfg={
    # DAC
    "qubit_ch": [GEN_CH0]*6, # Qubit Channel Port, Full-speed DAC
    "res_ch": [GEN_CH1]*6, # Single Tone Readout Port, Full-speed DAC
    "qubit_ch_ef": [GEN_CH0]*6, # Qubit ef Channel, Full-speed DAC
    "mux_ch": GEN_CH4,

    "nqz_qubit": 1,
    "nqz_res": 2,
    # ADC
    "ro_ch": [RO_CH0] * 6, # tproc configured readout channel
    "mux_ro_chs": [RO_CH2, RO_CH3, RO_CH4, RO_CH5, RO_CH6, RO_CH7]
       }

# Readout Configuration
readout_cfg={
    "trig_time": 0.50, #0.75, # [Clock ticks] - get this value from TOF experiment
    "ro_length": 6.0, # [us]
    "mixer_freq": 7000, # [MHz] - used for mux_ch and interpolated_ch

    # Changes related to the resonator output channel
    "res_freq_ge": [7148.975, 7171.378, 7204.201, 7228.902, 7264.342, 7287.472], #7149.42 #[7149.387, 7171.327, 7204.104, 7228.775, 7264.330, 7287.479], #7149.72 [MHz]
    "res_gain_ge": [0.095, 0.05, 0.05, 0.05, 0.05, 0.05], # [DAC units]
    "res_freq_ef": [7149.096, 7171.048, 7203.846, 7228.515, 7263.947, 7287.106], # [MHz]
    "res_gain_ef": [0.05, 0.05, 0.05, 0.05, 0.05, 0.05], # [DAC units]
    "res_length": 6.0, # [us] (1.0 for res spec)
    "res_phase": [0, 0, 0, 0, 0, 0], # Rotation Angle From QICK Function 
    "threshold": [0, 0, 0, 0, 0, 0], # Threshold for Distinguish g/e, from QICK Function
    }

# Qubit Configuration
qubit_cfg={
    "qubit_freq_ge": [2963.28+0.015, 3156.624-0.12, 3098.97-0.172, 3285.28-0.084, 3255.62-0.12, 3294.95-0.191+0.023], # Freqs of Qubit g/e Transition
    "qubit_gain_ge": [0.344206, 0.427149, 0.286935, 0.300428, 0.151342, 0.221503],#[0.3172, 0.3414, 0.3256, 0.3324, 0.1598, 0.2343], [0.001, 0.001, 0.001, 0.01, 0.001, 0.0001] # Pulse Gain (qubit spec: [0.001, 0.001, 0.001, 0.01, 0.001, 0.0001])
    "qubit_length_ge": 25.0, #0.248, #15.0, # [us] for Constant Pulse
    "qubit_freq_ef": [2794.775-0.036, 2990.251-0.67, 2932.602, 3120.384, 3091.112, 0], # [MHz] Freqs of Qubit e/f Transition
    "qubit_gain_ef": [0.0891, 0.086, 0.03, 0.03, 0.03, 0.1],#[0.01, 0.05, 0.05, 0.05, 0.01, 0.5], # [DAC units] Pulse Gain
    "qubit_length_ef": 25.0, # [us] for Constant Pulse
    "qubit_phase": 0, # [deg]
    "sigma": [0.2/5, 0.2/5, 0.2/5, 1.0/5, 0.8/5, 1.0/5], # [us] for Gaussian Pulse
    "sigma_ef": [1.0/5, 1.0/5, 0.2/5, 0.2/5, 0.2/5, 1.0/5], # [us] for Gaussian Pulse
}
