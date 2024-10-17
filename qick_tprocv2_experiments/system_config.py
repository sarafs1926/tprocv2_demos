from qick import *
from qick.pyro import make_proxy

soc, soccfg = make_proxy(ns_host="192.168.1.144", ns_port=8000, proxy_name="rfsoc")
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

# ADC Readout Channels
RO_CH0 = 0
RO_CH1 = 1
RO_CH2 = 2
RO_CH3 = 3
RO_CH4 = 4
RO_CH5 = 5

# Qubit you want to work with
QUBIT_INDEX = 0

# Hardware Configuration
hw_cfg={
    # DAC
    "qubit_ch": [GEN_CH4]*6, # Qubit Channel Port, Full-speed DAC
    "res_ch": [GEN_CH6]*6, # Single Tone Readout Port, Full-speed DAC
    "qubit_ch_ef": [GEN_CH5]*6, # Qubit ef Channel, Full-speed DAC
    "nqz_qubit": 1,
    "nqz_res": 2,
    # ADC
    "ro_ch": [RO_CH1] * 6 # tproc configured readout channel
       }

# Readout Configuration
readout_cfg={
    "trig_time": 0.65, # [Clock ticks] - get this value from TOF experiment
    "ro_length": 10.0, # [us]

    # Changes related to the resonator output channel
    "res_freq_ge": [7142, 0, 7204.48, 0, 0, 0], #7149.72 [MHz]
    "res_gain_ge": [0.7, 0, 1.0, 0, 0, 0], # [DAC units]
    "res_freq_ef": [7149.44, 0, 0, 0, 0, 0], # [MHz]
    "res_gain_ef": [0.6, 0, 0, 0, 0, 0], # [DAC units]
    "res_length": 6.0, # [us] (1.0 for res spec)
    "res_phase": [0, 0, 0, 0, 0, 0], # Rotation Angle From QICK Function 
    "threshold": [0, 0, 0, 0, 0, 0], # Threshold for Distinguish g/e, from QICK Function
    }

# Qubit Configuration
qubit_cfg={
    "qubit_freq_ge": [2964.38 - 0.08, 0,3102.75,0,0,0], # Freqs of Qubit g/e Transition
    "qubit_gain_ge": [0.251, 0, 0.27, 0, 0, 0], # 0.251 0.405Time Rabi Constant Pulse Gain (qubit spec: 0.01)
    "qubit_length_ge": 15.0, #0.248, #15.0, # [us] for Constant Pulse
    "qubit_freq_ef": [2795.876,0,0,0,0,0], # Freqs of Qubit e/f Transition
    "qubit_gain_ef": [0.323, 0, 0, 0, 0, 0], # Qubit Spec and Time Rabi Constant Pulse Gain
    "qubit_length_ef": 15.0, # [us] for Constant Pulse
    "qubit_phase": 0, # [deg]
    "sigma": 2.0 / 5, # [us] for Gaussian Pulse
    "sigma_ef": 3.0 / 5, # [us] for Gaussian Pulse
}
