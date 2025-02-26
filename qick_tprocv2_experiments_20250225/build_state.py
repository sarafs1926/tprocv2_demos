from qick import *
import copy
import json

# open configuration file
with open('system_config.json', 'r') as f:
    configs = json.load(f)['config']

NUM_QUBITS = configs['num_qubits']

# Add DAC and ADC Channels
def add_qubit_channel(system_config, QubitIndex):
    hw_config = copy.deepcopy(system_config['hw_cfg'])
    var = ["qubit_ch", "res_ch", "ro_ch", "qubit_ch_ef"]
    for Index in var:
        value = hw_config[Index][QubitIndex]
        hw_config.update([(Index,value)])
    return hw_config

# Add Readout Parameters
def add_readout_cfg(system_config, QubitIndex):
    readout_config = copy.deepcopy(system_config['readout_cfg'])
    var = ["res_freq_ge", "res_gain_ge",
           "res_freq_ef", "res_gain_ef",
           "res_phase", "threshold"]
    for Index in var:
        value = readout_config[Index][QubitIndex]
        readout_config.update([(Index,value)])
    return readout_config

# Add Qubit Frequency and Constant Drive Gain Parameters
def add_qubit_cfg(system_config, QubitIndex):
    qubit_config = copy.deepcopy(system_config['qubit_cfg'])
    var = ["qubit_freq_ge", "qubit_freq_ef", 
           "qubit_gain_ge", "qubit_gain_ef",
           "sigma_ge", "sigma_ef"]
    for Index in var:
        value = qubit_config[Index][QubitIndex]
        qubit_config.update([(Index,value)])
    return qubit_config

# Build a Single Qubit State Dictionary
def qubit_state(system_config, QubitIndex):
    hw_cfg = add_qubit_channel(system_config, QubitIndex)
    readout_cfg = add_readout_cfg(system_config, QubitIndex)
    qubit_cfg = add_qubit_cfg(system_config, QubitIndex)
    return {**hw_cfg, **readout_cfg, **qubit_cfg}

def all_qubit_state(system_config):
    state = {}
    for QubitIndex in range(NUM_QUBITS):
        Qi_state = copy.deepcopy(qubit_state(system_config, QubitIndex))
        state.update([("Q"+str(QubitIndex),Qi_state)])
    return state

# print(qubit_state(system_config,1))

#print(all_qubit_state(system_config))

# state = {"Q0": qubit_state(system_config,0),"Q1": qubit_state(system_config,1),
          #"Q2": qubit_state(system_config,2),"Q3": qubit_state(system_config,3)}
