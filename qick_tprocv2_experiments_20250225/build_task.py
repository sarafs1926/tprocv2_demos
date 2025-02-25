from qick import *
import numpy as np
import time
import json

# the main program class
from qick.asm_v2 import AveragerProgramV2
# for defining sweeps
from qick.asm_v2 import QickSpan, QickSweep1D


def add_single_qubit_experiment(expt_cfg, expt_name, QubitIndex):
    # In most of cases, we have single parameter experiments and two parameter experiment.
    # For single parameter experiments, using "start", "stop", and "step" to calculate "expts".
    # For two parameter experiment, we follow the order "length" or "gain" (inner, RAverage) -- "freq" (outer, Python)
    expt_cfg = expt_cfg[expt_name]

    # # Simple Single Parameter Experiments
    if "stop" in expt_cfg:
        try:
            start = expt_cfg["start"][QubitIndex]
        except TypeError as e:
            start = expt_cfg["start"]
        
        try:
            stop = expt_cfg["stop"][QubitIndex]
        except TypeError as e:
            stop = expt_cfg["stop"]

        
    
        # np.arrange only support "start" > "stop"
        if start >= stop:
            print("Warning: Start value is smaller than Stop value, and it will cause 'expts' = 0.")
        
        expt_cfg.update([("start",start), ("stop", stop)])

    elif "start" in expt_cfg: # for time rabi
        start = expt_cfg["start"][QubitIndex]
        expts = expt_cfg['expts'][QubitIndex]

    # Decide what parameter we are changing.
    if expt_name == 'res_spec_ge':
        expt_cfg.update([('res_freq_ge', QickSweep1D('freqloop', start, stop))])
    elif expt_name == 'qubit_spec_ge':
        expt_cfg.update([('qubit_freq_ge', QickSweep1D('freqloop', start, stop))])
    elif expt_name == 'time_rabi_ge' or expt_name == 'qubit_temp':
        expt_cfg.update([('expts', expts), ('start', start)])
    elif expt_name == 'power_rabi_ge':
        expt_cfg.update([('qubit_gain_ge', QickSweep1D('gainloop', start, stop))])
    elif expt_name == 'Ramsey_ge' or expt_name == 'SpinEcho_ge' or expt_name == 'T1_ge' or expt_name == 'Ramsey_ef':
        expt_cfg.update([('wait_time', QickSweep1D('waitloop', start, stop))])
    elif expt_name == 'res_spec_ef':
        expt_cfg.update([('res_freq_ef', QickSweep1D('freqloop', start, stop))])
    elif expt_name == 'qubit_spec_ef':
        expt_cfg.update([('qubit_freq_ef', QickSweep1D('freqloop', start, stop))])
    elif expt_name == 'power_rabi_ef':
        expt_cfg.update([('qubit_gain_ef', QickSweep1D('gainloop', start, stop))])
    elif expt_name == "Readout_Optimization_ge" or expt_name == "Readout_Optimization_gef": # for readout optimization experiments
        freq_start = expt_cfg["freq_start"][QubitIndex]
        freq_step_size = expt_cfg['freq_step_size'][QubitIndex]
        gain_start = expt_cfg["gain_start"][QubitIndex]
        gain_step_size = expt_cfg['gain_step_size'][QubitIndex]
        expt_cfg.update([('freq_start', freq_start), ('freq_step_size', freq_step_size),
                         ('gain_start', gain_start), ('gain_step_size', gain_step_size)])
    return expt_cfg


def add_multi_qubit_experiment(expt_cfg, expt_name, QubitIndex):
    # Build Multi-qubit Experiment Configuration 
    pass

def add_qubit_experiment(expt_cfg, expt_name, Qubit_list):
    # Classify Single Qubit Experiment or Many Qubit Experiment

    if isinstance(Qubit_list,int) == True:
        QubitIndex = Qubit_list
        expt_config = add_single_qubit_experiment(expt_cfg, expt_name, QubitIndex)
    elif isinstance(Qubit_list,list) == True:
        QubitIndex = Qubit_list
        expt_config = add_multi_qubit_experiment(expt_cfg, expt_name, QubitIndex)
    return expt_config