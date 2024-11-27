from qick import *
from qick.pyro import make_proxy
import numpy as np
# from expt_config import expt_cfg
import time

# the main program class
from qick.asm_v2 import AveragerProgramV2
# for defining sweeps
from qick.asm_v2 import QickSpan, QickSweep1D

# soc, soccfg = make_proxy(ns_host="192.168.1.144", ns_port=8000, proxy_name="rfsoc")
# print(soccfg)

def add_single_qubit_experiment(expt_cfg, expt_name, QubitIndex):
    # In most of cases, we have single parameter experiments and two parameter experiment.
    # For single parameter experiments, using "start", "stop", and "step" to calculate "expts".
    # For two parameter experiment, we follow the order "length" or "gain" (inner, RAverage) -- "freq" (outer, Python)

    #exp_cfg = expt_cfg[expt_name]
    import copy
    expt_cfg_deep_copy = copy.deepcopy(expt_cfg)
    exp_cfg = expt_cfg_deep_copy[expt_name]

    # # Simple Single Parameter Experiments
    if "stop" in exp_cfg:
        start = exp_cfg["start"][QubitIndex]
        stop = exp_cfg["stop"][QubitIndex]
    
        # np.arrange only support "start" > "stop"
        if start >= stop:
            print("Warning: Start value is smaller than Stop value, and it will cause 'expts' = 0.")

        exp_cfg.update([("start",start), ("stop", stop)]) #this line is perminantly updating tprocv2_demos.qick_tprocv2_experiments_mux.expt_config import expt_cfg


    elif "start" in exp_cfg: # for time rabi
        start = exp_cfg["start"][QubitIndex]
        expts = exp_cfg['expts'][QubitIndex]


    # else: # for single shot IQ plot

    # Decide what parameter we are changing.
    if expt_name == 'res_spec_ge':
        exp_cfg.update([('res_freq_ge', QickSweep1D('freqloop', start, stop))])
    elif expt_name == 'qubit_spec_ge':
        exp_cfg.update([('qubit_freq_ge', QickSweep1D('freqloop', start, stop))])
    elif expt_name == 'time_rabi_ge' or expt_name == 'qubit_temp':
        exp_cfg.update([('expts', expts), ('start', start)])
    elif expt_name == 'power_rabi_ge':
        exp_cfg.update([('qubit_gain_ge', QickSweep1D('gainloop', start, stop))])
    elif expt_name == 'Ramsey_ge' or expt_name == 'SpinEcho_ge' or expt_name == 'T1_ge' or expt_name == 'Ramsey_ef':
        exp_cfg.update([('wait_time', QickSweep1D('waitloop', start, stop))])
    elif expt_name == 'res_spec_ef':
        exp_cfg.update([('res_freq_ef', QickSweep1D('freqloop', start, stop))])
    elif expt_name == 'qubit_spec_ef':
        exp_cfg.update([('qubit_freq_ef', QickSweep1D('freqloop', start, stop))])
    elif expt_name == 'power_rabi_ef':
        exp_cfg.update([('qubit_gain_ef', QickSweep1D('gainloop', start, stop))])
    
    return exp_cfg


def add_multi_qubit_experiment(expt_cfg, expt_name, QubitIndex):
    # Build Multi-qubit Experiment Configuration 
    pass

def add_qubit_experiment(expt_cfg, expt_name, Qubit_list):
    # Classify Single Qubit Experiment or Many Qubit Experiment

    if isinstance(Qubit_list,int) == True:
        QubitIndex = Qubit_list
        exp_config = add_single_qubit_experiment(expt_cfg, expt_name, QubitIndex)
    elif isinstance(Qubit_list,list) == True:
        QubitIndex = Qubit_list
        exp_config = add_multi_qubit_experiment(expt_cfg, expt_name, QubitIndex)

    return exp_config
