from qick import *

expt_cfg = {
    "tof": {
        "reps": 1,
        "soft_avgs": 1000,
        "relax_delay": 1000,  # [us] 
    },

    "res_spec_ge": {
        "reps": 1000,
        "py_avg": 10,
        "start": [7147, 0, 7202, 0, 0, 0], # [MHz]
        "stop":  [7152, 0, 7207, 0, 0, 0], # [MHz]
        "steps": 1000,
        "relax_delay": 0.5, # [us]
    },

    "qubit_spec_ge": {
        "reps": 1000,
        "py_avg": 10,
        "start": [2900, 0, 3060, 0, 0, 0],# [2920, 0, 3060, 0, 0, 0], # [MHz]
        "stop": [3000, 0, 3140, 0, 0, 0], # [MHz]
        "steps": 1000,
        "relax_delay": 1000, # [us]
    },

    "time_rabi_ge": {
        "reps": 100,
        "py_avg": 10,
        "start": [0.02] * 6, # [us]
        "expts":  [50] * 6,
        "step": 0.02, # [us]
        "relax_delay": 1000, # [us]
    },

    "power_rabi_ge": {
        "reps": 100,
        "py_avg": 10,
        "start": [0.0] * 6, # [DAC units]
        "stop":  [1.0] * 6, # [DAC units]
        "steps": 100,
        "relax_delay": 1000, # [us]
    },

    "Ramsey_ge": {
        "reps": 100,
        "py_avg": 10,
        "start": [0.0] * 6, # [us]
        "stop":  [100] * 6, # [us]
        "steps": 100,
        "ramsey_freq": 0.05,  # [MHz]
        "relax_delay": 1000, # [us]
        "wait_time": 0.0, # [us]
    },

    "SpinEcho_ge": {
        "reps": 100,
        "py_avg": 10,
        "start": [0.0] * 6, # [us]
        "stop":  [100] * 6, # [us]
        "steps": 100,
        "ramsey_freq": 0.05,  # [MHz]
        "relax_delay": 1000, # [us]
        "wait_time": 0.0, # [us]
    },

    "T1_ge":{
        "reps": 100,
        "py_avg": 10,
        "start": [0.0] * 6, # [us]
        "stop": [500.0] * 6, # [us]
        "steps": 50,
        "relax_delay": 1000, # [us]
        "wait_time": 0.0, # [us]
    },

    "res_spec_ef": {
        "reps": 100,
        "py_avg": 10,
        "start": [7148, 0, 7202, 0, 0, 0], # [MHz]
        "stop":  [7151, 0, 7207, 0, 0, 0], # [MHz]
        "steps": 200,
        "relax_delay": 1000, # [us]
    },

    "qubit_spec_ef": {
        "reps": 100,
        "py_avg": 10,
        "start": [2750, 0, 0, 0, 0, 0], # [MHz]
        "stop":  [2850, 0, 0, 0, 0, 0], # [MHz]
        "steps": 500,
        "relax_delay": 1000, # [us]
    },

    "qubit_temp": {
        "reps": 100,
        "py_avg": 10,
        "start": [0.02] * 6, # [us]
        "expts":  [200] * 6,
        "step": 0.02, # [us]
        "relax_delay": 1000, # [us]
    },

    "power_rabi_ef": {
        "reps": 1000,
        "py_avg": 10,
        "start": [0.0] * 6, # [DAC units]
        "stop":  [1.0] * 6, # [DAC units]
        "steps": 100,
        "relax_delay": 1000, # [us]
    },

    "Ramsey_ef": {
        "reps": 100,
        "py_avg": 10,
        "start": [0.0] * 6, # [us]
        "stop":  [100] * 6, # [us]
        "steps": 100,
        "ramsey_freq": 0.05,  # [MHz]
        "relax_delay": 1000, # [us]
        "wait_time": 0.0, # [us]
    },

    "IQ_plot":{
        "steps": 5000, # shots
        "py_avg": 1,
        "reps": 1,
        "relax_delay": 1000, # [us]
        "SS_ONLY": False,
    },

    "Readout_Optimization":{
        "steps": 20000, # shots
        "py_avg": 1,
        "gain_start" : [0, 0, 0, 0],
        "gain_stop" : [1, 0, 0, 0],
        "gain_step" : 0.1,
        "freq_start" : [6176.0, 0, 0, 0],
        "freq_stop" : [6178.0, 0, 0, 0],
        "freq_step" : 0.1,
        "relax_delay": 1000, # [us]
    },

}