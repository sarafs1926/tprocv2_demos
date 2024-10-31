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
        "start": [7147, 7169, 7202, 7226, 7262, 7285], # [MHz]
        "stop":  [7152, 7175, 7207, 7232, 7268, 7291], # [MHz]
        "steps": 1000,
        "relax_delay": 0.5, # [us]
    },

    "res_spec_mux": {
        "reps": 1000,
        "py_avg": 10,
        "start": -2.5, # [MHz]
        "expts":  500,
        "step": 0.01, # [us]
        "relax_delay": 0.5, # [us]
    },

    "qubit_spec_ge": {
        "reps": 1000,
        "py_avg": 10,
        "start": [2920, 3110, 3060, 3240, 3210, 3250], # [MHz]
        "stop": [3000, 3190, 3140, 3320, 3290, 3330], # [MHz]
        "steps": 500,
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
        "stop":  [300] * 6, # [us]
        "steps": 150,
        "ramsey_freq": 0.05,  # [MHz]
        "relax_delay": 1000, # [us]
        "wait_time": 0.0, # [us]
    },

    "SpinEcho_ge": {
        "reps": 100,
        "py_avg": 10,
        "start": [0.0] * 6, # [us]
        "stop":  [300] * 6, # [us]
        "steps": 150,
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
        "start": [7148, 7168, 7202, 7225, 7261, 7284], # [MHz]
        "stop":  [7151, 7174, 7207, 7231, 7267, 7290], # [MHz]
        "steps": 200,
        "relax_delay": 1000, # [us]
    },

    "qubit_spec_ef": {
        "reps": 100,
        "py_avg": 10,
        "start": [2750, 2950, 2850, 3050, 3050, 3000], # [MHz]
        "stop":  [2850, 3050, 3000, 3150, 3150, 3200], # [MHz]
        "steps": 500,
        "relax_delay": 1000, # [us]
    },

    "qubit_temp": {
        "reps": 1000,
        "py_avg": 10,
        "start": [0.02] * 6, # [us]
        "expts":  [100] * 6,
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
        "ge_meas": True, # Do we want to ro using ge res freq
    },

    "IQ_plot":{
        "steps": 10000, # shots
        "py_avg": 1,
        "reps": 1,
        "relax_delay": 1000, # [us]
        "SS_ONLY": True,
    },

    "Readout_Optimization_ge":{
        "steps": 10000, # shots
        "py_avg": 1,
        "gain_start" : [0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
        "gain_step_size" : [0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
        "gain_steps" : 10,
        "freq_start" : [7147.8, 7171.2, 7204.0, 7228.6, 7264.1, 7287.2],
        "freq_step_size" : [0.05, 0.025, 0.025, 0.025, 0.025, 0.025],
        "freq_steps" : 30,
        "relax_delay": 1000, # [us]
    },
    
}