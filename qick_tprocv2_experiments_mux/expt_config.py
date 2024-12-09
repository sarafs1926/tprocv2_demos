from qick import *
import numpy as np

VNA_res = np.array([6.20905, 6.26145, 6.321265, 6.401472, 6.467723, 6.5209414])*1000  # run 5
# VNA_res = np.array([6191.519, 6216, 6292.321, 6405.85, 6432.959, 6468.441,]) # run 4a

VNA_qubit = np.array([4184.10, 3821.43, 4156.88, 4459.12, 4471.18, 4998.04])  # Freqs of Qubit g/e Transition

expt_cfg = {
    "tof": {
        "reps": 1, #reps doesnt make a difference here, leave it at 1
        "soft_avgs": 300,
        "relax_delay": 0,  # [us]
    },

    "res_spec": {
        "reps": 100,
        "rounds": 1,
        "start": -3.5,  # [MHz]
        "step_size": 0.12,  # [MHz]
        "steps": 101,
        "relax_delay": 20,  # [us]
    },

    "qubit_spec_ge": {
        "reps": 100, #100
        "rounds": 1, #10
        "start": list(VNA_qubit-5), # [MHz]
        "stop":  list(VNA_qubit+5), # [MHz]
        "steps": 300,
        "relax_delay": 1000, # [us]
    },

    "power_rabi_ge": {
        "reps": 300, #100
        "rounds": 1, #5
        "start": [0.0] * 6, # [DAC units]
        "stop":  [1.0] * 6, # [DAC units]
        "steps": 100,
        "relax_delay": 1000, # [us]
    },

    "T1_ge": {
        "reps": 300, #300
        "rounds": 2, #1
        "start": [0.0] * 6,  # [us]
        "stop": [200,200,200,300,200,200], #[250.0] * 6,  # [us] ### Should be ~10x T1! Should change this per qubit.
        "steps": 80,
        "relax_delay": 1000,  # [us] ### Should be >10x T1!
        "wait_time": 0.0,  # [us]
    },

    "Ramsey_ge": {
        "reps": 300, #300
        "rounds": 2,#10
        "start": [0.0] * 6, # [us]
        "stop":  [60] * 6, # [us]
        "steps": 100,
        "ramsey_freq": 0.12,  # [MHz]
        "relax_delay": 1000, # [us] the time to wait to let the qubit to relax to gnd again after exciting it (make it way above T1)
        "wait_time": 0.0, # [us]
    },

    "SpinEcho_ge": {
        "reps": 300,
        "rounds": 2,
        "start": [0.0] * 6, # [us]
        "stop":  [60] * 6, # [us]
        "steps": 100,
        "ramsey_freq": 0.12,  # [MHz]
        "relax_delay": 1000, # [us]
        "wait_time": 0.0, # [us]
    },
#

#
#     "res_spec_ef": {
#         "reps": 100,
#         "py_avg": 10,
#         "start": [7148, 0, 7202, 0, 0, 0], # [MHz]
#         "stop":  [7151, 0, 7207, 0, 0, 0], # [MHz]
#         "steps": 200,
#         "relax_delay": 1000, # [us]
#     },
#
#     "qubit_spec_ef": {
#         "reps": 100,
#         "py_avg": 10,
#         "start": [2750, 0, 0, 0, 0, 0], # [MHz]
#         "stop":  [2850, 0, 0, 0, 0, 0], # [MHz]
#         "steps": 500,
#         "relax_delay": 1000, # [us]
#     },
#
#     "qubit_temp": {
#         "reps": 100,
#         "py_avg": 10,
#         "start": [0.02] * 6, # [us]
#         "expts":  [200] * 6,
#         "step": 0.02, # [us]
#         "relax_delay": 1000, # [us]
#     },
#
#     "power_rabi_ef": {
#         "reps": 1000,
#         "py_avg": 10,
#         "start": [0.0] * 6, # [DAC units]
#         "stop":  [1.0] * 6, # [DAC units]
#         "steps": 100,
#         "relax_delay": 1000, # [us]
#     },
#
#     "Ramsey_ef": {
#         "reps": 100,
#         "py_avg": 10,
#         "start": [0.0] * 6, # [us]
#         "stop":  [100] * 6, # [us]
#         "steps": 100,
#         "ramsey_freq": 0.05,  # [MHz]
#         "relax_delay": 1000, # [us]
#         "wait_time": 0.0, # [us]
#     },
#
#     "IQ_plot":{
#         "steps": 5000, # shots
#         "py_avg": 1,
#         "reps": 1,
#         "relax_delay": 1000, # [us]
#         "SS_ONLY": False,
#     },
#
    "Readout_Optimization":{
        "steps": 3000, # shots
        "py_avg": 1,
        "gain_start" : [0, 0, 0, 0],
        "gain_stop" : [1, 0, 0, 0],
        "gain_step" : 0.1,
        "freq_start" : [6176.0, 0, 0, 0],
        "freq_stop" : [6178.0, 0, 0, 0],
        "freq_step" : 0.1,
        "relax_delay": 1000, # [us]
    },
#
}