from numpy.ma.core import outer

from qick import *

# the main program class
from qick.asm_v2 import AveragerProgramV2
# for defining sweeps
from qick.asm_v2 import QickSpan, QickSweep1D

import json
import datetime
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import pickle

# Used for live plotting, need to run "python -m visdom.server" in the terminal and open the IP address in browser
# import visdom

from build_task import *
from build_state import *
from expt_config import *
from system_config import *

import pickle
import pandas as pd

# expt_name = "T1_ge"
# outerFolder_expt = outerFolder + expt_name + "/"
# file_name = outerFolder_expt + "T1_benchmark_Q1_2024-10-23_06-50-49_40000x.pkl"

#file_name = "/data/QICK_data/6transmon_run4a/2024-10-22/T1_ge/T1_benchmark_Q1_2024-10-22_23-40-29_40000x.pkl"
file_name = "/data/QICK_data/6transmon_run4a/2024-10-23/T1_ge/T1_benchmark_Q1_2024-10-23_15-09-25_40000x.pkl"

with open(file_name, "rb") as f:
    T1_est = pickle.load(f)
    T1_err = pickle.load(f)
    dates = pickle.load(f)
    I = pickle.load(f)
    Q = pickle.load(f)
    delay_times = pickle.load(f)

print(T1_est, T1_err, dates, I, Q, delay_times)
print("stop")
