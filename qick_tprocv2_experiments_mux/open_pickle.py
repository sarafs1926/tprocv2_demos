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

expt_name = "T1_ge"
outerFolder_expt = outerFolder + "/" + expt_name + "/"
file_name = "T1_benchmark_Q1_40000x.pkl"

with open(outerFolder_expt+file_name, "rb") as f:
    T1_est = pickle.load(f)
    T1_err = pickle.load(f)
    dates = pickle.load(f)

print(T1_est, T1_err, dates)

