import os
folder = os.getcwd()
os.chdir(folder + '/qick_tprocv2_experiments')

from qick import *
from qick.pyro import make_proxy

from malab import SlabFile, get_next_filename
# for now, all the tProc v2 classes need to be individually imported (can't use qick.*)

# the main program class
from qick.asm_v2 import AveragerProgramV2
# for defining sweeps
from qick.asm_v2 import QickSpan, QickSweep1D

import json
import datetime

# Used for live plotting, need to run "python -m visdom.server" in the terminal and open the IP address in browser
import visdom

from build_state import *
from qick_setup_funcs import *

# connect to the rfsoc and print soccfg
from rfsoc_connect import *


# ----- Experiment configurations ----- #
expt_name = "tof"
config, expt_name = initialize_configs(expt_name)
print(expt_name + '\n')
print(config)

# Update parameters to see TOF pulse with your setup
if MUX == 'True': # MUX case
    config.update(
        [('mixer_freq', 6000),
        ('res_freq', [6000 + 25*(i+1) for i in range(NUM_QUBITS)]), 
        ('res_gain', [0.8] * NUM_QUBITS), 
        ('res_length', 0.5), #[0.5] * NUM_QUBITS), 
        ('ro_length', 1.5)])
else: # non-MUX case
    config.update(
        [('mixer_freq', 6000),
        ('res_freq', 6025),
        ('res_gain', 0.8), 
        ('res_length', 0.5), 
        ('ro_length', 1.5)])

print(config)

##################
# Define Program #
##################

class LoopbackProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        # initialize the readout channels and pulses for DAC and ADC
        initialize_ro_chs(self, cfg, suffix='')
    
    def _body(self, cfg):
        # readout using proper configurations for MUX and non-MUX
        readout(self, cfg)


###################
# Run the Program
###################

prog =LoopbackProgram(soccfg, reps=config['reps'], final_delay=config['relax_delay'], cfg=config)
iq_list = prog.acquire_decimated(soc, soft_avgs=config['soft_avgs'])
t = prog.get_time_axis(ro_index=0)



#####################################
# ----- Saves data to a file ----- #
#####################################

prefix = str(datetime.date.today())
exp_name = expt_name + '_Q' + str(QUBIT_INDEX) + '_' + prefix
print('Experiment name: ' + exp_name)

data_path = DATA_PATH

fname = get_next_filename(data_path, exp_name, suffix='.h5')
print('Current data file: ' + fname)

with SlabFile(data_path + '\\' + fname, 'a') as f:
    # 'a': read/write/create

    # - Adds data to the file - #
    f.append('t', t)

    f.append('avgi', iq_list[0].T[0])
    f.append('avgq', iq_list[0].T[1])

    f.attrs['config'] = json.dumps(config) # cant save configs yet with QickParams

data = data_path + '\\' + fname
print('data path:', data)


