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
from expt_config import *
from system_config import *

# ----- Experiment configurations ----- #
expt_name = "tof"
QubitIndex = QUBIT_INDEX
Qubit = 'Q' + str(QubitIndex)

exp_cfg = expt_cfg[expt_name]
q_config = all_qubit_state(system_config)
config = {**q_config['Q' + str(QubitIndex)], **exp_cfg}

# Update parameters to see TOF pulse with your setup
config.update([('res_freq', 7100), ('res_gain', 0.8), ('res_length', 0.5), ('ro_length', 1.5)])
print(config)

##################
# Define Program #
##################

class LoopbackProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        ro_ch = cfg['ro_ch']
        res_ch = cfg['res_ch']
        
        self.declare_gen(ch=res_ch, nqz=cfg['nqz_res'])
        # pynq configured
        # self.declare_readout(ch=ro_ch, length=cfg['ro_len'], freq=cfg['freq'], gen_ch=gen_ch)
        
        # tproc configured
        self.declare_readout(ch=ro_ch, length=cfg['ro_length'])
        self.add_readoutconfig(ch=ro_ch, name="myro", freq=cfg['res_freq'], gen_ch=res_ch)

        self.add_pulse(ch=res_ch, name="myconst", ro_ch=ro_ch, 
                       style="const", 
                       length=cfg['res_length'], 
                       freq=cfg['res_freq'], 
                       phase=cfg['res_phase'],
                       gain=cfg['res_gain'],
                      )
    
    def _body(self, cfg):
        self.send_readoutconfig(ch=cfg['ro_ch'], name="myro", t=0)
        self.pulse(ch=cfg['res_ch'], name="myconst", t=0)
        self.trigger(ros=[cfg['ro_ch']], pins=[0], t=0)


###################
# Run the Program
###################

prog =LoopbackProgram(soccfg, reps=1, final_delay=config['relax_delay'], cfg=config)
iq_list = prog.acquire_decimated(soc, soft_avgs=config['soft_avgs'])
t = prog.get_time_axis(ro_index=0)



#####################################
# ----- Saves data to a file ----- #
#####################################

prefix = str(datetime.date.today())
exp_name = expt_name + '_Q' + str(QubitIndex) + '_' + prefix
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


