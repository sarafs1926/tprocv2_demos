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

from build_task import *
from build_state import *
from expt_config import *
from system_config import *

from qualang_tools.plot import Fit
import pprint as pp

from tqdm import tqdm

# ----- Experiment configurations ----- #
expt_name = "res_spec_mux"
# QubitIndex = QUBIT_INDEX
# Qubit = 'Q' + str(QubitIndex)

exp_cfg = expt_cfg[expt_name]
# q_config = all_qubit_state(system_config)
config = {**hw_cfg, **readout_cfg, **exp_cfg}
print(config)

##################
# Define Program #
##################

class SingleToneSpectroscopyProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        ro_chs = cfg['mux_ro_chs']
        gen_ch = cfg['mux_ch']
        
        self.declare_gen(ch=gen_ch, nqz=cfg['nqz_res'], ro_ch=ro_chs[0], 
                         mux_freqs=cfg['res_freq_ge'], 
                         mux_gains=cfg['res_gain_ge'],
                         mux_phases=cfg['res_phase'],
                         mixer_freq=cfg['mixer_freq'])
        for ch, f, ph in zip(ro_chs, cfg['res_freq_ge'], cfg['res_phase']):
            self.declare_readout(ch=ch, length=cfg['res_length'], freq=f, phase=ph, gen_ch=gen_ch)

        self.add_pulse(ch=gen_ch, name="res_pulse", 
                       style="const", 
                       length=cfg["res_length"],
                       mask=[0,1,2,3,4,5],
                      )
        
    def _body(self, cfg):
        self.pulse(ch=cfg['mux_ch'], name="res_pulse", t=0)
        self.trigger(ros=cfg['mux_ro_chs'], pins=[0], t=cfg['trig_time'])

###################
# Run the Program
###################

fpts=[config["start"] + ii*config["step"] for ii in range(config["expts"])]
fcenter = np.array(config['res_freq_ge'])

avgi = np.zeros((len(fcenter), len(fpts)))
avgq = np.zeros((len(fcenter), len(fpts)))
amps = np.zeros((len(fcenter), len(fpts)))
for index, f in enumerate(tqdm(fpts)):
    config["res_freq_ge"] = fcenter + f
    prog = SingleToneSpectroscopyProgram(soccfg, reps=config["reps"], final_delay=0.5, cfg=config)    
    iq_list = prog.acquire(soc, soft_avgs = config["py_avg"], progress=False)
    for i in range(len(fcenter)):
        avgi[i][index] = iq_list[i][:,0]
        avgq[i][index] = iq_list[i][:,1]
        amps[i][index] = np.abs(iq_list[i][:,0]+1j*iq_list[i][:,1])
amps = np.array(amps)
avgi = np.array(avgi)
avgq = np.array(avgq)

#####################################
# ----- Saves data to a file ----- #
#####################################

prefix = str(datetime.date.today())
exp_name = expt_name + '_' + prefix
print('Experiment name: ' + exp_name)

data_path = DATA_PATH

fname = get_next_filename(data_path, exp_name, suffix='.h5')
print('Current data file: ' + fname)

with SlabFile(data_path + '\\' + fname, 'a') as f:
    # 'a': read/write/create

    # - Adds data to the file - #
    f.append('fpts', fpts)
    f.append('amps', amps)
    f.append('fcenter', fcenter)
    f.append('avgi', avgi)
    f.append('avgq', avgq)

    del config['res_freq_ge'] # dont need to save this...
    # formats config into file as a single line
    f.attrs['config'] = json.dumps(config)

data = data_path + '\\' + fname