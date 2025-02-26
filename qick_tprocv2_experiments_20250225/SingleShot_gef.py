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
from qick_setup_funcs import *

from qualang_tools.plot import Fit
import pprint as pp

# connect to the rfsoc and print soccfg
from rfsoc_connect import *


# ----- Experiment configurations ----- #
expt_name = "IQ_plot"
config, expt_name = initialize_configs(expt_name)
print(expt_name + '\n')
print(config)

##################
# Define Program #
##################

# ----- single shot g experiment ----- #
class SingleShotProgram_g(AveragerProgramV2):
    def _initialize(self, cfg):
        # initialize the readout channels and pulses for DAC and ADC
        initialize_ro_chs(self, cfg, suffix='_ge')

        self.add_loop("shotloop", cfg["steps"]) # number of total shots

    def _body(self, cfg):
        # readout using proper configurations for MUX and non-MUX
        readout(self, cfg)


# ----- single shot e experiment ----- #
class SingleShotProgram_e(AveragerProgramV2):
    def _initialize(self, cfg):
        # initialize the readout channels and pulses for DAC and ADC
        initialize_ro_chs(self, cfg, suffix='_ge')

        self.add_loop("shotloop", cfg["steps"]) # number of total shots

        # initialize the qubit generator channel
        declare_gen_ch(self, cfg, cfg['qubit_ch'], usage='qubit', suffix='_ge')
        # initialize qubit pulse
        declare_pulse(self, cfg, cfg['qubit_ch'], usage='qubit', 
                      pulse_style='gauss', pulse_name='qubit_pulse', suffix='_ge')

    def _body(self, cfg):
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse", t=0)  #play pulse
        self.delay_auto(0.01, tag='wait')
        # readout using proper configurations for MUX and non-MUX
        readout(self, cfg)


# ----- single shot f experiment ----- #
class SingleShotProgram_f(AveragerProgramV2):
    def _initialize(self, cfg):
        qubit_ch = cfg['qubit_ch']
        qubit_ch_ef = cfg['qubit_ch_ef']
        
        # initialize the readout channels and pulses for DAC and ADC
        initialize_ro_chs(self, cfg, suffix='_ge')
        # initialize the qubit generator channel for pi pulse
        declare_gen_ch(self, cfg, qubit_ch, usage='qubit', suffix='_ge')
        # initialize the qubit generator channel for ef pulse
        declare_gen_ch(self, cfg, qubit_ch_ef, usage='qubit', suffix='_ef')

        self.add_loop("shotloop", cfg["steps"])
        
        # initialize qubit pi pulse
        declare_pulse(self, cfg, qubit_ch, usage='qubit', 
                      pulse_style='gauss', pulse_name='qubit_pi_pulse', suffix='_ge')

        # initialize qubit ef pulse
        declare_pulse(self, cfg, qubit_ch_ef, usage='qubit', 
                      pulse_style='gauss', pulse_name='qubit_pulse_ef', suffix='_ef')
    
    def _body(self, cfg):
        self.send_readoutconfig(ch=cfg['ro_ch'], name="myro", t=0)
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pi_pulse", t=0)  #play pulse
        self.delay_auto(0.01, tag='wait1')
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse_ef", t=0)  #play pulse
        self.delay_auto(0.01)
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pi_pulse", t=0)  #play pulse
        self.delay_auto(0.01)
        # readout using proper configurations for MUX and non-MUX
        readout(self, cfg)

###################
# Run the Program
###################
if config['SS_ONLY'] == 'True':
    ssp_g = SingleShotProgram_g(soccfg, reps=1, final_delay=config['relax_delay'], cfg=config)
    iq_list_g = ssp_g.acquire(soc, soft_avgs=1, progress=True)

    ssp_e = SingleShotProgram_e(soccfg, reps=1, final_delay=config['relax_delay'], cfg=config)
    iq_list_e = ssp_e.acquire(soc, soft_avgs=1, progress=True)

    ssp_f = SingleShotProgram_f(soccfg, reps=1, final_delay=config['relax_delay'], cfg=config)
    iq_list_f = ssp_f.acquire(soc, soft_avgs=1, progress=True)

    I_g = iq_list_g[0][0].T[0]
    Q_g = iq_list_g[0][0].T[1]
    I_e = iq_list_e[0][0].T[0]
    Q_e = iq_list_e[0][0].T[1]
    I_f = iq_list_f[0][0].T[0]
    Q_f = iq_list_f[0][0].T[1]

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
        f.append('I_g', I_g)
        f.append('Q_g', Q_g)
        f.append('I_e', I_e)
        f.append('Q_e', Q_e)    
        f.append('I_f', I_f)
        f.append('Q_f', Q_f)  

        # formats config into file as a single line
        f.attrs['config'] = json.dumps(config)

    data = data_path + '\\' + fname