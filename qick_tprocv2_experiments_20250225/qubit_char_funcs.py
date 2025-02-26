# ----- Imports ----- #
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

# from expt_config import *
# from system_config import *

# open configuration file
with open('system_config.json', 'r') as f:
    configs = json.load(f)

# ---- time of flight experiemnt ---- #
class LoopbackProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        # initialize the readout channels and pulses for DAC and ADC
        initialize_ro_chs(self, cfg, suffix='')
    
    def _body(self, cfg):
        # readout using proper configurations for MUX and non-MUX
        readout(self, cfg)

# ----- resonator spectrosocpy experiment ----- #
class SingleToneSpectroscopyProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        if configs['config']['MUX'] == 'False': # non-MUX case
            ro_ch = cfg['ro_ch']
            res_ch = cfg['res_ch']
            
            # declare the proper readout channels using
            declare_gen_ch(self, cfg, res_ch, usage = 'res', suffix = '_ge')
            self.declare_readout(ch=ro_ch, length=cfg['ro_length'])

            # add a loop for frequency sweep of readout pulse frequency
            self.add_loop("freqloop", cfg["steps"])
            self.add_readoutconfig(ch=ro_ch, name="myro", freq=cfg['res_freq' + '_ge'], gen_ch=res_ch)

            # declare a pulse for the non-MUX channel
            declare_pulse(self, cfg, res_ch, pulse_name='res_pulse', suffix='_ge')
            
        else: # MUX case
            # initialize the readout channels and pulses for DAC and ADC
            initialize_ro_chs(self, cfg, suffix='_ge')

    def _body(self, cfg):
        # readout using proper configurations for MUX and non-MUX
        readout(self, cfg)

# ----- qubit spectrosocpy experiment ----- #
class PulseProbeSpectroscopyProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        # initialize the readout channels and pulses for DAC and ADC
        initialize_ro_chs(self, cfg, suffix='_ge')

        # add a loop for frequency sweep of qubit pulse frequency
        self.add_loop("freqloop", cfg["steps"])

        # initialize the qubit generator channel
        declare_gen_ch(self, cfg, cfg['qubit_ch'], usage='qubit', suffix='_ge')
        # initialize qubit pulse
        declare_pulse(self, cfg, cfg['qubit_ch'], usage='qubit', 
                      pulse_style='const', pulse_name='qubit_pulse', suffix='_ge')

    def _body(self, cfg):
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse", t=0)  #play probe pulse

        # readout using proper configurations for MUX and non-MUX
        readout(self, cfg)

# ----- time rabi experiment ----- #
class LengthRabiProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        # initialize the readout channels and pulses for DAC and ADC
        initialize_ro_chs(self, cfg, suffix='_ge')

        # initialize the qubit generator channel
        declare_gen_ch(self, cfg, cfg['qubit_ch'], usage='qubit', suffix='_ge')
        # initialize qubit pulse
        declare_pulse(self, cfg, cfg['qubit_ch'], usage='qubit', 
                      pulse_style='const', pulse_name='qubit_pulse', suffix='_ge')
        

    def _body(self, cfg):
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse", t=0)  #play probe pulse
        
        self.delay_auto(t=0.02, tag='waiting')
        
        # readout using proper configurations for MUX and non-MUX
        readout(self, cfg)

# ----- power rabi experiment ----- #
class AmplitudeRabiProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        # initialize the readout channels and pulses for DAC and ADC
        initialize_ro_chs(self, cfg, suffix='_ge')

        self.add_loop("gainloop", cfg["steps"])
        
        # initialize the qubit generator channel
        declare_gen_ch(self, cfg, cfg['qubit_ch'], usage='qubit', suffix='_ge')
        # initialize qubit pulse
        declare_pulse(self, cfg, cfg['qubit_ch'], usage='qubit', 
                      pulse_style='gauss', pulse_name='qubit_pulse', suffix='_ge')
    
    def _body(self, cfg):
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse", t=0)  #play probe pulse
        
        self.delay_auto(t=0.01, tag='waiting')
        
        # readout using proper configurations for MUX and non-MUX
        readout(self, cfg)

# ----- ramsey experiment ----- #
class RamseyProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        # initialize the readout channels and pulses for DAC and ADC
        initialize_ro_chs(self, cfg, suffix='_ge')

        self.add_loop("waitloop", cfg["steps"])
        
       # initialize the qubit generator channel
        declare_gen_ch(self, cfg, cfg['qubit_ch'], usage='qubit', suffix='_ge')
        # initialize qubit pulse
        declare_pulse(self, cfg, cfg['qubit_ch'], usage='qubit', 
                      pulse_style='gauss', pulse_name='qubit_pulse1', suffix='_ge')
        
        # initialize qubit pulse
        declare_pulse(self, cfg, cfg['qubit_ch'], usage='qubit', 
                      pulse_style='gauss', pulse_name='qubit_pulse2', suffix='_ge')
    
    def _body(self, cfg):
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse1", t=0)  #play probe pulse
        
        self.delay_auto(cfg['wait_time']+0.01, tag='wait') # wait_time after last pulse
        
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse2", t=0)  #play pulse

        self.delay_auto(0.01) # wait_time after last pulse
        
        # readout using proper configurations for MUX and non-MUX
        readout(self, cfg)


# ----- spin echo experiment ----- #
class SpinEchoProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        # initialize the readout channels and pulses for DAC and ADC
        initialize_ro_chs(self, cfg, suffix='_ge')

        self.add_loop("waitloop", cfg["steps"])
        
       # initialize the qubit generator channel
        declare_gen_ch(self, cfg, cfg['qubit_ch'], usage='qubit', suffix='_ge')

        # initialize qubit pulse 1
        declare_pulse(self, cfg, cfg['qubit_ch'], usage='qubit', 
                      pulse_style='gauss', pulse_name='qubit_pulse1', suffix='_ge')
        
        # initialize qubit pi-pulse
        declare_pulse(self, cfg, cfg['qubit_ch'], usage='qubit', 
                      pulse_style='gauss', pulse_name='qubit_pulse2', suffix='_ge')
        
        # initialize qubit pulse 2
        declare_pulse(self, cfg, cfg['qubit_ch'], usage='qubit', 
                      pulse_style='gauss', pulse_name='qubit_pulse2', suffix='_ge')
    
        
    def _body(self, cfg):
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse1", t=0)  #play probe pulse
        
        self.delay_auto((cfg['wait_time']/2)+0.01, tag = 'wait1') # wait_time after last pulse (wait / 2)
        
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse_pi", t=0)  #play pulse
        
        self.delay_auto((cfg['wait_time']/2)+0.01, tag = 'wait2') # wait_time after last pulse (wait / 2)

        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse2", t=0)  #play pulse

        self.delay_auto(0.01) # wait_time after last pulse
        
        # readout using proper configurations for MUX and non-MUX
        readout(self, cfg)


# ----- t1 experiment ----- #
class T1Program(AveragerProgramV2):
    def _initialize(self, cfg):
        # initialize the readout channels and pulses for DAC and ADC
        initialize_ro_chs(self, cfg, suffix='_ge')

        self.add_loop("waitloop", cfg["steps"])
        
        # initialize the qubit generator channel
        declare_gen_ch(self, cfg, cfg['qubit_ch'], usage='qubit', suffix='_ge')
        
        # initialize qubit pulse
        declare_pulse(self, cfg, cfg['qubit_ch'], usage='qubit', 
                      pulse_style='gauss', pulse_name='qubit_pulse', suffix='_ge')
    
    def _body(self, cfg):
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse", t=0)  #play pulse
        
        self.delay_auto(cfg['wait_time']+0.01, tag='wait') # wait_time after last pulse

        # readout using proper configurations for MUX and non-MUX
        readout(self, cfg)

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
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_ge_pulse", t=0)  #play pulse
        self.delay_auto(0.01, tag='wait1')
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_ef_pulse", t=0)  #play pulse
        self.delay_auto(0.01)
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_ge_pulse", t=0)  #play pulse
        self.delay_auto(0.01)
        # readout using proper configurations for MUX and non-MUX
        readout(self, cfg)