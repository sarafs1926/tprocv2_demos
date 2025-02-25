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

if SS == 'True':
    # import single shot information for g-e calibration
    from SingleShot import SingleShotProgram_g, SingleShotProgram_e, SingleShotProgram_f
    from SingleShot import config as ss_config

# ----- Experiment configurations ----- #
expt_name = "Ramsey_ef"
config, expt_name = initialize_configs(expt_name)
print(expt_name + '\n')
print(config)

##################
# Define Program #
##################

class RamseyProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        qubit_ch = cfg['qubit_ch']
        qubit_ch_ef = cfg['qubit_ch_ef']
        
        # initialize the readout channels and pulses for DAC and ADC
        initialize_ro_chs(self, cfg, suffix='_ef')
        # initialize the qubit generator channel for pi pulse
        declare_gen_ch(self, cfg, qubit_ch, usage='qubit', suffix='_ge')
        # initialize the qubit generator channel for ef pulse
        declare_gen_ch(self, cfg, qubit_ch_ef, usage='qubit', suffix='_ef')

        self.add_loop("waitloop", cfg["steps"])
        
        # initialize qubit pi pulse
        declare_pulse(self, cfg, qubit_ch, usage='qubit', 
                    pulse_style='gauss', pulse_name='qubit_pi_pulse', suffix='_ge')

        # initialize qubit ef pulse 1
        self.add_gauss(ch=qubit_ch_ef, name="ramp1", sigma=cfg['sigma_ef'], length=cfg['sigma_ef']*5, even_length=True)
        self.add_pulse(ch=qubit_ch_ef, name="qubit_pulse1",
                    style="arb",
                    envelope="ramp1", 
                    freq=cfg['qubit_freq_ef'], 
                    phase=cfg['qubit_phase'],
                    gain=cfg['qubit_gain_ef'] / 2, 
                    )
        
        # initialize qubit ef pulse 2
        self.add_pulse(ch=qubit_ch_ef, name="qubit_pulse2",
                    style="arb", 
                    envelope="ramp1", 
                    freq=cfg['qubit_freq_ef'], 
                    phase=cfg['qubit_phase'] + cfg['wait_time']*360*cfg['ramsey_freq'], # current phase + time * 2pi * ramsey freq
                    gain=cfg['qubit_gain_ef'] / 2, 
                    )
    
    def _body(self, cfg):
        if config['ge_meas'] == False:
            self.send_readoutconfig(ch=cfg['ro_ch'], name="myro", t=0)

            self.pulse(ch=cfg['qubit_ch'], name="qubit_pi_pulse", t=0)
            self.delay_auto(0.01)

            self.pulse(ch=self.cfg["qubit_ch_ef"], name="qubit_pulse1", t=0)  #play pulse
            
            self.delay_auto(cfg['wait_time']+0.01, tag='wait') # wait_time after last pulse
            
            self.pulse(ch=self.cfg["qubit_ch_ef"], name="qubit_pulse2", t=0)  #play pulse

            self.delay_auto(0.01) # wait_time after last pulse
            
            # readout using proper configurations for MUX and non-MUX
            readout(self, cfg)
        else:
            self.send_readoutconfig(ch=cfg['ro_ch'], name="myro", t=0)

            self.pulse(ch=cfg['qubit_ch'], name="qubit_pi_pulse", t=0)
            self.delay_auto(0.01)

            self.pulse(ch=self.cfg["qubit_ch_ef"], name="qubit_pulse1", t=0)  #play pulse
            
            self.delay_auto(cfg['wait_time']+0.01, tag='wait') # wait_time after last pulse
            
            self.pulse(ch=self.cfg["qubit_ch_ef"], name="qubit_pulse2", t=0)  #play pulse

            self.delay_auto(0.01) # wait_time after last pulse
            
            self.pulse(ch=cfg['qubit_ch'], name="qubit_pi_pulse", t=0)
            self.delay_auto(0.01)

            # readout using proper configurations for MUX and non-MUX
            readout(self, cfg)

###################
# Run the Program
###################

ramsey=RamseyProgram(soccfg, reps=config['reps'], final_delay=config['relax_delay'], cfg=config)
py_avg = config['py_avg']

# for live plotting
IS_VISDOM = True
if IS_VISDOM:
    expt_I = expt_Q = expt_mags = expt_phases = expt_pop = None
    viz = visdom.Visdom()
    assert viz.check_connection(timeout_seconds=5), "Visdom server not connected!"
    viz.close(win=None) # close previous plots
    win1 = viz.line( X=np.arange(0, 1), Y=np.arange(0, 1),
        opts=dict(height=400, width=700, title='T2 Ramsey Experiment', showlegend=True, xlabel='expt_pts'))

    for ii in range(py_avg):
        iq_list = ramsey.acquire(soc, soft_avgs=1, progress=False)
        delay_times = ramsey.get_time_param('wait', "t", as_array=True)

        iq_list = iq_list[0][0].T
        # what is the correct shape/index?
        this_I = (iq_list[0])
        this_Q = (iq_list[1])

        if expt_I is None: # ii == 0
            expt_I, expt_Q = this_I, this_Q
        else:
            expt_I = (expt_I * ii + this_I) / (ii + 1.0)
            expt_Q = (expt_Q * ii + this_Q) / (ii + 1.0)

        expt_mags = np.abs(expt_I + 1j * expt_Q)  # magnitude
        expt_phases = np.angle(expt_I + 1j * expt_Q)  # phase

        viz.line(X = delay_times, Y = expt_mags, win=win1, name='I',
        opts=dict(height=400, width=700, title='T2 Ramsey Experiment', showlegend=True, xlabel='expt_pts'))


    amps = expt_mags
    
else:
    iq_list = ramsey.acquire(soc, soft_avgs=py_avg, progress=True)
    delay_times = ramsey.get_time_param('wait', "t", as_array=True)
    amps = np.abs(iq_list[0][0].dot([1,1j]))

# ### Fit ###
# fit = Fit()

# # Choose the suitable fitting function
# fit_result = fit.ramsey(delay_times, amps)

# fit_result = {
#     "f": fit_result['f'],
#     "phase": fit_result['phase'],
#     "T2": fit_result['T2'],
#     "amp": fit_result['amp'],
#     "initial_offset": fit_result['initial_offset'],
#     "final_offset": fit_result['final_offset']
#     }

# pp.pprint(fit_result)

if SS == 'True':
    print('performing single shot for g-e-f calibration')

    ssp_g = SingleShotProgram_g(soccfg, reps=1, final_delay=ss_config['relax_delay'], cfg=ss_config)
    iq_list_g = ssp_g.acquire(soc, soft_avgs=1, progress=True)

    ssp_e = SingleShotProgram_e(soccfg, reps=1, final_delay=ss_config['relax_delay'], cfg=ss_config)
    iq_list_e = ssp_e.acquire(soc, soft_avgs=1, progress=True)

    ssp_f = SingleShotProgram_f(soccfg, reps=1, final_delay=ss_config['relax_delay'], cfg=ss_config)
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
    f.append('delay_times', delay_times)
    f.append('amps', amps)
    if IS_VISDOM:
        f.append('avgi', expt_I)
        f.append('avgq', expt_Q)
    else:
        f.append('avgi', iq_list[0][0].T[0])
        f.append('avgq', iq_list[0][0].T[1])

    del config['wait_time']  # cant save QickParam
    # formats config into file as a single line
    f.attrs['config'] = json.dumps(config)
    # f.attrs['fit_result'] = json.dumps(fit_result)

    if SS == 'True':
        # - Adds ss data to the file - #
        f.append('I_g', I_g)
        f.append('Q_g', Q_g)
        f.append('I_e', I_e)
        f.append('Q_e', Q_e)    
        f.append('I_f', I_f)
        f.append('Q_f', Q_f)

        # formats config into file as a single line
        f.attrs['ss_config'] = json.dumps(ss_config)

data = data_path + '\\' + fname
