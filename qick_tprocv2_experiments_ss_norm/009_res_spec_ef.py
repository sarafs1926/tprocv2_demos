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


# ----- Experiment configurations ----- #
expt_name = "res_spec_ef"
QubitIndex = QUBIT_INDEX
Qubit = 'Q' + str(QubitIndex)

exp_cfg = add_qubit_experiment(expt_cfg, expt_name, QubitIndex)
q_config = all_qubit_state(system_config)
config = {**q_config['Q' + str(QubitIndex)], **exp_cfg}
print(config)

##################
# Define Program #
##################

class SingleToneSpectroscopyProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        ro_ch = cfg['ro_ch']
        res_ch = cfg['res_ch']
        qubit_ch = cfg['qubit_ch']

        self.declare_gen(ch=res_ch, nqz=cfg['nqz_res'])
        self.declare_gen(ch=qubit_ch, nqz=cfg['nqz_qubit'])

        self.declare_readout(ch=ro_ch, length=cfg['ro_length'])

        self.add_loop("freqloop", cfg["steps"])
        self.add_readoutconfig(ch=ro_ch, name="myro", freq=cfg['res_freq_ef'], gen_ch=res_ch)
        
        self.add_pulse(ch=res_ch, name="res_pulse", ro_ch=ro_ch, 
                       style="const", 
                       length=cfg['res_length'], 
                       freq=cfg['res_freq_ef'], 
                       phase=cfg['res_phase'],
                       gain=cfg['res_gain_ef'],
                      )
        
        self.add_gauss(ch=qubit_ch, name="ramp", sigma=cfg['sigma'], length=cfg['sigma']*5, even_length=True)
        self.add_pulse(ch=qubit_ch, name="qubit_pi_pulse", 
                       style="arb", 
                       envelope="ramp", 
                       freq=cfg['qubit_freq_ge'], 
                       phase=cfg['qubit_phase'],
                       gain=cfg['qubit_gain_ge'], 
                      )
        
    def _body(self, cfg):
        self.send_readoutconfig(ch=cfg['ro_ch'], name="myro", t=0)
        self.pulse(ch=cfg['qubit_ch'], name="qubit_pi_pulse", t=0)
        self.delay_auto(0.01)
        self.pulse(ch=cfg['res_ch'], name="res_pulse", t=0)
        self.trigger(ros=[cfg['ro_ch']], pins=[0], t=cfg['trig_time'])


###################
# Run the Program
###################

prog = SingleToneSpectroscopyProgram(soccfg, reps=config['reps'], final_delay=config['relax_delay'], cfg=config)
py_avg = config['py_avg']

# for live plotting
IS_VISDOM = True
if IS_VISDOM:
    expt_I = expt_Q = expt_mags = expt_phases = expt_pop = None
    viz = visdom.Visdom()
    assert viz.check_connection(timeout_seconds=5), "Visdom server not connected!"
    viz.close(win=None) # close previous plots
    win1 = viz.line( X=np.arange(0, 1), Y=np.arange(0, 1),
                    opts=dict(height=400, width=700, title='Resonator Spectroscopy', showlegend=True, xlabel='expt_pts'))

    for ii in range(py_avg):
        iq_list = prog.acquire(soc, soft_avgs=1, progress=False)
        freqs = prog.get_pulse_param("res_pulse", "freq", as_array=True)

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

        viz.line(X = freqs, Y = expt_mags, win=win1, name='I',
                 opts=dict(height=400, width=700, title='Resonator Spectroscopy', showlegend=True, xlabel='expt_pts'))

    amps = np.abs(expt_I + 1j*expt_Q)
    
else:
    iq_list = prog.acquire(soc, soft_avgs = py_avg, progress=True)
    freqs = prog.get_pulse_param("res_pulse", "freq", as_array=True)
    amps = np.abs(iq_list[0][0].dot([1,1j]))


# print("Fitting data...")
# ### Fit ###
# fit = Fit()
# # Choose the suitable fitting function
# fit_result = fit.reflection_resonator_spectroscopy(freqs, amps)

# fit_result = {
#         "f": fit_result['f'],
#         "kc": fit_result['kc'],
#         "ki": fit_result['ki'],
#         "k": fit_result['k'],
#         "offset": fit_result['offset'],
#         "slope": fit_result['slope']
#     }

# pp.pprint(fit_result)

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
    f.append('fpts', freqs)
    f.append('amps', amps)
    if IS_VISDOM:
        f.append('avgi', expt_I)
        f.append('avgq', expt_Q)
    else:
        f.append('avgi', iq_list[0][0].T[0])
        f.append('avgq', iq_list[0][0].T[1])

    del config['res_freq_ef'] # cant save QickParam
    # formats config into file as a single line
    f.attrs['config'] = json.dumps(config)
    # f.attrs['fit_result'] = json.dumps(fit_result)

data = data_path + '\\' + fname