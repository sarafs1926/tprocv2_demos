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

# import single shot information for g-e calibration
from SingleShot import SingleShotProgram_g, SingleShotProgram_e
from SingleShot import config as ss_config

# ----- Experiment configurations ----- #
expt_name = "qubit_spec_ge"
QubitIndex = QUBIT_INDEX
Qubit = 'Q' + str(QubitIndex)

exp_cfg = add_qubit_experiment(expt_cfg, expt_name, QubitIndex)
q_config = all_qubit_state(system_config)
config = {**q_config['Q' + str(QubitIndex)], **exp_cfg}
print(config)

##################
# Define Program #
##################

class PulseProbeSpectroscopyProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        ro_ch = cfg['ro_ch']
        res_ch = cfg['res_ch']
        qubit_ch = cfg['qubit_ch']
        
        self.declare_gen(ch=res_ch, nqz=cfg['nqz_res'])
        self.declare_gen(ch=qubit_ch, nqz=cfg['nqz_qubit'])
        self.declare_readout(ch=ro_ch, length=cfg['ro_length'])

        self.add_loop("freqloop", cfg["steps"])
        self.add_readoutconfig(ch=ro_ch, name="myro", freq=cfg['res_freq_ge'], gen_ch=res_ch)
        
        self.add_pulse(ch=res_ch, name="res_pulse", ro_ch=ro_ch, 
                       style="const", 
                       length=cfg['res_length'], 
                       freq=cfg['res_freq_ge'], 
                       phase=cfg['res_phase'],
                       gain=cfg['res_gain_ge'],
                      )
        
        self.add_pulse(ch=qubit_ch, name="qubit_pulse", ro_ch=ro_ch, 
                       style="const", 
                       length=cfg['qubit_length_ge'], 
                       freq=cfg['qubit_freq_ge'], 
                       phase= 0,
                       gain=cfg['qubit_gain_ge'],
                      )
    
    def _body(self, cfg):
        self.send_readoutconfig(ch=cfg['ro_ch'], name="myro", t=0)
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse", t=0)  #play probe pulse

        self.pulse(ch=cfg['res_ch'], name="res_pulse", t=0)
        self.trigger(ros=[cfg['ro_ch']], pins=[0], t=cfg['trig_time'])

###################
# Run the Program
###################

qspec=PulseProbeSpectroscopyProgram(soccfg, reps=1000, final_delay= 0.5, cfg=config)
py_avg = config['py_avg']

# for live plotting
IS_VISDOM = True
if IS_VISDOM:
    expt_I = expt_Q = expt_mags = expt_phases = expt_pop = None
    viz = visdom.Visdom()
    assert viz.check_connection(timeout_seconds=5), "Visdom server not connected!"
    viz.close(win=None) # close previous plots
    win1 = viz.line( X=np.arange(0, 1), Y=np.arange(0, 1),
        opts=dict(height=400, width=700, title='Qubit Spectroscopy', showlegend=True, xlabel='expt_pts'))

    for ii in range(py_avg):
        iq_list = qspec.acquire(soc, soft_avgs=1, progress=False)
        freqs = qspec.get_pulse_param('qubit_pulse', "freq", as_array=True)

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

        viz.line(X = freqs, Y = expt_mags, win=win1, name='I')

    amps = np.abs(expt_I + 1j*expt_Q)
    
else:
    iq_list = qspec.acquire(soc, soft_avgs=py_avg, progress=True)
    freqs = qspec.get_pulse_param('qubit_pulse', "freq", as_array=True)
    amps = np.abs(iq_list[0][0].dot([1,1j]))

### Fit ###
fit = Fit()
# Choose the suitable fitting function
fit_result = fit.transmission_resonator_spectroscopy(freqs, amps)

# data = fit_result

fit_result = {
        "f": fit_result['f'],
        "kc": fit_result['kc'],
        "ki": fit_result['ki'],
        "k": fit_result['k'],
        "offset": fit_result['offset']
    }

pp.pprint(fit_result)

print('performing single shot for g-e calibration')

ssp_g = SingleShotProgram_g(soccfg, reps=1, final_delay=ss_config['relax_delay'], cfg=ss_config)
iq_list_g = ssp_g.acquire(soc, soft_avgs=1, progress=True)

ssp_e = SingleShotProgram_e(soccfg, reps=1, final_delay=ss_config['relax_delay'], cfg=ss_config)
iq_list_e = ssp_e.acquire(soc, soft_avgs=1, progress=True)

I_g = iq_list_g[0][0].T[0]
Q_g = iq_list_g[0][0].T[1]
I_e = iq_list_e[0][0].T[0]
Q_e = iq_list_e[0][0].T[1]

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

    del config['qubit_freq_ge']  # cant save QickParam
    # formats config into file as a single line
    f.attrs['config'] = json.dumps(config) # cant save configs yet with QickParams
    f.attrs['fit_result'] = json.dumps(fit_result)

    # - Adds ss data to the file - #
    f.append('I_g', I_g)
    f.append('Q_g', Q_g)
    f.append('I_e', I_e)
    f.append('Q_e', Q_e)    

    # formats config into file as a single line
    f.attrs['ss_config'] = json.dumps(ss_config)

data = data_path + '\\' + fname