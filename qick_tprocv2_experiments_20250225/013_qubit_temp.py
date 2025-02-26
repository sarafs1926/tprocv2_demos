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
expt_name = "qubit_temp"
config, expt_name = initialize_configs(expt_name)
print(expt_name + '\n')
print(config)

##################
# Define Program #
##################

class QubitTemperatureProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        qubit_ch = cfg['qubit_ch']
        qubit_ch_ef = cfg['qubit_ch_ef']
        
        # initialize the readout channels and pulses for DAC and ADC
        initialize_ro_chs(self, cfg, suffix='_ge')
        # initialize the qubit generator channel for pi pulse
        declare_gen_ch(self, cfg, qubit_ch, usage='qubit', suffix='_ge')
        # initialize the qubit generator channel for ef pulse
        declare_gen_ch(self, cfg, qubit_ch_ef, usage='qubit', suffix='_ef')
        
        # initialize qubit pi pulse
        declare_pulse(self, cfg, qubit_ch, usage='qubit', 
                      pulse_style='gauss', pulse_name='qubit_pi_pulse', suffix='_ge')

        # initialize qubit ef pulse
        declare_pulse(self, cfg, qubit_ch_ef, usage='qubit', 
                      pulse_style='const', pulse_name='qubit_pulse_ef', suffix='_ef')
        
    def _body(self, cfg):
        self.send_readoutconfig(ch=cfg['ro_ch'], name="myro", t=0)
        # self.pulse(ch=cfg['qubit_ch'], name="qubit_pi_pulse", t=0)
        # self.delay_auto(0.01)
        self.pulse(ch=self.cfg["qubit_ch_ef"], name="qubit_pulse_ef", t=0)  #play probe pulse
        self.delay_auto(0.02)
        self.pulse(ch=cfg['qubit_ch'], name="qubit_pi_pulse", t=0)
        self.delay_auto(0.01)
        # readout using proper configurations for MUX and non-MUX
        readout(self, cfg)


class QubitTemperatureRefProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        qubit_ch = cfg['qubit_ch']
        qubit_ch_ef = cfg['qubit_ch_ef']
        
        # initialize the readout channels and pulses for DAC and ADC
        initialize_ro_chs(self, cfg, suffix='_ef')
        # initialize the qubit generator channel for pi pulse
        declare_gen_ch(self, cfg, qubit_ch, usage='qubit', suffix='_ge')
        # initialize the qubit generator channel for ef pulse
        declare_gen_ch(self, cfg, qubit_ch_ef, usage='qubit', suffix='_ef')
        
        # initialize qubit pi pulse
        declare_pulse(self, cfg, qubit_ch, usage='qubit', 
                      pulse_style='gauss', pulse_name='qubit_pi_pulse', suffix='_ge')

        # initialize qubit ef pulse
        declare_pulse(self, cfg, qubit_ch_ef, usage='qubit', 
                      pulse_style='gauss', pulse_name='qubit_pulse_ef', suffix='_ef')

    def _body(self, cfg):
        self.send_readoutconfig(ch=cfg['ro_ch'], name="myro", t=0)
        self.pulse(ch=cfg['qubit_ch'], name="qubit_pi_pulse", t=0)
        self.delay_auto(0.01)
        self.pulse(ch=self.cfg["qubit_ch_ef"], name="qubit_pulse_ef", t=0)  #play probe pulse
        self.delay_auto(0.02)
        self.pulse(ch=cfg['qubit_ch'], name="qubit_pi_pulse", t=0)
        self.delay_auto(0.01)
        # readout using proper configurations for MUX and non-MUX
        readout(self, cfg)

###################
# Run the Program
###################

expt_pts=[config["start"] + ii*config["step"] for ii in range(config["expts"])]

# for live plotting
IS_VISDOM = True
if IS_VISDOM:
    expt_I = expt_Q = expt_mags = expt_phases = expt_pop = None
    expt_I_ref = expt_Q_ref = expt_mags_ref = expt_phases_ref = expt_pop_ref = None
    viz = visdom.Visdom()
    assert viz.check_connection(timeout_seconds=5), "Visdom server not connected!"
    viz.close(win=None) # close previous plots
    win1 = viz.line( X=np.arange(0, 1), Y=np.arange(0, 1),
        opts=dict(height=400, width=700, title='Length Rabi', showlegend=True, xlabel='expt_pts'))

    
py_avg = config['py_avg']
results_tot=[]
results_tot_ref=[]
for ii in range(py_avg):
    results=[]
    pulse_lengths=[]
    results_ref=[]
    for length in expt_pts:
        config["qubit_length_ef"]=length

        rabi=QubitTemperatureProgram(soccfg,
                            reps=config['reps'],
                            final_delay=config['relax_delay'],
                            cfg=config)
        data = rabi.acquire(soc, soft_avgs = 1, progress=False)
        
        results.append(data[0][0])
        pulse_lengths.append(rabi.get_pulse_param('qubit_pulse_ef', 'length', as_array=True))

        rabi_ref=QubitTemperatureRefProgram(soccfg,
                            reps=config['reps'],
                            final_delay=config['relax_delay'],
                            cfg=config)
        data_ref = rabi_ref.acquire(soc, soft_avgs = 1, progress=False)
        
        results_ref.append(data_ref[0][0])
        # wait_times.append(rabi.get_time_param('waiting', 't', as_array=True))

    iq_list = np.array(results).T
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

    iq_list_ref = np.array(results_ref).T
    # what is the correct shape/index?
    this_I_ref = (iq_list_ref[0])
    this_Q_ref = (iq_list_ref[1])

    if expt_I_ref is None: # ii == 0
        expt_I_ref, expt_Q_ref = this_I_ref, this_Q_ref
    else:
        expt_I_ref = (expt_I_ref * ii + this_I_ref) / (ii + 1.0)
        expt_Q_ref = (expt_Q_ref * ii + this_Q_ref) / (ii + 1.0)

    expt_mags_ref = np.abs(expt_I_ref + 1j * expt_Q_ref)  # magnitude
    expt_phases_ref = np.angle(expt_I_ref + 1j * expt_Q_ref)  # phase

    if IS_VISDOM:
        viz.line(X = pulse_lengths, Y = expt_mags_ref, win=win1, name='ref')
        viz.line(X = pulse_lengths, Y = expt_mags, win=win1, name='I', update='append')
amps = expt_mags
amps_ref = expt_mags_ref

# # ### Fit ###
# fit = Fit()
# # Choose the suitable fitting function
# fit_result = fit.rabi(expt_pts, amps)

# fit_result = {
#         "f": fit_result['f'],
#         "phase": fit_result['phase'],
#         "T": fit_result['T'],
#         "amp": fit_result['amp'],
#         "offset": fit_result['offset']
#     }

# fit_result_ref = fit.rabi(expt_pts, amps_ref)

# fit_result_ref = {
#         "f": fit_result_ref['f'],
#         "phase": fit_result_ref['phase'],
#         "T": fit_result_ref['T'],
#         "amp": fit_result_ref['amp'],
#         "offset": fit_result_ref['offset']
#     }

# pp.pprint(fit_result_ref)


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
    f.append('lengths', expt_pts)
    f.append('amps', amps)
    f.append('amps_ref', amps_ref)
    if IS_VISDOM:
        f.append('avgi', expt_I)
        f.append('avgq', expt_Q)
        f.append('avgi_ref', expt_I_ref)
        f.append('avgq_ref', expt_Q_ref)
    else:
        f.append('avgi', iq_list[0][0].T[0])
        f.append('avgq', iq_list[0][0].T[1])
        f.append('avgi_ref', iq_list_ref[0][0].T[0])
        f.append('avgq_ref', iq_list_ref[0][0].T[1])

    # formats config into file as a single line
    f.attrs['config'] = json.dumps(config)
    # f.attrs['fit_result'] = json.dumps(fit_result)
    # f.attrs['fit_result_ref'] = json.dumps(fit_result_ref)

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


