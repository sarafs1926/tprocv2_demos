import os
os.chdir(os.getcwd() + '/qick_tprocv2_experiments')

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

from tqdm import tqdm

# connect to the rfsoc and print soccfg
from rfsoc_connect import *

# ----- Experiment configurations ----- #
expt_name = "res_spec_ge"
config, expt_name = initialize_configs(expt_name)
print(expt_name + '\n')
print(config)

##################
# Define Program #
##################

class SingleToneSpectroscopyProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        if MUX == 'False': # non-MUX case
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


if MUX == 'False': # non-MUX program
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
            iq_list = prog.acquire(soc, soft_avgs=1, progress=False, remove_offset = False)
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
        f.append('fpts', freqs)
        f.append('amps', amps)
        if IS_VISDOM:
            f.append('avgi', expt_I)
            f.append('avgq', expt_Q)
        else:
            f.append('avgi', iq_list[0][0].T[0])
            f.append('avgq', iq_list[0][0].T[1])

        del config['res_freq_ge'] # cant save QickParam
        # formats config into file as a single line
        f.attrs['config'] = json.dumps(config)
        # f.attrs['fit_result'] = json.dumps(fit_result)

    data = data_path + '\\' + fname

else: # MUX program
    ###################
    # Run the Program
    ###################

    fpts=[config["start"] + ii*config["step"] for ii in range(config["expts"])]
    fcenter = np.array(config['freqs'])

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