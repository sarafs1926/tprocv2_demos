from build_task import *
from build_state import *
from expt_config import *
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import datetime
import copy
import visdom

from NetDrivers import E36300

class BiasQubitSpectroscopy:
    def __init__(self, QubitIndex, outerFolder, round_num, signal, save_figs, experiment=None, live_plot=None):
        self.QubitIndex = QubitIndex
        self.outerFolder = outerFolder
        self.expt_name = "bias_qubit_spec_ge"
        self.signal = signal
        self.save_figs = save_figs
        self.experiment = experiment
        self.Qubit = 'Q' + str(self.QubitIndex)
        self.exp_cfg = expt_cfg[self.expt_name]
        self.round_num = round_num
        if experiment is not None:
            self.q_config = all_qubit_state(self.experiment)
            self.live_plot = live_plot
            self.exp_cfg = add_qubit_experiment(expt_cfg, self.expt_name, self.QubitIndex)
            self.config = {**self.q_config[self.Qubit], **self.exp_cfg}

            print(f'Q {self.QubitIndex + 1} Round {self.round_num} Qubit Spec configuration: ')

        ## Some of the above can probably be gotten rid of since we don't necessarily want to do this during RR

    def sweep_bias(self, soccfg, soc, start_volt, stop_volt, volt_pts):
        ## Is there a better way to define and read in voltage sweep params?
        vsweep = np.linspace(start_volt, stop_volt, volt_pts, endpoint = True)

        Bias_PS_ip = ['192.168.0.44', '192.168.0.44', '192.168.0.44', '192.168.0.41'] #IP address of bias PS (qubits 1-3 are the same PS)
        Bias_ch = [1, 2, 3, 1] #Channel number of qubit 1-4 on associated PS
        qubit_index = int(self.QubitIndex)
        BiasPS = E36300(Bias_PS_ip[qubit_index], server_port = 5025)

        BiasPS.setVoltage(Bias_ch[qubit_index], 0)
        BiasPS.setOutputState(Bias_ch[qubit_index], enable=True) #Why is this weird

        I_arr = []
        Q_arr = []
        freq_arr = []

        for index, v in enumerate(tqdm(vsweep)):
            BiasPS.setVoltage(Bias_ch[qubit_index], v)

            #add wait time? Probably

            qspec = PulseProbeSpectroscopyProgram(soccfg, reps=self.config['reps'], final_delay = self.exp_cfg['relax_delay'], cfg=self.config)
            iq_list = qspec.acquire(soc, soft_avgs = self.exp_cfg["rounds"], progress=True)
            I = iq_list[self.QubitIndex][0, :, 0]
            Q = iq_list[self.QubitIndex][0, :, 1]
            freq = qspec.get_pulse_param('qubit_pulse', "freq", as_array=True)
            I_arr.append(I)
            Q_arr.append(Q)
            freq_arr.append(freq)
        BiasPS.setOutputState(Bias_ch[qubit_index], enable-False)
        return I_arr, Q_arr, freq_arr




class PulseProbeSpectroscopyProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        ro_ch = cfg['ro_ch']
        res_ch = cfg['res_ch']
        qubit_ch = cfg['qubit_ch']

        self.declare_gen(ch=res_ch, nqz=cfg['nqz_res'], ro_ch=ro_ch[0],
                         mux_freqs=cfg['res_freq_ge'],
                         mux_gains=cfg['res_gain_ge'],
                         mux_phases=cfg['res_phase'],
                         mixer_freq=cfg['mixer_freq'])
        for ch, f, ph in zip(cfg['ro_ch'], cfg['res_freq_ge'], cfg['ro_phase']):
            self.declare_readout(ch=ch, length=cfg['res_length'], freq=f, phase=ph, gen_ch=res_ch)

        self.add_pulse(ch=res_ch, name="res_pulse",
                       style="const",
                       length=cfg["res_length"],
                       mask=[0, 1, 2, 3],
                       )

        self.declare_gen(ch=qubit_ch, nqz=cfg['nqz_qubit'], mixer_freq=cfg['qubit_mixer_freq'])
        self.add_pulse(ch=qubit_ch, name="qubit_pulse", ro_ch=ro_ch[0],
                       style="const",
                       length=cfg['qubit_length_ge'],
                       freq=cfg['qubit_freq_ge'],
                       phase=0,
                       gain=cfg['qubit_gain_ge'],
                       )

        self.add_loop("freqloop", cfg["steps"])

    def _body(self, cfg):
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse", t=0)  # play probe pulse
        self.delay_auto(t=0.01, tag='waiting')  # Wait til qubit pulse is done before proceeding
        self.pulse(ch=cfg['res_ch'], name="res_pulse", t=0)
        self.trigger(ros=cfg['ro_ch'], pins=[0], t=cfg['trig_time'])