import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from datetime import datetime
from build_task import *
from build_state import *
from expt_config import *
from system_config import *
import copy


class AmplitudeRabiExperiment:
    def __init__(self, QubitIndex, outerFolder, round_num, qubit_freq, signal, save_figs, experiment):
        self.QubitIndex = QubitIndex
        self.outerFolder = outerFolder
        self.expt_name = "power_rabi_ge"
        self.Qubit = 'Q' + str(self.QubitIndex)
        self.exp_cfg = expt_cfg[self.expt_name]
        self.experiment = experiment
        self.q_config = all_qubit_state(self.experiment)
        self.round_num = round_num
        self.qubit_freq = qubit_freq
        self.signal = signal
        self.save_figs = save_figs

        self.exp_cfg = add_qubit_experiment(expt_cfg, self.expt_name, self.QubitIndex)
        self.config = {**self.q_config[self.Qubit], **self.exp_cfg}

    def run(self, soccfg, soc):
        # defaults to 5, just make it to only look at this qubit
        res_gains = self.set_res_gain_ge(self.QubitIndex)
        self.config.update([('res_gain_ge', res_gains)])

        # now update for qubit frequency
        self.config.update([('qubit_freq_ge', self.qubit_freq)])

        #look at the config before we do the experiment
        print(f'Q {self.QubitIndex + 1} Round {self.round_num} Rabi configuration: ', self.config)

        amp_rabi = AmplitudeRabiProgram(soccfg, reps=self.exp_cfg['reps'], final_delay=self.exp_cfg['relax_delay'], cfg=self.config)
        iq_list = amp_rabi.acquire(soc, soft_avgs=self.exp_cfg["rounds"], progress=True)
        gains = amp_rabi.get_pulse_param('qubit_pulse', "gain", as_array=True)

        self.plot_results( iq_list, gains)
        return

    def set_res_gain_ge(self, QUBIT_INDEX, num_qubits=6):
        """Sets the gain for the selected qubit to 1, others to 0."""
        res_gain_ge = [0] * num_qubits  # Initialize all gains to 0
        if 0 <= QUBIT_INDEX < num_qubits:  # makes sure you are within the range of options
            res_gain_ge[QUBIT_INDEX] = 1  # Set the gain for the selected qubit
        return res_gain_ge

    def cosine(self, x, a, b, c, d):
        return a * np.cos(2. * np.pi * b * x - c * 2 * np.pi) + d

    def create_folder_if_not_exists(self, folder_path):
        import os
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

    def plot_results(self, iq_list, gains):

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        plt.rcParams.update({'font.size': 18})

        I = iq_list[self.QubitIndex][0, :, 0]
        Q = iq_list[self.QubitIndex][0, :, 1]

        if 'Q' in self.signal:
            q1_amp = Q
        else:
            q1_amp = I
        q1_a_guess = (np.max(q1_amp) - np.min(q1_amp)) / 2
        q1_b_guess = 1 / gains[-1]
        q1_c_guess = 0
        q1_d_guess = np.mean(q1_amp)

        q1_guess = [q1_a_guess, q1_b_guess, q1_c_guess, q1_d_guess]
        q1_popt, q1_pcov = curve_fit(self.cosine, gains, q1_amp, maxfev=100000, p0=q1_guess)
        q1_fit_cosine = self.cosine(gains, *q1_popt)

        first_three_avg = np.mean(q1_fit_cosine[:3])
        last_three_avg = np.mean(q1_fit_cosine[-3:])

        if last_three_avg > first_three_avg:
            pi_amp = gains[np.argmax(q1_fit_cosine)]
        else:
            pi_amp = gains[np.argmin(q1_fit_cosine)]

        ax1.plot(gains, I, label="Gain (a.u.)", linewidth=2)
        ax1.set_ylabel("I Amplitude (a.u.)", fontsize=20)
        ax1.tick_params(axis='both', which='major', labelsize=16)

        ax2.plot(gains, Q, label="Q", linewidth=2)
        ax2.set_xlabel("Gain (a.u.)", fontsize=20)
        ax2.set_ylabel("Q Amplitude (a.u.)", fontsize=20)
        ax2.tick_params(axis='both', which='major', labelsize=16)

        if 'Q' in self.signal:
            ax2.plot(gains, q1_fit_cosine, '-', color='red', linewidth=3, label="Fit")
        else:
            ax1.plot(gains, q1_fit_cosine, '-', color='red', linewidth=3, label="Fit")


        plt.tight_layout()

        plot_middle = (ax1.get_position().x0 + ax1.get_position().x1) / 2
        fig.text(plot_middle, 0.98,
                 f"Rabi Q{self.QubitIndex + 1}, pi gain %.2f" % pi_amp + f", {self.config['sigma'] * 1000} ns sigma" + f", {self.config['reps']} avgs",
                 fontsize=24, ha='center', va='top')
        plt.subplots_adjust(top=0.93)
        self.experiment.outerFolder_expt = self.outerFolder + "/" + self.expt_name + "/"
        self.experiment.create_folder_if_not_exists = self.outerFolder + "/" + self.expt_name + "/"
        self.create_folder_if_not_exists(self.experiment.outerFolder_expt)
        now = datetime.datetime.now()
        formatted_datetime = now.strftime("%Y-%m-%d_%H-%M-%S")
        file_name = self.experiment.outerFolder_expt + f"R_{self.round_num}_" + f"Q_{self.QubitIndex+1}_" + f"{formatted_datetime}_" + self.expt_name + f"_q{self.QubitIndex+1}.png"

        if self.save_figs:
            fig.savefig(file_name, dpi=300, bbox_inches='tight')
        plt.close(fig)


class AmplitudeRabiProgram(AveragerProgramV2):
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
                       mask=[0, 1, 2, 3, 4, 5],
                       )

        self.declare_gen(ch=qubit_ch, nqz=cfg['nqz_qubit'])
        self.add_gauss(ch=qubit_ch, name="ramp", sigma=cfg['sigma'], length=cfg['sigma'] * 5, even_length=True)
        self.add_pulse(ch=qubit_ch, name="qubit_pulse",
                       style="arb",
                       envelope="ramp",
                       freq=cfg['qubit_freq_ge'],
                       phase=cfg['qubit_phase'],
                       gain=cfg['qubit_gain_ge'],
                       )

        self.add_loop("gainloop", cfg["steps"])

    def _body(self, cfg):
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse", t=0)
        self.delay_auto(t=0.01, tag='waiting')
        self.pulse(ch=cfg['res_ch'], name="res_pulse", t=0)
        self.trigger(ros=cfg['ro_ch'], pins=[0], t=cfg['trig_time'])
