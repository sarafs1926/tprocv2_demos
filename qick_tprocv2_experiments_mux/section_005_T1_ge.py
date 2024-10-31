
import datetime
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import h5py
# Assuming these are defined elsewhere and importable
from build_task import *
from build_state import *
from expt_config import *
from system_config import *
class T1Program(AveragerProgramV2):
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
                       gain=cfg['pi_amp'],
                       )

        self.add_loop("waitloop", cfg["steps"])

    def _body(self, cfg):
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse", t=0)  # play probe pulse
        self.delay_auto(cfg['wait_time'] + 0.01, tag='wait')  # wait_time after last pulse
        self.pulse(ch=cfg['res_ch'], name="res_pulse", t=0)
        self.trigger(ros=cfg['ro_ch'], pins=[0], t=cfg['trig_time'])


class T1Measurement:
    def __init__(self, QubitIndex, outerFolder, round_num):
        self.QubitIndex = QubitIndex
        self.outerFolder = outerFolder
        self.expt_name = "T1_ge"
        self.Qubit = 'Q' + str(self.QubitIndex)
        self.exp_cfg = expt_cfg[self.expt_name]
        self.q_config = all_qubit_state(system_config)
        self.round_num = round_num

        self.exp_cfg = add_qubit_experiment(expt_cfg, self.expt_name, self.QubitIndex)
        self.config_orig = {**self.q_config[self.Qubit], **self.exp_cfg}

        import copy
        self.config = copy.deepcopy(self.config_orig)

        print(f'Q {self.QubitIndex + 1} Round {round_num} T1 configuration: ', self.config)

        self.q1_t1 = []
        self.q1_t1_err = []
        self.dates = []

    def run(self, soccfg, soc):
        # defaults to 5, just make it to only look at this qubit
        res_gains = self.set_res_gain_ge(self.QubitIndex)
        self.config.update([('res_gain_ge', res_gains)])

        now = datetime.datetime.now()
        self.dates.append(now)

        t1 = T1Program(soccfg, reps=self.exp_cfg['reps'], final_delay=self.exp_cfg['relax_delay'], cfg=self.config)
        iq_list = t1.acquire(soc, soft_avgs=self.exp_cfg['rounds'], progress=True)
        delay_times = t1.get_time_param('wait', "t", as_array=True)

        self.plot_results(iq_list, delay_times, now, self.config, self.QubitIndex)

    def set_res_gain_ge(self, QUBIT_INDEX, num_qubits=6):
        """Sets the gain for the selected qubit to 1, others to 0."""
        res_gain_ge = [0] * num_qubits  # Initialize all gains to 0
        if 0 <= QUBIT_INDEX < num_qubits:  # makes sure you are within the range of options
            res_gain_ge[QUBIT_INDEX] = 1  # Set the gain for the selected qubit
        return res_gain_ge

    def exponential(self, x, a, b, c, d):
        return a * np.exp(-(x - b) / c) + d

    def create_folder_if_not_exists(self, folder_path):
        import os
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

    def plot_results(self, iq_list, delay_times, now, config, QubitIndex):
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        plt.rcParams.update({'font.size': 18})

        I = iq_list[self.QubitIndex][0, :, 0]
        Q = iq_list[self.QubitIndex][0, :, 1]

        signal = Q
        if signal[-1] > signal[0]:
            q1_a_guess = -(np.max(signal) - np.min(signal))
        else:
            q1_a_guess = (np.max(signal) - np.min(signal))

        q1_b_guess = 0
        q1_c_guess = delay_times[-1] / 5
        q1_d_guess = np.min(signal)

        q1_guess = [q1_a_guess, q1_b_guess, q1_c_guess, q1_d_guess]
        q1_popt, q1_pcov = curve_fit(self.exponential, delay_times, signal, maxfev=100000, p0=q1_guess)
        q1_fit_exponential = self.exponential(delay_times, *q1_popt)

        T1_est = q1_popt[2]
        T1_err = np.sqrt(q1_pcov[2][2])

        # I subplot
        ax1.plot(delay_times, I, label="Gain (a.u.)", linewidth=2)
        ax1.set_ylabel("I Amplitude (a.u.)", fontsize=20)
        ax1.tick_params(axis='both', which='major', labelsize=16)
        # ax1.axvline(freq_q, color='orange', linestyle='--', linewidth=2)

        # Q subplot
        ax2.plot(delay_times, Q, label="Q", linewidth=2)
        ax2.plot(delay_times, q1_fit_exponential, '-', color='red', linewidth=3, label="Fit")
        ax2.set_xlabel("Delay time (us)", fontsize=20)
        ax2.set_ylabel("Q Amplitude (a.u.)", fontsize=20)
        ax2.tick_params(axis='both', which='major', labelsize=16)
        # ax2.axvline(freq_q, color='orange', linestyle='--', linewidth=2)

        plt.tight_layout()

        plot_middle = (ax1.get_position().x0 + ax1.get_position().x1) / 2

        fig.text(plot_middle, 0.98,
                 f"T1 Q{QubitIndex + 1}, pi gain %.2f" % config[
                     'pi_amp'] + f", {config['sigma'] * 1000} ns sigma" + f", {config['reps']}*{config['rounds']} avgs," + f" T1 = {T1_est:.3f} ± {T1_err:.3f} µs",
                 fontsize=24, ha='center', va='top')

        # Adjust the top margin to make room for the title
        plt.subplots_adjust(top=0.93)
        outerFolder_expt = outerFolder + "/" + self.expt_name + "/"
        self.create_folder_if_not_exists(outerFolder_expt)
        formatted_datetime = now.strftime("%Y-%m-%d_%H-%M-%S")
        file_name = outerFolder_expt + f"R_{self.round_num}" + f"Q_{self.QubitIndex+1}" + f"{formatted_datetime}_" + self.expt_name + f"_q{QubitIndex + 1}.png"
        fig.savefig(file_name, dpi=300, bbox_inches='tight')  # , facecolor='white'
        plt.close(fig)

        self.q1_t1.append(T1_est)
        self.q1_t1_err.append(T1_err)


        # Create the HDF5 file and save data
        h5_filename = outerFolder_expt + f"R_{self.round_num}" + f"{formatted_datetime}_" + self.expt_name + f"_q{QubitIndex + 1}.h5"
        with h5py.File(h5_filename, 'w') as f:
            # Save T1 and error
            f.create_dataset("T1_estimates", data=np.array(self.q1_t1))
            f.create_dataset("T1_errors", data=np.array(self.q1_t1_err))

            # Save dates
            f.create_dataset("dates", data=np.array(self.dates, dtype='S'))

            # Save I, Q, and delay times
            f.create_dataset("I", data=I)
            f.create_dataset("Q", data=Q)
            f.create_dataset("delay_times", data=delay_times)

            # Save q1_fit_exponential so we can plot it later
            f.create_dataset("q1_fit_exponential", data=q1_fit_exponential)

            # Save the parameters needed for plotting the text
            f.create_dataset("pi_amp", data=np.array([config['pi_amp']]))
            f.create_dataset("sigma", data=np.array([config['sigma']]))
            f.create_dataset("reps", data=np.array([config['reps']]))
            f.create_dataset("rounds", data=np.array([config['rounds']]))
            f.create_dataset("T1_est", data=np.array([T1_est]))
            f.create_dataset("T1_err", data=np.array([T1_err]))

            # Save plot_middle value
            f.create_dataset("plot_middle", data=np.array([plot_middle]))

            # Save Qubit_num as QubitIndex + 1
            f.create_dataset("Qubit_num", data=np.array([QubitIndex + 1]))
