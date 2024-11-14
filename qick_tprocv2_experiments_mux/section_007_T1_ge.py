from prompt_toolkit.key_binding.bindings.named_commands import self_insert
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from build_task import *
from build_state import *
from expt_config import *
from system_config import *
import copy
import visdom

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

        self.declare_gen(ch=qubit_ch, nqz=cfg['nqz_qubit'], mixer_freq=4000)
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
    def __init__(self, QubitIndex, outerFolder, round_num, signal, save_figs, experiment, live_plot):
        self.QubitIndex = QubitIndex
        self.outerFolder = outerFolder
        self.expt_name = "T1_ge"
        self.Qubit = 'Q' + str(self.QubitIndex)
        self.experiment = experiment
        self.exp_cfg = expt_cfg[self.expt_name]
        self.q_config = all_qubit_state(self.experiment)
        self.round_num = round_num
        self.live_plot = live_plot
        self.signal = signal
        self.save_figs = save_figs

        self.exp_cfg = add_qubit_experiment(expt_cfg, self.expt_name, self.QubitIndex)
        self.config = {**self.q_config[self.Qubit], **self.exp_cfg}
        print(f'Q {self.QubitIndex + 1} Round {self.round_num} T1 configuration: ', self.config)

    def run(self, soccfg, soc):
        now = datetime.datetime.now()
        t1 = T1Program(soccfg, reps=self.exp_cfg['reps'], final_delay=self.exp_cfg['relax_delay'], cfg=self.config)

        if self.live_plot:
            I, Q, delay_times = self.live_plotting(t1, soc)
        else:
            iq_list = t1.acquire(soc, soft_avgs=self.exp_cfg['rounds'], progress=True)
            I = iq_list[self.QubitIndex][0, :, 0]
            Q = iq_list[self.QubitIndex][0, :, 1]
            delay_times = t1.get_time_param('wait', "t", as_array=True)

        T1_est, T1_err, I, Q, q1_fit_exponential = self.plot_results(I, Q, delay_times, now, self.config, self.QubitIndex)
        return  T1_est, T1_err, I, Q, delay_times, q1_fit_exponential

    def live_plotting(self, t1, soc):
        I = Q = expt_mags = expt_phases = expt_pop = None
        viz = visdom.Visdom()
        assert viz.check_connection(timeout_seconds=5), "Visdom server not connected!"
        viz.close(win=None)  # close previous plots
        win1 = viz.line(X=np.arange(0, 1), Y=np.arange(0, 1),
                        opts=dict(height=400, width=700, title='T2 Ramsey Experiment', showlegend=True,
                                  xlabel='expt_pts'))

        for ii in range(self.config["rounds"]):
            iq_list = t1.acquire(soc, soft_avgs=1, progress=True)
            delay_times = t1.get_time_param('wait', "t", as_array=True)

            this_I = iq_list[self.QubitIndex][0, :, 0]
            this_Q = iq_list[self.QubitIndex][0, :, 1]

            if I is None:  # ii == 0
                I, Q = this_I, this_Q
            else:
                I = (I * ii + this_I) / (ii + 1.0)
                Q = (Q * ii + this_Q) / (ii + 1.0)

            if 'I' in self.signal:
                signal = I
                plot_sig = 'I'
            elif 'Q' in self.signal:
                signal = Q
                plot_sig = 'Q'
            else:
                if abs(I[-1] - I[0]) > abs(Q[-1] - Q[0]):
                    signal = I
                    plot_sig = 'I'
                else:
                    signal = Q
                    plot_sig = 'Q'

            viz.line(X=delay_times, Y=signal, win=win1, name=plot_sig)
        return I, Q, delay_times

    def exponential(self, x, a, b, c, d):
        return a * np.exp(-(x - b) / c) + d

    def create_folder_if_not_exists(self, folder_path):
        import os
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

    def exponential(self, x, a, b, c, d):
        return a * np.exp(- (x - b) / c) + d

    def plot_results(self, I, Q, delay_times, now, config, QubitIndex):
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        plt.rcParams.update({'font.size': 18})

        if 'I' in self.signal:
            signal = I
            plot_sig='I'
        elif 'Q' in self.signal:
            signal = Q
            plot_sig = 'Q'
        else:
            if abs(I[-1]-I[0]) > abs(Q[-1]-Q[0]):
                signal = I
                plot_sig = 'I'
            else:
                signal = Q
                plot_sig = 'Q'

        # Initial guess for parameters
        q1_a_guess = np.max(signal) - np.min(signal)  # Initial guess for amplitude (a)
        q1_b_guess = 0  # Initial guess for time shift (b)
        q1_c_guess = (delay_times[-1] - delay_times[0]) / 5  # Initial guess for decay constant (T1)
        q1_d_guess = np.min(signal)  # Initial guess for baseline (d)

        # Form the guess array
        q1_guess = [q1_a_guess, q1_b_guess, q1_c_guess, q1_d_guess]

        # Define bounds to constrain T1 (c) to be positive, but allow amplitude (a) to be negative
        lower_bounds = [-np.inf, -np.inf, 0, -np.inf]  # Amplitude (a) can be negative/positive, but T1 (c) > 0
        upper_bounds = [np.inf, np.inf, np.inf, np.inf]  # No upper bound on parameters

        # Perform the fit using the 'trf' method with bounds
        q1_popt, q1_pcov = curve_fit(self.exponential, delay_times, signal,
                                     p0=q1_guess, bounds=(lower_bounds, upper_bounds),
                                     method='trf', maxfev=10000)

        # Generate the fitted exponential curve
        q1_fit_exponential = self.exponential(delay_times, *q1_popt)

        # Extract T1 and its error
        T1_est = q1_popt[2]  # Decay constant T1
        T1_err = np.sqrt(q1_pcov[2][2]) if q1_pcov[2][2] >= 0 else float('inf')  # Ensure error is valid

        # I subplot
        ax1.plot(delay_times, I, label="Gain (a.u.)", linewidth=2)
        ax1.set_ylabel("I Amplitude (a.u.)", fontsize=20)
        ax1.tick_params(axis='both', which='major', labelsize=16)
        # ax1.axvline(freq_q, color='orange', linestyle='--', linewidth=2)

        # Q subplot
        ax2.plot(delay_times, Q, label="Q", linewidth=2)
        ax2.set_xlabel("Delay time (us)", fontsize=20)
        ax2.set_ylabel("Q Amplitude (a.u.)", fontsize=20)
        ax2.tick_params(axis='both', which='major', labelsize=16)
        # ax2.axvline(freq_q, color='orange', linestyle='--', linewidth=2)

        if 'I' in plot_sig:
            ax1.plot(delay_times, q1_fit_exponential, '-', color='red', linewidth=3, label="Fit")
        else:
            ax2.plot(delay_times, q1_fit_exponential, '-', color='red', linewidth=3, label="Fit")

        # Adjust spacing
        plt.tight_layout()

        # Calculate the middle of the plot area
        plot_middle = (ax1.get_position().x0 + ax1.get_position().x1) / 2

        # Add title, centered on the plot area
        fig.text(plot_middle, 0.98,
                 f"T1 Q{QubitIndex + 1}, pi gain %.2f" % config[
                     'pi_amp'] + f", {config['sigma'] * 1000} ns sigma" + f", {config['reps']}*{config['rounds']} avgs," + f" T1 = {T1_est:.3f} ± {T1_err:.3f} µs",
                 fontsize=24, ha='center', va='top')

        # Adjust the top margin to make room for the title
        plt.subplots_adjust(top=0.93)
        self.experiment.outerFolder_expt = self.outerFolder + "/" + self.expt_name + "/"
        self.create_folder_if_not_exists(self.experiment.outerFolder_expt)
        formatted_datetime = now.strftime("%Y-%m-%d_%H-%M-%S")
        file_name = self.experiment.outerFolder_expt + f"R_{self.round_num}_" + f"Q_{self.QubitIndex+1}_" + f"{formatted_datetime}_" + self.expt_name + f"_q{QubitIndex + 1}.png"
        if self.save_figs:
            fig.savefig(file_name, dpi=300, bbox_inches='tight')  # , facecolor='white'
        plt.close(fig)
        return T1_est, T1_err, I, Q, q1_fit_exponential

