from prompt_toolkit.key_binding.bindings.named_commands import self_insert
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from build_task import *
from build_state import *
from expt_config import *
from system_config import *
import copy
import visdom
from scipy import optimize
from sklearn import preprocessing
import matplotlib.pyplot as plt
from typing import List, Union
import itertools
import json
import numpy as np
import warnings
from scipy.optimize import OptimizeWarning
class Fit:
    """
    This class takes care of the fitting to the measured data.
    It includes:
        - Fitting to: linear line
                      T1 experiment
                      Ramsey experiment
                      transmission resonator spectroscopy
                      reflection resonator spectroscopy
        - Printing the initial guess and fitting results
        - Plotting the data and the fitting function
        - Saving the data
    """

    # Remove optimize warnings
    warnings.simplefilter("ignore", RuntimeWarning)
    warnings.simplefilter("ignore", OptimizeWarning)

    @staticmethod
    def linear(
        x_data: Union[np.ndarray, List[float]],
        y_data: Union[np.ndarray, List[float]],
        guess=None,
        verbose=False,
        plot=False,
        save=False,
    ) -> dict:
        """
        Create a linear fit of the form

        .. math::
        f(x) = a * x + b

        for unknown parameters :
             a - The slope of the function
             b - The free parameter of the function

         :param x_data: The data on the x-axis
         :param y_data: The data on the y-axis
         :param dict guess: Dictionary containing the initial guess for the fitting parameters (guess=dict(a=20))
         :param verbose: if True prints the initial guess and fitting results
         :param plot: if True plots the data and the fitting function
         :param save: if not False saves the data into a json file
                      The id of the file is save='id'. The name of the json file is `id.json`
         :return: A dictionary of (fit_func, a, b)

        """

        # Normalizing the vectors
        xn = preprocessing.normalize([x_data], return_norm=True)
        yn = preprocessing.normalize([y_data], return_norm=True)
        x = xn[0][0]
        y = yn[0][0]
        x_normal = xn[1][0]
        y_normal = yn[1][0]

        # Finding an initial guess to the slope
        a0 = (y[-1] - y[0]) / (x[-1] - x[0])

        # Finding an initial guess to the free parameter
        b0 = y[0]

        # Check user guess
        if guess is not None:
            for key in guess.keys():
                if key == "a":
                    a0 = float(guess[key]) * x_normal / y_normal
                elif key == "b":
                    b0 = float(guess[key]) / y_normal
                else:
                    raise Exception(
                        f"The key '{key}' specified in 'guess' does not match a fitting parameters for this function."
                    )
        # Print the initial guess if verbose=True
        if verbose:
            print(f"Initial guess:\n" f" a = {a0 * y_normal / x_normal:.3f}, \n" f" b = {b0 * y_normal:.3f}")

        # Fitting function
        def func(x_var, c0, c1):
            return a0 * c0 * x_var + b0 * c1

        def fit_type(x_var, a):
            return func(x_var, a[0], a[1])

        popt, pcov = optimize.curve_fit(func, x, y, p0=[1, 1])
        perr = np.sqrt(np.diag(pcov))

        # Output the fitting function and its parameters
        out = {
            "fit_func": lambda x_var: fit_type(x_var / x_normal, popt) * y_normal,
            "a": [
                popt[0] * a0 * y_normal / x_normal,
                perr[0] * a0 * y_normal / x_normal,
            ],
            "b": [popt[1] * b0 * y_normal, perr[1] * b0 * y_normal],
        }
        # Print the fitting results if verbose=True
        if verbose:
            print(
                f"Fitting results:\n"
                f" a = {out['a'][0]:.3f} +/- {out['a'][1]:.3f}, \n"
                f" b = {out['b'][0]:.3f} +/- {out['b'][1]:.3f}"
            )
        # Plot the data and the fitting function if plot=True
        if plot:
            plt.plot(x_data, fit_type(x, popt) * y_normal)
            plt.plot(
                x_data,
                y_data,
                ".",
                label=f"a  = {out['a'][0]:.1f} +/- {out['a'][1]:.1f} \n b  = {out['b'][0]:.1f} +/- {out['b'][1]:.1f}",
            )
            plt.legend(loc="upper right")
        # Save the data in a json file named 'id.json' if save=id
        if save:
            fit_params = dict(itertools.islice(out.items(), 1, len(out)))
            fit_params["x_data"] = x_data.tolist()
            fit_params["y_data"] = y_data.tolist()
            fit_params["y_fit"] = (func(x, popt[0], popt[1]) * y_normal).tolist()
            json_object = json.dumps(fit_params)
            if save[-5:] == ".json":
                save = save[:-5]
            with open(f"{save}.json", "w") as outfile:
                outfile.write(json_object)

        return out

class T2RProgram(AveragerProgramV2):
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
        self.add_pulse(ch=qubit_ch, name="qubit_pulse1",
                       style="arb",
                       envelope="ramp",
                       freq=cfg['qubit_freq_ge'] ,
                       phase=cfg['qubit_phase'],
                       gain=cfg['pi_amp'] / 2,
                       )

        self.add_pulse(ch=qubit_ch, name="qubit_pulse2",
                       style="arb",
                       envelope="ramp",
                       freq=cfg['qubit_freq_ge'],
                       phase=cfg['qubit_phase'] + cfg['wait_time']*360*cfg['ramsey_freq'], # current phase + time * 2pi * ramsey freq
                       gain=cfg['pi_amp'] / 2,
                      )

        self.add_loop("waitloop", cfg["steps"])

    def _body(self, cfg):
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse1", t=0)  # play probe pulse
        self.delay_auto(cfg['wait_time'] + 0.01, tag='wait')  # wait_time after last pulse
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse2", t=0)  # play probe pulse
        self.delay_auto(0.01)  # wait_time after last pulse
        self.pulse(ch=cfg['res_ch'], name="res_pulse", t=0)
        self.trigger(ros=cfg['ro_ch'], pins=[0], t=cfg['trig_time'])


class T2RMeasurement:
    def __init__(self, QubitIndex, outerFolder, round_num, signal, save_figs, experiment, live_plot):
        self.QubitIndex = QubitIndex
        self.outerFolder = outerFolder
        self.expt_name = "Ramsey_ge"
        self.Qubit = 'Q' + str(self.QubitIndex)
        self.experiment = experiment
        self.exp_cfg = expt_cfg[self.expt_name]
        self.q_config = all_qubit_state(self.experiment)
        self.round_num = round_num
        # self.qubit_freq = qubit_freq
        self.signal = signal
        self.save_figs = save_figs
        self.live_plot = live_plot
        self.exp_cfg = add_qubit_experiment(expt_cfg, self.expt_name, self.QubitIndex)
        self.config = {**self.q_config[self.Qubit], **self.exp_cfg}

    def t2_fit(self, x_data, y_data, verbose = False, guess=None, plot=False):
        #convert us to ns

        # Normalizing the vectors
        xn = preprocessing.normalize([x_data], return_norm=True)
        yn = preprocessing.normalize([y_data], return_norm=True)
        x = xn[0][0]
        y = yn[0][0]
        x_normal = xn[1][0]
        y_normal = yn[1][0]

        # Compute the FFT for guessing the frequency
        fft = np.fft.fft(y)
        f = np.fft.fftfreq(len(x))
        # Take the positive part only
        fft = fft[1: len(f) // 2]
        f = f[1: len(f) // 2]
        # Remove the DC peak if there is one
        if (np.abs(fft)[1:] - np.abs(fft)[:-1] > 0).any():
            first_read_data_ind = np.where(np.abs(fft)[1:] - np.abs(fft)[:-1] > 0)[0][0]  # away from the DC peak
            fft = fft[first_read_data_ind:]
            f = f[first_read_data_ind:]

        # Finding a guess for the frequency
        out_freq = f[np.argmax(np.abs(fft))]
        guess_freq = out_freq / (x[1] - x[0])

        # The period is 1 / guess_freq --> number of oscillations --> peaks decay to get guess_T2
        period = int(np.ceil(1 / out_freq))
        peaks = (
                np.array([np.std(y[i * period: (i + 1) * period]) for i in range(round(len(y) / period))]) * np.sqrt(
            2) * 2
        )

        # Finding a guess for the decay (slope of log(peaks))
        if len(peaks) > 1:
            guess_T2 = -1 / ((np.log(peaks)[-1] - np.log(peaks)[0]) / (period * (len(peaks) - 1))) * (x[1] - x[0])
        else:
            guess_T2 = 100 / x_normal

        # Finding a guess for the offsets
        initial_offset = np.mean(y[:period])
        final_offset = np.mean(y[-period:])

        # Finding a guess for the phase
        guess_phase = np.angle(fft[np.argmax(np.abs(fft))]) - guess_freq * 2 * np.pi * x[0]

        # Check user guess
        if guess is not None:
            for key in guess.keys():
                if key == "f":
                    guess_freq = float(guess[key]) * x_normal
                elif key == "phase":
                    guess_phase = float(guess[key])
                elif key == "T2":
                    guess_T2 = float(guess[key]) * x_normal
                elif key == "amp":
                    peaks[0] = float(guess[key]) / y_normal
                elif key == "initial_offset":
                    initial_offset = float(guess[key]) / y_normal
                elif key == "final_offset":
                    final_offset = float(guess[key]) / y_normal
                else:
                    raise Exception(
                        f"The key '{key}' specified in 'guess' does not match a fitting parameters for this function."
                    )

        # Print the initial guess if verbose=True
        if verbose:
            print(
                f"Initial guess:\n"
                f" f = {guess_freq / x_normal:.3f}, \n"
                f" phase = {guess_phase:.3f}, \n"
                f" T2 = {guess_T2 * x_normal:.3f}, \n"
                f" amp = {peaks[0] * y_normal:.3f}, \n"
                f" initial offset = {initial_offset * y_normal:.3f}, \n"
                f" final_offset = {final_offset * y_normal:.3f}"
            )

        # Fitting function
        def func(x_var, a0, a1, a2, a3, a4, a5):
            return final_offset * a4 * (1 - np.exp(-x_var / (guess_T2 * a1))) + peaks[0] / 2 * a2 * (
                    np.exp(-x_var / (guess_T2 * a1))
                    * (a5 * initial_offset / peaks[0] * 2 + np.cos(2 * np.pi * a0 * guess_freq * x + a3))
            )

        def fit_type(x_var, a):
            return func(x_var, a[0], a[1], a[2], a[3], a[4], a[5])

        popt, pcov = optimize.curve_fit(
            func,
            x,
            y,
            p0=[1, 1, 1, guess_phase, 1, 1],
        )

        perr = np.sqrt(np.diag(pcov))

        # Output the fitting function and its parameters
        out = {
            "fit_func": lambda x_var: fit_type(x_var / x_normal, popt) * y_normal,
            "f": [popt[0] * guess_freq / x_normal, perr[0] * guess_freq / x_normal],
            "phase": [popt[3] % (2 * np.pi), perr[3] % (2 * np.pi)],
            "T2": [(guess_T2 * popt[1]) * x_normal, perr[1] * guess_T2 * x_normal],
            "amp": [peaks[0] * popt[2] * y_normal, perr[2] * peaks[0] * y_normal],
            "initial_offset": [
                popt[5] * initial_offset * y_normal,
                perr[5] * initial_offset * y_normal,
            ],
            "final_offset": [
                final_offset * popt[4] * y_normal,
                perr[4] * final_offset * y_normal,
            ],
        }
        # Print the fitting results if verbose=True
        if verbose:
            print(
                f"Fitting results:\n"
                f" f = {out['f'][0] * 1000:.3f} +/- {out['f'][1] * 1000:.3f} MHz, \n"
                f" phase = {out['phase'][0]:.3f} +/- {out['phase'][1]:.3f} rad, \n"
                f" T2 = {out['T2'][0]:.2f} +/- {out['T2'][1]:.3f} ns, \n"
                f" amp = {out['amp'][0]:.2f} +/- {out['amp'][1]:.3f} a.u., \n"
                f" initial offset = {out['initial_offset'][0]:.2f} +/- {out['initial_offset'][1]:.3f}, \n"
                f" final_offset = {out['final_offset'][0]:.2f} +/- {out['final_offset'][1]:.3f} a.u."
            )
        # Plot the data and the fitting function if plot=True
        if plot:
            plt.plot(x_data, fit_type(x, popt) * y_normal)
            plt.plot(
                x_data,
                y_data,
                ".",
                label=f"T2  = {out['T2'][0]:.1f} +/- {out['T2'][1]:.1f}ns \n f = {out['f'][0] * 1000:.3f} +/- {out['f'][1] * 1000:.3f} MHz",
            )
            plt.legend(loc="upper right")
        t2r_est = out['T2'][0] #in ns
        t2r_err = out['T2'][1] #in ns
        return fit_type(x, popt) * y_normal, t2r_est, t2r_err

    def run(self, soccfg, soc):
        # defaults to 5, just make it to only look at this qubit
        res_gains = self.set_res_gain_ge(self.QubitIndex)
        self.config.update([('res_gain_ge', res_gains)])

        # now update for qubit frequency
        #current_freqs = self.config['qubit_freq_ge']
        #current_freqs[self.QubitIndex] = self.qubit_freq  # update with found freq from Qubit Spec
        # self.config.update([('qubit_freq_ge', self.qubit_freq)])

        # look at the config before we do the experiment
        print(f'Q {self.QubitIndex + 1} Round {self.round_num} T2R configuration: ', self.config)

        now = datetime.datetime.now()

        ramsey = T2RProgram(soccfg, reps=self.exp_cfg['reps'], final_delay=self.exp_cfg['relax_delay'],
                         cfg=self.config)
        # for live plotting
        if self.live_plot:
            I = Q = expt_mags = expt_phases = expt_pop = None
            viz = visdom.Visdom()
            assert viz.check_connection(timeout_seconds=5), "Visdom server not connected!"
            viz.close(win=None)  # close previous plots
            win1 = viz.line(X=np.arange(0, 1), Y=np.arange(0, 1),
                            opts=dict(height=400, width=700, title='T2 Ramsey Experiment', showlegend=True,
                                      xlabel='expt_pts'))

            for ii in range(self.config["rounds"]):
                iq_list = ramsey.acquire(soc, soft_avgs=1, progress=True)
                delay_times = ramsey.get_time_param('wait', "t", as_array=True)

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

        else:
            iq_list = ramsey.acquire(soc, soft_avgs=self.exp_cfg['rounds'], progress=True)
            I = iq_list[self.QubitIndex][0, :, 0]
            Q = iq_list[self.QubitIndex][0, :, 1]
            delay_times = ramsey.get_time_param('wait', "t", as_array=True)

        fit,t2r_est, t2r_err = self.t2_fit(delay_times, I)
        I, Q = self.plot_results(I, Q, delay_times, now, self.config, self.QubitIndex, fit, t2r_est, t2r_err)
        return  t2r_est, t2r_err, I, Q, delay_times, fit


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

    def plot_results(self, I, Q, delay_times, now, config, QubitIndex, fit, t2r_est, t2r_err):
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
            ax1.plot(delay_times, fit, '-', color='red', linewidth=3, label="Fit")
        else:
            ax2.plot(delay_times, fit, '-', color='red', linewidth=3, label="Fit")

        # Adjust spacing
        plt.tight_layout()

        # Calculate the middle of the plot area
        plot_middle = (ax1.get_position().x0 + ax1.get_position().x1) / 2

        # Add title, centered on the plot area
        fig.text(plot_middle, 0.98,
                 f"T1 Q{QubitIndex + 1}, pi gain %.2f" % config[
                     'pi_amp'] + f", {config['sigma'] * 1000} ns sigma" + f", {config['reps']}*{config['rounds']} avgs," + f" T1 = {t2r_est:.3f} ± {t2r_err:.3f} µs",
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
        return I, Q

