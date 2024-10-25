from build_task import *
from build_state import *
from expt_config import *
from system_config import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import datetime


def lorentzian(f, f0, gamma, A, B):
    return A * gamma ** 2 / ((f - f0) ** 2 + gamma ** 2) + B


def max_offset_difference_with_x(x_values, y_values, offset):
    max_average_difference = -1
    corresponding_x = None

    # average all 3 to avoid noise spikes
    for i in range(len(y_values) - 2):
        # group 3 vals
        y_triplet = y_values[i:i + 3]

        # avg differences for these 3 vals
        average_difference = sum(abs(y - offset) for y in y_triplet) / 3

        # see if this is the highest difference yet
        if average_difference > max_average_difference:
            max_average_difference = average_difference
            # x value for the middle y value in the 3 vals
            corresponding_x = x_values[i + 1]

    return corresponding_x, max_average_difference


def fit_lorenzian(I, Q, freqs, freq_q):
    # guesses
    initial_guess_I = [freq_q, 1, np.max(I), np.min(I)]  # x guess (which is very off here), amplitude guess, offset
    initial_guess_Q = [freq_q, 1, np.max(Q), np.min(Q)]

    # fitting the Lorentzian
    params_I, _ = curve_fit(lorentzian, freqs, I, p0=initial_guess_I)
    params_Q, _ = curve_fit(lorentzian, freqs, Q, p0=initial_guess_Q)

    x_max_diff_I, max_diff_I = max_offset_difference_with_x(freqs, I, params_I[3])
    x_max_diff_Q, max_diff_Q = max_offset_difference_with_x(freqs, Q, params_Q[3])

    # guesses
    initial_guess_I = [x_max_diff_I, 1, np.max(I),
                       np.min(I)]  # x guess (which is now accurate), amplitude guess, offset
    initial_guess_Q = [x_max_diff_Q, 1, np.max(Q), np.min(Q)]

    # fitting the Lorentzian
    params_I, _ = curve_fit(lorentzian, freqs, I, p0=initial_guess_I)
    params_Q, _ = curve_fit(lorentzian, freqs, Q, p0=initial_guess_Q)

    # make line from the fits
    I_fit = lorentzian(freqs, *params_I)
    Q_fit = lorentzian(freqs, *params_Q)

    mean_I = params_I[0]  # the mean which is from f0 from the fitted parameters for I
    mean_Q = params_Q[0]  # the mean which is from f0 from the fitted parameters for Q

    # find which fit has the widest curve, becasue data is so noisy im going to assume a thin curve is fitting to a noise peak if there is nothing there
    fwhm_I = 2 * params_I[1]
    fwhm_Q = 2 * params_Q[1]

    # Determine which fit has the widest curve
    if fwhm_I > fwhm_Q:
        widest_fit = "I"
        widest_curve_mean = mean_I
        widest_fwhm = fwhm_I
    else:
        widest_fit = "Q"
        widest_curve_mean = mean_Q
        widest_fwhm = fwhm_Q

    # Print the FWHM for the fit with the widest curve
    print(f"The widest FWHM is for the {widest_fit} data: {widest_fwhm}")

    return mean_I, mean_Q, I_fit, Q_fit, widest_curve_mean, widest_fwhm


def create_folder_if_not_exists(folder_path):
    import os
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)


class QubitSpectroscopy:
    def __init__(self, QubitIndex, outerFolder):
        self.QubitIndex = QubitIndex
        self.outerFolder = outerFolder
        self.expt_name = "qubit_spec_ge"
        self.Qubit = 'Q' + str(self.QubitIndex)
        self.exp_cfg = add_qubit_experiment(expt_cfg, self.expt_name, self.QubitIndex)
        self.q_config = all_qubit_state(system_config)
        self.config = {**self.q_config['Q' + str(self.QubitIndex)], **self.exp_cfg}

    def run(self, soccfg, soc):
        qspec = PulseProbeSpectroscopyProgram(soccfg, reps=self.config['reps'], final_delay=0.5, cfg=self.config)
        iq_list = qspec.acquire(soc, soft_avgs=self.exp_cfg["rounds"], progress=True)
        freqs = qspec.get_pulse_param('qubit_pulse', "freq", as_array=True)

        # Prepare your data
        I = iq_list[self.QubitIndex][0, :, 0]
        Q = iq_list[self.QubitIndex][0, :, 1]
        freqs = np.array(freqs)
        freq_q = freqs[np.argmax(I)]

        mean_I, mean_Q, I_fit, Q_fit, widest_curve_mean, widest_fwhm = fit_lorenzian(I, Q, freqs, freq_q)

        # Plotting and saving (same as before)
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        plt.rcParams.update({'font.size': 18})

        # I subplot
        ax1.plot(freqs, I, label='I', linewidth=2)
        ax1.plot(freqs, I_fit, 'r--', label='Lorentzian Fit')
        ax1.set_ylabel("I Amplitude (a.u.)", fontsize=20)
        ax1.tick_params(axis='both', which='major', labelsize=16)
        ax1.axvline(widest_curve_mean, color='orange', linestyle='--', linewidth=2)
        ax1.legend()

        # Q subplot
        ax2.plot(freqs, Q, label='Q', linewidth=2)
        ax2.plot(freqs, Q_fit, 'r--', label='Lorentzian Fit')
        ax2.set_xlabel("Qubit Frequency (MHz)", fontsize=20)
        ax2.set_ylabel("Q Amplitude (a.u.)", fontsize=20)
        ax2.tick_params(axis='both', which='major', labelsize=16)
        ax2.axvline(widest_curve_mean, color='orange', linestyle='--', linewidth=2)
        ax2.legend()

        # Adjust spacing
        plt.tight_layout()

        # Calculate the middle of the plot area
        plot_middle = (ax1.get_position().x0 + ax1.get_position().x1) / 2

        # Add title, centered on the plot area
        fig.text(plot_middle, 0.98,
                 f"Qubit Spectroscopy Q{self.QubitIndex + 1}, %.2f MHz" % widest_curve_mean + f" FWHM: {round(widest_fwhm, 1)}" + f", {self.config['reps']} avgs",
                 fontsize=24, ha='center', va='top')

        # Adjust the top margin to make room for the title
        plt.subplots_adjust(top=0.93)

        ### Save figure
        outerFolder_expt = self.outerFolder + "/" + self.expt_name + "/"
        create_folder_if_not_exists(outerFolder_expt)
        now = datetime.datetime.now()
        formatted_datetime = now.strftime("%Y-%m-%d_%H-%M-%S")
        file_name = outerFolder_expt + f"{formatted_datetime}_" + self.expt_name + f"_q{self.QubitIndex + 1}.png"

        fig.savefig(file_name, dpi=300, bbox_inches='tight')  # , facecolor='white'
        print("Saving to: ", file_name)
        plt.close(fig)


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
                       mask=[0, 1, 2, 3, 4, 5],
                       )

        self.declare_gen(ch=qubit_ch, nqz=cfg['nqz_qubit'])
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

