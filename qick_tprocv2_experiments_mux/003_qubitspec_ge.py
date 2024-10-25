from qick import *

# the main program class
from qick.asm_v2 import AveragerProgramV2
# for defining sweeps
from qick.asm_v2 import QickSpan, QickSweep1D

import json
import datetime

# Used for live plotting, need to run "python -m visdom.server" in the terminal and open the IP address in browser
# import visdom

from build_task import *
from build_state import *
from expt_config import *
from system_config import *

# from qualang_tools.plot import Fit
import pprint as pp
import matplotlib.pyplot as plt
import numpy as np

# # import single shot information for g-e calibration
# from SingleShot import SingleShotProgram_g, SingleShotProgram_e
# from SingleShot import config as ss_config
from scipy.optimize import curve_fit


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

    # print("Offset for I data:", params_I[3])  #B in lorenzian def
    # print("Offset for Q data:", params_Q[3])  #B in lorenzian def

    x_max_diff_I, max_diff_I = max_offset_difference_with_x(freqs, I, params_I[3])
    x_max_diff_Q, max_diff_Q = max_offset_difference_with_x(freqs, Q, params_Q[3])

    # print(f"Max difference for I data: {max_diff_I} at x = {x_max_diff_I}")
    # print(f"Max difference for Q data: {max_diff_Q} at x = {x_max_diff_Q}")

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
        # relax delay ...


###################
# Run the Program
###################

qspec = PulseProbeSpectroscopyProgram(soccfg, reps=config['reps'], final_delay=0.5, cfg=config)
iq_list = qspec.acquire(soc, soft_avgs=exp_cfg["rounds"], progress=True)
freqs = qspec.get_pulse_param('qubit_pulse', "freq", as_array=True)

# Prepare your data
I = iq_list[QubitIndex][0, :, 0]
Q = iq_list[QubitIndex][0, :, 1]
freqs = np.array(freqs)
freq_q = freqs[np.argmax(I)]

mean_I, mean_Q, I_fit, Q_fit, widest_curve_mean, widest_fwhm = fit_lorenzian(I, Q, freqs, freq_q)

# Plot the data and fits
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
         f"Qubit Spectroscopy Q{QubitIndex + 1}, %.2f MHz" % widest_curve_mean + f" FWHM: {round(widest_fwhm, 1)}" + f", {config['reps']}*{config['rounds']} avgs",
         fontsize=24, ha='center', va='top')

# Adjust the top margin to make room for the title
plt.subplots_adjust(top=0.93)

### Save figure
outerFolder_expt = outerFolder + "/" + expt_name + "/"
create_folder_if_not_exists(outerFolder_expt)
now = datetime.datetime.now()
formatted_datetime = now.strftime("%Y-%m-%d_%H-%M-%S")
file_name = outerFolder_expt + f"{formatted_datetime}_" + expt_name + f"_q{QubitIndex + 1}.png"

fig.savefig(file_name, dpi=300, bbox_inches='tight')  # , facecolor='white'
plt.close(fig)

# #####################################
# # ----- Saves data to a file ----- #
# #####################################
#
# prefix = str(datetime.date.today())
# exp_name = expt_name + '_Q' + str(QubitIndex) + '_' + prefix
# print('Experiment name: ' + exp_name)
#
# data_path = DATA_PATH
#
# fname = get_next_filename(data_path, exp_name, suffix='.h5')
# print('Current data file: ' + fname)
#
# with SlabFile(data_path + '\\' + fname, 'a') as f:
#     # 'a': read/write/create
#
#     # - Adds data to the file - #
#     f.append('fpts', freqs)
#     f.append('amps', amps)
#     if IS_VISDOM:
#         f.append('avgi', expt_I)
#         f.append('avgq', expt_Q)
#     else:
#         f.append('avgi', iq_list[0][0].T[0])
#         f.append('avgq', iq_list[0][0].T[1])
#
#     del config['qubit_freq_ge']  # cant save QickParam
#     # formats config into file as a single line
#     f.attrs['config'] = json.dumps(config) # cant save configs yet with QickParams
#     f.attrs['fit_result'] = json.dumps(fit_result)
#
#     # - Adds ss data to the file - #
#     f.append('I_g', I_g)
#     f.append('Q_g', Q_g)
#     f.append('I_e', I_e)
#     f.append('Q_e', Q_e)
#
#     # formats config into file as a single line
#     f.attrs['ss_config'] = json.dumps(ss_config)
#
# data = data_path + '\\' + fname