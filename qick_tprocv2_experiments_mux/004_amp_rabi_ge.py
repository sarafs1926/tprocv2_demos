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

# ----- Experiment configurations ----- #
expt_name = "power_rabi_ge"
QubitIndex = QUBIT_INDEX
Qubit = 'Q' + str(QubitIndex)

exp_cfg = add_qubit_experiment(expt_cfg, expt_name, QubitIndex)
q_config = all_qubit_state(system_config)
config = {**q_config['Q' + str(QubitIndex)], **exp_cfg}
print(config)


##################
# Define Program #
##################

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
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse", t=0)  # play probe pulse
        self.delay_auto(t=0.01, tag='waiting')  # Wait til qubit pulse is done before proceeding
        self.pulse(ch=cfg['res_ch'], name="res_pulse", t=0)
        self.trigger(ros=cfg['ro_ch'], pins=[0], t=cfg['trig_time'])
        # relax delay ...


###################
# Run the Program
###################

amp_rabi = AmplitudeRabiProgram(soccfg, reps=exp_cfg['reps'], final_delay=exp_cfg['relax_delay'], cfg=config)
iq_list = amp_rabi.acquire(soc, soft_avgs=exp_cfg["rounds"], progress=True)
gains = amp_rabi.get_pulse_param('qubit_pulse', "gain", as_array=True)


def cosine(x, a, b, c, d):
    return a * np.cos(2. * np.pi * b * x - c * 2 * np.pi) + d


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
plt.rcParams.update({'font.size': 18})  # Set base font size

I = iq_list[QubitIndex][0, :, 0]
Q = iq_list[QubitIndex][0, :, 1]

## Fit
q1_amp = I
q1_a_guess = (np.max(q1_amp) - np.min(q1_amp)) / 2
q1_b_guess = 1 / gains[-1]
q1_c_guess = 0
q1_d_guess = np.mean(q1_amp)

q1_guess = [q1_a_guess, q1_b_guess, q1_c_guess, q1_d_guess]
q1_popt, q1_pcov = curve_fit(cosine, gains, q1_amp, maxfev=100000, p0=q1_guess)
q1_fit_cosine = cosine(gains, *q1_popt)

# figure out whether to use argmax or argmin based on height of first and last values
first_three_avg = np.mean(q1_fit_cosine[:3])
last_three_avg = np.mean(q1_fit_cosine[-3:])

if last_three_avg > first_three_avg:
    pi_amp = gains[np.argmax(q1_fit_cosine)]
else:
    pi_amp = gains[np.argmin(q1_fit_cosine)]

# pi_amp = gains[np.argmax(q1_fit_cosine)] #change to argmax/min depending on shape of curve

# I subplot
ax1.plot(gains, I, label="Gain (a.u.)", linewidth=2)
ax1.set_ylabel("I Amplitude (a.u.)", fontsize=20)
ax1.plot(gains, q1_fit_cosine, '-', color='red', linewidth=3, label="Fit")
ax1.tick_params(axis='both', which='major', labelsize=16)
# ax1.axvline(freq_q, color='orange', linestyle='--', linewidth=2)

# Q subplot
ax2.plot(gains, Q, label="Q", linewidth=2)
ax2.set_xlabel("Gain (a.u.)", fontsize=20)
ax2.set_ylabel("Q Amplitude (a.u.)", fontsize=20)
ax2.tick_params(axis='both', which='major', labelsize=16)
# ax2.axvline(freq_q, color='orange', linestyle='--', linewidth=2)

# Adjust spacing
plt.tight_layout()

# Calculate the middle of the plot area
plot_middle = (ax1.get_position().x0 + ax1.get_position().x1) / 2

# Add title, centered on the plot area
fig.text(plot_middle, 0.98,
         f"Rabi Q{QubitIndex + 1}, pi gain %.2f" % pi_amp + f", {config['sigma'] * 1000} ns sigma" +  f", {config['reps']}*{config['rounds']} avgs",
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