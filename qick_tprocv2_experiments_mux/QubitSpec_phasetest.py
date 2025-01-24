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
expt_name = "qubit_spec_ge"
QubitIndex = 3
Qubit = 'Q' + str(QubitIndex)

exp_cfg = add_qubit_experiment(expt_cfg, expt_name, QubitIndex)
q_config = all_qubit_state(system_config)
config = {**q_config['Q' + str(QubitIndex)], **exp_cfg}
print(config)


##################
# Define Program #
##################

class PulseProbeSpectroscopyProgram(AveragerProgramV2):
    def __init__(self, cfg, list_of_all_qubits, **kwargs):
        super().__init__(cfg, **kwargs)
        self.list_of_all_qubits = list_of_all_qubits

    def _initialize(self, cfg):
        super()._initialize(cfg)
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
                       mask=self.list_of_all_qubits,
                       )

        self.declare_gen(ch=qubit_ch, nqz=cfg['nqz_qubit'], mixer_freq=cfg['mixer_freq_q'])
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

# mean_I, mean_Q, I_fit, Q_fit, widest_curve_mean, widest_fwhm = fit_lorenzian(I, Q, freqs, freq_q)

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
         f"Qubit Spectroscopy Q{QubitIndex + 1}, %.2f MHz" % widest_curve_mean + f" FWHM: {round(widest_fwhm, 1)}" + f", {config['reps']} avgs",
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
