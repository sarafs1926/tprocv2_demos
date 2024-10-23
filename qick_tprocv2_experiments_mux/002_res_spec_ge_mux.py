from qick import *

# the main program class
from qick.asm_v2 import AveragerProgramV2
# for defining sweeps
from qick.asm_v2 import QickSpan, QickSweep1D

import json
import datetime
import matplotlib.pyplot as plt
import numpy as np

# Used for live plotting, need to run "python -m visdom.server" in the terminal and open the IP address in browser
# import visdom
import datetime
from tqdm import tqdm

from build_state import *
from expt_config import *
from system_config import * # Proxy is loaded in here


# ----- Experiment configurations ----- #
expt_name = "res_spec"
QubitIndex = QUBIT_INDEX
Qubit = 'Q' + str(QubitIndex)

##################
# Define Program #
##################

class SingleToneSpectroscopyProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        ro_chs = cfg['ro_ch']
        res_ch = cfg['res_ch']

        self.declare_gen(ch=res_ch, nqz=cfg['nqz_res'], ro_ch=ro_chs[0],
                         mux_freqs=cfg['res_freq_ge'],
                         mux_gains=cfg['res_gain_ge'],
                         mux_phases=cfg['res_phase'],
                         mixer_freq=cfg['mixer_freq'])
        for ch, f, ph in zip(cfg['ro_ch'], cfg['res_freq_ge'], cfg['ro_phase']):
            self.declare_readout(ch=ch, length=cfg['res_length'], freq=f, phase=ph, gen_ch=res_ch)

        self.add_pulse(ch=res_ch, name="mymux",
                       style="const",
                       length=cfg["res_length"],
                       mask=[0, 1, 2, 3, 4, 5],
                       )

    def _body(self, cfg):
        self.trigger(ros=cfg['ro_ch'], pins=[0], t=cfg['trig_time'], ddr4=True)
        self.pulse(ch=cfg['res_ch'], name="mymux", t=0)
        # relax delay ...


###################
# Run the Program
###################

# N benchmark
n = 1
j=0
while j < n:
    j+=1

    exp_cfg = expt_cfg[expt_name]
    q_config = all_qubit_state(system_config)
    config = {**q_config['Q' + str(QubitIndex)], **exp_cfg}
    print(config)

    fpts = exp_cfg["start"] + exp_cfg["step_size"] * np.arange(exp_cfg["steps"])
    fcenter = config['res_freq_ge']
    amps = np.zeros((len(fcenter), len(fpts)))
    for index, f in enumerate(tqdm(fpts)):
        config["res_freq_ge"] = fcenter + f  # array of 4 pulse freq
        prog = SingleToneSpectroscopyProgram(soccfg, reps=exp_cfg["reps"], final_delay=0.5, cfg=config)
        iq_list = prog.acquire(soc, soft_avgs=exp_cfg["rounds"], progress=False)
        for i in range(len(config['res_freq_ge'])):
            amps[i][index] = np.abs(iq_list[i][:, 0] + 1j * iq_list[i][:, 1])
    amps = np.array(amps)

    ###################
    # Plot
    ###################
    res_freqs = []

    # Increase figure size
    plt.figure(figsize=(12, 8))

    # Set larger font sizes
    plt.rcParams.update({
        'font.size': 14,  # Base font size
        'axes.titlesize': 18,  # Title font size
        'axes.labelsize': 16,  # Axis label font size
        'xtick.labelsize': 14,  # X-axis tick label size
        'ytick.labelsize': 14,  # Y-axis tick label size
        'legend.fontsize': 14,  # Legend font size
    })

    for i in range(6):
        plt.subplot(2, 3, i + 1)
        plt.plot(fpts + config['res_freq_ge'][i], amps[i], '-', linewidth=1.5)
        freq_r = fpts[np.argmin(amps[i])] + config['res_freq_ge'][i]
        res_freqs.append(freq_r)
        plt.axvline(freq_r, linestyle='--', color='orange', linewidth=1.5)
        plt.xlabel("Frequency (MHz)", fontweight='normal')
        plt.ylabel("Amplitude (a.u.)", fontweight='normal')
        plt.title(f"Resonator {i + 1} {freq_r:.3f} MHz", pad=10)

        # Adjust y-axis limits to show more of the dip
        plt.ylim(plt.ylim()[0] - 0.05 * (plt.ylim()[1] - plt.ylim()[0]), plt.ylim()[1])

    # Add a main title to the figure
    plt.suptitle(f"MUXed resonator spectroscopy {exp_cfg['reps']} avgs", fontsize=24, y=0.95)

    plt.tight_layout(pad=2.0)

    ### Save figure
    outerFolder_expt = outerFolder + "/" + expt_name + "/"
    create_folder_if_not_exists(outerFolder_expt)
    now = datetime.datetime.now()
    formatted_datetime = now.strftime("%Y-%m-%d_%H-%M-%S")
    file_name = outerFolder_expt + f"{formatted_datetime}_" + expt_name + ".png"
    plt.savefig(file_name, dpi=300);

    ### Save params
    res_gains = config["res_gain_ge"]
    res_freqs = [round(x, 3) for x in res_freqs]
    print("Resonator gains:", res_gains)
    print("Resonator freqs:", res_freqs)

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
#     f.append('t', t)
#
#     f.append('avgi', iq_list[0].T[0])
#     f.append('avgq', iq_list[0].T[1])
#
#     f.attrs['config'] = json.dumps(config) # cant save configs yet with QickParams
#
# data = data_path + '\\' + fname
# print('data path:', data)
#
#
