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

from build_state import *
from expt_config import *
from system_config import * # Proxy is loaded in here


# ----- Experiment configurations ----- #
expt_name = "tof"
QubitIndex = QUBIT_INDEX
Qubit = 'Q' + str(QubitIndex)

exp_cfg = expt_cfg[expt_name]
q_config = all_qubit_state(system_config)
config = {**q_config['Q' + str(QubitIndex)], **exp_cfg}

# Update parameters to see TOF pulse with your setup
config.update([('trig_time', 0.0)]) #('res_gain', 0.8), ('res_length', 0.5), ('ro_length', 1.5)])
print(config)

##################
# Define Program #
##################

class MuxProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        ro_chs = cfg['ro_ch']
        gen_ch = cfg['res_ch']

        self.declare_gen(ch=gen_ch, nqz=cfg['nqz_res'], ro_ch=ro_chs[0],
                         mux_freqs=cfg['res_freq_ge'],
                         mux_gains=cfg['res_gain_ge'],
                         mux_phases=cfg['res_phase'],
                         mixer_freq=cfg['mixer_freq'])
        for ch, f, ph in zip(cfg['ro_ch'], cfg['res_freq_ge'], cfg['ro_phase']):
            self.declare_readout(ch=ch, length=cfg['res_length'], freq=f, phase=ph, gen_ch=gen_ch)

        self.add_pulse(ch=gen_ch, name="mymux",
                       style="const",
                       length=cfg["res_length"],
                       mask=[0, 1, 2, 3, 4, 5],
                       )

    def _body(self, cfg):
        self.trigger(ros=cfg['ro_ch'], pins=[0], t=cfg['trig_time'], ddr4=True)
        self.pulse(ch=cfg['res_ch'], name="mymux", t=0)


###################
# Run the Program
###################

# N benchmark
n = 1 #2
j=0
while j < n:
    j+=1
    prog = MuxProgram(soccfg, reps=1, final_delay=0.5, cfg=config)
    iq_list = prog.acquire_decimated(soc, soft_avgs=config['soft_avgs'])

    ###################
    # Plot
    ###################

    t = prog.get_time_axis(ro_index=0)
    fig, axes = plt.subplots(len(config['ro_ch']), 1, figsize=(12,12))
    for i, ch in enumerate(config['ro_ch']):
        plot = axes[i]
        plot.plot(t, iq_list[i][:,0], label="I value")
        plot.plot(t, iq_list[i][:,1], label="Q value")
        plot.plot(t, np.abs(iq_list[i].dot([1,1j])), label="magnitude")
        plot.legend()
        plot.set_ylabel("a.u.")
        plot.set_xlabel("us");
        ######## SET THIS VALUE AS YOUR config['trig_time'] #########
        plot.axvline(0.75, c='r')

    # Save figure
    outerFolder_expt = outerFolder + "/" + expt_name + "/"
    create_folder_if_not_exists(outerFolder_expt)
    now = datetime.datetime.now()
    formatted_datetime = now.strftime("%Y-%m-%d_%H-%M-%S")
    file_name = outerFolder_expt + f"{formatted_datetime}_" + expt_name + ".png"
    plt.savefig(file_name, dpi=300);


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
