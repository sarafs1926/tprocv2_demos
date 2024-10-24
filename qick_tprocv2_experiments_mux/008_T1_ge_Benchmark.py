from orca.orca import start

from qick import *

# the main program class
from qick.asm_v2 import AveragerProgramV2
# for defining sweeps
from qick.asm_v2 import QickSpan, QickSweep1D

import json
import datetime
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import pickle

# Used for live plotting, need to run "python -m visdom.server" in the terminal and open the IP address in browser
# import visdom

from build_task import *
from build_state import *
from expt_config import *
from system_config import *

# from qualang_tools.plot import Fit
import pprint as pp

# import single shot information for g-e calibration
# from SingleShot import SingleShotProgram_g, SingleShotProgram_e
# from SingleShot import config as ss_config

# ----- Experiment configurations ----- #
expt_name = "T1_ge"
QubitIndex = QUBIT_INDEX
Qubit = 'Q' + str(QubitIndex)

exp_cfg = add_qubit_experiment(expt_cfg, expt_name, QubitIndex)
q_config = all_qubit_state(system_config)
config = {**q_config['Q' + str(QubitIndex)], **exp_cfg}
print(config)

##################
# Define Program #
##################
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
        self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse", t=0)  #play probe pulse
        self.delay_auto(cfg['wait_time']+0.01, tag='wait') # wait_time after last pulse
        self.pulse(ch=cfg['res_ch'], name="res_pulse", t=0)
        self.trigger(ros=cfg['ro_ch'], pins=[0], t=cfg['trig_time'])
        # relax delay ...

###################
# Run the Program
###################

# N benchmark
n = 40000 # 2

q1_t1 = [];  q1_t1_err = []; dates = [];

j = 0
start_time = datetime.datetime.now()
formatted_starttime = start_time.strftime("%Y-%m-%d_%H-%M-%S")
while j < n:
    j += 1
    now = datetime.datetime.now()
    dates.append(now)
    print("iteration is",j)

    t1=T1Program(soccfg, reps=exp_cfg['reps'], final_delay=exp_cfg['relax_delay'], cfg=config)
    iq_list = t1.acquire(soc, soft_avgs=exp_cfg['rounds'], progress=True)
    delay_times = t1.get_time_param('wait', "t", as_array=True)


    def exponential(x, a, b, c, d):
        return a * np.exp(- (x - b) / c) + d


    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    plt.rcParams.update({'font.size': 18})  # Set base font size

    I = iq_list[QubitIndex][0, :, 0]
    Q = iq_list[QubitIndex][0, :, 1]

    # ## Fit
    # q1_a_guess = np.max(Q)-np.min(Q)   #change to np.max(Q)-np.min(Q) depending on direction

    # figure out which way the curve is oriented based on the averages of first and last 3 points
    signal = Q

    ### Olivia's way
    # average_first_three = np.mean(signal[:3])
    # average_last_three = np.mean(signal[-3:])
    #
    # if average_first_three < average_last_three:
    #     q1_a_guess = np.max(signal) - np.min(signal)
    # else:
    #     q1_a_guess = np.min(signal) - np.max(signal)
    #
    # Sara's way
    if signal[-1] > signal[0]:
        q1_a_guess = -(np.max(signal)-np.min(signal))
    else:
        q1_a_guess = (np.max(signal)-np.min(signal))

    q1_b_guess = 0
    q1_c_guess = delay_times[-1] / 5 # if not fitting change /6 or /5
    q1_d_guess = np.min(signal)

    q1_guess = [q1_a_guess, q1_b_guess, q1_c_guess, q1_d_guess]
    q1_popt, q1_pcov = curve_fit(exponential, delay_times, signal, maxfev=100000, p0=q1_guess)
    q1_fit_exponential = exponential(delay_times, *q1_popt)

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

    # Adjust spacing
    plt.tight_layout()

    # Calculate the middle of the plot area
    plot_middle = (ax1.get_position().x0 + ax1.get_position().x1) / 2

    # Add title, centered on the plot area
    fig.text(plot_middle, 0.98,
             f"T1 Q{QubitIndex + 1}, pi gain %.2f" % config['pi_amp'] + f", {config['sigma'] * 1000} ns sigma" + f", {config['reps']}*{config['rounds']} avgs," + f" T1 = {T1_est:.3f} ± {T1_err:.3f} µs",
             fontsize=24, ha='center', va='top')

    # Adjust the top margin to make room for the title
    plt.subplots_adjust(top=0.93)

    ### Save figure
    outerFolder_expt = outerFolder + "/" + expt_name + "/"
    create_folder_if_not_exists(outerFolder_expt)
    formatted_datetime = now.strftime("%Y-%m-%d_%H-%M-%S")
    file_name = outerFolder_expt + f"{formatted_datetime}_" + expt_name + f"_q{QubitIndex+1}.png"
    fig.savefig(file_name, dpi=300, bbox_inches='tight')  # , facecolor='white'
    plt.close(fig)

    ### Save the T1 and error into arrays and pickle them
    q1_t1.append(T1_est)
    q1_t1_err.append(T1_err)
    print(q1_t1, q1_t1_err, dates)

    with open(outerFolder_expt + "T1_benchmark_Q1_" + f"{formatted_starttime}_"+ str(n) + "x.pkl", 'wb') as f:
        pickle.dump(q1_t1, f)
        pickle.dump(q1_t1_err, f)
        pickle.dump(dates,f)
        pickle.dump(I, f)
        pickle.dump(Q, f)
        pickle.dump(delay_times, f)


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
#     f.append('delay_times', delay_times)
#     f.append('amps', amps)
#     if IS_VISDOM:
#         f.append('avgi', expt_I)
#         f.append('avgq', expt_Q)
#     else:
#         f.append('avgi', iq_list[0][0].T[0])
#         f.append('avgq', iq_list[0][0].T[1])
#
#     del config['wait_time']  # cant save QickParam
#     # formats config into file as a single line
#     f.attrs['config'] = json.dumps(config)
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
