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
import h5py
import os

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

        # Print to verify the configuration
        print(f"Initializing with res_gain_ge: {cfg['res_gain_ge']}")

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
        # relax delay ...


###################
# Run the Program
###################

# N benchmark
n = 400000  # 2

q1_t1 = [];
q1_t1_err = [];
dates = [];

j = 0
start_time = datetime.datetime.now()
formatted_starttime = start_time.strftime("%Y-%m-%d_%H-%M-%S")
while j < n:
    j += 1
    now = datetime.datetime.now()
    dates.append(now)
    print("iteration is", j)

    t1 = T1Program(soccfg, reps=exp_cfg['reps'], final_delay=exp_cfg['relax_delay'], cfg=config)
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
    signal = I

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
    """"
    if signal[-1] > signal[0]:
        q1_a_guess = -(np.max(signal) - np.min(signal))
    else:
        q1_a_guess = (np.max(signal) - np.min(signal))

    q1_b_guess = 0
    q1_c_guess = delay_times[-1] / 6  # if not fitting change /6 or /5
    q1_d_guess = np.min(signal)

    q1_guess = [q1_a_guess, q1_b_guess, q1_c_guess, q1_d_guess]
    q1_popt, q1_pcov = curve_fit(exponential, delay_times, signal, maxfev=100000, p0=q1_guess)
    q1_fit_exponential = exponential(delay_times, *q1_popt)

    T1_est = q1_popt[2]
    T1_err = np.sqrt(q1_pcov[2][2])

    """

    #Trying new fitting
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
    q1_popt, q1_pcov = curve_fit(exponential, delay_times, signal,
                                 p0=q1_guess, bounds=(lower_bounds, upper_bounds),
                                 method='trf', maxfev=10000)

    # Generate the fitted exponential curve
    q1_fit_exponential = exponential(delay_times, *q1_popt)

    # Extract T1 and its error
    T1_est = q1_popt[2]  # Decay constant T1
    T1_err = np.sqrt(q1_pcov[2][2]) if q1_pcov[2][2] >= 0 else float('inf')  # Ensure error is valid

    # I subplot
    ax1.plot(delay_times, I, label="Gain (a.u.)", linewidth=2)
    ax1.set_ylabel("I Amplitude (a.u.)", fontsize=20)
    ax1.plot(delay_times, q1_fit_exponential, '-', color='red', linewidth=3, label="Fit")
    ax1.tick_params(axis='both', which='major', labelsize=16)
    # ax1.axvline(freq_q, color='orange', linestyle='--', linewidth=2)

    # Q subplot
    ax2.plot(delay_times, Q, label="Q", linewidth=2)
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
             f"T1 Q{QubitIndex + 1}, pi gain %.2f" % config[
                 'pi_amp'] + f", {config['sigma'] * 1000} ns sigma" + f", {config['reps']}*{config['rounds']} avgs," + f" T1 = {T1_est:.3f} ± {T1_err:.3f} µs",
             fontsize=24, ha='center', va='top')

    # Adjust the top margin to make room for the title
    plt.subplots_adjust(top=0.93)

    ### Save figure
    outerFolder_expt = outerFolder + "/" + expt_name + "/"
    create_folder_if_not_exists(outerFolder_expt)
    formatted_datetime = now.strftime("%Y-%m-%d_%H-%M-%S")
    file_name = outerFolder_expt + f"{formatted_datetime}_" + expt_name + f"_q{QubitIndex + 1}.png"
    fig.savefig(file_name, dpi=300, bbox_inches='tight')  # , facecolor='white'
    plt.close(fig)

    ### Save the T1 and error into arrays
    q1_t1.append(T1_est)
    q1_t1_err.append(T1_err)
    print(q1_t1, q1_t1_err, dates)

    # Create the HDF5 file and save data
    h5_filename = outerFolder_expt + f"{formatted_datetime}_" + expt_name + f"_q{QubitIndex + 1}.h5"
    with h5py.File(h5_filename, 'w') as f:
        # Save T1 and error
        f.create_dataset("T1_estimates", data=np.array(q1_t1))
        f.create_dataset("T1_errors", data=np.array(q1_t1_err))

        # Save dates
        f.create_dataset("dates", data=np.array(dates, dtype='S'))

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