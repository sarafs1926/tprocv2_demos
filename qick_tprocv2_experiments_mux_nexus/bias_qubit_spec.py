from build_task import *
from build_state import *
from expt_config import *
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import csv
import datetime
import time

from NetDrivers import E36300

class BiasQubitSpectroscopy:
    def __init__(self, QubitIndex, outerFolder, experiment):
        self.QubitIndex = QubitIndex
        self.outerFolder = outerFolder
        self.expt_name = "bias_qubit_spec_ge"
        self.experiment = experiment
        self.Qubit = 'Q' + str(self.QubitIndex)
        self.exp_cfg = expt_cfg[self.expt_name]
        self.q_config = all_qubit_state(self.experiment)
        self.exp_cfg = add_qubit_experiment(expt_cfg, self.expt_name, self.QubitIndex)
        self.config = {**self.q_config[self.Qubit], **self.exp_cfg}

        print(f'Q {self.QubitIndex + 1} Qubit Spec configuration: ', self.config)

    def run(self, soccfg, soc, start_volt, stop_volt, volt_pts, plot_sweeps=True, plot_2d=True):

        vsweep = np.linspace(start_volt, stop_volt, volt_pts, endpoint=True)
        Is, Qs, amps, freqs = self.sweep_bias(soccfg, soc, vsweep)

        if plot_sweeps:
            self.plot_sweeps(vsweep, Is, Qs, freqs)

        if plot_2d:
            self.plot2d(vsweep, Is, Qs, amps, freqs)

        return

    def sweep_bias(self, soccfg, soc, vsweep):

        Bias_PS_ip = ['192.168.0.44', '192.168.0.44', '192.168.0.44', '192.168.0.41'] #IP address of bias PS (qubits 1-3 are the same PS)
        Bias_ch = [1, 2, 3, 1] #Channel number of qubit 1-4 on associated PS
        qubit_index = int(self.QubitIndex)

        print(f"Qubit_index {qubit_index}")
        print(f"PSip {Bias_PS_ip[qubit_index]}")
        print(f"PSch {Bias_ch[qubit_index]}")
        BiasPS = E36300(Bias_PS_ip[qubit_index], server_port = 5025)

        BiasPS.setVoltage(0, Bias_ch[qubit_index])
        BiasPS.enable(Bias_ch[qubit_index])

        I_arr = []
        Q_arr = []
        amps_arr = []
        freq_arr = []

        for index, v in enumerate(vsweep):
            print(f"Setting bias to {v}V")
            BiasPS.setVoltage(v, Bias_ch[qubit_index])
            time.sleep(5)


            qspec = PulseProbeSpectroscopyProgram(soccfg, reps=self.config['reps'], final_delay = self.exp_cfg['relax_delay'], cfg=self.config)
            iq_list = qspec.acquire(soc, soft_avgs = self.exp_cfg["rounds"], progress=True)
            I = iq_list[self.QubitIndex][0, :, 0]
            Q = iq_list[self.QubitIndex][0, :, 1]
            amps = np.abs(I + 1j * Q)
            freqs = qspec.get_pulse_param('qubit_pulse', "freq", as_array=True)
            I_arr.append(I)
            Q_arr.append(Q)
            amps_arr.append(amps)
            freq_arr.append(freqs)
        BiasPS.disable(Bias_ch[qubit_index])
        print(freq_arr[0])

        '''outerFolder_expt = os.path.join(self.outerFolder, 'bias_spec')
        self.experiment.create_folder_if_not_exists(outerFolder_expt)
        now = datetime.datetime.now()
        formatted_datetime = now.strftime("%Y-%m-%d_%H-%M-%S")
        file_name_Iarr = os.path.join(outerFolder_expt, f"{formatted_datetime}_BiasSpec_Q{self.QubitIndex + 1}_Iarr")
        file_name_Qarr = os.path.join(outerFolder_expt, f"{formatted_datetime}_BiasSpec_Q{self.QubitIndex + 1}_Qarr")
        file_name_Amparr = os.path.join(outerFolder_expt, f"{formatted_datetime}_BiasSpec_Q{self.QubitIndex + 1}_Amparr")
        #file_name_freqarr = os.path.join(outerFolder_expt, f"{formatted_datetime}_BiasSpec_Q{self.QubitIndex + 1}_freqarr")
        with open(f"{file_name_Iarr}.csv", 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(I_arr)
        with open(f"{file_name_Qarr}.csv", 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(Q_arr)
        with open(f"{file_name_Amparr}.csv", 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(amps_arr)
        # with open(f"{file_name_freqarr}.csv", 'w', newline='') as f:
        #     writer = csv.writer(f)
        #     writer.writerows(freq_arr)'''

        return I_arr, Q_arr, amps_arr, freq_arr

    def plot_sweeps(self, vsweep, I_arr, Q_arr, freq_arr):
        #plt.figure(figsize=(12, 8))

        # Set larger font sizes
        plt.rcParams.update({
            'font.size': 14,  # Base font size
            'axes.titlesize': 18,  # Title font size
            'axes.labelsize': 16,  # Axis label font size
            'xtick.labelsize': 14,  # X-axis tick label size
            'ytick.labelsize': 14,  # Y-axis tick label size
            'legend.fontsize': 14,  # Legend font size
        })

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex='all')
        #plt.rcParams.update({'font.size': 18})
        ax1.set_ylabel("I Amplitude (a.u.)", fontsize=20)
        ax1.tick_params(axis='both', which='major', labelsize=16)
        ax2.set_ylabel("Q Amplitude (a.u.)", fontsize=20)
        ax2.tick_params(axis='both', which='major', labelsize=16)
        ax2.set_xlabel("Qubit Frequency (MHz)", fontsize=20)

        for volt_index in range(len(vsweep)):
            ax1.plot(freq_arr[volt_index], I_arr[volt_index], linewidth=2, label=round(vsweep[volt_index],3))
            ax2.plot(freq_arr[volt_index], Q_arr[volt_index], line_width=2, label = round(vsweep[volt_index],3))

        ax1.legend(fontsize='6', title='Voltage')
        ax2.legend(fontsize='6', title='Voltage')
        fig.suptitle(f"Qubit Spectroscopy Q{self.QubitIndex+1} at Voltage Bias Points", fontsize=24)

        plt.tight_layout()

        # Adjust the top margin to make room for the title
        plt.subplots_adjust(top=0.93)

        outerFolder_expt = os.path.join(self.outerFolder, 'bias_spec')
        self.experiment.create_folder_if_not_exists(outerFolder_expt)
        now = datetime.datetime.now()
        formatted_datetime = now.strftime("%Y-%m-%d_%H-%M-%S")
        file_name = os.path.join(outerFolder_expt, f"{formatted_datetime}_BiasSpec_Q{self.QubitIndex+1}_sweeps.png")
        fig.savefig(file_name, dpi=300, bbox_inches='tight')
        plt.close(fig)
        return

    def plot2d(self, vsweep, I_arr, Q_arr, amps_arr, freq_arr):
        plt.imshow(amps_arr, aspect='auto', origin='lower',
                   extent=[freq_arr[0], freq_arr[-1], vsweep[0], vsweep[-1]])
        plt.colorbar(label="Amplitude (a.u.)")
        plt.xlabel("Qubit Frequency (MHz)")
        plt.ylabel("Voltage Bias (V)")
        plt.title(f"Bias Spectrocopy for Q{self.QubitIndex+1}")

        # Save the plot
        outerFolder_expt = os.path.join(self.outerFolder, 'bias_spec')
        self.experiment.create_folder_if_not_exists(outerFolder_expt)
        now = datetime.datetime.now()
        formatted_datetime = now.strftime("%Y-%m-%d_%H-%M-%S")
        file_name = os.path.join(outerFolder_expt, f"{formatted_datetime}_BiasSpec_Q{self.QubitIndex+1}_2d.png")
        plt.savefig(file_name, dpi=300, bbox_inches='tight')
        plt.close()
        return

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
                       mask=[0, 1, 2, 3],
                       )

        self.declare_gen(ch=qubit_ch, nqz=cfg['nqz_qubit'], mixer_freq=cfg['qubit_mixer_freq'])
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

# Biasclass = bias qubit spec with init stuff then can call function