import matplotlib.pyplot as plt
from qick.asm_v2 import AveragerProgramV2
from tqdm import tqdm
from build_state import *
from expt_config import *
import copy
import datetime

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

class ResonanceSpectroscopy:
    def __init__(self, QubitIndex, outerFolder, round_num, save_figs, experiment):
        self.QubitIndex = QubitIndex
        self.outerFolder = outerFolder
        self.expt_name = "res_spec"
        self.Qubit = 'Q' + str(self.QubitIndex)
        self.round_num = round_num
        self.save_figs = save_figs
        self.experiment = experiment

        self.exp_cfg = expt_cfg[self.expt_name]
        self.q_config = all_qubit_state(self.experiment)
        self.config = {**self.q_config[self.Qubit], **self.exp_cfg}


    def run(self, soccfg, soc):
        #defaults to 5, just make it to only look at this qubit
        res_gains = self.set_res_gain_ge(self.QubitIndex)
        self.config.update([('res_gain_ge', res_gains)])

        print(f'Q {self.QubitIndex + 1} Round {self.round_num} Res Spec configuration: ', self.config)

        fpts = self.exp_cfg["start"] + self.exp_cfg["step_size"] * np.arange(self.exp_cfg["steps"])
        fcenter = self.config['res_freq_ge']
        amps = np.zeros((len(fcenter), len(fpts)))

        for index, f in enumerate(tqdm(fpts)):
            self.config["res_freq_ge"] = fcenter + f

            prog = SingleToneSpectroscopyProgram(soccfg, reps=self.exp_cfg["reps"], final_delay=0.5, cfg=self.config)

            iq_list = prog.acquire(soc, soft_avgs=self.exp_cfg["rounds"], progress=False)
            for i in range(len(self.config['res_freq_ge'])):
                amps[i][index] = np.abs(iq_list[i][:, 0] + 1j * iq_list[i][:, 1])
        amps = np.array(amps)
        res_freqs = self.plot_results(fpts, fcenter, amps) #return freqs from plotting loop

        return res_freqs

    def set_res_gain_ge(self, QUBIT_INDEX, num_qubits=6):
        """Sets the gain for the selected qubit to 1, others to 0."""
        res_gain_ge = [0] * num_qubits  # Initialize all gains to 0
        if 0 <= QUBIT_INDEX < num_qubits:  # makes sure you are within the range of options
            res_gain_ge[QUBIT_INDEX] = 1  # Set the gain for the selected qubit
        return res_gain_ge

    def plot_results(self, fpts, fcenter, amps):
        res_freqs = []
        plt.figure(figsize=(12, 8))
        plt.rcParams.update({
            'font.size': 14,
            'axes.titlesize': 18,
            'axes.labelsize': 16,
            'xtick.labelsize': 14,
            'ytick.labelsize': 14,
            'legend.fontsize': 14,
        })

        for i in range(6):
            plt.subplot(2, 3, i + 1)
            plt.plot(fpts + fcenter[i], amps[i], '-', linewidth=1.5)
            freq_r = fpts[np.argmin(amps[i])] + fcenter[i]
            res_freqs.append(freq_r)
            plt.axvline(freq_r, linestyle='--', color='orange', linewidth=1.5)
            plt.xlabel("Frequency (MHz)")
            plt.ylabel("Amplitude (a.u.)")
            plt.title(f"Resonator {i + 1} {freq_r:.3f} MHz", pad=10)
            plt.ylim(plt.ylim()[0] - 0.05 * (plt.ylim()[1] - plt.ylim()[0]), plt.ylim()[1])

        plt.suptitle(f"MUXed resonator spectroscopy {self.exp_cfg['reps']} avgs", fontsize=24, y=0.95)
        plt.tight_layout(pad=2.0)

        outerFolder_expt = self.outerFolder + "/" + self.expt_name + "/"
        self.experiment.create_folder_if_not_exists(outerFolder_expt)
        now = datetime.datetime.now()
        formatted_datetime = now.strftime("%Y-%m-%d_%H-%M-%S")
        file_name = outerFolder_expt + f"R_{self.round_num}_" + f"Q_{self.QubitIndex+1}_" + f"{formatted_datetime}_" + self.expt_name + ".png"
        if self.save_figs:
            plt.savefig(file_name, dpi=300)
        plt.close()

        res_freqs = [round(x, 3) for x in res_freqs]
        print("Resonator freqs:", res_freqs)

        return res_freqs


