import matplotlib.pyplot as plt
from qick.asm_v2 import AveragerProgramV2
from tqdm import tqdm
from build_state import *
from expt_config import *
from system_config import *


class ResonanceSpectroscopy:
    def __init__(self, QubitIndex, outerFolder):
        self.QubitIndex = QubitIndex
        self.outerFolder = outerFolder
        self.expt_name = "res_spec"
        self.Qubit = 'Q' + str(self.QubitIndex)

    def _create_program(self, config):
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
                               length=config["res_length"],
                               mask=[0, 1, 2, 3, 4, 5],
                               )

            def _body(self, cfg):
                self.trigger(ros=cfg['ro_ch'], pins=[0], t=cfg['trig_time'], ddr4=True)
                self.pulse(ch=cfg['res_ch'], name="mymux", t=0)

        return SingleToneSpectroscopyProgram

    def run(self, socfg, soc, n=1):
        for j in range(n):
            exp_cfg = expt_cfg[self.expt_name]
            q_config = all_qubit_state(system_config)
            config = {**q_config['Q' + str(self.QubitIndex)], **exp_cfg}
            print(config)

            fpts = exp_cfg["start"] + exp_cfg["step_size"] * np.arange(exp_cfg["steps"])
            fcenter = config['res_freq_ge']
            amps = np.zeros((len(fcenter), len(fpts)))

            SingleToneSpectroscopyProgram = self._create_program(config)

            for index, f in enumerate(tqdm(fpts)):
                config["res_freq_ge"] = fcenter + f
                prog = SingleToneSpectroscopyProgram(soccfg, reps=exp_cfg["reps"], final_delay=0.5, cfg=config)
                iq_list = prog.acquire(soc, soft_avgs=exp_cfg["rounds"], progress=False)
                for i in range(len(config['res_freq_ge'])):
                    amps[i][index] = np.abs(iq_list[i][:, 0] + 1j * iq_list[i][:, 1])
            amps = np.array(amps)

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

            plt.suptitle(f"MUXed resonator spectroscopy {exp_cfg['reps']} avgs", fontsize=24, y=0.95)
            plt.tight_layout(pad=2.0)

            outerFolder_expt = self.outerFolder + "/" + self.expt_name + "/"
            create_folder_if_not_exists(outerFolder_expt)
            now = datetime.datetime.now()
            formatted_datetime = now.strftime("%Y-%m-%d_%H-%M-%S")
            file_name = outerFolder_expt + f"{formatted_datetime}_" + self.expt_name + ".png"
            plt.savefig(file_name, dpi=300)
            plt.close()

            res_gains = config["res_gain_ge"]
            res_freqs = [round(x, 3) for x in res_freqs]
            print("Resonator gains:", res_gains)
            print("Resonator freqs:", res_freqs)

