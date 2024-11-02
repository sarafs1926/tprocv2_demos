from qick.asm_v2 import AveragerProgramV2
import matplotlib.pyplot as plt
from build_state import *
from expt_config import *
from system_config import *

class TOFExperiment:
    def __init__(self, QubitIndex, outerFolder, round_num, trigger_time=0):
        # every time a class instance is created, these definitions are set
        self.expt_name = "tof"
        self.QubitIndex = QubitIndex
        self.outerFolder = outerFolder
        self.Qubit = 'Q' + str(QubitIndex)

        self.exp_cfg = expt_cfg[self.expt_name]
        self.q_config = all_qubit_state(system_config)
        self.config = {**self.q_config[self.Qubit], **self.exp_cfg}
        self.round_num = round_num

        # Update parameters to see TOF pulse with your setup
        # self.config.update([('trig_time', trigger_time)])  #starting off with TOF = 0, we will change this later to be the found TOF
        # print("initial readout phase is")
        # print(self.config['ro_phase'])
        # new_phases = [np.float64(70.40572044377751), np.float64(-162.07287664422063), np.float64(-91.29642515281115), np.float64(-53.58106709833659), np.float64(13.776008961640695), np.float64(-7.341398813031154)]
        # self.config.update([('ro_phase', new_phases)])  #starting off with TOF = 0, we will change this later to be the found TOF
        # print("new readout phase is")
        print(self.config['ro_phase'])
        print(f'Q {self.QubitIndex + 1} Round {round_num} TOF configuration: ',self.config)


    def run(self, soccfg, soc):
        class MuxProgram(AveragerProgramV2):
            def _initialize(self, cfg):
                ro_chs = cfg['ro_ch']
                gen_ch = cfg['res_ch']

                self.declare_gen(
                    ch=gen_ch, nqz=cfg['nqz_res'], ro_ch=ro_chs[0],
                    mux_freqs=cfg['res_freq_ge'],
                    mux_gains=cfg['res_gain_ge'],
                    mux_phases=cfg['res_phase'],
                    mixer_freq=cfg['mixer_freq']
                )
                for ch, f, ph in zip(cfg['ro_ch'], cfg['res_freq_ge'], cfg['ro_phase']):
                    self.declare_readout(
                        ch=ch, length=cfg['res_length'], freq=f, phase=ph, gen_ch=gen_ch
                    )

                self.add_pulse(
                    ch=gen_ch, name="mymux",
                    style="const",
                    length=cfg["res_length"],
                    mask=[0, 1, 2, 3, 4, 5]
                )

            def _body(self, cfg):
                self.trigger(ros=cfg['ro_ch'], pins=[0], t=cfg['trig_time'], ddr4=True)
                self.pulse(ch=cfg['res_ch'], name="mymux", t=0)


        prog = MuxProgram(soccfg, reps=1, final_delay=0.5, cfg=self.config)
        iq_list = prog.acquire_decimated(soc, soft_avgs=self.config['soft_avgs'])
        self.plot_results(prog, iq_list)

    def plot_results(self, prog, iq_list):
        t = prog.get_time_axis(ro_index=0)
        fig, axes = plt.subplots(len(self.config['ro_ch']), 1, figsize=(12, 12))
        phase_offsets=[]
        for i, ch in enumerate(self.config['ro_ch']):
            plot = axes[i]
            plot.plot(t, iq_list[i][:, 0], label="I value")
            plot.plot(t, iq_list[i][:, 1], label="Q value")
            plot.plot(t, np.abs(iq_list[i].dot([1, 1j])), label="magnitude")
            plot.legend()
            plot.set_ylabel("a.u.")
            plot.set_xlabel("us")
            plot.axvline(0.75, c='r')

            phase_offset = np.angle(iq_list[i].dot([1, 1j]).sum(), deg=True)
            # print("measured phase %f degrees" % (phase_offset))
            phase_offsets.append(phase_offset)
        #
        print(phase_offsets)
        # Save
        outerFolder_expt = self.outerFolder + "/" + self.expt_name + "/"
        create_folder_if_not_exists(outerFolder_expt)
        now = datetime.datetime.now()
        formatted_datetime = now.strftime("%Y-%m-%d_%H-%M-%S")
        file_name = outerFolder_expt + f"R_{self.round_num}" + f"Q_{self.QubitIndex+1}" + f"{formatted_datetime}_" + self.expt_name + ".png"
        plt.savefig(file_name, dpi=300)
        plt.close(fig)

