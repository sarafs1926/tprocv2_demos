from qick.asm_v2 import AveragerProgramV2
import matplotlib.pyplot as plt
from build_state import *
from expt_config import *
from system_config import *

class TOFExperiment:
    def __init__(self, QubitIndex, outerFolder, experiment, round_num = 1, save_figs = True, title = False):
        # every time a class instance is created, these definitions are set
        self.expt_name = "tof"
        self.QubitIndex = QubitIndex
        self.outerFolder = outerFolder
        self.Qubit = 'Q' + str(QubitIndex)
        self.exp_cfg = expt_cfg[self.expt_name]
        self.experiment = experiment
        self.save_figs = save_figs
        self.title = title

        self.q_config = all_qubit_state(self.experiment)
        self.round_num = round_num
        if 'All' in self.QubitIndex:
            self.config = {**self.q_config['Q0'], **self.exp_cfg}
            print(f'Q {self.QubitIndex} Round {round_num} TOF configuration: ', self.config)
        else:
            self.config = {**self.q_config[self.Qubit], **self.exp_cfg}
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
                self.trigger(ros=cfg['ro_ch'], pins=[0], t=0, ddr4=True)
                self.pulse(ch=cfg['res_ch'], name="mymux", t=0)

        prog = MuxProgram(soccfg, reps=1, final_delay=0.5, cfg=self.config)
        iq_list = prog.acquire_decimated(soc, soft_avgs=self.config['soft_avgs'])
        if self.save_figs:
            (average_y_mag_values_last, average_y_mag_values_mid, average_y_mag_values_oct, DAC_attenuator1, DAC_attenuator2, ADC_attenuator) = self.plot_results(prog, iq_list)
        else:
            (average_y_mag_values_last, average_y_mag_values_mid, average_y_mag_values_oct, DAC_attenuator1, DAC_attenuator2, ADC_attenuator) = None, None, None, None, None, None,

        return (average_y_mag_values_last, average_y_mag_values_mid, average_y_mag_values_oct, DAC_attenuator1, DAC_attenuator2, ADC_attenuator)


    def plot_results(self, prog, iq_list):
        t = prog.get_time_axis(ro_index=0)
        fig, axes = plt.subplots(len(self.config['ro_ch']), 1, figsize=(12, 12))
        phase_offsets=[]
        average_y_mag_values_mid = []
        average_y_I_values_mid = []
        average_y_Q_values_mid = []
        average_y_mag_values_oct = []
        average_y_I_values_oct = []
        average_y_Q_values_oct = []
        average_y_mag_values_last = []
        average_y_I_values_last = []
        average_y_Q_values_last = []
        for i, ch in enumerate(self.config['ro_ch']):
            plot = axes[i]
            plot.plot(t, iq_list[i][:, 0], label="I value")
            plot.plot(t, iq_list[i][:, 1], label="Q value")
            magnitude = np.abs(iq_list[i].dot([1, 1j]))
            plot.plot(t, magnitude, label="magnitude")
            plot.legend()
            plot.set_ylabel("a.u.")
            plot.set_xlabel("us")
            plot.axvline(0.75, c='r')

            phase_offset = np.angle(iq_list[i].dot([1, 1j]).sum(), deg=True)
            # print("measured phase %f degrees" % (phase_offset))
            phase_offsets.append(phase_offset)


            # Find indices of the middle three x-values
            mid_index = len(t) // 2
            indices_mid = [mid_index - 15, mid_index, mid_index + 15] #average 7 values

            one_eighth_index = len(t) // 15
            indices_oct = [one_eighth_index - 1, one_eighth_index, one_eighth_index + 1]

            indices_last = slice(-15, None)  # this will grab the last 7 elements

            # Calculate average y-values for I, Q, and magnitude
            avg_i_mid = np.mean(iq_list[i][indices_mid, 0])
            avg_q_mid = np.mean(iq_list[i][indices_mid, 1])
            avg_mag_mid = np.mean(magnitude[indices_mid])

            # Calculate average y-values for I, Q, and magnitude
            avg_i_oct = np.mean(iq_list[i][indices_oct, 0])
            avg_q_oct = np.mean(iq_list[i][indices_oct, 1])
            avg_mag_oct = np.mean(magnitude[indices_oct])

            # Calculate average y-values for I, Q, and magnitude
            avg_i_last = np.mean(iq_list[i][indices_last, 0])
            avg_q_last = np.mean(iq_list[i][indices_last, 1])
            avg_mag_last = np.mean(magnitude[indices_last])

            # Append the average magnitude to the list, you can change this to average I or Q.
            average_y_mag_values_mid.append(avg_mag_mid)
            average_y_I_values_mid.append(avg_i_mid)
            average_y_Q_values_mid.append(avg_q_mid)

            average_y_mag_values_oct.append(avg_mag_oct)
            average_y_I_values_oct.append(avg_i_oct)
            average_y_Q_values_oct.append(avg_q_oct)

            average_y_mag_values_last.append(avg_mag_last)
            average_y_I_values_last.append(avg_i_last)
            average_y_Q_values_last.append(avg_q_last)
        if self.title:
            plt.suptitle(f"TOF DAC_Att_1:{self.experiment.DAC_attenuator1} DAC_Att_2:{self.experiment.DAC_attenuator2} ADC_Att:{self.experiment.ADC_attenuator}", fontsize=24, y=0.95)


        # Save
        if self.save_figs:
            outerFolder_expt = self.outerFolder + "/" + self.expt_name + "/"
            self.experiment.create_folder_if_not_exists(outerFolder_expt)
            now = datetime.datetime.now()
            formatted_datetime = now.strftime("%Y-%m-%d_%H-%M-%S")
            if 'All' in self.QubitIndex:
                file_name = outerFolder_expt + f"R_{self.round_num}" + f"Q_{self.QubitIndex}" + f"{formatted_datetime}_" + self.expt_name + ".png"
            else:
                file_name = outerFolder_expt + f"R_{self.round_num}" + f"Q_{self.QubitIndex+1}" + f"{formatted_datetime}_" + self.expt_name + ".png"
            plt.savefig(file_name, dpi=50)
            plt.close(fig)

        return average_y_mag_values_last, average_y_mag_values_mid, average_y_mag_values_oct, self.experiment.DAC_attenuator1, self.experiment.DAC_attenuator2, self.experiment.ADC_attenuator



