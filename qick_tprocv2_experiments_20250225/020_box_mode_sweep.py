# perform box mode sweep for T2 Ramsey experiment

# reomote control digital attenuator
import Pyro5.api
signal_core_01 = Pyro5.api.Proxy("PYRONAME:signal_core_zcu216")    # use name server object lookup uri shortcut

import os

# from malab.instruments import rfSignalCore
# ## parameter for probe LO
# address_red = "1000299D"
folder = os.getcwd()
os.chdir(folder + '/qick_tprocv2_experiments')

from qick import *
from qick.pyro import make_proxy

from malab import SlabFile, get_next_filename
# for now, all the tProc v2 classes need to be individually imported (can't use qick.*)

# the main program class
from qick.asm_v2 import AveragerProgramV2
# for defining sweeps
from qick.asm_v2 import QickSpan, QickSweep1D

import json
import datetime

# Used for live plotting, need to run "python -m visdom.server" in the terminal and open the IP address in browser
import visdom

from build_task import *
from build_state import *
from qick_setup_funcs import *

from qualang_tools.plot import Fit
import pprint as pp

# connect to the rfsoc and print soccfg
from rfsoc_connect import *

if SS == 'True':
    # import single shot information for g-e calibration
    from SingleShot import SingleShotProgram_g, SingleShotProgram_e
    from SingleShot import config as ss_config

##################
# Define Program #
##################

class RamseyProgram(AveragerProgramV2):
    def _initialize(self, cfg):
        # initialize the readout channels and pulses for DAC and ADC
        initialize_ro_chs(self, cfg, suffix='_ge')

        self.add_loop("waitloop", cfg["steps"])
        
       # initialize the qubit generator channel
        declare_gen_ch(self, cfg, cfg['qubit_ch'], usage='qubit', suffix='_ge')
        # initialize qubit pulse
        declare_pulse(self, cfg, cfg['qubit_ch'], usage='qubit', 
                      pulse_style='gauss', pulse_name='ramsey1', suffix='_ge')
        
        # initialize qubit pulse
        declare_pulse(self, cfg, cfg['qubit_ch'], usage='qubit', 
                      pulse_style='gauss', pulse_name='ramsey2', suffix='_ge')
    
    def _body(self, cfg):
        self.pulse(ch=self.cfg["qubit_ch"], name="ramsey1", t=0)  #play probe pulse
        
        self.delay_auto(cfg['wait_time']+0.01, tag='wait') # wait_time after last pulse
        
        self.pulse(ch=self.cfg["qubit_ch"], name="ramsey2", t=0)  #play pulse

        self.delay_auto(0.01) # wait_time after last pulse
        
        # readout using proper configurations for MUX and non-MUX
        readout(self, cfg)

freq_list = np.linspace(0,10,400) * 1e9
power_LO_list = [5]

tot_start_time = time.time()
for index_power, power_LO in enumerate(power_LO_list):
    for index, LO_freq in enumerate(freq_list):
        start_time = time.time()
        LO_freq_int = int(LO_freq)
        # print(power_LO)
        ## update the hardware parameters
        # print('starting signalcores for red sideband: ' + address_red)
        # sc1 = rfSignalCore.SignalCore(name="SignalCore", address=address_red)
        # sc1.set_power(power_LO)
        # sc1.set_frequency(LO_freq_int)
        # sc1.close_device()

        # os.chdir('C:/Users/G41Lab/Dropbox/People/Santi/code/pyro_test/basic_test')
        # read_power = signal_core_01.set_power(power_LO)
        # read_freq = signal_core_01.set_freq(LO_freq_int)
        while True:
            try:
                read_power = signal_core_01.set_power(power_LO)
                read_freq = signal_core_01.set_freq(LO_freq_int)
                break
            except:
                print('Error: signal core not found')

        print('Reading freq =', read_freq)
        print('Reading power =', read_power)
        # os.chdir(folder + '/qick_tprocv2_experiments')

        # ----- Experiment configurations ----- #
        expt_name = "Ramsey_ge"
        config, expt_name = initialize_configs(expt_name)
        print(expt_name + '\n')
        print(config)
        ramsey=RamseyProgram(soccfg, reps=config['reps'], final_delay=config['relax_delay'], cfg=config)
        py_avg = config['py_avg']

        # for live plotting
        IS_VISDOM = True
        if IS_VISDOM:
            expt_I = expt_Q = expt_mags = expt_phases = expt_pop = None
            viz = visdom.Visdom()
            assert viz.check_connection(timeout_seconds=5), "Visdom server not connected!"
            viz.close(win=None) # close previous plots
            win1 = viz.line( X=np.arange(0, 1), Y=np.arange(0, 1),
                opts=dict(height=400, width=700, title='T2 Ramsey Experiment', showlegend=True, xlabel='expt_pts'))

            for ii in range(py_avg):
                iq_list = ramsey.acquire(soc, soft_avgs=1, progress=False)
                delay_times = ramsey.get_time_param('wait', "t", as_array=True)
                delay_times -= delay_times[0]
                iq_list = iq_list[0][0].T
                # what is the correct shape/index?
                this_I = (iq_list[0])
                this_Q = (iq_list[1])

                if expt_I is None: # ii == 0
                    expt_I, expt_Q = this_I, this_Q
                else:
                    expt_I = (expt_I * ii + this_I) / (ii + 1.0)
                    expt_Q = (expt_Q * ii + this_Q) / (ii + 1.0)

                expt_mags = np.abs(expt_I + 1j * expt_Q)  # magnitude
                expt_phases = np.angle(expt_I + 1j * expt_Q)  # phase

                viz.line(X = delay_times, Y = expt_mags, win=win1, name='I',
                opts=dict(height=400, width=700, title='T2 Ramsey Experiment', showlegend=True, xlabel='expt_pts'))

            amps = expt_mags

        if SS == 'True':
            print('performing single shot for g-e calibration')

            ssp_g = SingleShotProgram_g(soccfg, reps=1, final_delay=ss_config['relax_delay'], cfg=ss_config)
            iq_list_g = ssp_g.acquire(soc, soft_avgs=1, progress=True)

            ssp_e = SingleShotProgram_e(soccfg, reps=1, final_delay=ss_config['relax_delay'], cfg=ss_config)
            iq_list_e = ssp_e.acquire(soc, soft_avgs=1, progress=True)

            I_g = iq_list_g[0][0].T[0]
            Q_g = iq_list_g[0][0].T[1]
            I_e = iq_list_e[0][0].T[0]
            Q_e = iq_list_e[0][0].T[1]
            
        #####################################
        # ----- Saves data to a file ----- #
        #####################################

        prefix = str(datetime.date.today())
        exp_name = expt_name + '_Q' + str(QUBIT_INDEX) + '_' + prefix
        print('Experiment name: ' + exp_name)

        data_path = DATA_PATH

        fname = get_next_filename(data_path, exp_name, suffix='.h5')
        print('Current data file: ' + fname)

        with SlabFile(data_path + '\\' + fname, 'a') as f:
            # 'a': read/write/create

            # - Adds data to the file - #
            f.append('delay_times', delay_times)
            f.append('amps', amps)
            f.append('LO_freq', [LO_freq])
            f.append('power_LO', [power_LO])
            if IS_VISDOM:
                f.append('avgi', expt_I)
                f.append('avgq', expt_Q)
            else:
                f.append('avgi', iq_list[0][0].T[0])
                f.append('avgq', iq_list[0][0].T[1])

            del config['wait_time']  # cant save QickParam
            # formats config into file as a single line
            f.attrs['config'] = json.dumps(config)
            # f.attrs['fit_result'] = json.dumps(fit_result)

            if SS == 'True':
                # - Adds ss data to the file - #
                f.append('I_g', I_g)
                f.append('Q_g', Q_g)
                f.append('I_e', I_e)
                f.append('Q_e', Q_e)    

                # formats config into file as a single line
                f.attrs['ss_config'] = json.dumps(ss_config)

        data = data_path + '\\' + fname

        end_time = time.time()
        print('Time taken for iteration =', end_time - start_time)

tot_end_time = time.time()
print('Total time taken =', tot_end_time - tot_start_time)