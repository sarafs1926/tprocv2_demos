# scan T1 over time (or as a function of freq) for multiple qubits using mux ro

import os
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
import time

# Used for live plotting, need to run "python -m visdom.server" in the terminal and open the IP address in browser
import visdom

# from build_task import *
# from build_state import *
from qick_setup_funcs import *

from qualang_tools.plot import Fit
import pprint as pp

from tqdm import tqdm

# connect to the rfsoc and print soccfg
from rfsoc_connect import *

if SS == 'True':
    # import single shot information for g-e calibration
    from SingleShot import SingleShotProgram_g, SingleShotProgram_e, SingleShotProgram_f
    from SingleShot import config as ss_config


from malab import *
from malab.datamanagement import SlabFile
from numpy import *
import os
import datetime
import os.path
from malab.instruments.PNAX import N5242A
from malab.instruments.voltsource import *
from malab.dsfit import fithanger_new_withQc

##################
# expt Funcitons #
##################

def T1():
    '''
    Perform t1 measurement, save data, and return fitted t1 time.
    '''
    start_time = time.time()
    # ----- Experiment configurations ----- #
    expt_name = "T1_mux_ge"
    with open('system_config.json', 'r') as f:
        configs = json.load(f)
    system_config = configs['config']
    expt_cfg = configs['expt_cfg']
    readout_cfg = system_config['readout_cfg']
    hw_cfg = add_qubit_channel(system_config, QUBIT_INDEX)
    qubit_cfg = system_config['qubit_cfg']
    # add the parameters for the specific experiment we want to do
    exp_cfg = add_qubit_experiment(expt_cfg, 'T1_ge', QUBIT_INDEX) # qubit index doesnt matter here

    config = {**hw_cfg, **readout_cfg, **qubit_cfg, **exp_cfg}

    print(expt_name + '\n')
    print(config)

    ##################
    # Define Program #
    ##################
    class T1Program(AveragerProgramV2):
        def _initialize(self, cfg):
            # initialize the readout channels and pulses for DAC and ADC
            initialize_ro_chs(self, cfg, suffix='_ge')

            self.add_loop("waitloop", cfg["steps"])
            
            # initialize the qubit generator channel
            declare_gen_ch(self, cfg, cfg['qubit_ch'], usage='qubit', suffix='_ge')
            
            # initialize qubit pulse
            suffix = '_ge'
            usage = 'qubit'
            ch = cfg['qubit_ch']
            for qubit in range(num_qubits):
                self.add_gauss(ch=ch, name="ramp" + suffix, sigma=cfg['sigma' + suffix][qubit], length=cfg['sigma' + suffix][qubit]*5, even_length=True)
                self.add_pulse(ch=ch, name='qubit_pulse' + str(qubit), ro_ch=cfg['mux_ro_chs'][qubit], 
                        style="arb", 
                        envelope="ramp" + suffix, 
                        freq=cfg[usage + '_freq' + suffix][qubit], 
                        phase=cfg[usage + '_phase'],
                        gain=cfg[usage + '_gain' + suffix][qubit], 
                        )
                # declare_pulse(self, cfg, cfg['qubit_ch'], usage='qubit', 
                #             pulse_style='gauss', pulse_name='qubit_pulse' + str(qubit), suffix='_ge')
            
        def _body(self, cfg):
            # for qubit in range(num_qubits):
            #     self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse" + str(qubit))  #play pulse
            #     self.delay_auto(0.01)

            self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse" + str(3))
            self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse" + str(7))
            self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse" + str(5))
            self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse" + str(2))
            self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse" + str(4))
            self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse" + str(0))
            self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse" + str(1))
            self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse" + str(6))
            
            self.delay_auto(cfg['wait_time']+0.01, tag='wait') # wait_time after last pulse

            # readout using proper configurations for MUX and non-MUX
            readout(self, cfg)

    ###################
    # Run the Program
    ###################

    t1=T1Program(soccfg, reps=config['reps'], final_delay=config['relax_delay'], cfg=config)
    py_avg = config['py_avg']

    # for live plotting
    IS_VISDOM = True
    if IS_VISDOM:
        expt_I = expt_Q = expt_mags = expt_phases = expt_pop = None
        expt_mags = []
        viz = visdom.Visdom()
        assert viz.check_connection(timeout_seconds=5), "Visdom server not connected!"
        viz.close(win=None) # close previous plots
        win1 = viz.line( X=np.arange(0, 1), Y=np.arange(0, 1),
            opts=dict(height=400, width=700, title='T1 Experiment', showlegend=True, xlabel='expt_pts'))

        for ii in range(py_avg):
            iq_list = t1.acquire(soc, soft_avgs=1, progress=False)
            delay_times = t1.get_time_param('wait', "t", as_array=True)
            delay_times -= delay_times[0]

            this_I = []
            this_Q = []
            for qubit in range(num_qubits):
                this_I.append([iq_list[qubit].T[0]])
                this_Q.append([iq_list[qubit].T[1]])
            
            this_I = np.array(this_I)
            this_Q = np.array(this_Q)

            if expt_I is None: # ii == 0
                expt_I, expt_Q = this_I, this_Q
            else:
                expt_I = (expt_I * ii + this_I) / (ii + 1.0)
                expt_Q = (expt_Q * ii + this_Q) / (ii + 1.0)

            expt_mags = []
            for qubit in range(num_qubits):
                expt_mags.append(np.abs(expt_I[qubit] + 1j * expt_Q[qubit])[0])  # magnitude
            expt_mags = np.array(expt_mags)
            # expt_mags = np.abs(expt_I + 1j * expt_Q)  # magnitude
            expt_phases = np.angle(expt_I + 1j * expt_Q)  # phase
            # print(expt_mags.shape)
            for qubit in range(num_qubits):
                if qubit == 0:
                    viz.line(X = delay_times, Y = expt_mags[qubit], win=win1, name='Qubit' + str(qubit),
                        opts=dict(height=400, width=700, title='T1 Experiment', showlegend=True, xlabel='expt_pts'))
                else:
                    viz.line(X = delay_times, Y = expt_mags[qubit], win=win1, name='Qubit' + str(qubit), update='append',
                        opts=dict(height=400, width=700, title='T1 Experiment', showlegend=True, xlabel='expt_pts'))

        amps = expt_mags
        
    else:
        iq_list = t1.acquire(soc, soft_avgs=py_avg, progress=True)
        delay_times = t1.get_time_param('wait', "t", as_array=True)
        amps = np.abs(iq_list[0][0].dot([1,1j]))

    t1_array = []
    for qubit in range(num_qubits):
        # ### Fit ###
        fit = Fit()
        # Choose the suitable fitting function
        fit_result = fit.T1(delay_times, amps[qubit].T[0])
        fit_result = {
                "T1": fit_result['T1'],
                "amp": fit_result['amp'],
                "final_offset": fit_result['final_offset']
            }
        # pp.pprint(fit_result)
        t1_array.append(fit_result['T1'][0])

    # pp.pprint(fit_result)

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
    exp_name = expt_name + '_' + prefix
    print('Experiment name: ' + exp_name)

    data_path = DATA_PATH

    fname = get_next_filename(data_path, exp_name, suffix='.h5')
    print('Current data file: ' + fname)
    end_time = time.time()
    tot_time = end_time - start_time
    with SlabFile(data_path + '\\' + fname, 'a') as f:
        # 'a': read/write/create
        f.append('expt_time', [int(tot_time)])
        # - Adds data to the file - #
        f.append('delay_times', delay_times)
        f.append('amps', amps)
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

    return t1_array

#################
# dac functions #
#################

### YOKO setup
YOKO_4 = ""
YOKO_2 = ""
YOKO_5 = ""

dcflux2 = YokogawaGS200(address=YOKO_5) # dc coil
dcflux2.recv_length = 1024
print(dcflux2.recv_length)
# print(dcflux2.get_id())
# print(dcflux2.get_id())
dcflux2.set_mode('current')
dcflux2.set_range(0.01)  # 10mA range
dcflux2.set_output(True)

# initial NWA configuration values
num_qubits = int(8)

### DACs setup
'''
labelled by flux lines
C5: Q0, C6: Q1, C7: Q2, C8: Q3
D1: Q4, D2: Q5, D3: Q6, D4: D4
'''
dac_C5toC8 = AD5780(address='')
dac_D1toD4 = AD5780(address='')
# dacs[qubit index][dac address or channel number]
dacs = [[dac_C5toC8, 1], [dac_C5toC8, 2],[dac_C5toC8, 3],[dac_C5toC8, 4],
        [dac_D1toD4, 1], [dac_D1toD4, 2],[dac_D1toD4, 3],[dac_D1toD4, 4]]
time.sleep(1)

#V/s or mA/s
dac_rate = 0.05

# converts DAC units to current value
def digit_to_curr(x):
    return 20*x/(2**18)-10

# initial offset of DAC channels, measured with volt meter
dac_offset = array([-0.002, -0.004, -0.003, -0.003, -0.004, -0.009, -0.02, -0.024])
diag_offset = dac_offset # set so background is zero, not needed here

def drive_YOKO(pt):

    print("Driving YOKO at (%.3f) mA" % (pt[0]))
    dcflux2.ramp_current(pt[0]*1e-3, 0.0005) # second argument is rate?
    time.sleep(0.2)

def drive_DACs(pt):

    # I'm leaving this hardcoded for 8 qubits for now - Hebah 2024 12 16
    print("Driving DAC at ( %.3f, %.3f, %.3f, %.3f,%.3f, %.3f, %.3f, %.3f ) mA" %(pt[1], pt[2], pt[3], pt[4], pt[5], pt[6], pt[7], pt[8]))

    for ii in range(len(dacs)):
        print('\n')
        print('Flux line for qubit ', ii)
        print(dacs[ii][0])
        print("Driving DAC at (%.3f) mA" % (pt[ii+1]-dac_offset[ii]))
        dacs[ii][0].get_voltage(dacs[ii][1])
        dacs[ii][0].get_voltage(dacs[ii][1])
        dacs[ii][0].get_voltage(dacs[ii][1])
        dacs[ii][0].get_voltage(dacs[ii][1])
        dac_old = digit_to_curr(int(dacs[ii][0].get_voltage(dacs[ii][1])[:-2]))
        print('old reading = ', dac_old)
        dacs[ii][0].ramp(dacs[ii][1], pt[ii+1], dac_rate)
        time.sleep(1.5 * abs(dac_old - pt[ii+1]) / dac_rate + 5)
        time.sleep(5.0)
        dacs[ii][0].get_voltage(dacs[ii][1])
        dacs[ii][0].get_voltage(dacs[ii][1])
        dacs[ii][0].get_voltage(dacs[ii][1])
        dacs[ii][0].get_voltage(dacs[ii][1])  # to clear the readout queue

def drive_DAC(pt, index):
    ii = index
    print('Driving DAC ' + str(index) + ' at ( %.3f ) mA' %(pt[ii+1]))

    
    print('\n')
    print('Flux line for qubit ', ii)
    print(dacs[ii][0])
    print("Driving DAC at (%.3f) mA" % (pt[ii+1]-dac_offset[ii]))
    dacs[ii][0].get_voltage(dacs[ii][1])
    dacs[ii][0].get_voltage(dacs[ii][1])
    dacs[ii][0].get_voltage(dacs[ii][1])
    dacs[ii][0].get_voltage(dacs[ii][1])
    dac_old = digit_to_curr(int(dacs[ii][0].get_voltage(dacs[ii][1])[:-2]))
    print('old reading = ', dac_old)
    dacs[ii][0].ramp(dacs[ii][1], pt[ii+1], dac_rate)
    time.sleep(1.5 * abs(dac_old - pt[ii+1]) / dac_rate + 5)
    time.sleep(5.0)
    dacs[ii][0].get_voltage(dacs[ii][1])
    dacs[ii][0].get_voltage(dacs[ii][1])
    dacs[ii][0].get_voltage(dacs[ii][1])
    dacs[ii][0].get_voltage(dacs[ii][1])  # to clear the readout queue

def drive_Yoko_and_DACs(pt):
    drive_YOKO(pt)
    drive_DACs(pt)


t1_array = []

# pt = [0]*9
# drive_Yoko_and_DACs(pt)
iterations = 1
for i in range(iterations):
    start_time = time.time()
    
    expt_name = 'T1_ge_mux'

    # run t1 mux experiment
    print('Starting T1')
    t1_value = T1()
    t1_array.append([t1_value]) # save t1 value
    print('all qubit T1 =', t1_value, 'us')

    end_time = time.time()

    print('Time taken for iteration =', end_time - start_time)


# data_path = "M:/malab/_Data/20250210 - Santi - RFSoC tprocv2 - LL8qubit meas/20250211 - t1 data"
# lookup_file = 'lookup_file'
# with SlabFile(data_path + '\\' + lookup_file, 'a') as f:
#     # 'a': read/write/create

#     # - Adds data to the file - #
#     # f.append('res0_freqs', res_freq_array)
#     # f.append('qubit0_freqs', qubit_freq_array)
#     f.append('t1_times', t1_array)
#     # f.append('flux_loc', pt)

# # drive_Yoko_and_DACs([0, -0.94526638, -0.12800869, 0.04248557, 0.21889241, 0.2840921,
# # 0.44851674, 0.55921558, 0.87927871]) # zero all the dacs