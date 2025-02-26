# perform collective decay experiment between two qubits

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

def ResonatorSpectrosocpy():
    '''
    Perform resonator spectroscopy, save data, and return resonator frequency found by looking for minimum amplitude.
    '''

    start_time = time.time()

    # ----- Experiment configurations ----- #
    expt_name = "res_spec_ge"
    config, expt_name = initialize_configs(expt_name)
    print(expt_name + '\n')
    print(config)

    ##################
    # Define Program #
    ##################

    class SingleToneSpectroscopyProgram(AveragerProgramV2):
        def _initialize(self, cfg):
            if MUX == 'False': # non-MUX case
                ro_ch = cfg['ro_ch']
                res_ch = cfg['res_ch']
                
                # declare the proper readout channels using
                declare_gen_ch(self, cfg, res_ch, usage = 'res', suffix = '_ge')
                self.declare_readout(ch=ro_ch, length=cfg['ro_length'])

                # add a loop for frequency sweep of readout pulse frequency
                self.add_loop("freqloop", cfg["steps"])
                self.add_readoutconfig(ch=ro_ch, name="myro", freq=cfg['res_freq' + '_ge'], gen_ch=res_ch)

                # declare a pulse for the non-MUX channel
                declare_pulse(self, cfg, res_ch, pulse_name='res_pulse', suffix='_ge')
                
            else: # MUX case
                # initialize the readout channels and pulses for DAC and ADC
                initialize_ro_chs(self, cfg, suffix='_ge')

        def _body(self, cfg):
            # readout using proper configurations for MUX and non-MUX
            readout(self, cfg)


    if MUX == 'False': # non-MUX program
        ###################
        # Run the Program
        ###################

        prog = SingleToneSpectroscopyProgram(soccfg, reps=config['reps'], final_delay=config['relax_delay'], cfg=config)
        py_avg = config['py_avg']

        # for live plotting
        IS_VISDOM = True
        if IS_VISDOM:
            expt_I = expt_Q = expt_mags = expt_phases = expt_pop = None
            viz = visdom.Visdom()
            assert viz.check_connection(timeout_seconds=5), "Visdom server not connected!"
            viz.close(win=None) # close previous plots
            win1 = viz.line( X=np.arange(0, 1), Y=np.arange(0, 1),
                            opts=dict(height=400, width=700, title='Resonator Spectroscopy', showlegend=True, xlabel='expt_pts'))

            for ii in range(py_avg):
                iq_list = prog.acquire(soc, soft_avgs=1, progress=False, remove_offset = False)
                freqs = prog.get_pulse_param("res_pulse", "freq", as_array=True)

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

                viz.line(X = freqs, Y = expt_mags, win=win1, name='I',
                        opts=dict(height=400, width=700, title='Resonator Spectroscopy', showlegend=True, xlabel='expt_pts'))

            amps = np.abs(expt_I + 1j*expt_Q)
            
        else:
            iq_list = prog.acquire(soc, soft_avgs = py_avg, progress=True)
            freqs = prog.get_pulse_param("res_pulse", "freq", as_array=True)
            amps = np.abs(iq_list[0][0].dot([1,1j]))

        #####################################
        # ----- Saves data to a file ----- #
        #####################################

        prefix = str(datetime.date.today())
        exp_name = expt_name + '_Q' + str(QUBIT_INDEX) + '_' + prefix
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
            f.append('fpts', freqs)
            f.append('amps', amps)
            if IS_VISDOM:
                f.append('avgi', expt_I)
                f.append('avgq', expt_Q)
            else:
                f.append('avgi', iq_list[0][0].T[0])
                f.append('avgq', iq_list[0][0].T[1])

            del config['res_freq_ge'] # cant save QickParam
            # formats config into file as a single line
            f.attrs['config'] = json.dumps(config)
            # f.attrs['fit_result'] = json.dumps(fit_result)

        data = data_path + '\\' + fname

    else: # MUX program
        ###################
        # Run the Program
        ###################

        fpts=[config["start"] + ii*config["step"] for ii in range(config["expts"])]
        fcenter = np.array(config['freqs'])

        avgi = np.zeros((len(fcenter), len(fpts)))
        avgq = np.zeros((len(fcenter), len(fpts)))
        amps = np.zeros((len(fcenter), len(fpts)))
        for index, f in enumerate(tqdm(fpts)):
            config["res_freq_ge"] = fcenter + f
            prog = SingleToneSpectroscopyProgram(soccfg, reps=config["reps"], final_delay=0.5, cfg=config)    
            iq_list = prog.acquire(soc, soft_avgs = config["py_avg"], progress=False)
            for i in range(len(fcenter)):
                avgi[i][index] = iq_list[i][:,0]
                avgq[i][index] = iq_list[i][:,1]
                amps[i][index] = np.abs(iq_list[i][:,0]+1j*iq_list[i][:,1])
        amps = np.array(amps)
        avgi = np.array(avgi)
        avgq = np.array(avgq)

        #####################################
        # ----- Saves data to a file ----- #
        #####################################

        prefix = str(datetime.date.today())
        exp_name = expt_name + '_' + prefix
        print('Experiment name: ' + exp_name)

        data_path = DATA_PATH

        fname = get_next_filename(data_path, exp_name, suffix='.h5')
        print('Current data file: ' + fname)

        with SlabFile(data_path + '\\' + fname, 'a') as f:
            # 'a': read/write/create

            # - Adds data to the file - #
            f.append('fpts', fpts)
            f.append('amps', amps)
            f.append('fcenter', fcenter)
            f.append('avgi', avgi)
            f.append('avgq', avgq)

            del config['res_freq_ge'] # dont need to save this...
            # formats config into file as a single line
            f.attrs['config'] = json.dumps(config)

        data = data_path + '\\' + fname
    resonator_freq = freqs[np.argmin(amps)]
    return resonator_freq

def QubitSpectrosocpy():
    '''
    Perform qubit spectroscopy, save data, and return qubit frequency found by looking for maximum amplitude.
    '''
    start_time = time.time()
    # ----- Experiment configurations ----- #
    expt_name = "qubit_spec_ge"
    config, expt_name = initialize_configs(expt_name)
    print(expt_name + '\n')
    print(config)

    ##################
    # Define Program #
    ##################

    class PulseProbeSpectroscopyProgram(AveragerProgramV2):
        def _initialize(self, cfg):
            # initialize the readout channels and pulses for DAC and ADC
            initialize_ro_chs(self, cfg, suffix='_ge')

            # add a loop for frequency sweep of qubit pulse frequency
            self.add_loop("freqloop", cfg["steps"])

            # initialize the qubit generator channel
            declare_gen_ch(self, cfg, cfg['qubit_ch'], usage='qubit', suffix='_ge')
            # initialize qubit pulse
            declare_pulse(self, cfg, cfg['qubit_ch'], usage='qubit', 
                        pulse_style='const', pulse_name='qubit_pulse', suffix='_ge')

        def _body(self, cfg):
            self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse", t=0)  #play probe pulse
            self.delay_auto(0.01, tag='wait') # wait_time after last pulse
            # readout using proper configurations for MUX and non-MUX
            readout(self, cfg)

    ###################
    # Run the Program #
    ###################

    qspec=PulseProbeSpectroscopyProgram(soccfg, reps=config['reps'], final_delay=config['relax_delay'], cfg=config)
    py_avg = config['py_avg']

    # for live plotting
    IS_VISDOM = True
    if IS_VISDOM:
        expt_I = expt_Q = expt_mags = expt_phases = expt_pop = None
        viz = visdom.Visdom()
        assert viz.check_connection(timeout_seconds=5), "Visdom server not connected!"
        viz.close(win=None) # close previous plots
        win1 = viz.line( X=np.arange(0, 1), Y=np.arange(0, 1),
            opts=dict(height=400, width=700, title='Qubit Spectroscopy', showlegend=True, xlabel='expt_pts'))

        for ii in range(py_avg):
            iq_list = qspec.acquire(soc, soft_avgs=1, progress=False)
            freqs = qspec.get_pulse_param('qubit_pulse', "freq", as_array=True)

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

            viz.line(X = freqs, Y = expt_mags, win=win1, name='I',
            opts=dict(height=400, width=700, title='Qubit Spectroscopy', showlegend=True, xlabel='expt_pts'))

        amps = np.abs(expt_I + 1j*expt_Q)
        
    else:
        iq_list = qspec.acquire(soc, soft_avgs=py_avg, progress=True)
        freqs = qspec.get_pulse_param('qubit_pulse', "freq", as_array=True)
        amps = np.abs(iq_list[0][0].dot([1,1j]))

    # ### Fit ###
    # fit = Fit()
    # # Choose the suitable fitting function
    # fit_result = fit.transmission_resonator_spectroscopy(freqs, amps)

    # # data = fit_result

    # fit_result = {
    #         "f": fit_result['f'],
    #         "kc": fit_result['kc'],
    #         "ki": fit_result['ki'],
    #         "k": fit_result['k'],
    #         "offset": fit_result['offset']
    #     }

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
    exp_name = expt_name + '_Q' + str(QUBIT_INDEX) + '_' + prefix
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
        f.append('fpts', freqs)
        f.append('amps', amps)
        if IS_VISDOM:
            f.append('avgi', expt_I)
            f.append('avgq', expt_Q)
        else:
            f.append('avgi', iq_list[0][0].T[0])
            f.append('avgq', iq_list[0][0].T[1])

        del config['qubit_freq_ge']  # cant save QickParam
        # formats config into file as a single line
        f.attrs['config'] = json.dumps(config) # cant save configs yet with QickParams
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

    qubit_freq = freqs[np.argmax(amps)]
    return qubit_freq

def T1():
    '''
    Perform t1 measurement, save data, and return fitted t1 time.
    '''
    start_time = time.time()
    # ----- Experiment configurations ----- #
    expt_name = "T1_ge"
    config, expt_name = initialize_configs(expt_name)
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
            declare_pulse(self, cfg, cfg['qubit_ch'], usage='qubit', 
                        pulse_style='const', pulse_name='qubit_pulse', suffix='_ge')
        
        def _body(self, cfg):
            self.pulse(ch=self.cfg["qubit_ch"], name="qubit_pulse", t=0)  #play pulse
            
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
        viz = visdom.Visdom()
        assert viz.check_connection(timeout_seconds=5), "Visdom server not connected!"
        viz.close(win=None) # close previous plots
        win1 = viz.line( X=np.arange(0, 1), Y=np.arange(0, 1),
            opts=dict(height=400, width=700, title='T1 Experiment', showlegend=True, xlabel='expt_pts'))

        for ii in range(py_avg):
            iq_list = t1.acquire(soc, soft_avgs=1, progress=False)
            delay_times = t1.get_time_param('wait', "t", as_array=True)

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
            opts=dict(height=400, width=700, title='T1 Experiment', showlegend=True, xlabel='expt_pts'))


        amps = np.abs(expt_I + 1j*expt_Q)
        
    else:
        iq_list = t1.acquire(soc, soft_avgs=py_avg, progress=True)
        delay_times = t1.get_time_param('wait', "t", as_array=True)
        amps = np.abs(iq_list[0][0].dot([1,1j]))

    ### Fit ###
    fit = Fit()

    # Choose the suitable fitting function
    fit_result = fit.T1(delay_times, amps)

    fit_result = {
            "T1": fit_result['T1'],
            "amp": fit_result['amp'],
            "final_offset": fit_result['final_offset']
        }

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
    exp_name = expt_name + '_Q' + str(QUBIT_INDEX) + '_' + prefix
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

    return fit_result['T1'][0]

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

import q1_q2_avoided_crossing_flux_values
from q1_q2_avoided_crossing_flux_values import freqs_q1, freqs_q2, corrected_dac_values

flux_array = np.array(corrected_dac_values[33]) # in unit of mA
flux_array = [[arr.tolist() for arr in flux_array]]

# we will readout using [q0, q1] resonator in this case
res_freqs = np.array([[6.170] * len(flux_array), [6.305] * len(flux_array)]) * 1000 # MHz
res_start = res_freqs - np.array([[5] * len(res_freqs[0]), [5] * len(res_freqs[1])])
res_stop = res_freqs + np.array([[5] * len(res_freqs[0]), [5] * len(res_freqs[1])])

qubit_freqs =  np.array([[freqs_q1[33]],[freqs_q2[33]]]) * 1000 # MHz
qubit_start = qubit_freqs - np.array([[20] * len(qubit_freqs[0]), [10] * len(qubit_freqs[1])])
qubit_stop = qubit_freqs +  np.array([[20] * len(qubit_freqs[0]), [30] * len(qubit_freqs[1])])

res1_freq_array = []
qubit1_freq_array = []
res2_freq_array = []
qubit2_freq_array = []
qubit1_t1_array = []
qubit2_t1_array = []

# print(flux_array)
# print(res_freqs)
# print(qubit_freqs)

# for i, pt in enumerate(flux_array):
#     start_time = time.time()
    
#     print('With crosstalk correction!')
#     print('flux actual =', pt)
#     print('\n\n')

#     with open('system_config.json', 'r') as f:
#         config = json.load(f)

#     # adjust experiment configuration values
#     config['config']['dac'] = pt

#     with open('system_config.json', 'w') as f:
#         json.dump(config, f, indent=3)

#     # now drive the YOKO and DACs
#     # drive_Yoko_and_DACs(pt)

#     # do the experiments for qubit 1 then qubit 2
#     for qubit in range(2):
#         # qubit = 1
#         QUBIT_INDEX = qubit+1
#         with open('system_config.json', 'r') as f:
#             config = json.load(f)

#         # adjust experiment configuration values
#         config['config']['qubit_index'] = QUBIT_INDEX

#         with open('system_config.json', 'w') as f:
#             json.dump(config, f, indent=3)

#         expt_name = 'res_spec_ge'
#         with open('system_config.json', 'r') as f:
#             config = json.load(f)

#         # adjust experiment configuration values
#         config['expt_cfg'][expt_name]['start'][QUBIT_INDEX] = res_start[qubit][i]
#         config['expt_cfg'][expt_name]['stop'][QUBIT_INDEX] = res_stop[qubit][i]

#         with open('system_config.json', 'w') as f:
#             json.dump(config, f, indent=3)

#         # run resonator spectroscopy experiment
#         print('Starting Resonator Spectroscopy')
#         res_freq = ResonatorSpectrosocpy()
#         if qubit == 0:
#             res1_freq_array.append(res_freq) # save res freq value
#         else:
#             res2_freq_array.append(res_freq) # save res freq value

#         print('Resonator frequency =', res_freq, 'MHz')

#         expt_name = 'qubit_spec_ge'
#         with open('system_config.json', 'r') as f:
#             config = json.load(f)

#         # adjust experiment configuration values
#         config['config']['readout_cfg']['res_freq_ge'][QUBIT_INDEX] = res_freq
#         config['expt_cfg'][expt_name]['start'][QUBIT_INDEX] = qubit_start[qubit][i]
#         config['expt_cfg'][expt_name]['stop'][QUBIT_INDEX] = qubit_stop[qubit][i]

#         with open('system_config.json', 'w') as f:
#             json.dump(config, f, indent=3)
            
#         # run qubit spectroscopy experiment
#         print('Starting Qubit Spectroscopy')
#         qubit_freq = QubitSpectrosocpy()
#         if qubit == 0:
#             qubit1_freq_array.append(qubit_freq) # save res freq value
#         else:
#             qubit2_freq_array.append(qubit_freq) # save res freq value
        
#         print('Qubit frequency =', qubit_freq, 'MHz')

#         # expt_name = 'T1_ge'
#         # with open('system_config.json', 'r') as f:
#         #     config = json.load(f)

#         # # adjust experiment configuration values
#         # config['config']['readout_cfg']['res_freq_ge'][QUBIT_INDEX] = 6170.350609588622 ####res_freq
#         # if qubit == 0:
#         #     config['config']['qubit_cfg']['qubit_freq_ge'][QUBIT_INDEX] = 4182.745759963988
#         # else:
#         #     config['config']['qubit_cfg']['qubit_freq_ge'][QUBIT_INDEX] = 4202

#         # with open('system_config.json', 'w') as f:
#         #     json.dump(config, f, indent=3)

#         # # run t1 experiment
#         # print('Starting T1')
#         # t1_value = T1()
#         # if qubit == 0:
#         #     qubit1_t1_array.append(t1_value) # save res freq value
#         # else:
#         #     qubit2_t1_array.append(t1_value) # save res freq value
        
#         # print('T1 =', t1_value, 'us')

#     end_time = time.time()

#     print('Time taken for iteration =', end_time - start_time)

# # save all data info

# data_path = "M:/malab/_Data/20250130 - Santi - RFSoC tprocv2 - LL8qubit J meas Q1 Q2/data/collective_decay"
# lookup_file = 'lookup_file'
# with SlabFile(data_path + '\\' + lookup_file, 'a') as f:
#     # 'a': read/write/create

#     # - Adds data to the file - #
#     f.append('res1_freqs', res1_freq_array)
#     f.append('res2_freqs', res2_freq_array)
#     f.append('qubit1_freqs', qubit1_freq_array)
#     f.append('qubit2_freqs', qubit2_freq_array)
#     f.append('qubit1_t1', qubit1_t1_array)
#     f.append('qubit2_t1', qubit2_t1_array)
#     f.append('flux_array', flux_array)

drive_Yoko_and_DACs([0] * 9) # zero all the dacs