import sys
import os
import numpy as np
import datetime
sys.path.append(os.path.abspath("/home/quietuser/Documents/GitHub/tprocv2_demos/qick_tprocv2_experiments_mux/"))
from section_001_time_of_flight import TOFExperiment
from section_002_res_spec_ge_mux import ResonanceSpectroscopy
from section_004_qubit_spec_ge import QubitSpectroscopy
from section_006_amp_rabi_ge import AmplitudeRabiExperiment
from section_007_T1_ge import T1Measurement
from section_005_single_shot_ge import SingleShot
from section_008_save_T1_T2 import Save
from section_009_T2R_ge import T2RMeasurement
from system_config import QICK_experiment
from section_003_punch_out_ge_mux import PunchOut
from expt_config import *

n= 1
save_r = 1 #how many rounds to save after
signal = 'Q' #'None' # where the signal is (after ss/angle optimization). Put 'None' if no optimization has happened
save_figs = True
live_plot = True # for live plotting open http://localhost:8097/ on firefox
fit_data = False
outerFolder = "/data/QICK_data/6transmon_run4a/" + str(datetime.date.today()) + "/"


Qs = [0,1,2,3,4,5]
Qs=[1]
t1_data = {Q: {'T1': [None]*save_r, 'Errors': [None]*save_r, 'Dates': [None]*save_r,
               'I': [None]*save_r, 'Q':[None]*save_r, 'Delay Times': [None]*save_r,
               'Fit': [None]*save_r, 'Round Num': [None]*save_r, 'Batch Num': [None]*save_r} for Q in range(6)}

t2r_data = {Q: {'T2': [None]*save_r, 'Errors': [None]*save_r, 'Dates': [None]*save_r,
                'I': [None]*save_r, 'Q':[None]*save_r, 'Delay Times': [None]*save_r,
                'Fit': [None]*save_r, 'Round Num': [None]*save_r, 'Batch Num': [None]*save_r} for Q in range(6)}


#initalize variables to 0
batch_num=0
j = 0
angles=[]
#res_leng_vals = [, , , , , , ]
res_gain = [0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.9]
while j < n:
    j += 1
    for QubitIndex in Qs:
        #Get the config for this qubit
        experiment = QICK_experiment(outerFolder, DAC_attenuator1 = 5, DAC_attenuator2 = 10)

        #Mask out all other resonators except this one
        res_gains = experiment.set_gain_filter_ge(QubitIndex, IndexGain=res_gain[QubitIndex]) #change index gain after optimization
        experiment.readout_cfg['res_gain_ge'] = res_gains
        #experiment.readout_cfg['res_length'] = res_leng_vals[QubitIndex]
        # ---------------------TOF------------------------
        # tof        = TOFExperiment(QubitIndex, outerFolder, j, save_figs, experiment)
        # tof.run(experiment.soccfg, experiment.soc)
        # del tof

        #---------------------Res spec---------------------
        res_spec   = ResonanceSpectroscopy(QubitIndex, outerFolder, j, save_figs, experiment)
        res_freqs = res_spec.run(experiment.soccfg, experiment.soc)
        experiment.readout_cfg['res_freq_ge'] = res_freqs #change after optimization, add offset value to each of the freqs in this list [r + offset for r in res_freqs]
        print("res_freqs: ",res_freqs)

        del res_spec

        # #
        # # #-----------------Roll Signal into I---------------
        # # #get the average theta value, then use that to rotate the signal. Plug that value into system_config res_phase
        # # leng=4
        # # ss = SingleShot(QubitIndex, outerFolder, experiment, j, save_figs)
        # # fid, angle, iq_list_g, iq_list_e = ss.run(experiment.soccfg, experiment.soc)
        # # angles.append(angle)
        # # #print(angles)
        # # #print('avg theta: ', np.average(angles))
        # # del ss
        # #
        # #--------------------Qubit spec--------------------
        q_spec = QubitSpectroscopy(QubitIndex, outerFolder, j, signal, save_figs, experiment, live_plot)
        qubit_freq = q_spec.run(experiment.soccfg, experiment.soc)
        experiment.qubit_cfg['qubit_freq_ge'][QubitIndex] = float(qubit_freq)
        print('qubit freq: for qubit ', QubitIndex + 1 ,' is: ',float(qubit_freq))
        del q_spec

        #-----------------------Rabi-----------------------
        rabi = AmplitudeRabiExperiment(QubitIndex, outerFolder, j, signal, save_figs, experiment, live_plot, fit_data)
        rabi.run(experiment.soccfg, experiment.soc)
        del rabi
        #
        # #------------------------T1-------------------------
        # t1 = T1Measurement(QubitIndex, outerFolder, j, signal, save_figs, experiment, live_plot, fit_data)
        # t1_est, t1_err, I, Q, delay_times, q1_fit_exponential = t1.run(experiment.soccfg, experiment.soc)
        #
        # #---------------------Collect T1 Results----------------
        #
        # t1_data[QubitIndex]['T1'][j - batch_num*save_r - 1] = t1_est
        # t1_data[QubitIndex]['Errors'][j - batch_num*save_r - 1] = t1_err
        # t1_data[QubitIndex]['Dates'][j - batch_num*save_r - 1] = datetime.datetime.now()
        # t1_data[QubitIndex]['I'][j - batch_num*save_r - 1] = I
        # t1_data[QubitIndex]['Q'][j - batch_num*save_r - 1] = Q
        # t1_data[QubitIndex]['Delay Times'][j - batch_num*save_r - 1] = delay_times
        # t1_data[QubitIndex]['Fit'][j - batch_num*save_r - 1] = q1_fit_exponential
        # t1_data[QubitIndex]['Round Num'][j - batch_num * save_r - 1] = j
        # t1_data[QubitIndex]['Batch Num'][j - batch_num * save_r - 1] = batch_num
        #
        # #------------------------T2R-------------------------
        # signal = 'I'
        # t2r = T2RMeasurement(QubitIndex, outerFolder, j, signal, save_figs, experiment, live_plot, fit_data)
        # t2r_est, t2r_err, I, Q, delay_times, fit_ramsey = t2r.run(experiment.soccfg, experiment.soc)
        #
        # #---------------------Collect T2 Results----------------
        # t2r_data[QubitIndex]['T2'][j - batch_num*save_r - 1] = t2r_est
        # t2r_data[QubitIndex]['Errors'][j - batch_num*save_r - 1] = t2r_err
        # t2r_data[QubitIndex]['Dates'][j - batch_num*save_r - 1] = datetime.datetime.now()
        # t2r_data[QubitIndex]['I'][j - batch_num*save_r - 1] = I
        # t2r_data[QubitIndex]['Q'][j - batch_num*save_r - 1] = Q
        # t2r_data[QubitIndex]['Delay Times'][j - batch_num*save_r - 1] = delay_times
        # t2r_data[QubitIndex]['Fit'][j - batch_num*save_r - 1] = fit_ramsey
        # t2r_data[QubitIndex]['Round Num'][j - batch_num * save_r - 1] = j
        # t2r_data[QubitIndex]['Batch Num'][j - batch_num * save_r - 1] = batch_num


        del experiment

    #-----------------------Potentially Save---------------
    # Check if you are at the right round number. If so, then save all of the data and change the round num so you replace data starting next round
    if j % save_r == 0:
        batch_num+=1
        saver_t1 = Save(outerFolder,t1_data, batch_num, save_r)
        saver_t1.save_to_h5('T1')
        del saver_t1
        saver_t2r = Save(outerFolder, t2r_data, batch_num, save_r)
        saver_t2r.save_to_h5('T2')
        del saver_t2r





