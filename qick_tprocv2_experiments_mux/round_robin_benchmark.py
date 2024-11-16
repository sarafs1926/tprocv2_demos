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
from section_008_save_data_to_h5 import Data_H5
from section_009_T2R_ge import T2RMeasurement
from system_config import QICK_experiment
from section_003_punch_out_ge_mux import PunchOut
from expt_config import *


n= 1000
save_r = 1          # how many rounds to save after
signal = 'Q'        #'I', or 'Q' depending on where the signal is (after optimization). Put 'None' if no optimization has happened
save_figs = False   # save plots for everything as you go along the RR script?
live_plot = True    # for live plotting open http://localhost:8097/ on firefox
fit_data = False    # fit the data here and save or plot the fits?
save_data_h5 = True # save all of the data to h5 files?
outerFolder = "/data/QICK_data/6transmon_run4a/" + str(datetime.date.today()) + "/"

Qs = [0,1,2,3,4,5]
#Qs = [0,1]
res_data = {Q: {'Dates': [None]*save_r, 'freq_pts': [None]*save_r, 'freq_center': [None]*save_r, 'Amps': [None]*save_r,
                'Found Freqs': [None]*save_r, 'Round Num': [None]*save_r, 'Batch Num': [None]*save_r} for Q in range(6)}

qspec_data = {Q: {'Dates': [None]*save_r, 'I': [None]*save_r, 'Q':[None]*save_r, 'Frequencies': [None]*save_r,
                'I Fit': [None]*save_r, 'Q Fit': [None]*save_r, 'Round Num': [None]*save_r, 'Batch Num': [None]*save_r} for Q in range(6)}

rabi_data = {Q: {'Dates': [None]*save_r, 'I': [None]*save_r, 'Q':[None]*save_r, 'Gains': [None]*save_r,
                'Fit': [None]*save_r, 'Round Num': [None]*save_r, 'Batch Num': [None]*save_r} for Q in range(6)}

t1_data = {Q: {'T1': [None]*save_r, 'Errors': [None]*save_r, 'Dates': [None]*save_r,
               'I': [None]*save_r, 'Q':[None]*save_r, 'Delay Times': [None]*save_r,
               'Fit': [None]*save_r, 'Round Num': [None]*save_r, 'Batch Num': [None]*save_r} for Q in range(6)}

t2r_data = {Q: {'T2': [None]*save_r, 'Errors': [None]*save_r, 'Dates': [None]*save_r,
                'I': [None]*save_r, 'Q':[None]*save_r, 'Delay Times': [None]*save_r,
                'Fit': [None]*save_r, 'Round Num': [None]*save_r, 'Batch Num': [None]*save_r} for Q in range(6)}


batch_num=0
j = 0
angles=[]
#res_leng_vals = [, , , , , , ]
res_gain = [0.7, 0.9, 0.7, 0.7, 0.7, 0.9, 0.9]
while j < n:
    j += 1
    for QubitIndex in Qs:
        #Get the config for this qubit
        experiment = QICK_experiment(outerFolder, DAC_attenuator1 = 5, DAC_attenuator2 = 10)

        #Mask out all other resonators except this one
        res_gains = experiment.set_gain_filter_ge(QubitIndex, IndexGain=res_gain[QubitIndex])
        experiment.readout_cfg['res_gain_ge'] = res_gains
        #experiment.readout_cfg['res_length'] = res_leng_vals[QubitIndex]

        # ---------------------TOF------------------------
        # tof        = TOFExperiment(QubitIndex, outerFolder, j, save_figs, experiment)
        # tof.run(experiment.soccfg, experiment.soc)
        # del tof

        #---------------------Res spec---------------------
        res_spec   = ResonanceSpectroscopy(QubitIndex, outerFolder, j, save_figs, experiment)
        res_freqs, freq_pts, freq_center, amps = res_spec.run(experiment.soccfg, experiment.soc)
        experiment.readout_cfg['res_freq_ge'] = res_freqs #change after optimization, add offset value to each of the freqs in this list [r + offset for r in res_freqs]
        print("Resonator frequencies found: ",res_freqs)
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
        q_spec = QubitSpectroscopy(QubitIndex, outerFolder, j, 'I', save_figs, experiment, live_plot)
        qspec_I, qspec_Q, qspec_freqs, qspec_I_fit, qspec_Q_fit, qubit_freq = q_spec.run(experiment.soccfg, experiment.soc)
        experiment.qubit_cfg['qubit_freq_ge'][QubitIndex] = float(qubit_freq)
        print('Qubit freq for qubit ', QubitIndex + 1 ,' is: ',float(qubit_freq))
        del q_spec

        #-----------------------Rabi-----------------------
        rabi = AmplitudeRabiExperiment(QubitIndex, outerFolder, j, signal, save_figs, experiment, live_plot, fit_data)
        rabi_I, rabi_Q, rabi_gains, rabi_fit = rabi.run(experiment.soccfg, experiment.soc)
        del rabi

        #------------------------T1-------------------------
        t1 = T1Measurement(QubitIndex, outerFolder, j, signal, save_figs, experiment, live_plot, fit_data)
        t1_est, t1_err, t1_I, t1_Q, t1_delay_times, q1_fit_exponential = t1.run(experiment.soccfg, experiment.soc)

        # ------------------------T2R-------------------------
        signal = 'I'
        t2r = T2RMeasurement(QubitIndex, outerFolder, j, signal, save_figs, experiment, live_plot, fit_data)
        t2r_est, t2r_err, t2r_I, t2r_Q, t2r_delay_times, fit_ramsey, config_to_save = t2r.run(experiment.soccfg, experiment.soc)

        if save_data_h5:
            # ---------------------Collect Res Spec Results----------------
            res_data[QubitIndex]['Dates'][j - batch_num * save_r - 1] = datetime.datetime.now()
            res_data[QubitIndex]['freq_pts'][j - batch_num * save_r - 1] = freq_pts
            res_data[QubitIndex]['freq_center'][j - batch_num * save_r - 1] = freq_center
            res_data[QubitIndex]['Amps'][j - batch_num * save_r - 1] = amps
            res_data[QubitIndex]['Found Freqs'][j - batch_num * save_r - 1] = res_freqs
            res_data[QubitIndex]['Round Num'][j - batch_num * save_r - 1] = j
            res_data[QubitIndex]['Batch Num'][j - batch_num * save_r - 1] = batch_num

            # ---------------------Collect QSpec Results----------------
            qspec_data[QubitIndex]['Dates'][j - batch_num * save_r - 1] = datetime.datetime.now()
            qspec_data[QubitIndex]['I'][j - batch_num * save_r - 1] = qspec_I
            qspec_data[QubitIndex]['Q'][j - batch_num * save_r - 1] = qspec_Q
            qspec_data[QubitIndex]['Frequencies'][j - batch_num * save_r - 1] = qspec_freqs
            qspec_data[QubitIndex]['I Fit'][j - batch_num * save_r - 1] = qspec_I_fit
            qspec_data[QubitIndex]['Q Fit'][j - batch_num * save_r - 1] = qspec_Q_fit
            qspec_data[QubitIndex]['Round Num'][j - batch_num * save_r - 1] = j
            qspec_data[QubitIndex]['Batch Num'][j - batch_num * save_r - 1] = batch_num

            # ---------------------Collect Rabi Results----------------
            rabi_data[QubitIndex]['Dates'][j - batch_num * save_r - 1] = datetime.datetime.now()
            rabi_data[QubitIndex]['I'][j - batch_num * save_r - 1] = rabi_I
            rabi_data[QubitIndex]['Q'][j - batch_num * save_r - 1] = rabi_Q
            rabi_data[QubitIndex]['Gains'][j - batch_num * save_r - 1] = rabi_gains
            rabi_data[QubitIndex]['Fit'][j - batch_num * save_r - 1] = rabi_fit
            rabi_data[QubitIndex]['Round Num'][j - batch_num * save_r - 1] = j
            rabi_data[QubitIndex]['Batch Num'][j - batch_num * save_r - 1] = batch_num

            #---------------------Collect T1 Results----------------
            t1_data[QubitIndex]['T1'][j - batch_num*save_r - 1] = t1_est
            t1_data[QubitIndex]['Errors'][j - batch_num*save_r - 1] = t1_err
            t1_data[QubitIndex]['Dates'][j - batch_num*save_r - 1] = datetime.datetime.now()
            t1_data[QubitIndex]['I'][j - batch_num*save_r - 1] = t1_I
            t1_data[QubitIndex]['Q'][j - batch_num*save_r - 1] = t1_Q
            t1_data[QubitIndex]['Delay Times'][j - batch_num*save_r - 1] = t1_delay_times
            t1_data[QubitIndex]['Fit'][j - batch_num*save_r - 1] = q1_fit_exponential
            t1_data[QubitIndex]['Round Num'][j - batch_num * save_r - 1] = j
            t1_data[QubitIndex]['Batch Num'][j - batch_num * save_r - 1] = batch_num

            #---------------------Collect T2 Results----------------
            t2r_data[QubitIndex]['T2'][j - batch_num*save_r - 1] = t2r_est
            t2r_data[QubitIndex]['Errors'][j - batch_num*save_r - 1] = t2r_err
            t2r_data[QubitIndex]['Dates'][j - batch_num*save_r - 1] = datetime.datetime.now()
            t2r_data[QubitIndex]['I'][j - batch_num*save_r - 1] = t2r_I
            t2r_data[QubitIndex]['Q'][j - batch_num*save_r - 1] = t2r_Q
            t2r_data[QubitIndex]['Delay Times'][j - batch_num*save_r - 1] = t2r_delay_times
            t2r_data[QubitIndex]['Fit'][j - batch_num*save_r - 1] = fit_ramsey
            t2r_data[QubitIndex]['Round Num'][j - batch_num * save_r - 1] = j
            t2r_data[QubitIndex]['Batch Num'][j - batch_num * save_r - 1] = batch_num

            # ---------------------save last config---------------------
            saver_config = Data_H5(outerFolder)
            saver_config.save_config(config_to_save)
            del saver_config

        del experiment

    # --------------------------------Potentially Save---------------------------
    if save_data_h5:
        # Check if you are at the right round number. If so, then save all of the data and change the round num so you replace data starting next round
        if j % save_r == 0:
            batch_num+=1

            # --------------------------save Res Spec-----------------------
            saver_res = Data_H5(outerFolder, res_data, batch_num, save_r)
            saver_res.save_to_h5('Res')
            del saver_res

            # --------------------------save QSpec-----------------------
            saver_qspec = Data_H5(outerFolder, qspec_data, batch_num, save_r)
            saver_qspec.save_to_h5('QSpec')
            del saver_qspec

            # --------------------------save Rabi-----------------------
            saver_rabi = Data_H5(outerFolder, rabi_data, batch_num, save_r)
            saver_rabi.save_to_h5('Rabi')
            del saver_rabi

            # --------------------------save t1-----------------------
            saver_t1 = Data_H5(outerFolder, t1_data, batch_num, save_r)
            saver_t1.save_to_h5('T1')
            del saver_t1

            #--------------------------save t2-----------------------
            saver_t2r = Data_H5(outerFolder, t2r_data, batch_num, save_r)
            saver_t2r.save_to_h5('T2')
            del saver_t2r






