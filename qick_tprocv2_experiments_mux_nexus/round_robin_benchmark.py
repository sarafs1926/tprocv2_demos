import sys
import os
import numpy as np
import datetime
import time
import visdom

sys.path.append(os.path.abspath("/home/nexusadmin/Documents/GitHub/tprocv2_demos/qick_tprocv2_experiments_mux_nexus/"))
from section_001_time_of_flight import TOFExperiment
from section_002_res_spec_ge_mux import ResonanceSpectroscopy
from section_004_qubit_spec_ge import QubitSpectroscopy
from section_006_amp_rabi_ge import AmplitudeRabiExperiment
from section_007_T1_ge import T1Measurement
from section_005_single_shot_ge import SingleShot
from section_008_save_data_to_h5 import Data_H5
from section_009_T2R_ge import T2RMeasurement
from section_010_T2E_ge import T2EMeasurement
from system_config import QICK_experiment
from section_003_punch_out_ge_mux import PunchOut
from expt_config import expt_cfg

n= 1000
save_r = 1            # how many rounds to save after
signal = 'None'       #'I', or 'Q' depending on where the signal is (after optimization). Put 'None' if no optimization has happened
save_figs = True    # save plots for everything as you go along the RR script?
live_plot = False      # for live plotting open http://localhost:8097/ on firefox
fit_data = False      # fit the data here and save or plot the fits?, normally kept as False to avoid issues over night
save_data_h5 = False   # save all of the data to h5 files?
outerFolder = os.path.join("/home/nexusadmin/qick/NEXUS_sandbox/Data/", str(datetime.date.today()))
custom_Ramsey=True

# Note if these Qs are not in order and dont start w 0, data wont save!
Qs = [0,1,2,3] #this no longer specifies the qubits we want to look at, see line 68. We can make this better later.

# #optimization outputs
#res_leng_vals = [10.5, 3.55, 3.55, 3.55] # From NR25
res_leng_vals = [9, 5.5, 6.5, 9] # From 12/9 optimization
res_gain = [0.46, 0.38, 0.46, 0.42] # Note 0.15 could be better for all, not sure yet
# spec_gain= 0.85
# res_gain = [spec_gain, 1, spec_gain, spec_gain]
freq_offsets = [-0.08, 0.0, -0.24, -0.32]

def create_data_dict(keys, save_r, qs):
    return {Q: {key: np.empty(save_r, dtype=object) for key in keys} for Q in range(len(qs))}

# Define what to save to h5 files
res_keys = ['Dates', 'freq_pts', 'Amps', 'Found Freqs', 'Round Num', 'Batch Num']
qspec_keys = ['Dates', 'I', 'Q', 'Frequencies', 'I Fit', 'Q Fit', 'Round Num', 'Batch Num']
rabi_keys = ['Dates', 'I', 'Q', 'Gains', 'Fit', 'Round Num', 'Batch Num']
t1_keys = ['T1', 'Errors', 'Dates', 'I', 'Q', 'Delay Times', 'Fit', 'Round Num', 'Batch Num']
t2r_keys = ['T2', 'Errors', 'Dates', 'I', 'Q', 'Delay Times', 'Fit', 'Round Num', 'Batch Num']
t2e_keys = ['T2E', 'Errors', 'Dates', 'I', 'Q', 'Delay Times', 'Fit', 'Round Num', 'Batch Num']

#initialize a dictionary to store those values
res_data = create_data_dict(res_keys, save_r, Qs)
qspec_data = create_data_dict(qspec_keys, save_r, Qs)
rabi_data = create_data_dict(rabi_keys, save_r, Qs)
t1_data = create_data_dict(t1_keys, save_r, Qs)
t2r_data = create_data_dict(t2r_keys, save_r, Qs)
t2e_data = create_data_dict(t2e_keys, save_r, Qs)

batch_num=0
j = 0
angles=[]
while j < n:
    j += 1
    #for QubitIndex in Qs:
    for QubitIndex in [0,1,2,3]: #specify which qubits you want to look at inside the brackets
        #Get the config for this qubit
        experiment = QICK_experiment(outerFolder)

        #Mask out all other resonators except this one
        res_gains = experiment.mask_gain_res(QubitIndex, IndexGain=res_gain[QubitIndex])
        experiment.readout_cfg['res_gain_ge'] = res_gains
        experiment.readout_cfg['res_length'] = res_leng_vals[QubitIndex]

    #    #---------------------TOF------------------------
    #    tof        = TOFExperiment(QubitIndex, outerFolder, experiment, j, save_figs)
    #    tof.run(experiment.soccfg, experiment.soc)
    #    del tof
    #
        #---------------------Res spec---------------------
        res_spec   = ResonanceSpectroscopy(QubitIndex, outerFolder, j, save_figs, experiment)
        res_freqs, freq_pts, amps, sys_config_to_save = res_spec.run(experiment.soccfg, experiment.soc)
        #experiment.readout_cfg['res_freq_ge'] = res_freqs
        offset = freq_offsets[QubitIndex] #use optimized offset values
        offset_res_freqs = [r + offset for r in res_freqs]
        experiment.readout_cfg['res_freq_ge'] = offset_res_freqs
        # experiment.readout_cfg['res_freq_ge'] = [6187.151, 5827.778, 6074.035, 5958.633] # Punched out freqs
        del res_spec
    #
    #     # #-----------------Roll Signal into I---------------
    #     # #get the average theta value, then use that to rotate the signal. Plug that value into system_config res_phase
    #     # leng=4
    #     # ss = SingleShot(QubitIndex, outerFolder, experiment, j, save_figs)
    #     # fid, angle, iq_list_g, iq_list_e = ss.run(experiment.soccfg, experiment.soc)
    #     # angles.append(angle)
    #     # #print(angles)
    #     # #print('avg theta: ', np.average(angles))
    #     # del ss

        # #--------------------Qubit spec--------------------
        # q_spec = QubitSpectroscopy(QubitIndex, outerFolder, j, signal, save_figs, experiment, live_plot)
        # qspec_I, qspec_Q, qspec_freqs, qspec_I_fit, qspec_Q_fit, qubit_freq = q_spec.run(experiment.soccfg, experiment.soc)
        # experiment.qubit_cfg['qubit_freq_ge'][QubitIndex] = float(qubit_freq)
        # print('Qubit freq for qubit ', QubitIndex + 1 ,' is: ',float(qubit_freq))
        # del q_spec
        #
        # #-----------------------Rabi-----------------------
        # rabi = AmplitudeRabiExperiment(QubitIndex, outerFolder, j, signal, save_figs, experiment, live_plot)
        # rabi_I, rabi_Q, rabi_gains, rabi_fit, pi_amp  = rabi.run(experiment.soccfg, experiment.soc)
        # experiment.qubit_cfg['pi_amp'][QubitIndex] = float(pi_amp)
        # print('Pi amplitude for qubit ', QubitIndex + 1, ' is: ', float(pi_amp))
        # del rabi
        # #
        # #------------------Single Shot Measurements---------------
        # ss = SingleShot(QubitIndex, outerFolder, experiment, round_num=0, save_figs=True)
        # fid, angle, iq_list_g, iq_list_e = ss.run(experiment.soccfg, experiment.soc)
        #
        # I_g = iq_list_g[QubitIndex][0].T[0]
        # Q_g = iq_list_g[QubitIndex][0].T[1]
        # I_e = iq_list_e[QubitIndex][0].T[0]
        # Q_e = iq_list_e[QubitIndex][0].T[1]
        #
        # fid, threshold, rotation_angle, ig_new, ie_new = ss.hist_ssf(
        #     data=[I_g, Q_g, I_e, Q_e], cfg=ss.config, plot=True)
        #
    #     #------------------------T1-------------------------
    #     t1 = T1Measurement(QubitIndex, outerFolder, j, signal, save_figs, experiment, live_plot, fit_data)
    #     t1_est, t1_err, t1_I, t1_Q, t1_delay_times, q1_fit_exponential = t1.run(experiment.soccfg, experiment.soc)
    #     del t1
    #
    #     #------------------------T2R-------------------------
    #     t2r = T2RMeasurement(QubitIndex, outerFolder, j, signal, save_figs, experiment, live_plot, fit_data, custom_Ramsey)
    #     t2r_est, t2r_err, t2r_I, t2r_Q, t2r_delay_times, fit_ramsey = t2r.run(experiment.soccfg, experiment.soc)
    #     del t2r
    #
    #     # ------------------------T2E-------------------------
    #     t2e = T2EMeasurement(QubitIndex, outerFolder, j, signal, save_figs, experiment, live_plot, fit_data)
    #     t2e_est, t2e_err, t2e_I, t2e_Q, t2e_delay_times, fit_ramsey_t2e, sys_config_to_save = t2e.run(experiment.soccfg, experiment.soc)
    #     del t2e
    # # #
        if save_data_h5:
            # ---------------------Collect Res Spec Results----------------
            res_data[QubitIndex]['Dates'][j - batch_num * save_r - 1] = time.mktime(datetime.datetime.now().timetuple())
            res_data[QubitIndex]['freq_pts'][j - batch_num * save_r - 1] = freq_pts
            res_data[QubitIndex]['Amps'][j - batch_num * save_r - 1] = amps
            res_data[QubitIndex]['Found Freqs'][j - batch_num * save_r - 1] = res_freqs
            res_data[QubitIndex]['Round Num'][j - batch_num * save_r - 1] = j
            res_data[QubitIndex]['Batch Num'][j - batch_num * save_r - 1] = batch_num

            # ---------------------Collect QSpec Results----------------
            qspec_data[QubitIndex]['Dates'][j - batch_num * save_r - 1] = time.mktime(datetime.datetime.now().timetuple())
            qspec_data[QubitIndex]['I'][j - batch_num * save_r - 1] = qspec_I
            qspec_data[QubitIndex]['Q'][j - batch_num * save_r - 1] = qspec_Q
            qspec_data[QubitIndex]['Frequencies'][j - batch_num * save_r - 1] = qspec_freqs
            qspec_data[QubitIndex]['I Fit'][j - batch_num * save_r - 1] = qspec_I_fit
            qspec_data[QubitIndex]['Q Fit'][j - batch_num * save_r - 1] = qspec_Q_fit
            qspec_data[QubitIndex]['Round Num'][j - batch_num * save_r - 1] = j
            qspec_data[QubitIndex]['Batch Num'][j - batch_num * save_r - 1] = batch_num

            #---------------------Collect Rabi Results----------------
            rabi_data[QubitIndex]['Dates'][j - batch_num * save_r - 1] = time.mktime(datetime.datetime.now().timetuple())
            rabi_data[QubitIndex]['I'][j - batch_num * save_r - 1] = rabi_I
            rabi_data[QubitIndex]['Q'][j - batch_num * save_r - 1] = rabi_Q
            rabi_data[QubitIndex]['Gains'][j - batch_num * save_r - 1] = rabi_gains
            rabi_data[QubitIndex]['Fit'][j - batch_num * save_r - 1] = rabi_fit
            rabi_data[QubitIndex]['Round Num'][j - batch_num * save_r - 1] = j
            rabi_data[QubitIndex]['Batch Num'][j - batch_num * save_r - 1] = batch_num

            #---------------------Collect T1 Results----------------
            t1_data[QubitIndex]['T1'][j - batch_num*save_r - 1] = t1_est
            t1_data[QubitIndex]['Errors'][j - batch_num*save_r - 1] = t1_err
            t1_data[QubitIndex]['Dates'][j - batch_num*save_r - 1] = time.mktime(datetime.datetime.now().timetuple())
            t1_data[QubitIndex]['I'][j - batch_num*save_r - 1] = t1_I
            t1_data[QubitIndex]['Q'][j - batch_num*save_r - 1] = t1_Q
            t1_data[QubitIndex]['Delay Times'][j - batch_num*save_r - 1] = t1_delay_times
            t1_data[QubitIndex]['Fit'][j - batch_num*save_r - 1] = q1_fit_exponential
            t1_data[QubitIndex]['Round Num'][j - batch_num * save_r - 1] = j
            t1_data[QubitIndex]['Batch Num'][j - batch_num * save_r - 1] = batch_num

            #---------------------Collect T2 Results----------------
            t2r_data[QubitIndex]['T2'][j - batch_num*save_r - 1] = t2r_est
            t2r_data[QubitIndex]['Errors'][j - batch_num*save_r - 1] = t2r_err
            t2r_data[QubitIndex]['Dates'][j - batch_num*save_r - 1] = time.mktime(datetime.datetime.now().timetuple())
            t2r_data[QubitIndex]['I'][j - batch_num*save_r - 1] = t2r_I
            t2r_data[QubitIndex]['Q'][j - batch_num*save_r - 1] = t2r_Q
            t2r_data[QubitIndex]['Delay Times'][j - batch_num*save_r - 1] = t2r_delay_times
            t2r_data[QubitIndex]['Fit'][j - batch_num*save_r - 1] = fit_ramsey
            t2r_data[QubitIndex]['Round Num'][j - batch_num * save_r - 1] = j
            t2r_data[QubitIndex]['Batch Num'][j - batch_num * save_r - 1] = batch_num

            #---------------------Collect T2E Results----------------
            t2e_data[QubitIndex]['T2E'][j - batch_num*save_r - 1] = t2e_est
            t2e_data[QubitIndex]['Errors'][j - batch_num*save_r - 1] = t2e_err
            t2e_data[QubitIndex]['Dates'][j - batch_num*save_r - 1] = time.mktime(datetime.datetime.now().timetuple())
            t2e_data[QubitIndex]['I'][j - batch_num*save_r - 1] = t2e_I
            t2e_data[QubitIndex]['Q'][j - batch_num*save_r - 1] = t2e_Q
            t2e_data[QubitIndex]['Delay Times'][j - batch_num*save_r - 1] = t2e_delay_times
            t2e_data[QubitIndex]['Fit'][j - batch_num*save_r - 1] = fit_ramsey_t2e
            t2e_data[QubitIndex]['Round Num'][j - batch_num * save_r - 1] = j
            t2e_data[QubitIndex]['Batch Num'][j - batch_num * save_r - 1] = batch_num

            # ---------------------save last system config and expt_cfg---------------------
            saver_config = Data_H5(outerFolder)
            saver_config.save_config(sys_config_to_save, expt_cfg)
            del saver_config
    #
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
            del res_data

            # --------------------------save QSpec-----------------------
            saver_qspec = Data_H5(outerFolder, qspec_data, batch_num, save_r)
            saver_qspec.save_to_h5('QSpec')
            del saver_qspec
            del qspec_data

            # --------------------------save Rabi-----------------------
            saver_rabi = Data_H5(outerFolder, rabi_data, batch_num, save_r)
            saver_rabi.save_to_h5('Rabi')
            del saver_rabi
            del rabi_data

            # --------------------------save t1-----------------------
            saver_t1 = Data_H5(outerFolder, t1_data, batch_num, save_r)
            saver_t1.save_to_h5('T1')
            del saver_t1
            del t1_data

            #--------------------------save t2r-----------------------
            saver_t2r = Data_H5(outerFolder, t2r_data, batch_num, save_r)
            saver_t2r.save_to_h5('T2')
            del saver_t2r
            del t2r_data

            #--------------------------save t2e-----------------------
            saver_t2e = Data_H5(outerFolder, t2e_data, batch_num, save_r)
            saver_t2e.save_to_h5('T2E')
            del saver_t2e
            del t2e_data

            # reset all dictionaries to none for safety
            res_data = create_data_dict(res_keys, save_r, Qs)
            qspec_data = create_data_dict(qspec_keys, save_r, Qs)
            rabi_data = create_data_dict(rabi_keys, save_r, Qs)
            t1_data = create_data_dict(t1_keys, save_r, Qs)
            t2r_data = create_data_dict(t2r_keys, save_r, Qs)
            t2e_data = create_data_dict(t2e_keys, save_r, Qs)




