from section_002_res_spec_ge_mux import ResonanceSpectroscopy
from section_004_qubit_spec_ge import QubitSpectroscopy
from section_006_amp_rabi_ge import AmplitudeRabiExperiment
from section_007_T1_ge import T1Measurement
from section_008_save_data_to_h5 import Data_H5
from section_009_T2R_ge import T2RMeasurement
from expt_config import *
import glob

date = '2024-11-15'
outerFolder = "/data/QICK_data/6transmon_run4a/" + date + "/"

save_figs = True
fit_saved = False
signal = 'Q'

loader_config_instance = Data_H5(outerFolder)
config = loader_config_instance.load_config()
del loader_config_instance

# ------------------------------------------Load/Plot/Save Res Spec------------------------------------
outerFolder_expt = outerFolder + "/Data_h5/Res_ge/"
h5_files = glob.glob(os.path.join(outerFolder_expt, "*.h5"))

for h5_file in h5_files:
    save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]
    H5_class_instance = Data_H5(h5_file)
    #H5_class_instance.print_h5_contents(h5_file)
    load_data = H5_class_instance.load_from_h5(data_type=  'Res', save_r = int(save_round))
    for q_key in load_data['Res']:
        dates = load_data['Res'][q_key].get('Dates', [])
        freq_pts = load_data['Res'][q_key].get('freq_pts', [])
        freq_center = load_data['Res'][q_key].get('freq_center', [])
        freqs_found = load_data['Res'][q_key].get('Found Freqs', [])
        amps = load_data['Res'][q_key].get('Amps', [])
        round_num = load_data['Res'][q_key].get('Round Num', [])
        batch_num = load_data['Res'][q_key].get('Batch Num', [])

        if len(freq_pts[0]) > 0:
            dates = dates[0][0]
            freq_pts = freq_pts[0][0]
            freq_center = freq_center[0][0]
            freqs_found = freqs_found[0][0]
            amps = amps[0][0]
            round_num = str(round_num[0][0]).split('.')[0]
            batch_num = batch_num[0][0]

            res_class_instance = ResonanceSpectroscopy(q_key, outerFolder, round_num, save_figs)
            res_class_instance.plot_results(freq_pts, freq_center, amps)
            del res_class_instance

    del H5_class_instance

# ----------------------------------------------Load/Plot/Save QSpec------------------------------------
outerFolder_expt = outerFolder + "/Data_h5/QSpec_ge/"
h5_files = glob.glob(os.path.join(outerFolder_expt, "*.h5"))

for h5_file in h5_files:
    save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]
    H5_class_instance = Data_H5(h5_file)
    load_data = H5_class_instance.load_from_h5(data_type=  'QSpec', save_r = int(save_round))

    for q_key in load_data['QSpec']:
        dates = load_data['QSpec'][q_key].get('Dates', [])
        I = load_data['QSpec'][q_key].get('I', [])
        Q = load_data['QSpec'][q_key].get('Q', [])
        I_fit = load_data['QSpec'][q_key].get('I Fit', [])
        Q_fit = load_data['QSpec'][q_key].get('Q Fit', [])
        freqs = load_data['QSpec'][q_key].get('Frequencies', [])
        round_num = load_data['QSpec'][q_key].get('Round Num', [])
        batch_num = load_data['QSpec'][q_key].get('Batch Num', [])

        if len(I[0])>0:
            dates = dates[0][0]
            I = I[0][0]
            Q = Q[0][0]
            I_fit = I_fit[0][0]
            Q_fit = Q_fit[0][0]
            freqs = freqs[0][0]
            round_num = str(round_num[0][0]).split('.')[0]
            batch_num = batch_num[0][0]

            qspec_class_instance = QubitSpectroscopy(q_key, outerFolder, round_num, signal, save_figs)
            qspec_class_instance.plot_results(I, Q, freqs, config)
            del qspec_class_instance

    del H5_class_instance

# ------------------------------------------------Load/Plot/Save Rabi---------------------------------------
outerFolder_expt = outerFolder + "/Data_h5/Rabi_ge/"
h5_files = glob.glob(os.path.join(outerFolder_expt, "*.h5"))

for h5_file in h5_files:

    save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]
    H5_class_instance = Data_H5(h5_file)
    load_data = H5_class_instance.load_from_h5(data_type=  'Rabi', save_r = int(save_round))

    for q_key in load_data['Rabi']:
        dates= load_data['Rabi'][q_key].get('Dates', [])
        I = load_data['Rabi'][q_key].get('I', [])
        Q = load_data['Rabi'][q_key].get('Q', [])
        gains = load_data['Rabi'][q_key].get('Gains', [])
        fit = load_data['Rabi'][q_key].get('Fit', [])
        round_num = load_data['Rabi'][q_key].get('Round Num', [])
        batch_num = load_data['Rabi'][q_key].get('Batch Num', [])

        if len(I[0])>0:
            dates = dates[0][0]
            I = I[0][0]
            Q = Q[0][0]
            gains = gains[0][0]
            if fit_saved:
                fit = fit[0][0]
            round_num = str(round_num[0][0]).split('.')[0]
            batch_num = batch_num[0][0]

            rabi_class_instance = AmplitudeRabiExperiment(q_key, outerFolder, round_num, signal, save_figs)
            rabi_class_instance.plot_results(I, Q, gains, config)
            del rabi_class_instance

    del H5_class_instance

# ------------------------------------------------Load/Plot/Save T1----------------------------------------------
outerFolder_expt = outerFolder + "/Data_h5/T1_ge/"
h5_files = glob.glob(os.path.join(outerFolder_expt, "*.h5"))

for h5_file in h5_files:

    save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]
    H5_class_instance = Data_H5(h5_file)
    load_data = H5_class_instance.load_from_h5(data_type=  'T1', save_r = int(save_round))

    for q_key in load_data['T1']:
        T1 = load_data['T1'][q_key].get('T1', [])
        errors = load_data['T1'][q_key].get('Errors', [])
        dates= load_data['T1'][q_key].get('Dates', [])
        I = load_data['T1'][q_key].get('I', [])
        Q = load_data['T1'][q_key].get('Q', [])
        delay_times = load_data['T1'][q_key].get('Delay Times', [])
        fit = load_data['T1'][q_key].get('Fit', [])
        round_num = load_data['T1'][q_key].get('Round Num', [])
        batch_num = load_data['T1'][q_key].get('Batch Num', [])

        if len(I[0])>0:
            dates = dates[0][0]
            I = I[0][0]
            Q = Q[0][0]
            delay_times = delay_times[0][0]
            if fit_saved:
                T1 = T1[0][0]
                T1_err = T1_err[0][0]
                fit = fit[0][0]
            round_num = str(round_num[0][0]).split('.')[0]
            batch_num = batch_num[0][0]

            T1_class_instance = T1Measurement(q_key, outerFolder, round_num, signal, save_figs)
            T1_class_instance.plot_results(I, Q, delay_times, dates[0], config)
            del T1_class_instance

    del H5_class_instance

# -------------------------------------------------------Load/Plot/Save T2------------------------------------------
outerFolder_expt = outerFolder + "/Data_h5/T2_ge/"
h5_files = glob.glob(os.path.join(outerFolder_expt, "*.h5"))

for h5_file in h5_files:
    save_round = h5_file.split('Num_per_batch')[-1].split('.')[0]
    H5_class_instance = Data_H5(h5_file)
    load_data = H5_class_instance.load_from_h5(data_type=  'T2', save_r = int(save_round))

    for q_key in load_data['T2']:
        T2 = load_data['T2'][q_key].get('T2', [])
        errors = load_data['T2'][q_key].get('Errors', [])
        dates = load_data['T2'][q_key].get('Dates', [])
        I = load_data['T2'][q_key].get('I', [])
        Q = load_data['T2'][q_key].get('Q', [])
        delay_times = load_data['T2'][q_key].get('Delay Times', [])
        fit = load_data['T2'][q_key].get('Fit', [])
        round_num = load_data['T2'][q_key].get('Round Num', [])
        batch_num = load_data['T2'][q_key].get('Batch Num', [])

        if len(I[0]) > 0:
            dates = dates[0][0]
            I = I[0][0]
            Q = Q[0][0]
            delay_times = delay_times[0][0]
            if fit_saved:
                T1 = T1[0][0]
                T1_err = T1_err[0][0]
                fit = fit[0][0]
            round_num = str(round_num[0][0]).split('.')[0]
            batch_num = batch_num[0][0]

            T2_class_instance = T2RMeasurement(q_key, outerFolder, round_num, signal, save_figs)
            fitted, t2r_est, t2r_err = T2_class_instance.t2_fit(delay_times, I)
            T2_class_instance.plot_results(I, Q, delay_times, dates[0], fitted, t2r_est, t2r_err,  config)
            del T2_class_instance

    del H5_class_instance