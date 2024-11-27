import sys
import os
import h5py
import matplotlib.pyplot as plt
import numpy as np

sys.path.append(os.path.abspath("/home/quietuser/Documents/GitHub/tprocv2_demos/qick_tprocv2_experiments_mux/"))
from system_config import QICK_experiment
from section_003_punch_out_ge_mux import PunchOut
import datetime
from section_001_time_of_flight import TOFExperiment

#sweep_ADC_attenuator = [10,15,20,25,30] #np.linspace(5,20, 4)
sweep_ADC_attenuator = np.linspace(0,30, 100)
outerFolder = "/data/QICK_data/6transmon_run4a/" + str(datetime.date.today()) + "/TOF_ADC_test/"

avg_mag_y_mid_DAC_atts = [] #list of 2 lists for each DAC att, each with 6 sublists for each Q
avg_mag_y_oct_DAC_atts = []
avg_mag_y_last_DAC_atts = []

DAC_att_1 = [2,5]
for DAC_att1 in DAC_att_1:
    avg_mag_y_mid = []
    avg_mag_y_oct = []
    avg_mag_y_last = []
    DAC_attenuator1 = None
    DAC_attenuator2 = None
    ADC_attenuator = None
    for att_1 in sweep_ADC_attenuator:
        experiment = QICK_experiment(outerFolder,DAC_attenuator1 = DAC_att1, DAC_attenuator2 = 10, ADC_attenuator = round(float(att_1),3))
        #put Arianna style code here

        tof = TOFExperiment('All', outerFolder, experiment,save_figs = True, title=True)
        (average_y_mag_values_last, average_y_mag_values_mid, average_y_mag_values_oct, DAC_attenuator1, DAC_attenuator2, ADC_attenuator) = tof.run(experiment.soccfg, experiment.soc)

        avg_mag_y_mid.append(average_y_mag_values_mid)
        avg_mag_y_oct.append(average_y_mag_values_oct)
        avg_mag_y_last.append(average_y_mag_values_last)

        del tof
        del experiment

    avg_mag_y_mid_DAC_atts.append(avg_mag_y_mid)
    avg_mag_y_last_DAC_atts.append(avg_mag_y_last)
    avg_mag_y_oct_DAC_atts.append(avg_mag_y_oct)

    #--------------------------plot individually------------------------
    fig, axes = plt.subplots(nrows=len(avg_mag_y_mid[0]), ncols=1, figsize=(8, 12)) # Adjust figsize
    for i in range(len(avg_mag_y_mid[0])):
        axes[i].plot(sweep_ADC_attenuator, [avg_mag_y_mid[l][i] for l in range(len(avg_mag_y_mid))], label='Mag')
        axes[i].legend(loc='best', title='Signal')
        axes[i].set_title(f'Q {i+1}')
        axes[i].set_xlabel('ADC Attenuator value')
        axes[i].set_ylabel('a.u')

    plt.suptitle(f'Average Magnitude for Middle 15 Points in TOF DAC_Att_1:{DAC_attenuator1} DAC_Att_2:{DAC_attenuator2} ADC_Att:{ADC_attenuator}') #Overall title
    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) #Adjust layout to prevent overlap
    plt.savefig(outerFolder + f'ADC_att_TOF_Middle_magnitude_plots_DAC_Att1_{DAC_attenuator1}_DAC_Att2_{DAC_attenuator2}_ADC_Att_{ADC_attenuator}.png')

    filename = outerFolder + f'ADC_att_TOF_Middle_magnitude_data_DAC_Att1_{DAC_attenuator1}_DAC_Att2_{DAC_attenuator2}_ADC_Att_{ADC_attenuator}.h5'  # Use a descriptive filename
    with h5py.File(filename, 'w') as hf:
        hf.create_dataset('sweep_ADC_attenuator', data=sweep_ADC_attenuator)
        hf.create_dataset('avg_mag_y', data=avg_mag_y_mid)
        hf.create_dataset('DAC_attenuator1', data=DAC_attenuator1)
        hf.create_dataset('DAC_attenuator2', data=DAC_attenuator2)
        hf.create_dataset('ADC_attenuator', data=ADC_attenuator)

    fig, axes = plt.subplots(nrows=len(avg_mag_y_last[0]), ncols=1, figsize=(8, 12)) # Adjust figsize
    for i in range(len(avg_mag_y_last[0])):
        axes[i].plot(sweep_ADC_attenuator, [avg_mag_y_last[l][i] for l in range(len(avg_mag_y_last))], label='Mag')
        axes[i].legend(loc='best', title='Signal')
        axes[i].set_title(f'Q {i+1}')
        axes[i].set_xlabel('ADC Attenuator value')
        axes[i].set_ylabel('a.u')

    plt.suptitle(f'Average Magnitude for Last 15 Points in TOF  DAC_Att_1:{DAC_attenuator1} DAC_Att_2:{DAC_attenuator2} ADC_Att:{ADC_attenuator}') #Overall title
    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) #Adjust layout to prevent overlap
    plt.savefig(outerFolder + f'ADC_att_TOF_Last_magnitude_plots_DAC_Att1_{DAC_attenuator1}_DAC_Att2_{DAC_attenuator2}_ADC_Att_{ADC_attenuator}.png')

    filename = outerFolder + 'ADC_att_TOF_last_magnitude_data_DAC_Att1_{DAC_attenuator1}_DAC_Att2_{DAC_attenuator2}_ADC_Att_{ADC_attenuator}.h5'  # Use a descriptive filename
    with h5py.File(filename, 'w') as hf:
        hf.create_dataset('sweep_ADC_attenuator', data=sweep_ADC_attenuator)
        hf.create_dataset('avg_mag_y', data=avg_mag_y_last)
        hf.create_dataset('DAC_attenuator1', data=DAC_attenuator1)
        hf.create_dataset('DAC_attenuator2', data=DAC_attenuator2)
        hf.create_dataset('ADC_attenuator', data=ADC_attenuator)

    fig, axes = plt.subplots(nrows=len(avg_mag_y_oct[0]), ncols=1, figsize=(8, 12)) # Adjust figsize
    for i in range(len(avg_mag_y_oct[0])):
        axes[i].plot(sweep_ADC_attenuator, [avg_mag_y_oct[l][i] for l in range(len(avg_mag_y_oct))], label='Mag')
        axes[i].legend(loc='best', title='Signal')
        axes[i].set_title(f'Q {i+1}')
        axes[i].set_xlabel('ADC Attenuator value')
        axes[i].set_ylabel('a.u')

    plt.suptitle(f'Average Magnitude for First Octant 15 Points in TOF DAC_Att_1:{DAC_attenuator1} DAC_Att_2:{DAC_attenuator2} ADC_Att:{ADC_attenuator}') #Overall title
    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) #Adjust layout to prevent overlap
    plt.savefig(outerFolder + f'ADC_att_TOF_first_oct_magnitude_plots_DAC_Att1_{DAC_attenuator1}_DAC_Att2_{DAC_attenuator2}_ADC_Att_{ADC_attenuator}.png')

    filename = outerFolder + f'ADC_att_TOF_first_oct_magnitude_data_DAC_Att1_{DAC_attenuator1}_DAC_Att2_{DAC_attenuator2}_ADC_Att_{ADC_attenuator}.h5'  # Use a descriptive filename
    with (h5py.File(filename, 'w') as hf):
        hf.create_dataset('sweep_ADC_attenuator', data=sweep_ADC_attenuator)
        hf.create_dataset('avg_mag_y', data=avg_mag_y_oct)
        hf.create_dataset('DAC_attenuator1', data=DAC_attenuator1)
        hf.create_dataset('DAC_attenuator2', data=DAC_attenuator2)
        hf.create_dataset('ADC_attenuator', data=ADC_attenuator)

#---------------------------plot ratio----------------------------
#---------------------------middle 7 points avg----------------------------
fig, axes = plt.subplots(nrows=len(avg_mag_y_mid_DAC_atts[0][0]), ncols=1, figsize=(8, 12))
dacBefore_i = avg_mag_y_mid_DAC_atts[0]
dacAfter_i = avg_mag_y_mid_DAC_atts[1]
ratio = []
for i in range(len(dacAfter_i)):
    row_result = []
    for j in range(len(dacAfter_i[i])):
        row_result.append(dacBefore_i[i][j] / dacAfter_i[i][j])
    ratio.append(row_result)
for i in range(len(avg_mag_y_mid_DAC_atts[0][0])):
    axes[i].plot(sweep_ADC_attenuator, [ratio[l][i] for l in range(len(ratio))], label='Mag')
    axes[i].legend(loc='best', title='Signal')
    axes[i].set_title(f'Q {i + 1}')
    axes[i].set_xlabel('ADC Attenuator value')
    axes[i].set_ylabel(f'DAC ({DAC_att_1[0]}, 10)/({DAC_att_1[1]}, 10)')

plt.suptitle(f'Ratio of Magnitudes for Middle 15 Points in TOF')
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(outerFolder + f'ADC_att_TOF_Middle_magnitude_plots.png')

filename = outerFolder + f'ADC_att_TOF_middle_magnitude_data_ratio.h5'
with (h5py.File(filename, 'w') as hf):
    hf.create_dataset('sweep_ADC_attenuator', data=sweep_ADC_attenuator)
    hf.create_dataset('ratio', data=ratio)

#-------------------------First Octant 7 points avg-----------------------------
fig, axes = plt.subplots(nrows=len(avg_mag_y_oct_DAC_atts[0][0]), ncols=1, figsize=(8, 12))
dacBefore_i = avg_mag_y_oct_DAC_atts[0]
dacAfter_i = avg_mag_y_oct_DAC_atts[1]
ratio = []
for i in range(len(dacAfter_i)):
    row_result = []
    for j in range(len(dacAfter_i[i])):
        row_result.append(dacBefore_i[i][j] / dacAfter_i[i][j])
    ratio.append(row_result)
for i in range(len(avg_mag_y_oct_DAC_atts[0][0])):
    axes[i].plot(sweep_ADC_attenuator, [ratio[l][i] for l in range(len(ratio))], label='Mag')
    axes[i].legend(loc='best', title='Signal')
    axes[i].set_title(f'Q {i + 1}')
    axes[i].set_xlabel('ADC Attenuator value')
    axes[i].set_ylabel(f'DAC ({DAC_att_1[0]}, 10)/({DAC_att_1[1]}, 10)')

plt.suptitle(f'Ratio of Magnitudes for fist octant 15 Points in TOF')
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(outerFolder + f'ADC_att_TOF_first_oct_magnitude_plots.png')

filename = outerFolder + f'ADC_att_TOF_first_oct_magnitude_data_ratio.h5'
with (h5py.File(filename, 'w') as hf):
    hf.create_dataset('sweep_ADC_attenuator', data=sweep_ADC_attenuator)
    hf.create_dataset('ratio', data=ratio)

#----------------------------Last points 7 points avg---------------------------------
fig, axes = plt.subplots(nrows=len(avg_mag_y_last_DAC_atts[0][0]), ncols=1, figsize=(8, 12))
dacBefore_i = avg_mag_y_last_DAC_atts[0]
dacAfter_i = avg_mag_y_last_DAC_atts[1]
ratio = []
for i in range(len(dacAfter_i)):
    row_result = []
    for j in range(len(dacAfter_i[i])):
        row_result.append(dacBefore_i[i][j] / dacAfter_i[i][j])
    ratio.append(row_result)
for i in range(len(avg_mag_y_last_DAC_atts[0][0])):
    axes[i].plot(sweep_ADC_attenuator, [ratio[l][i] for l in range(len(ratio))], label='Mag')
    axes[i].legend(loc='best', title='Signal')
    axes[i].set_title(f'Q {i + 1}')
    axes[i].set_xlabel('ADC Attenuator value')
    axes[i].set_ylabel(f'DAC ({DAC_att_1[0]}, 10)/({DAC_att_1[1]}, 10)')

plt.suptitle(f'Ratio of Magnitudes for Last 15 Points in TOF')
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(outerFolder + f'ADC_att_TOF_last_magnitude_plots.png')

filename = outerFolder + f'ADC_att_TOF_last_magnitude_data_ratio.h5'
with (h5py.File(filename, 'w') as hf):
    hf.create_dataset('sweep_ADC_attenuator', data=sweep_ADC_attenuator)
    hf.create_dataset('ratio', data=ratio)