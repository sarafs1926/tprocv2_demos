from section_005_single_shot_ge import SingleShot
from system_config import *
import numpy as np
import matplotlib.pyplot as plt
import os

#output_folder = "/home/quietuser/Documents/Fit_Tests"
import datetime

def create_folder_if_not_exists(folder_path):
    """Creates a folder at the given path if it doesn't already exist."""
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

# Where do you want to save data
prefix = str(datetime.date.today())
output_folder = "/data/QICK_data/6transmon_run4a/" + prefix + "/SingleShot_Test/"
create_folder_if_not_exists(output_folder)

n = 2
j = 0
Qs=[0,1,2,3,4,5]
lengs = np.linspace(0.5,50,100)
max_fid = {i: [] for i in range(0, 6)}
while j < n:
    j += 1
    for QubitIndex in Qs:
        #QubitIndex = 0
        fids = []
        for leng in lengs:
            # ------------------------Single Shot-------------------------
            ss = SingleShot(QubitIndex, outerFolder, j, round(leng, 3))
            fid = ss.run(soccfg, soc)
            fids.append(fid)
            print(fids)
        plt.figure()
        maximum_index = np.argmax(fids)  # Find the index of the maximum value
        maximum_value = fids[maximum_index]  # Get the maximum value
        max_length = lengs[maximum_index]  # Corresponding x value

        max_fid[QubitIndex].append(maximum_value)  # Store the maximum value if needed
        plt.axvline(max_length, color='orange', linestyle='--', linewidth=2)  # Vertical line at max length
        plt.plot(lengs, fids)
        plt.xlabel('Readout and Pulse Length')
        plt.ylabel('Fidelity')
        plt.title(f'Maximum Fidelity at Length {max_length}: {maximum_value}')  # Update title
        plt.savefig(os.path.join(output_folder, 'fidelity' + f"{QubitIndex}" + '.png'), dpi=300)
        #plt.show()
        plt.close()
