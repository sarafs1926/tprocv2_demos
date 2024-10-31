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
from expt_config import *
from system_config import *

from qualang_tools.plot import Fit
import pprint as pp

# import single shot information for g-e calibration
from SingleShot_gef import SingleShotProgram_g, SingleShotProgram_e, SingleShotProgram_f

# ----- Experiment configurations ----- #
expt_name = "Readout_Optimization_gef"
QubitIndex = QUBIT_INDEX
Qubit = 'Q' + str(QubitIndex)

exp_cfg = add_qubit_experiment(expt_cfg, expt_name, QubitIndex)
q_config = all_qubit_state(system_config)
config = {**q_config['Q' + str(QubitIndex)], **exp_cfg}
print(config)

##################
# Define Program #
##################

I_g_array = []
Q_g_array = []
I_e_array = []
Q_e_array = []
I_f_array = []
Q_f_array = []
for freqs in range(config['freq_steps']):
    I_g_data = []
    Q_g_data = []
    I_e_data = []
    Q_e_data = []
    I_f_data = []
    Q_f_data = []
    for gains in range(config['gain_steps']):
        config['res_freq_ge'] = config['freq_start'] + config['freq_step_size']*freqs
        config['res_gain_ge'] = config['gain_start'] + config['gain_step_size']*gains

        ssp_g = SingleShotProgram_g(soccfg, reps=1, final_delay=config['relax_delay'], cfg=config)
        iq_list_g = ssp_g.acquire(soc, soft_avgs=1, progress=False)

        ssp_e = SingleShotProgram_e(soccfg, reps=1, final_delay=config['relax_delay'], cfg=config)
        iq_list_e = ssp_e.acquire(soc, soft_avgs=1, progress=False)

        ssp_f = SingleShotProgram_f(soccfg, reps=1, final_delay=config['relax_delay'], cfg=config)
        iq_list_f = ssp_f.acquire(soc, soft_avgs=1, progress=False)

        I_g = iq_list_g[0][0].T[0]
        Q_g = iq_list_g[0][0].T[1]
        I_e = iq_list_e[0][0].T[0]
        Q_e = iq_list_e[0][0].T[1]
        I_f = iq_list_f[0][0].T[0]
        Q_f = iq_list_f[0][0].T[1]

        I_g_data.append([I_g])
        Q_g_data.append([Q_g])
        I_e_data.append([I_e])
        Q_e_data.append([Q_e])
        I_f_data.append([I_f])
        Q_f_data.append([Q_f])

    I_g_array.append([I_g_data])
    Q_g_array.append([Q_g_data])
    I_e_array.append([I_e_data])
    Q_e_array.append([Q_e_data])
    I_f_array.append([I_f_data])
    Q_f_array.append([Q_f_data])

#####################################
# ----- Saves data to a file ----- #
#####################################

prefix = str(datetime.date.today())
exp_name = expt_name + '_Q' + str(QubitIndex) + '_' + prefix
print('Experiment name: ' + exp_name)

data_path = DATA_PATH

fname = get_next_filename(data_path, exp_name, suffix='.h5')
print('Current data file: ' + fname)

with SlabFile(data_path + '\\' + fname, 'a') as f:
    # 'a': read/write/create

    f.append('I_g_array', I_g_array)
    f.append('Q_g_array', Q_g_array)
    f.append('I_e_array', I_e_array)
    f.append('Q_e_array', Q_e_array)
    f.append('I_f_array', I_f_array)
    f.append('Q_f_array', Q_f_array)

    # formats config into file as a single line
    f.attrs['config'] = json.dumps(config)
    
data = data_path + '\\' + fname


