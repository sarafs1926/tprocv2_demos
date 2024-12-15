import h5py
import json
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import MaxNLocator

def create_folder_if_not_exists(folder):
    """Creates a folder at the given path if it doesn't already exist."""
    if not os.path.exists(folder):
        os.makedirs(folder)

def load_from_h5(filename):
    with h5py.File(filename, 'r') as hf:
        # Load attributes
        t1_vals = json.loads(hf.attrs['t1_vals'])
        t1_errs = json.loads(hf.attrs['t1_errs'])
        t1_std_values = json.loads(hf.attrs['t1_std_values'])
        t1_mean_values = json.loads(hf.attrs['t1_mean_values'])
        run_number = hf.attrs['run_number']
        run_notes = hf.attrs['run_notes']
        last_date = hf.attrs['last_date']

        t2r_vals = json.loads(hf.attrs['t2r_vals'])
        t2r_errs = json.loads(hf.attrs['t2r_errs'])
        t2r_std_values = json.loads(hf.attrs['t2r_std_values'])
        t2r_mean_values = json.loads(hf.attrs['t2r_mean_values'])

        t2e_vals = json.loads(hf.attrs['t2e_vals'])
        t2e_errs = json.loads(hf.attrs['t2e_errs'])
        t2e_std_values = json.loads(hf.attrs['t2e_std_values'])
        t2e_mean_values = json.loads(hf.attrs['t2e_mean_values'])

    return {
        't1_vals': t1_vals,
        't1_errs': t1_errs,
        't1_std_values': t1_std_values,
        't1_mean_values': t1_mean_values,
        'run_number': run_number,
        'run_notes': run_notes,
        'last_date': last_date,
        't2r_vals': t2r_vals,
        't2r_errs': t2r_errs,
        't2r_std_values': t2r_std_values,
        't2r_mean_values': t2r_mean_values,
        't2e_vals': t2e_vals,
        't2e_errs': t2e_errs,
        't2e_std_values': t2e_std_values,
        't2e_mean_values': t2e_mean_values,
    }

run_number_list = [1,2]  # Example: if you have more runs, e.g., [1, 2, 3]

if len(run_number_list) > 1:
    t1_data = {}
    t1_err = {}
    t2r_data = {}
    t2r_err = {}
    t2e_data = {}
    t2e_err = {}

    qubit_list = None

    #load data for each run
    for r in run_number_list:
        run_stats_folder = f"run_stats/run{r}/"
        filename = run_stats_folder + 'experiment_data.h5'
        loaded_data = load_from_h5(filename)

        t1_means = loaded_data['t1_mean_values']
        t1_stds = loaded_data['t1_std_values']
        t2r_means = loaded_data['t2r_mean_values']
        t2r_stds = loaded_data['t2r_std_values']
        t2e_means = loaded_data['t2e_mean_values']
        t2e_stds = loaded_data['t2e_std_values']
        run_notes = loaded_data['run_notes']

        #on the first run do this to make all of the lists and such
        if qubit_list is None:
            qubit_list = list(t1_means.keys())
            for qb in qubit_list:
                t1_data[qb] = []
                t1_err[qb] = []
                t2r_data[qb] = []
                t2r_err[qb] = []
                t2e_data[qb] = []
                t2e_err[qb] = []

        #get the data for this run and append to list for saving and plotting
        for qb in qubit_list:
            t1_data[qb].append(t1_means[qb])
            t1_err[qb].append(t1_stds[qb])
            t2r_data[qb].append(t2r_means[qb])
            t2r_err[qb].append(t2r_stds[qb])
            t2e_data[qb].append(t2e_means[qb])
            t2e_err[qb].append(t2e_stds[qb])

    #make one figure with subplots, one per qubit
    fig, axes = plt.subplots(nrows=len(qubit_list), ncols=1, sharex=True, figsize=(6, 4*len(qubit_list)))
    if len(qubit_list) == 1:
        #if there's only one qubit, axes isnt a list, so wrap it
        axes = [axes]

    x = run_number_list
    for i, qb in enumerate(qubit_list):
        ax = axes[i]
        ax.errorbar(x, t1_data[qb], yerr=t1_err[qb], fmt='o-', label='T1', capsize=3)
        ax.errorbar(x, t2r_data[qb], yerr=t2r_err[qb], fmt='o-', label='T2R', capsize=3)
        ax.errorbar(x, t2e_data[qb], yerr=t2e_err[qb], fmt='o-', label='T2E', capsize=3)

        ax.set_ylabel('Time (µs)')
        ax.set_title(qb)
        ax.legend()

        # -1 grabs the most recent run, because we only have one run and note right now. update this next run
        note_x = x[-1]


        #find the max err bar height at the last run, from (T1, T2R, T2E)
        highest_point = max(
            t1_data[qb][-1] + t1_err[qb][-1],
            t2r_data[qb][-1] + t2r_err[qb][-1],
            t2e_data[qb][-1] + t2e_err[qb][-1])

        #extend top of y-axis limit so the note is fully visible above the highest error bar.
        current_ylim = ax.get_ylim()
        new_upper_limit = max(current_ylim[1], highest_point * 1.3) #add 30% padding above, change if you want
        ax.set_ylim(current_ylim[0], new_upper_limit)

        #put the note above the
        ax.annotate(
            run_notes,
            xy=(note_x, highest_point),
            xytext=(0, 10), #offset by some amount, 10 right now, to go a little higher than the top of err bar
            textcoords='offset points',
            ha='center',
            va='bottom',
            bbox=dict(boxstyle='round', facecolor='lightblue', edgecolor='black', alpha=0.7)
            # arrowprops=dict(facecolor='black', shrink=0.05)  #if you want a fancy arrow
        )



    #set the bottom plot's x-axis label (shared for all) only do it for the bottom
    axes[-1].set_xlabel('Run Number')
    axes[-1].xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.tight_layout()
    analysis_folder = "/data/QICK_data/6transmon_run5/benchmark_analysis_plots/"
    create_folder_if_not_exists(analysis_folder)
    plt.savefig(analysis_folder + 'compare_runs.pdf', dpi=500)

elif len(run_number_list) == 1:
    #only 1 data point per qubit
    r = run_number_list[0]
    run_stats_folder = f"run_stats/run{r}/"
    filename = run_stats_folder + 'experiment_data.h5'
    loaded_data = load_from_h5(filename)

    t1_means = loaded_data['t1_mean_values']
    t1_stds = loaded_data['t1_std_values']
    t2r_means = loaded_data['t2r_mean_values']
    t2r_stds = loaded_data['t2r_std_values']
    t2e_means = loaded_data['t2e_mean_values']
    t2e_stds = loaded_data['t2e_std_values']
    run_notes = loaded_data['run_notes']

    qubit_list = list(t1_means.keys())
    x = [r]

    #one figure w all qubits
    fig, axes = plt.subplots(nrows=len(qubit_list), ncols=1, sharex=True, figsize=(6, 4*len(qubit_list)))
    if len(qubit_list) == 1:
        axes = [axes]

    for i, qb in enumerate(qubit_list):
        ax = axes[i]
        ax.errorbar(x, [t1_means[qb]], yerr=[t1_stds[qb]], fmt='o-', label='T1', capsize=3)
        ax.errorbar(x, [t2r_means[qb]], yerr=[t2r_stds[qb]], fmt='o-', label='T2R', capsize=3)
        ax.errorbar(x, [t2e_means[qb]], yerr=[t2e_stds[qb]], fmt='o-', label='T2E', capsize=3)

        ax.set_ylabel('Time (µs)')
        ax.set_title(qb)
        ax.legend()

        # -1 grabs the most recent run, because we only have one run and note right now. update this next run
        note_x = x[-1]
        # grab the T1 value and put the note above that
        note_y = t1_means[qb][-1]

        #add the note and make the background colored and pretty
        ax.annotate(
            run_notes,
            xy=(note_x, note_y),
            xytext=(0, 20),  # offset text 20 points above the data point
            textcoords='offset points',
            ha='center',
            va='bottom',
            bbox=dict(boxstyle='round', facecolor='lightblue', edgecolor='black', alpha=0.7),
            #arrowprops=dict(facecolor='black', shrink=0.05)
        )

    #set the bottom plot's x-axis label (shared for all) only do it for the bottom
    axes[-1].set_xlabel('Run Number')
    axes[-1].xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.tight_layout()
    analysis_folder = "/data/QICK_data/6transmon_run5/benchmark_analysis_plots/"
    create_folder_if_not_exists(analysis_folder)
    plt.savefig(analysis_folder + 'compare_runs.pdf', dpi=500)