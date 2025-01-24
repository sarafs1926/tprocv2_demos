import h5py
import os

class UpdateNote:
    def __init__(self, run_number, new_run_notes):
        self.run_number = run_number
        self.new_run_notes = new_run_notes

    def run(self):
        run_stats_folder = f"run_stats/run{self.run_number}/"
        filename = os.path.join(run_stats_folder, 'experiment_data.h5')

        with h5py.File(filename, 'r+') as f:

            #if there is no run notes it will throw err
            if 'run_notes' in f.attrs:
                f.attrs['run_notes'] = self.new_run_notes
            else:
                print("Error: 'run_notes' not found in the file.")