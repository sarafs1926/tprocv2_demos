import datetime
import numpy as np
import h5py

class Save:
    def __init__(self, outerFolder, data, batch_num, save_r):
        self.outerFolder_expt = outerFolder
        self.data = data
        self.batch_num = batch_num
        self.save_r = save_r

    def create_folder_if_not_exists(self, folder):
        """Creates a folder at the given path if it doesn't already exist."""
        import os
        if not os.path.exists(folder):
            os.makedirs(folder)

    def save_to_h5(self, t_type):
        self.outerFolder_expt = self.outerFolder_expt + "/" + f"{t_type}_ge" + "/"
        self.create_folder_if_not_exists(self.outerFolder_expt)
        formatted_datetime =  datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        h5_filename = self.outerFolder_expt + f"{formatted_datetime}_" + f"{t_type}_results_batch_{self.batch_num}_" + f"Num_per_batch{self.save_r}.h5"
        with h5py.File(h5_filename, 'w') as f:
            f.attrs['datestr']=formatted_datetime
            f.attrs[f'{t_type}_results_batch']=self.batch_num
            f.attrs['num_per_batch'] = self.save_r
            for QubitIndex, data in self.data.items():
                # Create a group for each qubit
                group = f.create_group(f'Q{QubitIndex + 1}')

                # Save estimates, errors, and dates
                if '1' in t_type:
                    group.create_dataset("T1", data=np.array(data['T1'], dtype=np.float32))
                if '2' in t_type:
                    group.create_dataset("T2", data=np.array(data['T2'], dtype=np.float32))
                group.create_dataset("Errors", data=np.array(data['Errors'], dtype=np.float32))
                group.create_dataset("Dates", data=np.array(data['Dates'], dtype='S'))