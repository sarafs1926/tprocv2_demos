import datetime
import numpy as np
import h5py

class Save:
    def __init__(self, outerFolder, t1_data, batch_num, save_r):
        self.outerFolder_expt = outerFolder + "/" + "T1_ge" + "/"
        self.t1_data = t1_data
        self.batch_num = batch_num
        self.save_r = save_r

    def save_to_h5(self):
        formatted_datetime =  datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        h5_filename = self.outerFolder_expt + f"{formatted_datetime}_" + f"T1_results_batch_{self.batch_num}_" + f"Num_per_batch{self.save_r}.h5"
        with h5py.File(h5_filename, 'w') as f:
            f.attrs['datestr']=formatted_datetime
            f.attrs['T1_results_batch']=self.batch_num
            f.attrs['num_per_batch'] = self.save_r
            for QubitIndex, data in self.t1_data.items():
                # Create a group for each qubit
                group = f.create_group(f'Q{QubitIndex + 1}')

                # Save T1 estimates, errors, and dates
                group.create_dataset("T1", data=data['T1'])
                group.create_dataset("Errors", data=data['Errors'])
                group.create_dataset("Dates", data=data['Dates'], dtype='S')