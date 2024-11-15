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

    def save_to_h5(self, data_type):
        self.outerFolder_expt = self.outerFolder_expt + "/Data_h5" + f"{data_type}_ge" + "/"
        self.create_folder_if_not_exists(self.outerFolder_expt)
        formatted_datetime =  datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        h5_filename = self.outerFolder_expt + f"{formatted_datetime}_" + f"{data_type}_results_batch_{self.batch_num}_" + f"Num_per_batch{self.save_r}.h5"
        with h5py.File(h5_filename, 'w') as f:
            f.attrs['datestr']=formatted_datetime
            f.attrs[f'{data_type}_results_batch']=self.batch_num
            f.attrs['num_per_batch'] = self.save_r
            for QubitIndex, data in self.data.items():
                # Create a group for each qubit
                group = f.create_group(f'Q{QubitIndex + 1}')

                # make the dataset even if none, if indeed it is none, just save an empty array
                def create_dataset(name, value):
                    if value is not None:
                        group.create_dataset(name, data=np.array(value, dtype='float32' if isinstance(value, (
                        list, np.ndarray)) and all(isinstance(i, (int, float)) for i in value) else 'S'))
                    else:
                        # Save an empty array or placeholder if the value is None
                        group.create_dataset(name, data=np.array([]))

                # Save everything with checks for 'None'. save the right datasets depending on "type"
                if 'Res' in data_type:
                    create_dataset("Dates", data.get('Dates'))
                    create_dataset("freq_pts", data.get('freq_pts'))
                    create_dataset("freq_center", data.get('freq_center'))
                    create_dataset("Amps", data.get('Amps'))
                    create_dataset("Found Freqs", data.get('Found Freqs'))
                    create_dataset("Round Num", data.get('Round Num'))
                    create_dataset("Batch Num", data.get('Batch Num'))

                if 'QSpec' in data_type:
                    create_dataset("Dates", data.get('Dates'))
                    create_dataset("I", data.get('I'))
                    create_dataset("Q", data.get('Q'))
                    create_dataset("Frequencies", data.get('Frequencies'))
                    create_dataset("I Fit", data.get('I Fit'))
                    create_dataset("Q Fit", data.get('Q Fit'))
                    create_dataset("Round Num", data.get('Round Num'))
                    create_dataset("Batch Num", data.get('Batch Num'))

                if 'Rabi' in data_type:
                    create_dataset("Dates", data.get('Dates'))
                    create_dataset("I", data.get('I'))
                    create_dataset("Q", data.get('Q'))
                    create_dataset("Gains", data.get('Gains'))
                    create_dataset("Fit", data.get('Fit'))
                    create_dataset("Round Num", data.get('Round Num'))
                    create_dataset("Batch Num", data.get('Batch Num'))

                if 'T1' in data_type:
                    create_dataset("T1", data.get('T1'))
                    create_dataset("Errors", data.get('Errors'))
                    create_dataset("Dates", data.get('Dates'))
                    create_dataset("I", data.get('I'))
                    create_dataset("Q", data.get('Q'))
                    create_dataset("Delay Times", data.get('Delay Times'))
                    create_dataset("Fit", data.get('Fit'))
                    create_dataset("Round Num", data.get('Round Num'))
                    create_dataset("Batch Num", data.get('Batch Num'))

                if 'T2' in data_type:
                    create_dataset("T2", data.get('T2'))
                    create_dataset("Errors", data.get('Errors'))
                    create_dataset("Dates", data.get('Dates'))
                    create_dataset("I", data.get('I'))
                    create_dataset("Q", data.get('Q'))
                    create_dataset("Delay Times", data.get('Delay Times'))
                    create_dataset("Fit", data.get('Fit'))
                    create_dataset("Round Num", data.get('Round Num'))
                    create_dataset("Batch Num", data.get('Batch Num'))
