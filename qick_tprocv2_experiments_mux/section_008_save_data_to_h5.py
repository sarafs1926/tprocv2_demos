import datetime
import numpy as np
import h5py

class Data_H5:
    def __init__(self, outerFolder, data = None, batch_num = 0, save_r = 0):
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
        self.outerFolder_expt = self.outerFolder_expt + "/Data_h5/" + f"{data_type}_ge" + "/"
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

                def create_dataset(name, value):
                    if value is not None:
                        if name == "Dates":
                            if isinstance(value, list) and 'None' in str(
                                    value[0]):  # Check for single 'None' string list
                                value = np.array([])
                            else:
                                try:
                                    value = np.array([dt.isoformat() for dt in value], dtype='S')
                                except (TypeError, AttributeError) as e:
                                    print(f"Warning: Could not convert dates for '{name}'. Error: {e}")
                                    value = np.array([str(dt) for dt in value], dtype='S')
                        elif isinstance(value, list) and 'None' in str(value[0]) :  # Check for single 'None' string list
                            value = np.array([])  # Save as empty array

                        else:
                            try:
                                value = np.array(value, dtype=np.float64)
                            except ValueError:
                                print(f"Warning: Could not convert data for '{name}' to float. Saving as string.")
                                value = np.array(value, dtype='S')
                        group.create_dataset(name, data=value)
                    else:
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

                elif 'QSpec' in data_type:
                    create_dataset("Dates", data.get('Dates'))
                    create_dataset("I", data.get('I'))
                    create_dataset("Q", data.get('Q'))
                    create_dataset("Frequencies", data.get('Frequencies'))
                    create_dataset("I Fit", data.get('I Fit'))
                    create_dataset("Q Fit", data.get('Q Fit'))
                    create_dataset("Round Num", data.get('Round Num'))
                    create_dataset("Batch Num", data.get('Batch Num'))

                elif 'Rabi' in data_type:
                    create_dataset("Dates", data.get('Dates'))
                    create_dataset("I", data.get('I'))
                    create_dataset("Q", data.get('Q'))
                    create_dataset("Gains", data.get('Gains'))
                    create_dataset("Fit", data.get('Fit'))
                    create_dataset("Round Num", data.get('Round Num'))
                    create_dataset("Batch Num", data.get('Batch Num'))

                elif 'T1' in data_type:
                    create_dataset("T1", data.get('T1'))
                    create_dataset("Errors", data.get('Errors'))
                    create_dataset("Dates", data.get('Dates'))
                    create_dataset("I", data.get('I'))
                    create_dataset("Q", data.get('Q'))
                    create_dataset("Delay Times", data.get('Delay Times'))
                    create_dataset("Fit", data.get('Fit'))
                    create_dataset("Round Num", data.get('Round Num'))
                    create_dataset("Batch Num", data.get('Batch Num'))

                elif 'T2' in data_type:
                    create_dataset("T2", data.get('T2'))
                    create_dataset("Errors", data.get('Errors'))
                    create_dataset("Dates", data.get('Dates'))
                    create_dataset("I", data.get('I'))
                    create_dataset("Q", data.get('Q'))
                    create_dataset("Delay Times", data.get('Delay Times'))
                    create_dataset("Fit", data.get('Fit'))
                    create_dataset("Round Num", data.get('Round Num'))
                    create_dataset("Batch Num", data.get('Batch Num'))

                elif 'T2E' in data_type:
                    create_dataset("T2E", data.get('T2E'))
                    create_dataset("Errors", data.get('Errors'))
                    create_dataset("Dates", data.get('Dates'))
                    create_dataset("I", data.get('I'))
                    create_dataset("Q", data.get('Q'))
                    create_dataset("Delay Times", data.get('Delay Times'))
                    create_dataset("Fit", data.get('Fit'))
                    create_dataset("Round Num", data.get('Round Num'))
                    create_dataset("Batch Num", data.get('Batch Num'))
                else:
                    print('Data type not supported, please pass one of the following: Res, QSpec, Rabi, T1, T2, T2E')

        #self.print_h5_contents(h5_filename)
        #print(h5_filename)

    def print_h5_contents(self, filename):
        try:
            with h5py.File(filename, 'r') as f:
                print("File Attributes:")
                for key, value in f.attrs.items():
                    print(f"  {key}: {value}")
                print("\nDataset Contents:")
                for key in f.keys():
                    group = f[key]
                    print(f"\nGroup: {key}")
                    for dataset_name in group.keys():
                        dataset = group[dataset_name]
                        print(f"  Dataset: {dataset_name}")
                        try:  # Handle various datatypes gracefully
                            data = dataset[()]
                            if isinstance(data,
                                          np.ndarray) and data.size > 100:  # Check array size for printing only sample
                                print(f"    Data (sample): {data[:50]} ... (truncated)")
                            else:
                                print(f"    Data: {data}")
                        except:
                            print("    Data: Could not print data")  # Catch potential errors printing complex objects

        except FileNotFoundError:
            print(f"Error: File '{filename}' not found.")
        except Exception as e:
            print(f"An error occurred: {e}")

    def load_from_h5(self, data_type,  save_r=1):  # Added save_r as parameter.
        """Loads data from an HDF5 file into specified dictionary format."""

        data = {data_type: {}}  # Initialize the main dictionary with the data_type.

        with h5py.File(self.outerFolder_expt, 'r') as f:
            for qubit_group in f.keys():
                qubit_index = int(qubit_group[1:]) - 1
                qubit_data = {}
                group = f[qubit_group]

                for dataset_name in group.keys():
                    # Attempt to map HDF5 keys to the target dictionaries' keys.
                    if data_type == 'Res':
                        target_keys = {'Dates': 'Dates', 'freq_pts': 'freq_pts', 'freq_center': 'freq_center',
                                       'Amps': 'Amps', 'Found Freqs': 'Found Freqs', 'Round Num': 'Round Num',
                                       'Batch Num': 'Batch Num'}
                    elif data_type == 'QSpec':
                        target_keys = {'Dates': 'Dates', 'I': 'I', 'Q': 'Q', 'Frequencies': 'Frequencies',
                                       'I Fit': 'I Fit', 'Q Fit': 'Q Fit', 'Round Num': 'Round Num',
                                       'Batch Num': 'Batch Num'}
                    elif data_type == 'Rabi':
                        target_keys = {'Dates': 'Dates', 'I': 'I', 'Q': 'Q', 'Gains': 'Gains', 'Fit': 'Fit',
                                       'Round Num': 'Round Num', 'Batch Num': 'Batch Num'}
                    elif data_type == 'T1':
                        target_keys = {'T1': 'T1', 'Errors': 'Errors', 'Dates': 'Dates', 'I': 'I', 'Q': 'Q',
                                       'Delay Times': 'Delay Times', 'Fit': 'Fit', 'Round Num': 'Round Num',
                                       'Batch Num': 'Batch Num'}
                    elif data_type == 'T2':
                        target_keys = {'T2': 'T2', 'Errors': 'Errors', 'Dates': 'Dates', 'I': 'I', 'Q': 'Q',
                                       'Delay Times': 'Delay Times', 'Fit': 'Fit', 'Round Num': 'Round Num',
                                       'Batch Num': 'Batch Num'}
                    elif data_type == 'T2E':
                        target_keys = {'T2E': 'T2E', 'Errors': 'Errors', 'Dates': 'Dates', 'I': 'I', 'Q': 'Q',
                                       'Delay Times': 'Delay Times', 'Fit': 'Fit', 'Round Num': 'Round Num',
                                       'Batch Num': 'Batch Num'}
                    else:
                        raise ValueError(f"Unsupported data_type: {data_type}")

                    try:
                        mapped_key = target_keys[dataset_name]  # Map HDF5 key to target key.
                        qubit_data[mapped_key] = [group[dataset_name][
                                                      ()]] * save_r  # Expand to match the desired length.

                    except KeyError:
                        print(
                            f"Warning: Key '{dataset_name}' not found in target dictionary for data_type '{data_type}'. Skipping.")
                        pass

                data[data_type][qubit_index] = qubit_data

        return data

    def convert_for_hdf5(self, data):
        if isinstance(data, dict):
            return {k: self.convert_for_hdf5(v) for k, v in data.items()}
        elif isinstance(data, list):
            return [self.convert_for_hdf5(x) for x in data]
        elif isinstance(data, np.ndarray):
            if data.dtype == 'O':
                try:
                    # Attempt to convert to string type. If it fails, raise an error
                    return np.array(data, dtype='U')
                except:
                    raise TypeError(f"Cannot convert NumPy object array to HDF5-compatible type: {data}")
            else:
                return data  # Already a suitable type
        elif isinstance(data, object):  # Added to handle objects
            try:
                return str(data)  # Convert to string
            except:
                raise TypeError(f"Cannot convert object to HDF5-compatible type: {data}")
        else:
            return data

    def save_config(self, sys_config, expt_cfg):
        sys_config = self.convert_for_hdf5(sys_config)
        expt_cfg = self.convert_for_hdf5(expt_cfg)

        outerFolder_expt = self.outerFolder_expt + "/Data_h5/config/"
        self.create_folder_if_not_exists(outerFolder_expt)

        with h5py.File(outerFolder_expt + "sys_config.h5" , "w") as hf:
            for key, value in sys_config.items():
                hf.create_dataset(key, data=value)

        with h5py.File(outerFolder_expt + "expt_cfg.h5" , "w") as hf:
            for key, value in expt_cfg.items():
                hf.create_dataset(key, data=value)

    def load_config(self):
        filename = self.outerFolder_expt + "/Data_h5/config/config.h5"
        loaded_config = {}
        with h5py.File(filename, "r") as hf:
            for key in hf.keys():
                loaded_config[key] = hf[key][()]  # The [()] gets the data as a NumPy array or scalar
        return loaded_config


