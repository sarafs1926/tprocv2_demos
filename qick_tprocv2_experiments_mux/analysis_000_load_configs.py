from section_008_save_data_to_h5 import Data_H5

class LoadConfigs:
    def __init__(self, outerFolder):
        self.outerFolder = outerFolder

    def run(self):
        loader_config_instance = Data_H5(self.outerFolder)
        sys_config = loader_config_instance.load_config('sys_config.h5')
        del loader_config_instance

        loader_config_instance = Data_H5(self.outerFolder)
        exp_config = loader_config_instance.load_config('expt_cfg.h5')
        del loader_config_instance

        return sys_config, exp_config