import sys
from qick.qick import *
from qick.rfboard import *

import time
import matplotlib.pyplot as plt

# Extended RFDC to allow setting/getting Nyquist Zone on ADC and DAC.
class eRFDC(RFDC):
    """
    Extends the xrfdc driver.
    Since operations on the RFdc tend to be slow (tens of ms), we cache the Nyquist zone and frequency.
    """
    bindto = ["xilinx.com:ip:usp_rf_data_converter:2.3",
              "xilinx.com:ip:usp_rf_data_converter:2.4",
              "xilinx.com:ip:usp_rf_data_converter:2.6"]

    def set_nyquist(self, blockid, nqz, blocktype='dac', force=False):
        # Check valid selection.
        if nqz not in [1,2]:
            raise RuntimeError("Nyquist zone must be 1 or 2")

        # Get tile and channel from id.
        tile, channel = [int(a) for a in blockid]

        # Need to update?
        #if not force and self.get_nyquist(blockid,blocktype) == nqz:
        #    return

        if blocktype == 'adc':
            self.adc_tiles[tile].blocks[channel].NyquistZone = nqz
            #self.dict['nqz'][blocktype][blockid] = nqz
        elif blocktype == 'dac':
            self.dac_tiles[tile].blocks[channel].NyquistZone = nqz
            #self.dict['nqz'][blocktype][blockid] = nqz
        else:
            raise RuntimeError("Blocktype %s not recognized" & blocktype)

    def get_nyquist(self, blockid, blocktype='dac'):
        # Get tile and channel from id.
        tile, channel = [int(a) for a in blockid]

        if blocktype == 'adc':
            return self.adc_tiles[tile].blocks[channel].NyquistZone
        elif blocktype == 'dac':
            return self.dac_tiles[tile].blocks[channel].NyquistZone
        else:
            raise RuntimeError("Blocktype %s not recognized" & blocktype)

class TestSoc(QickSoc):
    def __init__(self, bitfile=None, **kwargs):
        super().__init__(bitfile=bitfile, **kwargs)

        # Dictionary for converter to channel mapping.
        self.CH_MAP = {}
        self.CH_MAP['dac'] = {}
        for i,dac in enumerate(self.dacs):
            self.CH_MAP['dac'][i] = dac

        self.CH_MAP['adc'] = {}
        for i,adc in enumerate(self.adcs):
            self.CH_MAP['adc'][i] = adc
        
        # Map signal 
        self.map_signal_paths()

        # Signal Gnerator.
        self.generator = Generator(self, self.dacs['00']['fs'], self.axis_signal_gen_v6_c_0, self.axis_signal_gen_v6_0, self.axis_switch_v1_0)

        # Buffer.
        self.buffer = Buffer(self, self.adcs['10']['fs'], self.axis_switch_1, self.mr_buffer_et_0, self.axi_dma_0)

        # SPI used for Filter.
        if 'filter_spi' in self.ip_dict.keys():
            self.filter_spi.config(lsb="msb")

            # Programmable filter.
            self.filter = prog_filter(self.filter_spi, ch=0)

            # Program ADI_SPI_CONFIG_A register to 0x3C.
            self.filter.reg_wr(reg="ADI_SPI_CONFIG_A", value=0x3C)
        else:
            raise RuntimeError("%s: filter_spi for filter control not found." % self.__class__.__name__) 

        
    # Extend map_signal_paths.
    def map_signal_paths(self):
        # Run standard QickSoc map_signal_paths.        
        #super().map_signal_paths()
        pass
        
    def description(self):
        lines = []
        lines.append("\n\tBoard: " + self['board'])

        return "\nQICK configuration:\n"+"\n".join(lines)

    def set_nyquist(self, nqz=1, ch=0, btype='dac'):
        # Get converter id.
        block_id = self.CH_MAP[btype][ch]

        # Set nyquist zone.
        self.rf.set_nyquist(block_id, nqz, blocktype=btype)

    def get_nyquist(self, ch=0, btype='dac'):
        # Get converter id.
        block_id = self.CH_MAP[btype][ch]

        # Get nyquist zone.
        return self.rf.get_nyquist(block_id, blocktype=btype)

class Generator():
    def __init__(self, soc, fs, ctrl, gen, sw):
        self.soc = soc
        self.fs = fs
        self.ctrl = ctrl
        self.gen = gen
        self.sw = sw

        # Configure control block.
        self.ctrl.configure(self.fs, self.gen)

    def set(self, f, g=0.99, ch=0, debug=False):
        # Set generator parameters.
        self.ctrl.add(freq = f, gain = g, debug=debug)

        # Select channel.
        self.sw.sel(mst=ch)

    def set_nyquist(self, nqz):
        for dac in self.soc.dacs.keys():
            self.soc.rf.set_nyquist(dac, nqz)

class Buffer():
    def __init__(self, soc, fs, sw, buff, dma):
        self.soc = soc
        self.fs = fs
        self.sw = sw
        self.buff_ip = buff
        self.dma = dma

        # Pre-allocated buffer.
        self.buff = allocate(shape=self.buff_ip['maxlen'], dtype=np.int16)

    def get_data(self, ch=0):
        # Select channel.
        self.sw.sel(slv=ch)

        # Capture.
        self.buff_ip.enable()
        time.sleep(0.1)
        self.buff_ip.disable()
        
        # Transfer.
        return self.transfer()

    def transfer(self):
        # Start send data mode.
        self.buff_ip.dr_start_reg = 1

        # DMA data.
        self.dma.recvchannel.transfer(self.buff)
        self.dma.recvchannel.wait()

        # Stop send data mode.
        self.buff_ip.dr_start_reg = 0

        return self.buff

class AxisSwitchV1(SocIp):
    bindto = ['user.org:user:axis_switch_v1:1.0']
    REGISTERS = {'channel_reg': 0}

    def __init__(self, description):
        """
        Constructor method
        """
        super().__init__(description)

        # Number of bits.
        self.B = int(description['parameters']['B'])
        # Number of master interfaces.
        self.N = int(description['parameters']['N'])

    def sel(self, mst=0):
        if mst > self.N-1:
            print("%s: Master number %d does not exist in block." %
                  __class__.__name__)
            return

        # Select channel.
        self.channel_reg = mst

