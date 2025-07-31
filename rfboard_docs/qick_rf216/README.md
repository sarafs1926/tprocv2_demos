Some draft demo noteboks for the QICK box.

The notebooks you can look at are:
* For general use and some explanation, `demo_03_tprocv2-full.ipynb`
* For saturation testing specifically, `Saturation_Testing.ipynb`

Note that both notebooks use the September 28, 2024 firmware. If you are using newer firmware and software, 
you will need to modify the `set_filter()` function:

```
def set_filter(gen_ch, ro_ch, lpf, hpf, state=0):
    sw = 0xc0 + (hpf<<3) + lpf
    filt_bits = (state<<4) + state
    rfb_ch = soc.gens[gen_ch].rfb_ch
    with soc.board_sel.enable_context(rfb_ch.card_num):
        rfb_ch.filter.write_reg('WR0_SW', sw)
        rfb_ch.filter.write_reg('WR0_FILTER', filt_bits)

    rfb_ch = soc.avg_bufs[ro_ch].rfb_ch
    with soc.board_sel.enable_context(rfb_ch.card_num):
        rfb_ch.filter.write_reg('WR0_SW', sw)
        rfb_ch.filter.write_reg('WR0_FILTER', filt_bits)
```
