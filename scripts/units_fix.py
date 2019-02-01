import h5py
import numpy as np

units = "ns"

# input_filename = 'V20_ESSIntegration_2018-12-13T1620_agg.nxs'
input_filename = 'V20_ESS_example.nxs'
with h5py.File(input_filename, 'r+') as input_file:
    dataset = input_file['/entry/raw_event_data/event_time_offset']
    dataset.attrs.create("units", np.array(units).astype('|S' + str(len(units))))
