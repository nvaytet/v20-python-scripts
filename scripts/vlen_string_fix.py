import h5py
import numpy as np

path = '/entry/instrument'

input_filename = 'V20_ESSIntegration_2018-12-13_0942.nxs'
with h5py.File(input_filename, 'r+') as input_file:
    del input_file[path+'/name']
    dt = h5py.special_dtype(vlen=str)
    dataset = input_file[path].create_dataset('name', (1,), dtype=dt, data=u"V20")
    dataset.attrs.create("short_name", dtype=dt, data=u"v20")
