import h5py
import numpy as np


# input_filename = 'V20_ESSIntegration_2018-12-13T1620_agg.nxs'
input_filename = 'V20_ESSIntegration_2018-12-14T0925_agg.nxs'
with h5py.File(input_filename, 'r+') as input_file:
    del input_file['/entry/instrument/detector_1/transformations/orientation']
    attributes = {'units': 'deg',
                  'vector': [0.0, 1.0, 0.0],
                  'transformation_type': 'rotation',
                  'depends_on': '.',
                  'NX_class': 'NXtransformation'}
    dataset = input_file['/entry/instrument/detector_1/transformations/'].create_dataset('orientation', data=[90.0])
    for key in attributes:
        if isinstance(attributes[key], str):
            dataset.attrs.create(key, np.array(attributes[key]).astype('|S' + str(len(attributes[key]))))
        else:
            dataset.attrs.create(key, np.array(attributes[key]))
