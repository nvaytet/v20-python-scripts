import sys
sys.path.append("/home/nvaytet/work/code/mantid/builds/py3/bin")
from mantid.simpleapi import *
import numpy as np

theFile = "V20_ESSIntegration_2018-12-14T0925_stripped.nxs"
ws = LoadEventNexus(FileName=theFile, LoadLogs=True)

# Number and size of bins
nbins = 1000
# tofs_lims = ws.extractX()[0]
tof_min = 0
tof_max = 7.1e4
binsize = (tof_max-tof_min) / nbins
print("Binsize is",binsize)

# Rebin the data
rebinned = Rebin(ws, '{},{},{}'.format(tof_min,binsize,tof_max))
# Get Y data
y = rebinned.extractY()
print("Size of the rebinned Y data", np.shape(y))
nhist = rebinned.getNumberHistograms()

# Get detectors
detInfo = ws.detectorInfo()
ndets = detInfo.size()
detIds = detInfo.detectorIDs()
