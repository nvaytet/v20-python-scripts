import sys
sys.path.append("/home/nvaytet/work/code/mantid/builds/py2/bin")
from mantid.simpleapi import *
# theFile = "V20_ESSIntegration_2018-12-11_1743_agg_2.nxs"
theFile = "V20_ESSIntegration_2018-12-12_1209_agg.nxs"
ws = LoadEventNexus(FileName=theFile, LoadLogs=True)
# Number and size of bins
nbins = 1000
tofs_lims = ws.extractX()[0]
binsize = (tofs_lims[1]-tofs_lims[0]) / nbins
print("Binsize is",binsize)
rebinned = Rebin(ws, binsize)
