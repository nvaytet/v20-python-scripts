import sys
sys.path.append("/home/nvaytet/work/code/mantid/builds/py3/bin")
from mantid.simpleapi import *
import numpy as np
import matplotlib.pyplot as plt

################################################################################
# Make tof monitor plot
################################################################################

# theFile = "V20_ESSIntegration_2018-12-12_1209_agg.nxs"
# theFile = "V20_ESSIntegration_2018-12-14T0925_stripped.nxs"
theFile = "V20_ESSIntegration_2018-12-12_1209_stripped.nxs"

ws = LoadNexusMonitors(theFile, LoadOnly='Events')
print("Number of monitors = {}".format(ws.getNumberHistograms()))
print("Total number of events = {}".format(ws.getNumberEvents()))

s1 = ws.getSpectrum(0)
s2 = ws.getSpectrum(1)

print("Monitor 1: number of events:",s1.getNumberEvents())
print("Monitor 2: number of events:",s2.getNumberEvents())

print("Monitor 1: Min/Max TOFs:",s1.getTofMin(),s1.getTofMax())
print("Monitor 2: Min/Max TOFs:",s2.getTofMin(),s2.getTofMax())

print("Monitor 1: TOFs:",s1.getTofs())
print("Monitor 2: TOFs:",s1.getTofs())

# Number and size of bins
nbins = 1000
tof_min = 0
tof_max = 7.1e4
binsize = (tof_max-tof_min) / nbins
print("Binsize is",binsize)

# Rebin the data
rebinned = Rebin(ws, "{},{},{}".format(tof_min,binsize,tof_max))
# Get Y data
y = rebinned.extractY()
print("Size of the rebinned Y data", np.shape(y))
nhist = rebinned.getNumberHistograms()

# Make figure
fig = plt.figure()
ax1 = fig.add_subplot(111)

for j in range(nhist):
    # Get tof centre coordinates
    x_edges = rebinned.extractX()[j]
    x = np.zeros([nbins])
    for i in range(nbins):
        x[i] = 0.5*(x_edges[i] + x_edges[i+1])
    # And plot
    ax1.plot(x,y[j,:],label="Monitor {}".format(j+1))

# Finalise figure
ax1.set_xlabel("Time of flight (microseconds)")
ax1.set_ylabel("Monitor counts")
ax1.legend()
fig.savefig('monitors.pdf', bbox_inches='tight')
