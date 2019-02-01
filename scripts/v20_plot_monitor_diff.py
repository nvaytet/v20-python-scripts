import sys
sys.path.append("/home/nvaytet/work/code/mantid/builds/py2/bin")
from mantid.simpleapi import *
import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.colors import LogNorm

################################################################################
# Make a x-coordinate vs tof diagram
################################################################################

# Make figure
fig = plt.figure()
ax = []
ax.append(fig.add_subplot(311))
ax.append(fig.add_subplot(312))
ax.append(fig.add_subplot(313))

theFiles = ["/media/nvaytet/Samsung_T3/V20_ESSIntegration_2018-12-13T1620_agg.nxs",
            "/media/nvaytet/Samsung_T3/V20_ESSIntegration_2018-12-14T0925_agg.nxs"]

tof_min = 0
tof_max = 7.1e4

y_save = []

ax_counter = 0
for f in theFiles:

    # Load data file
    ws = LoadEventNexus(FileName=f, LoadLogs=False)

    # Number and size of bins
    nbins = 1000
    # tofs_lims = ws.extractX()[0]
    # binsize = (tofs_lims[1]-tofs_lims[0]) / nbins
    binsize = (tof_max-tof_min) / nbins
    print("Binsize is",binsize)

    # Rebin the data
    rebinned = Rebin(ws, binsize)
    # Get Y data
    y = rebinned.extractY()
    print("Size of the rebinned Y data", np.shape(y))
    nhist = rebinned.getNumberHistograms()

    # Get detectors
    detInfo = ws.detectorInfo()
    ndets = detInfo.size()
    detIds = detInfo.detectorIDs()

    # Find the monitor
    ispecs = []
    for i in range(ndets):
        if detInfo.isMonitor(i):
            ispecs.append(i)

    mon_counter = 0
    for j in ispecs:
        mon_counter += 1

        # Get tof centre coordinates
        x_edges = rebinned.extractX()[j]
        x_im = np.zeros([nbins])
        for i in range(nbins):
            x_im[i] = 0.5*(x_edges[i] + x_edges[i+1])

        x_save.append(x_im)
        y_save.append(y[j,:])
        ax[ax_counter].plot(x_im,y_save[-1],label="Monitor {}".format(mon_counter))

    ax_counter += 1


ydiff_mon1 = y_save[2] - y_save[0]
ydiff_mon2 = y_save[3] - y_save[1]

ax[2].plot(x_save[0],ydiff_mon1,label="Difference for monitor 1")
ax[2].plot(x_save[1],ydiff_mon2,label="Difference for monitor 2")

for i in range(len(ax)):
    ax[i].set_xlabel("Time of flight (microseconds)")
    ax[i].set_ylabel("Monitor counts")
fig.savefig('monitors.pdf', bbox_inches='tight')
