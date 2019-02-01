import sys
sys.path.append("/home/nvaytet/work/code/mantid/builds/py3/bin")
from mantid.simpleapi import *
import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.colors import LogNorm

################################################################################
# Make a x-coordinate vs tof diagram
################################################################################

files = ["V20_ESSIntegration_2018-12-14T0925_stripped.nxs", "V20_ESSIntegration_2018-12-12_1209_stripped.nxs"]

# Construct the figure axes
fig = plt.figure()
ratio = 0.7
sizex = 20.0
fig.set_size_inches(sizex,ratio*sizex)
ax = []
cbax = []
ncols = len(files)
nrows = 3
deltax = 0.07
deltay = 0.05
cbgap = 0.01
dx = (1.0 - deltax*(ncols-1)) / ncols
dy = (1.0 - deltay*(nrows-1)) / nrows
cbdx = 0.01
for i in range(ncols):
    ax.append(fig.add_axes([i*(dx+deltax), 2*(dy+deltay), dx, dy]))
    ax.append(fig.add_axes([i*(dx+deltax), 1*(dy+deltay), dx, dy]))
    ax.append(fig.add_axes([i*(dx+deltax), 0*(dy+deltay), dx, dy]))
    cbax.append(fig.add_axes([i*(dx+deltax)+dx+cbgap, 2*(dy+deltay), cbdx, dy]))

icol = 0

slice_indices = [45, 165, 173]

for f in files:

    # Read event data ==========================================================

    # Load data file
    ws = LoadEventNexus(FileName=f, LoadLogs=True)
    nevents = ws.getNumberEvents()

    # Number and size of bins
    nbins = 1000
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

    # Store all x coordinates in a dict
    event_dict = dict()
    for i in range(ndets):
        # Only process non-monitors
        if not detInfo.isMonitor(i):
            detPos = detInfo.position(i)
            # key = str(detPos[0]) # X position of detector
            key = str(detPos[1])  # X position is actually stored in Y
            if key not in event_dict.keys():
                event_dict[key] = np.zeros([nbins])
            thisID = detIds[i]
            idFromWs = rebinned.getDetector(i).getID()
            # Check that detector ids match
            if idFromWs != thisID:
                print("ids do not match:",idFromWs,thisID)
            event_dict[key] += y[i,:]

    # Number of x coordinates
    nx = len(event_dict.keys())

    # Allocate work arrays
    z_im = np.zeros([nx,nbins])
    y_im = np.zeros([nx])
    y_temp = np.zeros([nx])

    # Get list of keys
    keylist = list(event_dict.keys())

    # First sort the x positions
    j = 0
    for key in keylist:
        y_temp[j] = float(key)
        j += 1
    argsort = np.argsort(y_temp)

    # Now store into the y and z arrays
    j = 0
    for k in argsort:
        key = keylist[k]
        for i in range(nbins):
            z_im[j,i] = event_dict[key][i]
        y_im[j] = y_temp[k]
        j += 1

    print("Min/Max X coordinates:",y_im[0],y_im[-1])
    print("Min/Max in image:",np.amin(z_im),np.amax(z_im))

    # Get tof centre coordinates
    x_edges = rebinned.extractX()[0]
    x_im = np.zeros([nbins])
    for i in range(nbins):
        x_im[i] = 0.5*(x_edges[i] + x_edges[i+1])

    # Now load monitor data ====================================================

    ws_mon = LoadNexusMonitors(f, LoadOnly='Events')
    print("Number of monitors = {}".format(ws_mon.getNumberHistograms()))
    print("Total number of events = {}".format(ws_mon.getNumberEvents()))

    s1 = ws_mon.getSpectrum(0)
    s2 = ws_mon.getSpectrum(1)

    print("Monitor 1: number of events:",s1.getNumberEvents())
    print("Monitor 2: number of events:",s2.getNumberEvents())

    nevents_mon = [s1.getNumberEvents(), s2.getNumberEvents()]

    print("Monitor 1: Min/Max TOFs:",s1.getTofMin(),s1.getTofMax())
    print("Monitor 2: Min/Max TOFs:",s2.getTofMin(),s2.getTofMax())

    print("Monitor 1: TOFs:",s1.getTofs())
    print("Monitor 2: TOFs:",s1.getTofs())

    # Rebin the data
    rebinned_mon = Rebin(ws_mon, "{},{},{}".format(tof_min,binsize,tof_max))
    # Get Y data
    y_mon = rebinned_mon.extractY()
    print("Size of the rebinned Y data", np.shape(y_mon))
    nhist_mon = rebinned_mon.getNumberHistograms()

    # Now plot on axes =========================================================

    # Top panel: 2D image
    contf = ax[icol*nrows].imshow(z_im, origin='lower',extent=[tof_min,tof_max,y_im[0],y_im[-1]], aspect='auto')
    ax[icol*nrows].set_xlabel("Arrival time at detector (microseconds)")
    ax[icol*nrows].set_ylabel("Detector x position (cm)")
    ax[icol*nrows].set_title(f + " ({} events)".format(nevents))
    cb = plt.colorbar(contf,ax=ax[icol*nrows],cax=cbax[icol])
    cbax[icol].set_ylabel("Counts")
    cb.ax.yaxis.set_label_coords(-0.85,0.5)

    # Middle panel: slices from image
    width = 10
    for j in slice_indices:
        ysum = np.zeros([nbins])
        for i in range(nbins):
            ysum[i] = np.sum(z_im[j-width//2:j+width//2,i])
        ax[icol*nrows+1].plot(x_im,ysum)
        ax[icol*nrows].plot([tof_min,tof_max],[y_im[j],y_im[j]])
    ax[icol*nrows+1].set_xlim([tof_min,tof_max])
    ax[icol*nrows+1].set_xlabel("Arrival time at detector (microseconds)")
    # ax[icol*nrows+1].set_ylabel("Counts (integrated over x)")
    ax[icol*nrows+1].set_ylabel("Counts in slice (slice width = 10 pixels)")

    # Bottom panel: monitor + summed sample data
    colors = ['k','cyan','r']
    ysum = np.zeros([nbins])
    for i in range(nbins):
        ysum[i] = np.sum(z_im[:,i])
    ax[icol*nrows+2].plot(x_im,ysum, label="Sample ({} events)".format(nevents),color=colors[0])
    # Monitors
    for j in range(nhist_mon):
        # Get tof centre coordinates
        x_edges = rebinned_mon.extractX()[j]
        x_mon = np.zeros([nbins])
        for i in range(nbins):
            x_mon[i] = 0.5*(x_edges[i] + x_edges[i+1])
        # And plot
        ax[icol*nrows+2].plot(x_mon,y_mon[j,:], label="Monitor {} ({} events)".format(j+1,nevents_mon[j]),color=colors[j+1])

    ax[icol*nrows+2].set_xlim([tof_min,tof_max])
    ax[icol*nrows+2].set_xlabel("Arrival time at detector (microseconds)")
    ax[icol*nrows+2].set_ylabel("Counts (integrated over x)")
    ax[icol*nrows+2].legend()

    # Cleanup
    DeleteWorkspace(ws)
    DeleteWorkspace(ws_mon)
    DeleteWorkspace(rebinned)
    DeleteWorkspace(rebinned_mon)
    icol += 1

fig.savefig('monitor_and_sample.pdf',bbox_inches='tight')
