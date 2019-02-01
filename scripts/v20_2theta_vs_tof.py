import sys
sys.path.append("/home/nvaytet/work/code/mantid/builds/py3/bin")
from mantid.simpleapi import *
import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.colors import LogNorm

################################################################################
# Make a x-coordinate vs tof diagram
################################################################################

# Load data file
# theFile = "V20_ESSIntegration_2018-12-11_1743_agg_2.nxs"
# theFile = "V20_ESSIntegration_2018-12-13_0942_agg.nxs"
# theFile = "/media/nvaytet/Samsung_T3/V20_ESSIntegration_2018-12-12_1209_agg.nxs"
# theFile = "/media/nvaytet/Samsung_T3/V20_ESSIntegration_2018-12-13T1620_agg.nxs"
# theFile = "/media/nvaytet/Samsung_T3/V20_ESSIntegration_2018-12-14T0925_agg.nxs"
# theFile = "V20_ESSIntegration_2018-12-14T0925_agg.nxs"
# theFile = "/media/nvaytet/30c9d25c-0aba-427f-b8ea-3079e881dfce/experimental_data/V20/V20_ESSIntegration_2018-12-12_1209_agg.nxs"

# theFile = "V20_ESSIntegration_2018-12-14T0925_stripped.nxs"
# theFile = "V20_ESSIntegration_2018-12-12_1209_stripped.nxs"
theFile = "/media/nvaytet/30c9d25c-0aba-427f-b8ea-3079e881dfce/experimental_data/V20/v20-2018-12-14T16:12:26+0100/v20-2018-12-14T16:12:26+0100_stripped.nxs"

ws = LoadEventNexus(FileName=theFile, LoadLogs=True)


#ws2 = CreateSampleWorkspace();
#mon2 = LoadInstrument(ws2, FileName="V20_geometry.nxs", RewriteSpectraMap=True)
#inst2 = ws2.getInstrument()
#di2 = ws2.detectorInfo()
#ci2 = ws2.componentInfo()
#print("Workspace {0} has instrument: {1}".format(ws2.getName(), inst2.getName()))
#print("Instrument {0} has {1} components, including {2} monitors and {3} detectors".format(inst2.getName(), ci2.size(), len(mon2), di2.size()))

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


# Find min and max limits of 2Theta
min_2th = 100.0
max_2th = 0.0
for i in range(ndets):
    # Only process non-monitors
    if not detInfo.isMonitor(i):
        th = detInfo.twoTheta(i)
        min_2th = min(min_2th, th)
        max_2th = max(max_2th, th)

# Number of x coordinates
nx = 300
# make theta limits a bit wider to accomodate for all data points
dth = (max_2th - min_2th) / nx
min_2th -= 0.5*dth
max_2th += 0.5*dth
# Allocate work arrays
z_im = np.zeros([nx,nbins])
# y_im = np.zeros([nx])
rad2deg = 180.0/np.pi
y_im = np.linspace(min_2th*rad2deg, max_2th*rad2deg, nx)
# y_temp = np.zeros([nx])


# Now bin the data
dth = (max_2th - min_2th) / nx
for i in range(ndets):
    # Only process non-monitors
    if not detInfo.isMonitor(i):
        th = detInfo.twoTheta(i)
        binloc = (th - min_2th) / dth
        ibin = int(binloc)
        # print(binloc, ibin)
        z_im[ibin,:] += y[i,:]
#
#
#
# # Store all x coordinates in a dict
# event_dict = dict()
# for i in range(ndets):
#     # Only process non-monitors
#     if not detInfo.isMonitor(i):
#         # detPos = detInfo.position(i)
#         # # key = str(detPos[0]) # X position of detector
#         # key = str(detPos[1])  # X position is actually stored in Y
#         detPos = detInfo.twoTheta(i)
#         # key = str(detPos[0]) # X position of detector
#         key = str(detPos)  # X position is actually stored in Y
#         if key not in event_dict.keys():
#             event_dict[key] = np.zeros([nbins])
#         thisID = detIds[i]
#         idFromWs = rebinned.getDetector(i).getID()
#         # Check that detector ids match
#         if idFromWs != thisID:
#             print("ids do not match:",idFromWs,thisID)
#         event_dict[key] += y[i,:]
#
# # Number of x coordinates
# nx = len(event_dict.keys())
# print("Number of x coordinates",nx)
#
# # Allocate work arrays
# z_im = np.zeros([nx,nbins])
# y_im = np.zeros([nx])
# y_temp = np.zeros([nx])
#
# # Get list of keys
# keylist = list(event_dict.keys())
#
# # First sort the x positions
# j = 0
# for key in keylist:
#     y_temp[j] = float(key)
#     j += 1
# argsort = np.argsort(y_temp)
#
# # Now store into the y and z arrays
# j = 0
# for k in argsort:
#     key = keylist[k]
#     for i in range(nbins):
#         z_im[j,i] = event_dict[key][i]
#     y_im[j] = y_temp[k]
#     j += 1

print("Min/Max X coordinates:",y_im[0],y_im[-1])
print("Min/Max in image:",np.amin(z_im),np.amax(z_im))

# Get tof centre coordinates
x_edges = rebinned.extractX()[0]
x_im = np.zeros([nbins])
for i in range(nbins):
    x_im[i] = 0.5*(x_edges[i] + x_edges[i+1])



# # Crop to the region of interest
# xmin = -0.3
# xmax = 0.1
# tofmin = 0
# tofmax = 80000
#
# xmin_not_found = True
# xmax_not_found = True
# for i in range(nbins):
#     if xmin_not_found and (y_im[i] > xmin):
#         jmin = i
#         xmin_not_found = False










# Make figure
fig = plt.figure()
ratio = 1.0
sizex = 10.0
fig.set_size_inches(sizex,ratio*sizex)
max_tof = 80000
# ax1 = fig.add_subplot(211)
ax1 = fig.add_axes([0, 0.53, 1.0, 0.47])
# contf = ax1.imshow(z_im, origin='lower',extent=[x_im[0],x_im[-1],y_im[0],y_im[-1]], aspect='auto')
contf = ax1.imshow(z_im, origin='lower',extent=[tof_min,tof_max,y_im[0],y_im[-1]], aspect='auto')
# contf = ax1.imshow(z_im, origin='lower', aspect='auto')
# , norm=LogNorm())
# contf = ax1.contourf(x_im,y_im,z_im,nc=50)
# ax1.set_xlim([0,max_tof])
# ax1.set_ylim([-0.3,0.1])
ax1.set_xlabel("Arrival time at detector (microseconds)")
ax1.set_ylabel(r'$2\theta$ (degrees)')
ax3 = fig.add_axes([1.01, 0.53, 0.03, 0.47])
cb = plt.colorbar(contf,ax=ax1,cax=ax3)

ax2 = fig.add_axes([0, 0.0, 1.0, 0.47])

indices = [45, 165, 173]
width = 10
for j in indices:
    ysum = np.zeros([nbins])
    for i in range(nbins):
        ysum[i] = np.sum(z_im[j-width//2:j+width//2,i])
    ax2.plot(x_im,ysum)
    # yav = np.average(y_im[j-width//2:j+width//2])
    ax1.plot([tof_min,tof_max],[y_im[j],y_im[j]])
# ax2.set_xlim([0,max_tof])
ax2.set_xlim([tof_min,tof_max])
# ax2.set_xlabel("Time of flight (microseconds)")
ax2.set_xlabel("Arrival time at detector (microseconds)")
ax2.set_ylabel("Counts (integrated over x)")

fig.savefig('figure.pdf',bbox_inches='tight')
