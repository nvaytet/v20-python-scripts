import sys
sys.path.append("/home/nvaytet/work/code/mantid/builds/py2/bin")
from mantid.simpleapi import *
import numpy as np
import math

################################################################################
# Make a x-coordinate vs tof diagram
################################################################################

# Load data file
theFile = 'V20_ESSIntegration_2018-12-12_1209_agg.nxs'
ws = LoadEventNexus(FileName=theFile, LoadLogs=True)

tof_min = 0.0
tof_max = 7.1e4

# Number and size of bins
nbins = 1000
tofs_lims = ws.extractX()[0]
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

# Count slice dimensions
xcoords = dict()
ycoords = dict()
xc = []
yc = []
coords_indices = []
for i in range(ndets):
    # Only process non-monitors
    if not detInfo.isMonitor(i):
        detPos = detInfo.position(i)
        ykey = str(detPos[0]) # Y position of detector
        xkey = str(detPos[1])  # X position is actually stored in Y
        if xkey not in xcoords.keys():
            xcoords[xkey] = float(xkey)
            xc.append(float(xkey))
        if ykey not in ycoords.keys():
            ycoords[ykey] = float(ykey)
            yc.append(float(ykey))

nx = len(xcoords.keys())
ny = len(ycoords.keys())

print(nx,ny)

# Allocate space for master array
data_3d_array = np.zeros([nx,ny,nbins])

# Find limits
xmin = np.amin(xc)
xmax = np.amax(xc)
ymin = np.amin(yc)
ymax = np.amax(yc)

dx = (xmax-xmin)/(nx-1)
dy = (ymax-ymin)/(ny-1)

x0 = xmin - 0.5*dx
y0 = ymin - 0.5*dy

for i in range(ndets):
    # Only process non-monitors
    if not detInfo.isMonitor(i):
        detPos = detInfo.position(i)
        iy = int(math.floor((detPos[0]-y0)/dy))
        ix = int(math.floor((detPos[1]-x0)/dx))
        data_3d_array[ix,iy,:] = y[i,:]

# Write the array to disk.
# This has to be done in successive 2D slices as savetxt does not support 3D arrays
with file(theFile.replace("nxs","txt"), 'w') as outfile:
    # Write header information with array size and limits
    outfile.write('# Array shape: {0}\n'.format(data_3d_array.shape))
    outfile.write('# Limits: {} {} {} {} {} {}\n'.format(x0,xmax+0.5*dx,y0,ymax+0.5*dy,tof_min,tof_max))

    # Iterating through a ndimensional array produces slices along
    # the last axis. This is equivalent to data[i,:,:] in this case
    for data_slice in data_3d_array:
        np.savetxt(outfile, data_slice)
        # Writing out a break to indicate different slices (this is optional)
        outfile.write('# New slice\n')
