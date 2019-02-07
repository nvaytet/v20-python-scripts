import sys
sys.path.append("/home/nvaytet/work/code/mantid/builds/py2/bin")
# created by Owen Arnold (optional e-mail) on 2018-09-19
# last modified by Peter M. Kadletz (peter.kadletz@esss.se) on 2018-11-12
import fabio
import numpy as np
import glob
import csv
import math
import bisect
import os
from mantid.api import AlgorithmManager
from mantid.simpleapi import *
from wfm_stitching import WFMProcessor, get_wfm_windows, get_frame_shifts, make_frame_shifts

import matplotlib.pyplot as plt


dirpath = "/media/nvaytet/30c9d25c-0aba-427f-b8ea-3079e881dfce/experimental_data/V20/"


# Load data file

theFile = "V20_ESSIntegration_2018-12-12_1209_stripped.nxs"
ws = LoadEventNexus(FileName=theFile, LoadLogs=False)

vanFile = "v20-2018-12-14T16:12:26+0100"
van = LoadEventNexus(FileName=dirpath+vanFile+"/"+vanFile+"_stripped.nxs", LoadLogs=True)

# Number and size of bins
nbins = 1000
#tofs_lims = ws.extractX()[0]
#binsize = (tofs_lims[1]-tofs_lims[0]) / nbins

tof_min = 0
tof_max = 7.1e4

binsize = (tof_max-tof_min) / nbins
rebin_parameters = '{},{},{}'.format(tof_min,binsize,tof_max)
print("Binsize is",binsize)


# Get frame shifts using function in wfm_processing package
frame_shifts = make_frame_shifts(-6630, get_frame_shifts(frequency_wfm_psc=70))
print(get_frame_shifts(frequency_wfm_psc=70))
print("Frame shifts are:", frame_shifts)


# Get frame windows from sample
rebinned = Rebin(ws, rebin_parameters, PreserveEvents=False)
summed = SumSpectra(rebinned)
nz = summed.blocksize()
summed_data = np.zeros([nz,2])
summed_x = summed.extractX()[0]
summed_y = summed.extractY()[0]
for i in range(nz):
    summed_data[i,0] = 0.5*(summed_x[i]+summed_x[i+1])
    summed_data[i,1] = summed_y[i]

frame_parameters1 = get_wfm_windows(data=summed_data, output_format="tof", plot=True, win_threshold=0.30)
print(frame_parameters1)

# Get frame windows from background
rebinned_van = Rebin(van, rebin_parameters, PreserveEvents=False)
summed_van = SumSpectra(rebinned_van)
nz = summed_van.blocksize()
summed_data_van = np.zeros([nz,2])
summed_x_van = summed_van.extractX()[0]
summed_y_van = summed_van.extractY()[0]
for i in range(nz):
    summed_data_van[i,0] = 0.5*(summed_x_van[i]+summed_x_van[i+1])
    summed_data_van[i,1] = summed_y_van[i]

frame_parameters2 = get_wfm_windows(data=summed_data_van, output_format="tof", plot=True, win_threshold=0.30)
print(frame_parameters2)

# # Make a common set of frame parameters by averaging the ones from sample and background
# rebin_step = 64
# frame_parameters = []
# for i in range(0,len(frame_parameters1),2):
#     frame_parameters.append("{},{},{}".format(0.5 * (frame_parameters1[i] + frame_parameters2[i]),
#                                               rebin_step,
#                                               0.5 * (frame_parameters1[i+1] + frame_parameters2[i+1])))

frame_parameters = 0.5 * (np.array(frame_parameters1) + np.array(frame_parameters2))
print("The averaged frame parameters are:", frame_parameters)


# rebin region reduced to cut off frames that contain no data
# rebin_parameters = '8500,64,43000'
# rebin_parameters = '{},64,{}'.format(frame_parameters[0].split(',')[0],frame_parameters[-1].split(',')[2])
# rebin_parameters = '{},64,{}'.format(tofs_lims[0],tofs_lims[1])
# rebin_parameters = '{},64,{}'.format(tof_min,tof_max)
instrument_filename = "V20_geometry.nxs" # theFile # os.path.join(grand_parent_dir, 'IDF', 'V20_Definition_GP2.xml')
# mapping_file_to_write = os.path.join(this_dir, 'mapping_file.txt')


frames = np.reshape(frame_parameters, (len(frame_parameters)//2,2))

y1 = rebinned.extractY()
y2 = rebinned_van.extractY()


# for f in frames:
#     print(f)

# exit()

detInfo = ws.detectorInfo()
ndets = detInfo.size()

# Find min and max limits of 2Theta
rad2deg = 180.0/np.pi

ymin = 1000.0
ymax = 0.0
for i in range(ndets):
    th = detInfo.twoTheta(i) * rad2deg
    ymin = min(ymin, th)
    ymax = max(ymax, th)

# Number of x coordinates
nth = 300
# make theta limits a bit wider to accomodate for all data points
dth = (ymax - ymin) / nth
ymin -= 0.5*dth
ymax += 0.5*dth
# Allocate work arrays
z1_im = np.zeros([nth,nbins])
z2_im = np.zeros([nth,nbins])

# Now bin the data
dth = (ymax - ymin) / nth


# make raw sample image
for i in range(ndets):
    th = detInfo.twoTheta(i)  * rad2deg
    yloc = (th - ymin) / dth
    ybin = int(yloc)
    z1_im[ybin,:] += y1[i,:]
    z2_im[ybin,:] += y2[i,:]





nx = 300

z1 = np.zeros([nth,nx])
z2 = np.zeros([nth,nx])

# print("yminmax",np.amin(y),np.amax(y))

xmin = 0.0
xmax = 5.0e4

dx = (xmax-xmin) / nx


# Get tof centre coordinates
x_edges = rebinned.extractX()[0]
x_im = np.zeros([nbins])
for i in range(nbins):
    x_im[i] = 0.5*(x_edges[i] + x_edges[i+1])


# Now bin the data
counter = 0
for j in range(nbins):
    found_frame = False
    for f in range(len(frames)):
        # print(f,x_im[j],frames[f][0],frames[f][1])
        if (x_im[j] >= frames[f][0]) and (x_im[j] <= frames[f][1]):
            found_frame = True
            xloc = (x_im[j] + frame_shifts[f] - xmin) / float(dx)
            xbin = int(xloc)
            break
    if found_frame:
        z1[:,xbin] += z1_im[:,j]
        z2[:,xbin] += z2_im[:,j]



print(np.amin(z1), np.amax(z1))
print(np.amin(z2), np.amax(z2))



# Make figure
fig = plt.figure()
ax1 = fig.add_subplot(111)
contf = ax1.imshow(z1/z2, origin='lower',extent=[xmin,xmax,ymin,ymax], aspect='auto',vmax=20)
cb = plt.colorbar(contf,ax=ax1)
ax1.set_xlabel("Time of flight (microseconds)")
fig.savefig('stitched_image.pdf')
