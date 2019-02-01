import sys
sys.path.append("/home/nvaytet/work/code/mantid/builds/py2/bin")
# created by Owen Arnold (optional e-mail) on 2018-09-19
# last modified by Peter M. Kadletz (peter.kadletz@esss.se) on 2018-11-12
import fabio
import numpy as np
import glob
import csv
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
th_min = 78.0
th_max = 96.0

binsize = (tof_max-tof_min) / nbins
rebin_parameters = '{},{},{}'.format(tof_min,binsize,tof_max)
print("Binsize is",binsize)

# Define theta bins
nth = 50
th_bin = (th_max-th_min) / nth



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

# Make a common set of frame parameters by averaging the ones from sample and background
rebin_step = 64
frame_parameters = []
for i in range(0,len(frame_parameters1),2):
    frame_parameters.append("{},{},{}".format(0.5 * (frame_parameters1[i] + frame_parameters2[i]),
                                              rebin_step,
                                              0.5 * (frame_parameters1[i+1] + frame_parameters2[i+1])))

print("The averaged frame parameters are:", frame_parameters)

# rebin region reduced to cut off frames that contain no data
# rebin_parameters = '8500,64,43000'
# rebin_parameters = '{},64,{}'.format(frame_parameters[0].split(',')[0],frame_parameters[-1].split(',')[2])
# rebin_parameters = '{},64,{}'.format(tofs_lims[0],tofs_lims[1])
rebin_parameters = '{},64,{}'.format(tof_min,tof_max)
instrument_filename = "V20_geometry.nxs" # theFile # os.path.join(grand_parent_dir, 'IDF', 'V20_Definition_GP2.xml')
# mapping_file_to_write = os.path.join(this_dir, 'mapping_file.txt')


# start of processing
# processor_sample = WFMProcessor(frame_parameters1, frame_shifts)
# processor_reference = WFMProcessor(frame_parameters2, frame_shifts)
processor_sample = WFMProcessor(frame_parameters, frame_shifts)
processor_reference = WFMProcessor(frame_parameters, frame_shifts)


sample_stitched = processor_sample.process(summed, instrument_filename, rebin_parameters,scale=1, delete_temporary_workspaces=False)

# # print(sample_stitched)
# # exit()
#
# # Make figure
# fig = plt.figure()
# ax1 = fig.add_subplot(111)
#
#
# # for st in mtd.getObjectNames():
# #     # st = sample_stitched + "_%i"%(j+1)
# #     print(st)
# #     if st.startswith("summed_"):
# #
# #         stitched_ws = mtd[st]
# #         # Get tof centre coordinates
# #         x_edges = stitched_ws.extractX()[0]
# #         nbins = len(x_edges) - 1
# #         x_tof = np.zeros([nbins])
# #         for i in range(nbins):
# #             x_tof[i] = 0.5*(x_edges[i] + x_edges[i+1])
# #
# #         ax1.plot(x_tof,stitched_ws.extractY()[0],label=st)
#
#
# stitched_ws = mtd[sample_stitched]
# # Get tof centre coordinates
# x_edges = stitched_ws.extractX()[0]
# nbins = len(x_edges) - 1
# x_tof = np.zeros([nbins])
# for i in range(nbins):
#     x_tof[i] = 0.5*(x_edges[i] + x_edges[i+1])
# ax1.plot(x_tof,stitched_ws.extractY()[0])
#
# ax1.set_xlabel("Time of flight (microseconds)")
# # ax1.legend()
# fig.savefig('stitched.pdf',bbox_inches='tight')
#
#
# exit()

reference_stitched = processor_reference.process(summed_van, instrument_filename, rebin_parameters,scale=1,delete_temporary_workspaces=True)
normalized = Divide(sample_stitched, reference_stitched)
normalized_wav = ConvertUnits(normalized, Target='Wavelength')


# Make figure
fig = plt.figure()
ratio = 2.0
sizex = 10.0
fig.set_size_inches(sizex, ratio*sizex)
ax1 = fig.add_subplot(411)
ax2 = fig.add_subplot(412)
ax3 = fig.add_subplot(413)
ax4 = fig.add_subplot(414)

xmin = 0.0
xmax = 5.0e4

# Sample
stitched_ws = mtd[sample_stitched]
# Get tof centre coordinates
x_edges = stitched_ws.extractX()[0]
nbins = len(x_edges) - 1
x_tof = np.zeros([nbins])
for i in range(nbins):
    x_tof[i] = 0.5*(x_edges[i] + x_edges[i+1])
ax1.plot(x_tof,stitched_ws.extractY()[0], label="Sample", color='b')
ax1.set_xlim([xmin, xmax])
ax1.set_xlabel(r"Time of flight ($\mu s$)")
ax1.set_ylabel("Counts")
ax1.legend()

# Background
stitched_ws = mtd[reference_stitched]
# Get tof centre coordinates
x_edges = stitched_ws.extractX()[0]
nbins = len(x_edges) - 1
x_tof = np.zeros([nbins])
for i in range(nbins):
    x_tof[i] = 0.5*(x_edges[i] + x_edges[i+1])
ax2.plot(x_tof,stitched_ws.extractY()[0], label="Background", color='r')
ax2.set_xlim([xmin, xmax])
ax2.set_xlabel(r"Time of flight ($\mu s$)")
ax2.set_ylabel("Counts")
ax2.legend()

# Normalized sample
# stitched_ws = mtd[normalized]
# Get tof centre coordinates
x_edges = normalized.extractX()[0]
nbins = len(x_edges) - 1
x_tof = np.zeros([nbins])
for i in range(nbins):
    x_tof[i] = 0.5*(x_edges[i] + x_edges[i+1])
ax3.plot(x_tof,normalized.extractY()[0], label="Normalized sample", color='g')
ax3.set_xlim([xmin, xmax])
ax3.set_xlabel(r"Time of flight ($\mu s$)")
ax3.set_ylabel("Arbitrary units")
ax3.legend()

# Wavelength
x_edges = normalized_wav.extractX()[0]
nbins = len(x_edges) - 1
x_tof = np.zeros([nbins])
for i in range(nbins):
    x_tof[i] = 0.5*(x_edges[i] + x_edges[i+1])
ax4.plot(x_tof,normalized_wav.extractY()[0], label="Normalized sample", color='k')
# ax3.set_xlim([xmin, xmax])
ax4.set_xlabel(r"Wavelength ($\AA$)")
ax4.set_ylabel("Arbitrary units")
ax4.legend()

fig.savefig('stitched.pdf',bbox_inches='tight')


exit()



normalized=ReplaceSpecialValues(normalized, NaNValue=0, InfinityValue=0)
mapping_file = make_map_file(destination=mapping_file_to_write, nx_original=324, ny_original=324, nx_target = nx_target, ny_target = ny_target, start_spectrum_id=1)
normalized = GroupDetectors(normalized, MapFile=mapping_file)
normalized_wav = ConvertUnits(normalized, Target='Wavelength')
Fit(Function='name=LinearBackground,A0=5000,A1=0;name=UserFunction,Formula=h*erfc(a*(x-x0)),h=5000,a=-0.5,x0=4',
InputWorkspace=normalized_wav, Output='fitted', WorkspaceIndex=4, StartX=3.6, EndX=4.4)

# create a nice tile plot
#def tileplot(input_data, nx_target, ny_target):
tw = newTiledWindow()
tw.undock()
tilewindow_number = nx_target * ny_target
tilewindow_matrix = np.reshape(range(tilewindow_number), (nx_target, ny_target))
# transforming the array to match the order the order in the instrument view
# converting the array to an array containing integers
tilewindow_matrix = matrix_transformation(tilewindow_matrix).astype(int)
# plot array
for row in range(len(tilewindow_matrix[:,0])):
    for col in range(len(tilewindow_matrix[0,:])):
        plt = plotSpectrum(normalized_wav, tilewindow_matrix[row,col])
        l = plt.activeLayer()
        # Rescale the x-axis to a smaller region
        l.setAxisScale(Layer.Bottom, 1.8, 7.5)
        # Rescale the y-axis to a smaller region
        l.setAxisScale(Layer.Left, 2E3, 10E3)
        #l.setAxisScale(Layer.Left, 5E1, 7E3)
        tw.addWidget(plt, row, col)

tw2 = newTiledWindow()
tw2.undock()
plots = list()
spec_info = normalized_wav.spectrumInfo();
for i in range(len(spec_info)):
    pos = spec_info.position(i)
    # Append spectrum index and position
    plots.append((i,pos))

# Sort by Y (reversed) then by X as tiled rows go from top 0 to n bottom.
plots.sort(key=lambda t: (-t[1].Y(), t[1].X()))

# plot array
plots_counter = 0
for row in range(ny_target):
    for col in range(nx_target):
        plt = plotSpectrum(normalized_wav, plots[plots_counter][0])
        l = plt.activeLayer()
        # Rescale the x-axis to a smaller region
        l.setAxisScale(Layer.Bottom, 1.8, 7.5)
        # Rescale the y-axis to a smaller region
        l.setAxisScale(Layer.Left, 2E3, 10E3)
        #l.setAxisScale(Layer.Left, 5E1, 7E3)
        tw2.addWidget(plt, row, col)
        plots_counter+=1
