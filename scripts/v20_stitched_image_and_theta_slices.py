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
from wfm_stitching import WFMProcessor, get_wfm_windows, get_frame_shifts, make_frame_shifts, detect_peaks

import matplotlib.pyplot as plt

import math
import bisect


dirpath = "/media/nvaytet/30c9d25c-0aba-427f-b8ea-3079e881dfce/experimental_data/V20/"

file_list = [["V20_ESSIntegration_2018-12-12_1209", "Silicon", 4]
             ,["V20_ESSIntegration_2018-12-13T16:20:00", "Iron", 3.5]
             ]
             # ,["v20-2018-12-15T18:10:22+0100", "NAK", 2]]


# Make figure
fig = plt.figure()
ratio = 1.0
sizex = 20.0
fig.set_size_inches(sizex, ratio*sizex)
ax = []
ax.append(fig.add_subplot(321))
ax.append(fig.add_subplot(323))
ax.append(fig.add_subplot(325))
ax.append(fig.add_subplot(322))
ax.append(fig.add_subplot(324))
ax.append(fig.add_subplot(326))

colors = ['r', 'g', 'b', 'k', 'cyan', 'magenta', 'r', 'g', 'b', 'k', 'cyan', 'magenta', 'r', 'g', 'b', 'k', 'cyan', 'magenta']

vanFile = "v20-2018-12-14T16:12:26+0100"
van_raw = LoadEventNexus(FileName=dirpath+vanFile+"/"+vanFile+"_stripped.nxs", LoadLogs=True)

for ifile in range(len(file_list)):

    # ifile = 0

    sample_name = file_list[ifile][1]




    # Load data file

    # theFile = "V20_ESSIntegration_2018-12-12_1209_stripped.nxs" # Silicon
    # theFile = "V20_ESSIntegration_2018-12-13T16:20:00" # Iron
    theFile = file_list[ifile][0]
    ws_raw = LoadEventNexus(FileName=dirpath+theFile+"/"+theFile+"_stripped.nxs", LoadLogs=True)
    # ws_raw = LoadEventNexus(FileName=theFile, LoadLogs=True)



    # theFile = "V20_ESSIntegration_2018-12-12_1209_stripped.nxs"
    # ws = LoadEventNexus(FileName=theFile, LoadLogs=False)

    # vanFile = "v20-2018-12-14T16:12:26+0100"
    # van = LoadEventNexus(FileName=dirpath+vanFile+"/"+vanFile+"_stripped.nxs", LoadLogs=True)

    # Number and size of bins
    nbins = 1000
    tof_min = 0
    tof_max = 7.1e4
    binsize = (tof_max-tof_min) / nbins
    rebin_parameters = '{},{},{}'.format(tof_min,binsize,tof_max)

    # Get frame shifts using function in wfm_processing package
    frame_shifts = make_frame_shifts(-6630, get_frame_shifts(frequency_wfm_psc=70))
    print(get_frame_shifts(frequency_wfm_psc=70))
    print("Frame shifts are:", frame_shifts)

    # Get frame windows from sample
    rebinned = Rebin(ws_raw, rebin_parameters, PreserveEvents=False)
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
    rebinned_van = Rebin(van_raw, rebin_parameters, PreserveEvents=False)
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

    detInfo = ws_raw.detectorInfo()
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
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    contf = ax[3*ifile].imshow(z1/z2, origin='lower',extent=[xmin,xmax,ymin,ymax], aspect='auto',vmax=20)
    cb = plt.colorbar(contf,ax=ax[3*ifile])
    ax[3*ifile].set_xlabel("Time of flight (microseconds)")
    ax[3*ifile].set_ylabel(r"Scattering angle $2\theta$ (degrees)")


    # Number and size of bins
    nbins = 1000
    #tofs_lims = ws.extractX()[0]
    #binsize = (tofs_lims[1]-tofs_lims[0]) / nbins

    tof_min = 0
    tof_max = 7.1e4
    th_min = 78.0
    th_max = 96.0
    # th_min = 82.0
    # th_max = 90.0

    binsize = (tof_max-tof_min) / nbins
    rebin_parameters = '{},{},{}'.format(tof_min,binsize,tof_max)
    print("Binsize is",binsize)

    # Define theta bins
    nth = 10
    dth = (th_max-th_min) / nth

    # nrows = 5
    # ncols = 2



    # ax = []
    for k in range(nth):
    # for k in range(2):
        # ax.append(fig.add_subplot(nrows, ncols, k+1))

    # xmin = 0.0
    # xmax = 5.0e4

        th1 = k * dth + th_min
        th2 = th1 + dth

        # Crop the workspace
        limit1 = th1 * math.pi / 180
        limit2 = th2 * math.pi / 180

        wsInst = LoadEmptyInstrument(Filename="V20_geometry.nxs")
        ws1 = ConvertSpectrumAxis(wsInst, Target='theta')
        spec_info = ws1.spectrumInfo()
        axis_values = [spec.twoTheta for spec in spec_info]
        sorted(axis_values)
        lower = bisect.bisect_left(axis_values, limit1)
        upper = bisect.bisect_left(axis_values, limit2)
        ws = CropWorkspace(ws_raw, StartWorkspaceIndex=lower, EndWorkspaceIndex=upper-1)

        # ws2 = ConvertSpectrumAxis(van_raw, Target='theta')
        # spec_info = ws2.spectrumInfo()
        # axis_values = [spec.twoTheta for spec in spec_info]
        # sorted(axis_values)
        # lower = bisect.bisect_left(axis_values, limit1)
        # upper = bisect.bisect_left(axis_values, limit2)
        van = CropWorkspace(van_raw, StartWorkspaceIndex=lower, EndWorkspaceIndex=upper-1)







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

        frame_parameters1 = get_wfm_windows(data=summed_data, output_format="tof", plot=True, win_threshold=0.30, bg_threshold=1.2e-4)
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

        frame_parameters2 = get_wfm_windows(data=summed_data_van, output_format="tof", plot=True, win_threshold=0.30, bg_threshold=1.2e-4)
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

        reference_stitched = processor_reference.process(summed_van, instrument_filename, rebin_parameters,scale=1,delete_temporary_workspaces=True)
        normalized = Divide(sample_stitched, reference_stitched)
        normalized_wav = ConvertUnits(normalized, Target='Wavelength')


        # Wavelength
        x_edges = normalized_wav.extractX()[0]
        nbins = len(x_edges) - 1
        x_wav = np.zeros([nbins])
        for i in range(nbins):
            x_wav[i] = 0.5*(x_edges[i] + x_edges[i+1])
        y_wav = normalized_wav.extractY()[0]
        if k in [0, 5, 9]:
            ax[3*ifile+1].plot(x_wav,y_wav, label=r"$2\theta$"+" = {}".format(0.5*(th1+th2)))
        # # ax3.set_xlim([xmin, xmax])
        # ax[k].set_xlabel(r"Wavelength ($\AA$)")
        # ax[k].set_ylabel("Arbitrary units")
        # ax[k].legend()

        mpd = 20
        peaks = detect_peaks(y_wav, mpd=mpd, mph=file_list[ifile][2])
        # Now filter out peaks that are between start and end
        for p in range(len(peaks)):
            ax[3*ifile+2].plot(x_wav[peaks[p]],0.5*(th1+th2), 'o', label="peak {}".format(p+1), color=colors[p])
            if k == 0:
                ax[3*ifile+1].text(x_wav[peaks[p]], y_wav[peaks[p]], "peak {}".format(p+1), color='k',ha='left', va='bottom')

    ax[3*ifile].set_title(sample_name)
    ax[3*ifile+1].set_xlabel(r"Wavelength ($\AA$)")
    ax[3*ifile+1].set_ylabel("Intensity (arbitrary units)")
    ax[3*ifile+1].legend()

    ax[3*ifile+2].set_xlabel(r"Wavelength ($\AA$)")
    ax[3*ifile+2].set_ylabel(r"Scattering angle $2\theta$ (degrees)")
    ax[3*ifile+2].legend()

fig.savefig('v20_figure_ikon.pdf',bbox_inches='tight')
