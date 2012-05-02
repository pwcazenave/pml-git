#!/usr/bin/env python

# Stardard library
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
# Custom functions
import read_fvcom_results as rfvcom
import plot_unstruct_grid as gp
from readFVCOM import readFVCOM
from range_test_fit import calculateRegression


if __name__ == '__main__':

    # Be verbose?
    noisy = True

    getVars = ['x', 'y', 'xc', 'yc', 'zeta', 'art1', 'h', 'time', 'TCO2', 'PH', 'DYE',         'siglev']

    # Sigma layer
    layerIdx = 0

    # Start index for each input file
    #startIdx = 770 # Riqui's
    startIdx = 120 # Mine

    # Skip value for animation (not zero)
    skipIdx = 5

    base = '/data/medusa/pica/models/FVCOM/runCO2_leak'

    # Get a list of files
    fileNames = glob(base + '/output/rate_ranges/11days/co2_S5_*_0001.nc')

    # Coarse grid
    in2 = base + '/input/configs/inputV5/co2_grd.dat'

    # Output for calculated CO2
    maxCO2 = np.zeros(np.shape(fileNames))*np.nan
    inputRate = np.zeros(np.shape(fileNames))*np.nan

    for aa, in1 in enumerate(fileNames):

        FVCOM = readFVCOM(in1, getVars, noisy)
        [triangles, nodes, x, y, z] = gp.parseUnstructuredGridFVCOM(in2)

        # Set the leak index point (one less than MATLAB)
        if 'V5' in in1:
            leakIdx = 1315
        elif 'V7' in in1:
            leakIdx = 8186
        elif 'S5' in in1:
            leakIdx = 1315
        elif 'S7' in in1:
            leakIdx = 8186
        else:
            # Could be either...
            leakIdx = 1315
            #leakIdx = 8186

        # Get the input amount (use the first time step at the leak site)
        inputRate[aa] = FVCOM['DYE'][121,0,leakIdx]

        dt = (FVCOM['time'][1]-FVCOM['time'][0])*60*60*24

        # Do total CO2 analysis
        totalCO2inSystem = rfvcom.calculateTotalCO2(FVCOM, 'DYE', startIdx, layerIdx, leakIdx, dt, noisy)

        # Calculate the total CO2 in the system using Riqui's algorithm
        allVolumes = rfvcom.unstructuredGridVolume(FVCOM)
        startDay = (5*24)

        CO2, CO2Leak, maxCO2[aa] = rfvcom.CO2LeakBudget(FVCOM, leakIdx, startDay)

        # Get the concentration for the model
        concZ = FVCOM['DYE']/allVolumes
        # Get the total concentration at n=72 (i.e. 24 hours after DYE release)
        dayConcZ = np.sum(concZ[np.r_[0:25]+startDay,:,:])
        # Scale the DYE by the volumes
        scaledZ = FVCOM['DYE']*allVolumes

        sumZ = np.sum(scaledZ, axis=1)
        totalZ = np.sum(sumZ, axis=1)
        print 'Total DYE at day %i: %.2f' % (startDay, totalZ[startDay])
        plt.figure()
        plt.plot(FVCOM['time'], totalZ, '-x')


    # Reorder the results by the inputRate
    sortedData = np.transpose(np.sort([inputRate, maxCO2]))

    # Calculate a regression
    linX, linY = calculateRegression(sortedData[:,0], sortedData[:,1], 'lin')

    # Make a pretty picture
    plt.figure()
    plt.loglog(sortedData[:,0], sortedData[:,1],'g-x', label='Data')
    plt.loglog(linX, linY, 'r-+', label='Linear regression')
    plt.xlabel('Input rate')
    plt.ylabel('Total CO2 in domain')
    plt.legend(loc=2, frameon=False)



