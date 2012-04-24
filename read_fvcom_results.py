#!/usr/bin/env python

import numpy as np
import plot_unstruct_grid as gp
import matplotlib.pyplot as plt
from readFVCOM import readFVCOM
from sys import argv


if __name__ == '__main__':

    # Be verbose?
    noisy = True

    getVars = ['x', 'y', 'xc', 'yc', 'zeta', 'art1', 'h', 'time', 'TCO2', 'PH', 'DYE', 'u']

    # If running as a script:
    #in1 = argv[1]
    #in2 = argv[2]
    #FVCOM = readFVCOM(in1, getVars)
    #[triangles, nodes, x, y, z] = gp.parseUnstructuredGridFVCOM(argv[2])

    # Running from ipython
    # Riqui's
    #base = '/data/medusa/rito/models/FVCOM/runCO2_leak/'
    #in1 = base + '/output/co2_S5.1.1.2_0002.nc'
    #in1 = base + '/output/co2_S5.1.2.1_0002.nc'
    #in1 = base + '/output/co2_V5.1.2.1_0001.nc'
    #in1 = base + '/output/co2_V7.1.1.1_0001.nc'
    #in1 = base + '/output/co2_V7.1.2.1_0001.nc'
    #in1 = base + '/output/co2_V7.1.2.1_0002.nc'
    #in2 = base + '/input/inputV5/co2_grd.dat'
    #in2 = base + '/input/inputV7/co2_grd.dat'

    #base = '/data/medusa/pica/models/FVCOM/runCO2_leak'
    # Low Fine
    #in1 = base + '/output/low_rate/co2_S7_low_run_0001.nc'
    # High Fine
    #in1 = base + '/output/high_rate/co2_S7_high_run_0001.nc'
    # Fine grid
    #in2 = base + '/input/configs/inputV7/co2_grd.dat'

    base = '/data/medusa/pica/models/FVCOM/runCO2_leak'
    # Low Coarse
    #in1 = base + '/output/low_rate/co2_S5_low_run_0001.nc'
    # High Coarse
    #in1 = base + '/output/high_rate/co2_S5_high_run_0001.nc'
    in1 = base + '/output/high_rate/co2_S1_0001.nc'
    # Coarse grid
    in2 = base + '/input/configs/inputV5/co2_grd.dat'


    print 'Model result: ' + in1
    print 'Input grid: ' + in2

    # Read in the NetCDF file and the unstructured grid file
    FVCOM = readFVCOM(in1, getVars)
    [triangles, nodes, x, y, z] = gp.parseUnstructuredGridFVCOM(in2)

    # Set the leak index point (one less than MATLAB)
    if 'V5' in in1:
        leakIdx = 1315
    elif 'V7' in in1:
        leakIdx = 8186
    if 'S5' in in1:
        leakIdx = 1315
    elif 'S7' in in1:
        leakIdx = 8186

    # Sigma layer
    layerIdx = 0

    # Start index for each input file
    #startIdx = 770 # Riqui's
    startIdx = 50 # Mine

    # Get the relevant variable
    Z = FVCOM['zeta']

    # Calculate total CO2 input
    if False:
        TCO2 = np.zeros(FVCOM['time'].shape)

        print np.shape(TCO2)
        dt = 3600.0

        for i in xrange(startIdx, Z.shape[0]):
            if i > 0:
                #TCO2[i] = TCO2[i-1] + (FVCOM['zeta'][i,0,leakIdx].squeeze() * dt)
                if len(np.shape(Z)) == 3:
                    TCO2[i] = TCO2[i-1] + (Z[i,layerIdx,leakIdx].squeeze() * dt)
                else:
                    TCO2[i] = TCO2[i-1] + (Z[i,leakIdx].squeeze() * dt)

                print "Total DYE: " + str(TCO2[i]) + "\n\tDYE: " + str(Z[i,0,leakIdx].squeeze())

        # Scale to daily input. Input rate begins two days into model run
        TCO2 = TCO2/(FVCOM['time'].max()-FVCOM['time'].min()-2)

        print "\nTotal CO2: " + str(TCO2[-1])

        # Get the total CO2 in the system at the end of the simulation
        totalCO2inSystem = np.sum(Z)
        print totalCO2inSystem/float(60*60*24*14)

        # Make a pretty picture
        plt.figure(100)
        plt.clf()
        #plt.plot(FVCOM['time'],TCO2,'r-x')
        plt.plot(xrange(Z.shape[0]),TCO2,'r-x')
        plt.xlabel('Time')
        plt.ylabel('CO2 input')
        plt.show()


    # Static figure
    #gp.plotUnstructuredGrid(triangles, nodes, FVCOM['x'], FVCOM['y'], np.squeeze(Z[47,:]), '')

    # Animated output (use ipython here)
    if True:
        plt.figure(200)
        plt.clf()
        if len(np.shape(Z)) == 3:
            plt.tripcolor(
                FVCOM['x'],
                FVCOM['y'],
                triangles,
                np.squeeze(Z[0,layerIdx,:]),
                shading='interp')
        else:
            plt.tripcolor(
                FVCOM['x'],
                FVCOM['y'],
                triangles,
                np.squeeze(Z[0,:]),
                shading='interp')

        plt.axes().set_aspect('equal', 'datalim')
        plt.colorbar()
        #plt.clim(6, 8)
        plt.draw()
        for i in xrange(startIdx, len(FVCOM['time'])):
            if len(np.shape(Z)) == 3:
                plotZ = np.squeeze(Z[i,layerIdx,:]) # dim1=time, dim2=sigma, dim3=dye
            else:
                plotZ = np.squeeze(Z[i,:]) # dim1=time, dim2=dye

            if noisy:
                print str(i+1) + ' of ' + str(len(FVCOM['time'])) + ' (date: ' + str(FVCOM['time'][i]) + ')'
                print 'Min: %.2f Max: %.2f Range: %.2f Standard deviation: %.2f' % (plotZ.min(), plotZ.max(), plotZ.max()-plotZ.min(), plotZ.std())
            else:
                print
            # Update the plot
            plt.clf()
            plt.tripcolor(FVCOM['x'], FVCOM['y'], triangles, plotZ, shading='interp')
            plt.colorbar()
            #plt.clim(-1.5, 1.5)
            plt.axes().set_aspect('equal', 'datalim')
            plt.draw()
            plt.show()

