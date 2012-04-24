#!/usr/bin/env python

import numpy as np
import plot_unstruct_grid as gp
import matplotlib.pyplot as plt
from readFVCOM import readFVCOM
from sys import argv


def calculateTotalCO2(FVCOM, varPlot, startIdx, layerIdx, leakIdx, dt, noisy=False):
    # Calculate total CO2 input and plot accordingly

    Z = FVCOM[varPlot]

    TCO2 = np.zeros(FVCOM['time'].shape)

    for i in xrange(startIdx, Z.shape[0]):
        if i > 0:
            if len(np.shape(Z)) == 3:
                TCO2[i] = TCO2[i-1] + (Z[i,layerIdx,leakIdx].squeeze() * dt)
            else:
                TCO2[i] = TCO2[i-1] + (Z[i,leakIdx].squeeze() * dt)

            if noisy:
                print "Total " + varPlot + ": " + str(TCO2[i]) + "\n\t" + varPlot + ": " + str(Z[i,0,leakIdx].squeeze())

    # Scale to daily input. Input rate begins two days into model run
    TCO2Scaled = TCO2/(FVCOM['time'].max()-FVCOM['time'].min()-2)

    if noisy:
        print "\nTotal input per day: " + str(TCO2Scaled[-1])

    # Get the total CO2 in the system at the end of the simulation
    totalCO2inSystem = np.sum(Z)
    if noisy:
        print "Total in the system: " + str(totalCO2inSystem)

    # Make a pretty picture
    plt.figure(100)
    plt.clf()
    #plt.plot(FVCOM['time'],TCO2,'r-x')
    plt.plot(xrange(Z.shape[0]),np.squeeze(Z[:,layerIdx,leakIdx]),'r-x')
    plt.xlabel('Time')
    plt.ylabel(varPlot + ' input')
    plt.show()

def animateModelOutput(FVCOM, varPlot, startIdx, layerIdx, noisy=False):
    # Animated output (use ipython here)

    Z = FVCOM[varPlot]

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
        # Add the vectors
        plt.hold('on')
        UU = np.squeeze(FVCOM['u'][i,layerIdx,:])
        VV = np.squeeze(FVCOM['v'][i,layerIdx,:])
        CC = np.sqrt(np.power(UU,2)+np.power(VV,2))
        Q = plt.quiver(FVCOM['xc'], FVCOM['yc'], UU, VV, CC, edgecolors=('k'))
        plt.quiverkey(Q, 0.5, 0.92, 1, r'$2 ms^{-1}$', labelpos='W')
        plt.axes().set_aspect('equal', 'datalim')
        plt.draw()
        plt.show()


def addVectors(x, y, u, v, scale, noisy=False):
    # Little function to be called within another to add vectors to the plot
    plt.quiver(x, y, u, v)



if __name__ == '__main__':

    # Be verbose?
    noisy = True

    getVars = ['x', 'y', 'xc', 'yc', 'zeta', 'art1', 'h', 'time', 'TCO2', 'PH', 'DYE', 'u', 'v']

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
    #in1 = base + '/output/low_rate/co2_S7_run_0001.nc'
    # High Fine
    #in1 = base + '/output/high_rate/co2_S7_run_0001.nc'
    # Fine grid
    #in2 = base + '/input/configs/inputV7/co2_grd.dat'

    base = '/data/medusa/pica/models/FVCOM/runCO2_leak'
    # Low Coarse
    #in1 = base + '/output/low_rate/co2_S5_run_0001.nc'
    # High Coarse
    in1 = base + '/output/high_rate/co2_S5_run_0001.nc'
    # Coarse grid
    in2 = base + '/input/configs/inputV5/co2_grd.dat'

    # Currently running
    #base = '/data/medusa/pica/models/FVCOM/runCO2_leak'
    #in1 = base + '/output/high_rate/co2_S1_0001.nc'
    #in2 = base + '/input/configs/inputV5/co2_grd.dat'
    #in2 = base + '/input/configs/inputV7/co2_grd.dat'

    print 'Model result: ' + in1
    print 'Input grid: ' + in2

    # Read in the NetCDF file and the unstructured grid file
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

    # Sigma layer
    layerIdx = 0

    # Start index for each input file
    #startIdx = 770 # Riqui's
    startIdx = 48 # Mine

    # Do total CO2 analysis
    dt = 3600.0
    #calculateTotalCO2(FVCOM, 'DYE', startIdx, layerIdx, leakIdx, dt, noisy)

    # Animate some variable (ipython only)
    animateModelOutput(FVCOM, 'DYE', startIdx, layerIdx, noisy)

    # Static figure
    #gp.plotUnstructuredGrid(triangles, nodes, FVCOM['x'], FVCOM['y'], np.squeeze(Z[47,:]), '')


