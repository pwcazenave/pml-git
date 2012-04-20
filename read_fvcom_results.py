#!/usr/bin/env python

from netCDF4 import Dataset, MFDataset
import numpy as np
import plot_unstruct_grid as gp
import matplotlib.pyplot as plt

def readFVCOM(file, varList, noisy=False):
    """ 
    Read in the FVCOM results file and spit out numpy arrays for 
    each of the variables.
    """

    rootgrp = Dataset(file, 'r', format='NETCDF4')
    mfdata = MFDataset(file)

    FVCOM = {}
    for key, var in rootgrp.variables.items():
        if noisy:
            print 'Found ' + key,

        if key in varList:
            if noisy:
                print '(extracted)'
            FVCOM[key] = mfdata.variables[key][:]
        else:
            if noisy:
                print
    
    return FVCOM
        


if __name__ == '__main__':

    from sys import argv

    FVCOM = readFVCOM(argv[1], ['x', 'y', 'zeta', 'art1', 'h', 'zeta','time'])
    Z = FVCOM['zeta'] # dim1=time, dim2=sigma, dim3=element
    
    # Get the unstructured grid
    [triangles, nodes, x, y, z] = gp.parseUnstructuredGridFVCOM(argv[2])

    # Calculate total CO2 input
    #TCO2 = np.zeros(FVCOM['time'].shape)
    #print FVCOM['zeta'].shape
    #leakIdx = 1316
    ##leakIdx = 8187
    #dt = 3600.0
    #for i in xrange(FVCOM['zeta'].shape[0]):
    #    if i > 0:
    #        #print FVCOM['zeta'][i,0,leakIdx].squeeze()
    #        #TCO2[i] = TCO2[i-1] + (FVCOM['zeta'][i,0,leakIdx].squeeze() * dt)
    #        TCO2[i] = TCO2[i-1] + (FVCOM['zeta'][i,leakIdx].squeeze() * dt)

    # Make a pretty picture
#    plt.figure()
#    plt.plot(FVCOM['time'],TCO2,'r-x')
#    plt.xlabel('Time')
#    plt.xlabel('CO2 input')
#    plt.show()

    # Static figure
    gp.plotUnstructuredGrid(triangles, nodes, FVCOM['x'], FVCOM['y'], np.squeeze(Z[48,:]), '$CO_{2}$')

    # Animated output (use ipython here)
    #plt.figure()
    #plt.tripcolor(FVCOM['x'], FVCOM['y'], triangles, np.squeeze(Z[0,1,:]), shading='interp')
    #plt.draw()
    #for i in xrange(len(FVCOM['time'])):
    #    print i
    #    plotZ = np.squeeze(Z[i,1,:]) # dim1=time, dim2=sigma, dim3=dye
    #    # Update the plot
    #    plt.clf() 
    #    plt.tripcolor(FVCOM['x'], FVCOM['y'], triangles, plotZ, shading='interp')
    #    plt.draw()

