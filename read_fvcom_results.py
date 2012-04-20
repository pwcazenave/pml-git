#!/usr/bin/env python

from sys import argv
from netCDF4 import Dataset, MFDataset
import numpy as np
import plot_unstruct_grid as gp
import matplotlib.pyplot as plt

def readFVCOM(file, varList):
    """ 
    Read in the FVCOM results file and spit out numpy arrays for 
    each of the variables.
    """

    rootgrp = Dataset(file, 'r', format='NETCDF4')
    mfdata = MFDataset(file)

    dataOut = {}
    for key, var in rootgrp.variables.items():
        print 'Found ' + key
        if key in varList:
            print 'Extracting ' + key
            dataOut[key] = mfdata.variables[key][:]
    
    return dataOut
        


if __name__ == '__main__':

    dataOut = readFVCOM(argv[1], ['x', 'y', 'DYE', 'time'])
    # Get the DYE values
    Z = dataOut['DYE']
    
    # Get the unstructured grid
    [triangles, nodes, x, y, z] = gp.parseUnstructuredGridFVCOM(argv[2])

    # Static figure
    gp.plotUnstructuredGrid(triangles, nodes, dataOut['x'], dataOut['y'], np.squeeze(Z[50,1,:]), '$CO_{2}$')

    # Animated output (use ipython here)
    #plt.figure()
    #plt.tripcolor(dataOut['x'], dataOut['y'], triangles, np.squeeze(Z[0,1,:]), shading='interp')
    #plt.draw()
    #for i in xrange(len(dataOut['time'])):
    #    print i
    #    plotZ = np.squeeze(Z[i,1,:]) # dim1=time, dim2=sigma, dim3=dye
    #    # Update the plot
    #    plt.clf() 
    #    plt.tripcolor(dataOut['x'], dataOut['y'], triangles, plotZ, shading='interp')
    #    plt.draw()

