#!/usr/bin/env python

# Plot an unstructured grid.

import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import math

def parseUnstructuredGridSMS(mesh):
    """ Reads in the output of SMS Mesh generation """

    fileRead = open(mesh, 'r')
    lines = fileRead.readlines()
    fileRead.close()

    triangles = []
    x = []
    y = []
    z = []

    for line in lines:
        if 'E3T' in line:
            ttt = line.split()
            t1 = int(ttt[2])-1
            t2 = int(ttt[3])-1
            t3 = int(ttt[4])-1
            triangles.append([t1, t2, t3])
        elif 'ND ' in line:
            xy = line.split()
            x.append(float(xy[2]))
            y.append(float(xy[3]))
            z.append(float(xy[4]))

    # Convert to numpy arrays.
    triangle = np.asarray(triangles)
    X = np.asarray(x)
    Y = np.asarray(y)
    Z = np.asarray(z)

    return(triangle, X, Y, Z)

def parseUnstructuredGridFVCOM(mesh):
    """ Reads in the unstructured grid input for FVCOM """

    fileRead = open(mesh, 'r')
    # Skip the file header (two lines)
    lines = fileRead.readlines()[2:]
    fileRead.close()

    triangles = []
    x = []
    y = []
    z = []

    for line in lines:
        ttt = line.strip().split()
        if len(ttt) == 5:
            t1 = int(ttt[1])-1
            t2 = int(ttt[2])-1
            t3 = int(ttt[3])-1
            triangles.append([t1, t2, t3])
        elif len(ttt) == 4:
            x.append(float(ttt[1]))
            y.append(float(ttt[2]))
            z.append(float(ttt[3]))

    # Convert to numpy arrays.
    triangle = np.asarray(triangles)
    X = np.asarray(x)
    Y = np.asarray(y)
    Z = np.asarray(z)

    return(triangle, X, Y, Z)

def plotUnstructuredGrid(triangles, x, y, z):
    """ Takes the output of parseUnstructuredGrid() and plots it """

    plt.figure()
    if z.max()-z.min() != 0:
        plt.tripcolor(x, y, triangles, z, shading='interp')
        cb = plt.colorbar()
        cb.set_label('Depth (m)')

    plt.triplot(x, y, triangles, '-', color=[0.6, 0.6, 0.6])
    plt.gca().set_aspect('equal')
    plt.gca().axis('tight')
    plt.title('triplot of user-specified triangulation')
    plt.xlabel('x')
    plt.ylabel('y')

    plt.show()

# An SMS grid
#[triangles, x, y, z] = parseUnstructuredGridSMS('tamar_co2V4.2dm')
# The model input grid
[triangles, x, y, z] = parseUnstructuredGridFVCOM('co2_grd.dat')

# Let's have a look-see
plotUnstructuredGrid(triangles, x, y, z)


