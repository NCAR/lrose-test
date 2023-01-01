#!/usr/bin/env python

#===========================================================================
#
# Grid and plot irregular dataProduce plots for noise mon data by time
#
#===========================================================================

from __future__ import print_function

import os
import sys
import subprocess
from optparse import OptionParser
import numpy as np
from scipy.interpolate import griddata
from numpy import convolve
from numpy import linalg, array, ones
import matplotlib.pyplot as plt
from matplotlib import dates
import numpy.ma as ma
from numpy.random import uniform, seed
import h5py as h5
import math
import datetime
import contextlib
import cartopy
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.geodesic as cgds
from cartopy import feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import shapely

def main():

#   globals

    global options
    global debug

# parse the command line

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option('--debug',
                      dest='debug', default=False,
                      action="store_true",
                      help='Set debugging on')
    parser.add_option('--verbose',
                      dest='verbose', default=False,
                      action="store_true",
                      help='Set verbose debugging on')
    parser.add_option('--file',
                      dest='hdf5FilePath',
                      default='/data/dixon/caroline/West_Coast_wave_2005.h5',
                      help='HDF5 data file')
    parser.add_option('--field',
                      dest='fieldName',
                      default='significant_wave_height',
                      help='Name of field to be plotted. Options are: directionality_coefficient, energy_period, maximum_energy_direction, mean_absolute_period, mean_wave_direction, mean_zero-crossing_period, omni-directional_wave_power, peak_period, significant_wave_height, spectral_width, water_depth')
    parser.add_option('--miss',
                      dest='missingVal',
                      default=-999.0,
                      help='Missing value for field')
    parser.add_option('--width',
                      dest='figWidthMm',
                      default=200,
                      help='Width of figure in mm')
    parser.add_option('--height',
                      dest='figHeightMm',
                      default=200,
                      help='Height of figure in mm')
    parser.add_option('--nX',
                      dest='nX',
                      default=500,
                      help='Number of regular grid locations in X')
    parser.add_option('--nY',
                      dest='nY',
                      default=500,
                      help='Number of regular grid locations in Y')
    parser.add_option('--nPts',
                      dest='nPts',
                      default=500,
                      help='Number of points in data set')
    parser.add_option('--test',
                      dest='plotTest', default=False,
                      action="store_true",
                      help='Plot the test data')

    (options, args) = parser.parse_args()
    
    if (options.verbose):
        options.debug = True

    if (options.debug):
        print("Running %prog", file=sys.stderr)
        print("  hdf5FilePath: ", options.hdf5FilePath, file=sys.stderr)
        print("  fieldName: ", options.fieldName, file=sys.stderr)
        print("  missingVal: ", options.missingVal, file=sys.stderr)
        print("  nX: ", options.nX, file=sys.stderr)
        print("  nY: ", options.nY, file=sys.stderr)
        print("  nPts: ", options.nPts, file=sys.stderr)

    # render the test plot

    if (options.plotTest):
        doPlotTest()
        sys.exit(0)

    # open the HDF5 file

    h5File = h5.File(options.hdf5FilePath,'r')
    h5File.keys()

    # plot field data

    doPlotFieldData(h5File)
    
    # done
    
    sys.exit(0)
    
########################################################################
# Plot wave data

def doPlotFieldData(h5File):
    
    widthIn = float(options.figWidthMm) / 25.4
    htIn = float(options.figHeightMm) / 25.4
    fig1 = plt.figure(1, (widthIn, htIn))

    fieldName = options.fieldName
    fieldVals2D = h5File[fieldName][:]

    # field values from only one time

    if (len(fieldVals2D.shape) == 1):
        fVals = fieldVals2D
    else:
        fVals = fieldVals2D[0,:]

    # set missing to Nan
    
    fVals[fVals == float(options.missingVal)] = math.nan
    minVal = np.nanmin(fVals)
    maxVal = np.nanmax(fVals)

    # coordinates
    
    coords = h5File['coordinates'][:]
    lats = coords[:,0]
    lons = coords[:,1]

    if (options.debug):
        print("==========================================", file=sys.stderr)
        print("  fieldName: ", options.fieldName, file=sys.stderr)
        print("  missingVal: ", options.missingVal, file=sys.stderr)
        print("  field min, max: ", minVal, maxVal, file=sys.stderr)
        print("  nPts: ", len(fVals), file=sys.stderr)
        print("  fVals.shape: ",
              fVals.shape, file=sys.stderr)
        print("  len(fVals.shape): ",
              len(fVals.shape), file=sys.stderr)
        
    minLon = np.min(lons)
    maxLon = np.max(lons)
    minLat = np.min(lats)
    maxLat = np.max(lats)

    deltaLon = maxLon - minLon
    deltaLat = maxLat - minLat

    lowerLon = minLon - deltaLon / 50.0
    upperLon = maxLon + deltaLon / 50.0
    lowerLat = minLat - deltaLat / 50.0
    upperLat = maxLat + deltaLat / 50.0
    
    # create a grid with constant spacing
    
    xLons = np.linspace(lowerLon, upperLon, options.nX)
    yLats = np.linspace(lowerLat, upperLat, options.nY)

    # grid the data.
    gVals = griddata((lons, lats), fVals,
                           (xLons[None,:], yLats[:,None]), method='linear')

    if (options.debug):
        print("  xLons.shape: ", xLons.shape, file=sys.stderr)
        print("  yLats.shape: ", yLats.shape, file=sys.stderr)
        print("  gVals.shape: ", gVals.shape, file=sys.stderr)

    if (options.verbose):
        # print("  xLons: ", xLons, file=sys.stderr)
        # print("  yLats: ", yLats, file=sys.stderr)
        print("  fieldVals: ",
              fVals[fVals != math.nan],
              file=sys.stderr)
        print("  gVals: ",
              gVals[gVals != math.nan],
              file=sys.stderr)

    # set up axis for plotting, using Cartopy
    
    ax1 = newMap(fig1, minLon, maxLon, minLat, maxLat)

    # plot the gridded data as color-filled contours

    #CS = plt.contour(xLons,yLats,gVals,16,linewidths=0.5,colors='k')
    CS = ax1.contourf(xLons,yLats,gVals,64,cmap=plt.cm.jet,vmin=minVal,vmax=maxVal)
    cbar = fig1.colorbar(CS, ax=ax1, shrink=0.9) # draw colorbar
    
    # plot data points.
    # plt.scatter(lons,lats,marker='.',c='b',s=5)

    # plot coastlines and states
    
    ax1.coastlines('10m', 'orange', linewidth=1, zorder=3)
    ax1.add_feature(cfeature.STATES, linewidth=0.3, edgecolor='brown', zorder=3)

    # set plot limits
    
    plt.xlim(lowerLon, upperLon)
    plt.ylim(lowerLat, upperLat)

    # title
    
    plt.title(fieldName)

    # show it
    
    plt.show()

########################################################################
# Create map for plotting lat/lon grids

def newMap(fig, minLon, maxLon, minLat, maxLat):
    
    ## Create projection centered as the CWB radar image:
    proj = ccrs.PlateCarree()
    
    ## New axes with the specified projection:
    ax = fig.add_subplot(1, 1, 1, projection=proj)
    
    ## Set extent
    ax.set_extent([minLon, maxLon, minLat, maxLat])
    
    ## Add grid lines & labels:
    gl = ax.gridlines(crs=ccrs.PlateCarree(),
                      draw_labels=True,
                      linewidth=1,
                      color='lightgray',
                      alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.left_labels = True
    gl.right_labels = False
    gl.xlines = True
    gl.ylines = True
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 8, 'weight': 'bold'}
    gl.ylabel_style = {'size': 8, 'weight': 'bold'}
    return ax

########################################################################
# Plot test

def doPlotTest():
    
    nPts = int(options.nPts)
    minX = -2.0
    maxX = 2.0
    minY = -2.0
    maxY = 2.0
    xRange = maxX - minX
    yRange = maxY - minY

    # make up some randomly distributed data
    seed(1234)
    x = uniform(minX,maxX,nPts)
    y = uniform(minY,maxY,nPts)
    z = x*np.exp(-x**2-y**2)
    # define grid.
    xi = np.linspace(minX - xRange / 100.0,
                     maxX + xRange / 100.0,
                     options.nX)
    yi = np.linspace(minY - yRange / 100.0,
                     maxY + yRange / 100.0,
                     options.nY)
    # grid the data.
    zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
    # contour the gridded data, plotting dots at the randomly spaced data points.
    CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
    CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
    plt.colorbar() # draw colorbar
    # plot data points.
    plt.scatter(x,y,marker='o',c='b',s=5)
    plt.xlim(minX,maxX)
    plt.ylim(minY,maxY)
    plt.title('griddata test (%d points)' % nPts)
    plt.show()

########################################################################
# Run a command in a shell, wait for it to complete

def runCommand(cmd):

    if (options.debug):
        print("running cmd:",cmd, file=sys.stderr)
    
    try:
        retcode = subprocess.call(cmd, shell=True)
        if retcode < 0:
            print("Child was terminated by signal: ", -retcode, file=sys.stderr)
        else:
            if (options.debug):
                print("Child returned code: ", retcode, file=sys.stderr)
    except OSError as e:
        print("Execution failed:", e, file=sys.stderr)

########################################################################
# Run - entry point

if __name__ == "__main__":
   main()

