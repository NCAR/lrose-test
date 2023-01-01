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
    parser.add_option('--width',
                      dest='figWidthMm',
                      default=200,
                      help='Width of figure in mm')
    parser.add_option('--height',
                      dest='figHeightMm',
                      default=200,
                      help='Height of figure in mm')
    parser.add_option('--minX',
                      dest='minX',
                      default=-2.0,
                      help='Minimum X val')
    parser.add_option('--maxX',
                      dest='maxX',
                      default=2.0,
                      help='Maximum X val')
    parser.add_option('--minY',
                      dest='minY',
                      default=-2.0,
                      help='Minimum Y val')
    parser.add_option('--maxY',
                      dest='maxY',
                      default=2.0,
                      help='Maximum Y val')
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

    (options, args) = parser.parse_args()
    
    if (options.verbose):
        options.debug = True

    global nPts, minX, maxX, minY, maxY, xRange, yRange
    nPts = int(options.nPts)
    minX = float(options.minX)
    maxX = float(options.maxX)
    minY = float(options.minY)
    maxY = float(options.maxY)
    xRange = maxX - minX
    yRange = maxY - minY

    if (options.debug):
        print("Running %prog", file=sys.stderr)
        print("  hdf5FilePath: ", options.hdf5FilePath, file=sys.stderr)
        print("  minX: ", options.minX, file=sys.stderr)
        print("  maxX: ", options.maxX, file=sys.stderr)
        print("  minY: ", options.minY, file=sys.stderr)
        print("  maxY: ", options.maxY, file=sys.stderr)
        print("  nX: ", options.nX, file=sys.stderr)
        print("  nY: ", options.nY, file=sys.stderr)
        print("  nPts: ", options.nPts, file=sys.stderr)

    # render the test plot
    
    # doPlotTest()

    # open the HDF5 file

    h5File = h5.File(options.hdf5FilePath,'r')
    h5File.keys()

    # plot wave data

    doPlotWaveData(h5File)
    
    # done
    
    sys.exit(0)
    
########################################################################
# Plot wave data

def doPlotWaveData(h5File):
    
    widthIn = float(options.figWidthMm) / 25.4
    htIn = float(options.figHeightMm) / 25.4
    fig1 = plt.figure(1, (widthIn, htIn))

    wht = h5File['significant_wave_height']
    print("  wht: ", wht, file=sys.stderr)

    fieldName = 'mean_wave_direction'
    mwd = h5File[fieldName][:]
    print("  fieldName: ", mwd, file=sys.stderr)

    coords = h5File['coordinates'][:]

    lats = coords[:,0]
    lons = coords[:,1]

    # mean wave direction from only one time
    mwd_0 = mwd[0,:]
    
    # check dim sizes
    
    print(coords)
    print(coords.shape)

    print(lats.shape)
    print(lats)
    
    print(lons.shape)
    print(lons)
    
    print(mwd_0.shape)
    print(mwd_0)

    # set missing to Nan
    
    mwd_0[mwd_0 == -999] = math.nan
    print(mwd_0)

    minLon = np.min(lons)
    maxLon = np.max(lons)
    minLat = np.min(lats)
    maxLat = np.max(lats)

    deltaLon = maxLon - minLon
    deltaLat = maxLat - minLat

    ax1 = newMap(fig1, minLon, maxLon, minLat, maxLat)

    # create a grid with constant spacing
    
    xi = np.linspace(minLon - deltaLon / 100.0,
                     maxLon + deltaLon / 100.0,
                     options.nX)
    yi = np.linspace(minLat - deltaLat / 100.0,
                     maxLat - deltaLat / 100.0,
                     options.nY)

    # make up some randomly distributed data
    # seed(1234)
    # x = uniform(minX,maxX,nPts)
    # y = uniform(minY,maxY,nPts)
    # z = x*np.exp(-x**2-y**2)
    # # define grid.
    # xi = np.linspace(minX - xRange / 100.0,
    #                  maxX + xRange / 100.0,
    #                  options.nX)
    # yi = np.linspace(minY - yRange / 100.0,
    #                  maxY + yRange / 100.0,
    #                  options.nY)
    # grid the data.
    zi = griddata((lons, lats), mwd_0, (xi[None,:], yi[:,None]), method='cubic')
    # contour the gridded data, plotting dots at the randomly spaced data points.
    #CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
    #CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
    #CS = ax1.contour(xi,yi,zi)
    CS = ax1.contourf(xi,yi,zi)
    fig1.colorbar(CS) # draw colorbar
    # plot data points.
    #plt.scatter(lons,lats,marker='.',c='b',s=5)
    ax1.coastlines('10m', 'orange', linewidth=1, zorder=0)
    ax1.add_feature(cfeature.STATES, linewidth=0.3, edgecolor='brown')
    #ax1.scatter(lons,lats,marker='.',c='b')
    #ax1.coastlines('10m', 'darkgray', linewidth=1, zorder=0)
    plt.xlim(minLon,maxLon)
    plt.ylim(minLat,maxLat)
    plt.title(fieldName)
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

