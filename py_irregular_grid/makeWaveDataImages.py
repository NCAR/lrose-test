#!/usr/bin/env python

#===========================================================================
#
# Make plots of wave data, by calling plotWaveData
#
#===========================================================================

from __future__ import print_function

import os
import sys
import subprocess
import math
import datetime
import pytz
from optparse import OptionParser
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib import dates
import h5py as h5
import cartopy
import cartopy.crs as ccrs
from cartopy import feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

def main():

#   globals

    global options
    global debug

# parse the command line

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option('--debug',
                      dest='debug', default=True,
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
    parser.add_option('--time',
                      dest='searchTime',
                      default='2005 01 01 00 00 00',
                      help='Time for plot data')
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
    parser.add_option('--contours',
                      dest='plotContours', default=False,
                      action="store_true",
                      help='Plot line contours as well as filled contours')
    parser.add_option('--nContours',
                      dest='nContours',
                      default=16,
                      help='Number of line contours')
    parser.add_option('--nColors',
                      dest='nColors',
                      default=64,
                      help='Number of colors in filled contours')
    parser.add_option('--saveDir',
                      dest='saveDir',
                      default='/tmp/waveData/images',
                      help='Directory for saved images')

    (options, args) = parser.parse_args()
    
    if (options.verbose):
        options.debug = True

    # compute search time
    
    global searchTime
    tz = pytz.timezone('UTC')
    year, month, day, hour, minute, sec = options.searchTime.split()
    searchTime = datetime.datetime(int(year), int(month), int(day),
                                   int(hour), int(minute), int(sec), 0, tz)
    if (options.debug):
        print("Running %prog", file=sys.stderr)
        print("  hdf5FilePath: ", options.hdf5FilePath, file=sys.stderr)
        print("  searchTime: ", searchTime, file=sys.stderr)

    # make images for all fields

    fieldNames = ['directionality_coefficient',
                  'energy_period',
                  'maximum_energy_direction',
                  'mean_absolute_period',
                  'mean_wave_direction',
                  'mean_zero-crossing_period',
                  'omni-directional_wave_power',
                  'peak_period',
                  'significant_wave_height',
                  'spectral_width',
                  'water_depth']

    # loop through fields
    
    for fieldName in fieldNames:

        # create command to call
        
        cmd = "plotWaveData.py"

        if (options.verbose):
            cmd = cmd + " --verbose"
        elif (options.debug):
            cmd = cmd + " --debug"

        cmd = cmd + " --file " + options.hdf5FilePath
        cmd = cmd + " --field " + fieldName
        cmd = cmd + " --time \'" + options.searchTime + "\'"
        cmd = cmd + " --width " + str(options.figWidthMm)
        cmd = cmd + " --height " + str(options.figHeightMm)
        cmd = cmd + " --nX " + str(options.nX)
        cmd = cmd + " --nY " + str(options.nY)
        if (options.plotContours):
            cmd = cmd + " --contours "
            cmd = cmd + " --nContours " + str(options.nContours)
        cmd = cmd + " --nColors " + str(options.nColors)
        cmd = cmd + " --saveDir " + options.saveDir
        cmd = cmd + " --save "

        runCommand(cmd)

    # done
    
    sys.exit(0)
    
########################################################################
# Run a command in a shell, wait for it to complete

def runCommand(cmd):

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

