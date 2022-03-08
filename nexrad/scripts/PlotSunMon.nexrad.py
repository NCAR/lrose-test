#!/usr/bin/env python

#===========================================================================
#
# Produce plots for sun scan analysis for NEXRAD
#
#===========================================================================

from __future__ import print_function

import os
import sys
import subprocess
from optparse import OptionParser
import numpy as np
import numpy.ma as ma
from numpy import convolve
import matplotlib.pyplot as plt
from matplotlib import dates
import math
import datetime
import contextlib

def main():

#   globals

    global options
    global debug
    global meanLen
    
    suncalPathDefault = "/scr/cirrus3/rsfdata/projects/nexrad/sunscans/koun/text_table/sunscan.koun.trimmed.txt"
    
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
    parser.add_option('--suncal_file',
                      dest='suncalFilePath',
                      default=suncalPathDefault,
                      help='SunCal results text file path')
    parser.add_option('--title',
                      dest='title',
                      default='NEXRAD SUN SCAN ANALYSIS',
                      help='Title for plot')
    parser.add_option('--width',
                      dest='figWidthMm',
                      default=400,
                      help='Width of figure in mm')
    parser.add_option('--height',
                      dest='figHeightMm',
                      default=320,
                      help='Height of figure in mm')
    parser.add_option('--meanLen',
                      dest='meanLen',
                      default=15,
                      help='Len of moving mean filter')
    parser.add_option('--overlayDays',
                      dest='overlayDays', default=False,
                      action="store_true",
                      help='Overlay days on top of each other')
    
    (options, args) = parser.parse_args()
    
    if (options.verbose == True):
        options.debug = True

    if (options.debug == True):
        print("Running %prog", file=sys.stderr)
        print("  suncalFilePath: ", options.suncalFilePath, file=sys.stderr)

    # read in column headers for sunscan results

    iret, scHdrs, scData = readColumnHeaders(options.suncalFilePath)
    if (iret != 0):
        sys.exit(-1)

    # read in data for sunscan results

    (scTimes, scData) = readInputData(options.suncalFilePath, scHdrs, scData)

    # render the plot
    
    doPlot(scTimes, scData)

    sys.exit(0)
    
########################################################################
# Read columm headers for the data
# this is in the first line

def readColumnHeaders(filePath):

    colHeaders = []
    colData = {}

    fp = open(filePath, 'r')
    line = fp.readline()
    fp.close()
    
    commentIndex = line.find("#")
    if (commentIndex == 0):
        # header
        colHeaders = line.lstrip("# ").rstrip("\n").split()
        if (options.debug == True):
            print("colHeaders: ", colHeaders, file=sys.stderr)
    else:
        print("ERROR - readColumnHeaders", file=sys.stderr)
        print("  First line does not start with #", file=sys.stderr)
        return -1, colHeaders, colData
    
    for index, var in enumerate(colHeaders, start=0):
        colData[var] = []
        
    return 0, colHeaders, colData

########################################################################
# Read in the data

def readInputData(filePath, colHeaders, colData):

    # open file

    fp = open(filePath, 'r')
    lines = fp.readlines()
    fp.close()

    # read in a line at a time, set colData
    for line in lines:
        
        commentIndex = line.find("#")
        if (commentIndex >= 0):
            continue
            
        # data
        
        data = line.strip().split()

        for index, var in enumerate(colHeaders, start=0):
            if (var == 'count' or \
                var == 'year' or var == 'month' or var == 'day' or \
                var == 'hour' or var == 'min' or var == 'sec' or \
                var == 'nBeamsNoise'):
                colData[var].append(int(data[index]))
            else:
                if (isNumber(data[index])):
                    colData[var].append(float(data[index]))
                else:
                    colData[var].append(data[index])

    # load observation times array

    year = colData['year']
    month = colData['month']
    day = colData['day']
    hour = colData['hour']
    minute = colData['min']
    sec = colData['sec']

    startTime = datetime.datetime(year[0], month[0], day[0],
                                  hour[0], minute[0], sec[0])

    obsTimes = []
    for ii, var in enumerate(year, start=0):
        thisTime = datetime.datetime(year[ii], month[ii], day[ii],
                                     hour[ii], minute[ii], sec[ii])
        obsTimes.append(thisTime)

    return obsTimes, colData

########################################################################
# Check is a number

def isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

########################################################################
# Moving average filter

def movingAverage(values, filtLen):

    weights = np.repeat(1.0, filtLen)/filtLen
    sma = np.convolve(values, weights, 'valid')
    smaList = sma.tolist()
    halfLen = int(filtLen / 2)
    for ii in range(0, halfLen):
        smaList.insert(0, smaList[0])
        smaList.append(smaList[-1])
    return np.array(smaList).astype(np.double)

########################################################################
# Plot

def doPlot(scTimes, scData):

    # sunscan times
    
    stimes = np.array(scTimes).astype(datetime.datetime)

    fileName = options.suncalFilePath
    titleStr = "File: " + fileName
    hfmt = dates.DateFormatter('%y/%m/%d')

    # set up plot structure

    widthIn = float(options.figWidthMm) / 25.4
    htIn = float(options.figHeightMm) / 25.4

    fig1 = plt.figure(1, (widthIn, htIn))
    ax1 = fig1.add_subplot(4,1,1,xmargin=0.0)
    ax2 = fig1.add_subplot(4,1,2,xmargin=0.0)
    ax3 = fig1.add_subplot(4,1,3,xmargin=0.0)
    ax4 = fig1.add_subplot(4,1,4,xmargin=0.0)

    #oneDay = datetime.timedelta(1.0)
    #ax1.set_xlim([stimes[0] - oneDay, stimes[-1] + oneDay])
    #ax2.set_xlim([stimes[0] - oneDay, stimes[-1] + oneDay])
    #ax3.set_xlim([stimes[0] - oneDay, stimes[-1] + oneDay])
    #ax4.set_xlim([stimes[0] - oneDay, stimes[-1] + oneDay])

    startTime = scTimes[0]
    endTime = scTimes[-1]
    startOrdinal = startTime.toordinal()
    endOrdinal = endTime.toordinal()

    startTimeLimit = datetime.datetime(startTime.year, startTime.month, startTime.day,
                                       0, 0, 0)
    endTimeLimit = datetime.datetime(endTime.year, endTime.month, endTime.day,
                                     23, 59, 59)
    if (options.overlayDays):
        startTimeLimit = datetime.datetime(startTime.year, startTime.month, startTime.day,
                                           15, 0, 0)
        endTimeLimit = datetime.datetime(startTime.year, startTime.month, startTime.day,
                                         22, 0, 0)

    ax1.set_xlim([startTimeLimit, endTimeLimit])
    ax2.set_xlim([startTimeLimit, endTimeLimit])
    ax3.set_xlim([startTimeLimit, endTimeLimit])
    ax4.set_xlim([startTimeLimit, endTimeLimit])

    # loop through the days

    
    #if (options.debug):
    print("  startTime: ", startTime, file=sys.stderr)
    print("  startOrdinal: ", startOrdinal, file=sys.stderr)
    print("  endTime: ", endTime, file=sys.stderr)
    print("  endOrdinal: ", endOrdinal, file=sys.stderr)

    index = 0
    for dayOrd in range(startOrdinal, endOrdinal + 1):
        doPlotDay(scTimes, scData, dayOrd, index, ax1, ax2, ax3, ax4)
        index = index + 1

    # legends etc
    
    ax1.set_title("Elevation and Azimuth offset (deg)", fontsize=12)
    ax2.set_title("Mean power in solar disk (dBm)", fontsize=12)
    ax3.set_title("SS (dB)", fontsize=12)
    ax4.set_title("Mean correlation", fontsize=12)
    
    configureAxis(ax1, -9999.0, -9999.0, "Antenna errors (deg)", 'upper left')
    configureAxis(ax2, -9999.0, -9999.0, "Receiver power (dBm)", 'upper left')
    configureAxis(ax3, -9999.0, -9999.0, "SS (dB)", 'upper left')
    configureAxis(ax4, -9999.0, -9999.0, "Mean correlation", 'upper left')

    if (options.overlayDays):
        fig1.suptitle("DIURNAL ANALYSIS OF KOUN SUN SCANS, December 2012", fontsize=16)
    else:
        fig1.suptitle("5-DAY ANALYSIS OF KOUN SUN SCANS, December 2012", fontsize=16)

    fig1.autofmt_xdate()

    plt.tight_layout()
    fig1.subplots_adjust(bottom=0.10, left=0.06, right=0.97, top=0.92)
    plt.show()

########################################################################
# Plot data for a day

def doPlotDay(scTimes, scData, dayOrd, index, ax1, ax2, ax3, ax4):
    
    startTime = scTimes[0]
    plotTimes = scTimes
    
    isToday =[]
    for plotTime in plotTimes:
        if (plotTime.toordinal() == dayOrd):
            isToday.append(True)
        else:
            isToday.append(False)

    if (options.overlayDays):
        plotTimes = []
        for scTime in scTimes:
            plotTime = datetime.datetime(startTime.year, startTime.month, startTime.day,
                                         scTime.hour, scTime.minute, scTime.second)
            plotTimes.append(plotTime)
        
            # sunscan times
    
    stimes = np.array(plotTimes).astype(datetime.datetime)

    fileName = options.suncalFilePath
    titleStr = "File: " + fileName
    hfmt = dates.DateFormatter('%y/%m/%d')

    # sun angle offset - only use values < max valid

    elOffset = np.array(scData["centroidElOffset"]).astype(np.double)
    validElOffset = (np.isfinite(elOffset) & (np.absolute(elOffset) < 0.2) & isToday)
    meanElOffset = np.mean(elOffset[validElOffset])
    smoothedElOffset = movingAverage(elOffset[validElOffset], int(options.meanLen))

    azOffset = np.array(scData["centroidAzOffset"]).astype(np.double)
    validAzOffset = (np.isfinite(azOffset) & (np.absolute(azOffset) < 0.2) & isToday)
    meanAzOffset = np.mean(azOffset[validAzOffset])
    smoothedAzOffset = movingAverage(azOffset[validAzOffset], int(options.meanLen))
    
    np.set_printoptions(precision=3)
    np.set_printoptions(suppress=False)

    if (options.debug):
        print("  meanElOffset: ", meanElOffset, file=sys.stderr)
        print("  elOffset: ", elOffset, file=sys.stderr)
        print("  meanAzOffset: ", meanAzOffset, file=sys.stderr)
        print("  azOffset: ", azOffset, file=sys.stderr)

    # load up power
    
    maxValidPower = -55.0
    minValidPower = -70.0
    powerHc = np.array(scData["maxPowerDbmHc"]).astype(np.double)
    validPowerHc = (np.isfinite(powerHc) & isToday &\
                    (powerHc < maxValidPower) & \
                    (powerHc > minValidPower))
    if (options.debug):
        print("  powerHc: ", powerHc, file=sys.stderr)
        print("  validPowerHc: ", validPowerHc, file=sys.stderr)

    smoothedPowerHc = movingAverage(powerHc[validPowerHc], int(options.meanLen))

    powerVc = np.array(scData["maxPowerDbmVc"]).astype(np.double)
    validPowerVc = (np.isfinite(powerVc) & isToday & \
                    (powerVc < maxValidPower) & \
                    (powerVc > minValidPower))
    smoothedPowerVc = movingAverage(powerVc[validPowerVc], int(options.meanLen))
    
    # load up SS, xpol ratio
    
    SS = np.array(scData["SS"]).astype(np.double)
    validSS = (np.isfinite(SS) & (SS > -1.2) & (SS < -0.5) & isToday)
    if (options.debug):
        print("  SS: ", SS, file=sys.stderr)
        print("  validSS: ", validSS, file=sys.stderr)
    smoothedSS = movingAverage(SS[validSS], int(options.meanLen))
    
    corr00 = np.array(scData["corr00"]).astype(np.double)
    validCorr00 = (np.isfinite(corr00) & (corr00 < 0.05) & isToday)
    if (options.debug):
        print("  corr00: ", corr00, file=sys.stderr)
        print("  validCorr00: ", validCorr00, file=sys.stderr)
    smoothedCorr00 = movingAverage(corr00[validCorr00], int(options.meanLen))
    
    XpolR = np.array(scData["ratioDbmVcHc"]).astype(np.double)
    validXpolR = np.isfinite(XpolR) & isToday
    smoothedXpolR = movingAverage(XpolR[validXpolR], int(options.meanLen))
    
    goodTimesHc = stimes[validPowerHc]
    goodTimesVc = stimes[validPowerVc]
    powerHcGood = smoothedPowerHc
    powerVcGood = smoothedPowerVc

    # plot antenna errors - axis 1
    
    if (index == 0):
        ax1.plot(stimes[validElOffset], smoothedElOffset, \
                 label = 'El Offset (deg)', linewidth=1, color='red')
        ax1.plot(stimes[validAzOffset], smoothedAzOffset, \
                 label = 'Az Offset (deg)', linewidth=1, color='blue')
    else:
        ax1.plot(stimes[validElOffset], smoothedElOffset, \
                 linewidth=1, color='red')
        ax1.plot(stimes[validAzOffset], smoothedAzOffset, \
                 linewidth=1, color='blue')
        
    # plot mean power - axis 2

    if (index == 0):
        ax2.plot(stimes[validPowerHc], smoothedPowerHc, \
                 label = 'Power HC', linewidth=1, color='red')
        ax2.plot(stimes[validPowerVc], smoothedPowerVc, \
                 label = 'Power VC', linewidth=1, color='blue')
    else:
        ax2.plot(stimes[validPowerHc], smoothedPowerHc, \
                 linewidth=1, color='red')
        ax2.plot(stimes[validPowerVc], smoothedPowerVc, \
                 linewidth=1, color='blue')

    # plot SS - axis 3
    
    if (index == 0):
        ax3.plot(stimes[validSS], smoothedSS, \
                 label = 'SS', linewidth=1, color='red')
    else:
        ax3.plot(stimes[validSS], smoothedSS, \
                 linewidth=1, color='red')

    # plot correlation - axis 4

    if (index == 0):
        ax4.plot(stimes[validCorr00], smoothedCorr00, \
                 label = 'corr00', linewidth=1, color='blue')
    else:
        ax4.plot(stimes[validCorr00], smoothedCorr00, \
                 linewidth=1, color='blue')
        

########################################################################
# initialize legends etc

def configureAxis(ax, miny, maxy, ylabel, legendLoc):
    
    legend = ax.legend(loc=legendLoc, ncol=6)
    for label in legend.get_texts():
        label.set_fontsize('x-small')
    if (options.overlayDays):
        ax.set_xlabel("Time-of-day")
    else:
        ax.set_xlabel("Date-Hr")
    ax.set_ylabel(ylabel)
    ax.grid(True)
    if (miny > -9990 and maxy > -9990):
        ax.set_ylim([miny, maxy])
    if (options.overlayDays):
        hfmt = dates.DateFormatter('%H:%M')
        ax.xaxis.set_major_locator(dates.HourLocator(interval=1))
    else:
        hfmt = dates.DateFormatter('%y/%m/%d-%H')
        ax.xaxis.set_major_locator(dates.HourLocator(interval=6))
    ax.xaxis.set_major_formatter(hfmt)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(8) 

########################################################################
# get flux closest in time to the search time

def getClosestFlux(hcTime, fluxTimes, obsFlux):

    minDeltaTime = 1.0e99
    fltime = fluxTimes[0]
    flux = obsFlux[0]
    for ii, ftime in enumerate(fluxTimes, start=0):
        deltaTime = math.fabs((ftime - hcTime).total_seconds())
        if (deltaTime < minDeltaTime):
            minDeltaTime = deltaTime
            flux = obsFlux[ii]
            #if (flux > 150.0):
            #    flux = 150.0
            fltime = ftime
            
    return (fltime, flux)

########################################################################
# compute daily means for dbz bias

def computeDailyStats(times, vals):

    dailyTimes = []
    dailyMeans = []

    nptimes = np.array(times).astype(datetime.datetime)
    npvals = np.array(vals).astype(np.double)

    validFlag = np.isfinite(npvals)
    timesValid = nptimes[validFlag]
    valsValid = npvals[validFlag]
    
    startTime = nptimes[0]
    endTime = nptimes[-1]
    
    startDate = datetime.datetime(startTime.year, startTime.month, startTime.day, 0, 0, 0)
    endDate = datetime.datetime(endTime.year, endTime.month, endTime.day, 0, 0, 0)

    oneDay = datetime.timedelta(1)
    halfDay = datetime.timedelta(0.5)
    
    thisDate = startDate
    while (thisDate < endDate + oneDay):
        
        nextDate = thisDate + oneDay
        result = []
        
        sum = 0.0
        sumDeltaTime = datetime.timedelta(0)
        count = 0.0
        for ii, val in enumerate(valsValid, start=0):
            thisTime = timesValid[ii]
            if (thisTime >= thisDate and thisTime < nextDate):
                sum = sum + val
                deltaTime = thisTime - thisDate
                sumDeltaTime = sumDeltaTime + deltaTime
                count = count + 1
                result.append(val)
        if (count > 1):
            mean = sum / count
            meanDeltaTime = datetime.timedelta(0, sumDeltaTime.total_seconds() / count)
            dailyMeans.append(mean)
            dailyTimes.append(thisDate + meanDeltaTime)
            # print >>sys.stderr, " daily time, meanStrong: ", dailyTimes[-1], meanStrong
            result.sort()
            
        thisDate = thisDate + oneDay

    return (dailyTimes, dailyMeans)


########################################################################
# Run a command in a shell, wait for it to complete

def runCommand(cmd):

    if (options.debug == True):
        print("running cmd:",cmd, file=sys.stderr)
    
    try:
        retcode = subprocess.call(cmd, shell=True)
        if retcode < 0:
            print("Child was terminated by signal: ", -retcode, file=sys.stderr)
        else:
            if (options.debug == True):
                print("Child returned code: ", retcode, file=sys.stderr)
    except OSError as e:
        print("Execution failed:", e, file=sys.stderr)

########################################################################
# Run - entry point

if __name__ == "__main__":
   main()

