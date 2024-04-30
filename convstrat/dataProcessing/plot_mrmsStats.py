#!/usr/bin/env python3

#===========================================================================
#
# Process MRMS ECCO output
#
#===========================================================================

import sys
from optparse import OptionParser
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import geopandas
import netCDF4 as nc
from matplotlib import dates
import datetime
import pathlib
from matplotlib.dates import DateFormatter
import pandas as pd
import pickle
import pyUtils.createFileList
from mpl_toolkits.basemap import Basemap

def main():

    indir='/scr/cirrus2/rsfdata/projects/nexrad-mrms/statMats/'
    infile='mrmsStats_20220501_to_20220601.pickle'
    
    figdir='/scr/cirrus2/rsfdata/projects/nexrad-mrms/figures/eccoStats/'
    
    stringDate=infile[-26:-7]
    
    with open(indir+infile, 'rb') as handle:
        echoType2D=pickle.load(handle)
    
    xlims=[min(echoType2D['lon']),max(echoType2D['lon'])]
    ylims=[min(echoType2D['lat']),max(echoType2D['lat'])]
    
   # bounds = geopandas.read_file('/scr/cirrus2/rsfdata/projects/nexrad-mrms/PoliticalBoundaries_Shapefile/NA_PoliticalDivisions/data/bound_p/boundaries_p_2021_v3.shp')
    
    fig,ax= plt.subplots(4, 3, figsize=(10,13))
    
    for key,value in echoType2D.items():
        m = Basemap(projection='gnom', lat_0=57.3, lon_0=-6.2,
               width=90000, height=120000, resolution='i', ax=ax[1])
        m.fillcontinents(color="#FFDDCC", lake_color='#DDEEFF')
        m.drawmapboundary(fill_color="#DDEEFF")
        m.drawcoastlines()
        ax[key].set_title(key);
            
           # bounds.boundary.plot()
            
########################################################################
# Run - entry point

if __name__ == "__main__":
   main()

