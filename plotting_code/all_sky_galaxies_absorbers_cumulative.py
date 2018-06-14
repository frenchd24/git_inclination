#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: all_sky_galaxies_absorbers_cumulative.py, v1 06/14/18

Plot an all sky map of galaxies and absorbers, split into 10 bins of velocity

Cumulative = start with 0-1000, then 0-2000, etc


Based on:
$Id: GTupdate2_allSkyPlot2.py, v2 06/14/18

v2: Make 10 separate plots, one for each 1000 km/s velocity step

Makes all sky plots of a bunch of stuff 

MUST BE RUN INSIDE THE PYSALT CONDA ENVIRONMENT:
>>> source activate pysalt
>>> python GTupdate2_allSkyPlot.py
>>> ...sweet success...

'''

# Import all required packages.
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib.pylab as lab
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable



import sys
import os
import csv
# import string
import warnings
import numpy
# import atpy
import getpass
from utilities import *
from pylab import *
import pickle

import math

import matplotlib.pyplot as plt

from astropy import units as u
from astropy.coordinates import SkyCoord

from astropy import wcs
from astropy.io import fits

# from astropy.io.votable import parse,tree

# from vo.table import parse
# import vo.tree

###########################################################################

def median_low(l):
    # returns the closest element in a list to the median, rounding down

    # E.g.,
    # list = [1, 2, 3, 4, 5]
    # median_low(list) = 3
    #
    # list = [1, 2, 3, 4]
    # median_low(list) = 2
    
#     l.sort()
    l = np.array(l)
    med = np.median(l)
    
    diff = abs(np.subtract(med,l))
    
    diff = list(diff)
    l = list(l)
    indexMin = diff.index(min(diff))
    
    return l[indexMin]

    
def main():
    # check which computer we're on, and grab the appropriate file
    
    user = getpass.getuser()

    if user == 'frenchd':        
        inFilename = '/Users/frenchd/Research/gt/FinalGalaxyTable13_filtered.csv'
        all_filename = '/Users/frenchd/Research/inclination/git_inclination/all4.p'

    elif user =='David':
        pass
#         inFilename = '/Users/David/Research_Documents/GT_update2/FinalGalaxyTable10_filtered.csv'

    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()
        
        
        
    # open the files
        # open the files
    csv.field_size_limit(sys.maxsize)
    
    # read in the tables
    inFile = open(inFilename,'rU')
    reader = csv.DictReader(inFile)
    
    # absorbers
    all_file = open(all_filename,'r')
    all = pickle.load(all_file)
    all_file.close()
    
    all_Lya_vs = all['Lya_vs']
    all_Lya_Ws = all['Lya_Ws']
    all_RA_targets = all['RA_targets']
    all_Dec_targets = all['Dec_targets']
    
    
    # what are the null values equal to?
    nullFloat = -99.99
    nullInt = -99
    nullStr = 'x'
    
    # the header/column names
    fieldnames = (\
    'Name',\
    'NEDname',\
    'z',\
    'RAdeg',\
    'DEdeg',\
    'RAh',\
    'RAm',\
    'RAs',\
    'DE-',\
    'DEd',\
    'DEm',\
    'DEs',\
    'GLON',\
    'GLAT',\
    'Vhel',\
    'vcorr',\
    'distvcorr',\
    'RID_mean',\
    'RID_median',\
    'RID_std',\
    'RID_min',\
    'RID_max',\
    'bestDist',\
    'e_bestDist',\
    'MajDiam_ang',\
    'MinDiam_ang',\
    'e_MajDiam_ang',\
    'e_MinDiam_ang',\
    'MajDiam',\
    'MinDiam',\
    'e_MajDiam',\
    'e_MinDiam',\
    'R_vir',\
    'inc',\
    'adjustedInc',\
    'e_inc',\
    'PA',\
    'diameterKey',\
    'ratioKey',\
    'paKey',\
    'RC3_type',\
    'RC3_d25',\
    'RC3_r25',\
    'RC3_pa',\
    'group_num',\
    'group_mem',\
    'group_dist',\
    'MType',\
    'flag',\
    'distIndicator',\
    'lumClass',\
    'E(B-V)',\
    'Bmag',\
    'Bmag_key',\
    'Bmag_max',\
    'Bmag_max_key',\
    'Bmag_min',\
    'Bmag_min_key',\
    'Bmag_sdss',\
    'gmag_sdss',\
    'rmag_sdss',
    'zmag_sdss',\
    'Lstar_med',\
    'e_Lstar_med',\
    'Lstar_max',\
    'e_Lstar_max',\
    'Lstar_min',\
    'e_Lstar_min',\
    'Lstar_sdss',\
    'e_Lstar_sdss',\
    'altNames')
    
    
##########################################################################################
##########################################################################################    

    hubbleC = 71.0
    
    vcorrs = []
    vhels = []
    ras = []
    decs = []
    
    
    v1 = []
    ra1 = []
    dec1 = []
    Lya_v1 = []
    Lya_W1 = []
    Lya_ra1 = []
    Lya_dec1 = []

    v2 = []
    ra2 = []
    dec2 = []
    Lya_v2 = []
    Lya_W2 = []
    Lya_ra2 = []
    Lya_dec2 = []

    v3 = []
    ra3 = []
    dec3 = []
    Lya_v3 = []
    Lya_W3 = []
    Lya_ra3 = []
    Lya_dec3 = []
        
    v4 = []
    ra4 = []
    dec4 = []
    Lya_v4 = []
    Lya_W4 = []
    Lya_ra4 = []
    Lya_dec4 = []

    v5 = []
    ra5 = []
    dec5 = []
    Lya_v5 = []
    Lya_W5 = []
    Lya_ra5 = []
    Lya_dec5 = []
        
    v6 = []
    ra6 = []
    dec6 = []
    Lya_v6 = []
    Lya_W6 = []
    Lya_ra6 = []
    Lya_dec6 = []

    v7 = []
    ra7 = []
    dec7 = []
    Lya_v7 = []
    Lya_W7 = []
    Lya_ra7 = []
    Lya_dec7 = []

    v8 = []
    ra8 = []
    dec8 = []
    Lya_v8 = []
    Lya_W8 = []
    Lya_ra8 = []
    Lya_dec8 = []
    
    v9 = []
    ra9 = []
    dec9 = []
    Lya_v9 = []
    Lya_W9 = []
    Lya_ra9 = []
    Lya_dec9 = []

    v10 = []
    ra10 = []
    dec10 = []
    Lya_v10 = []
    Lya_W10 = []
    Lya_ra10 = []
    Lya_dec10 = []
            
    # do the work
    count = 0
    for l in reader:
        count +=1
        
        Name = l['Name']
        NEDname = l['NEDname']
        vcorr = eval(l['vcorr'])
        RID_mean = eval(l['RID_mean'])
        RID_median = eval(l['RID_median'])
        Vhel = eval(l['Vhel'])
        flag = eval(l['flag'])

        ra = eval(l['RAdeg'])
        dec = eval(l['DEdeg'])
        
#         if float(ra) <10. and float(ra) > -10:
        vcorrs.append(vcorr)
        vhels.append(Vhel)
        ras.append(ra)
        decs.append(dec)
        
        if flag !=1:
            go = True
        else:
            go = False
            
        if go:
            if Vhel <= 1000.:
                v1.append(Vhel)
                ra1.append(ra)
                dec1.append(dec)
        
            if Vhel <= 2000.:
                v2.append(Vhel)
                ra2.append(ra)
                dec2.append(dec)

            if Vhel <= 3000.:
                v3.append(Vhel)
                ra3.append(ra)
                dec3.append(dec)

            if Vhel <= 4000.:
                v4.append(Vhel)
                ra4.append(ra)
                dec4.append(dec)
            
            if Vhel <= 5000.:
                v5.append(Vhel)
                ra5.append(ra)
                dec5.append(dec)
            
            if Vhel <= 6000.:
                v6.append(Vhel)
                ra6.append(ra)
                dec6.append(dec)
            
            if Vhel <= 7000.:
                v7.append(Vhel)
                ra7.append(ra)
                dec7.append(dec)
            
            if Vhel <= 8000.:
                v8.append(Vhel)
                ra8.append(ra)
                dec8.append(dec)
            
            if Vhel <= 9000.:
                v9.append(Vhel)
                ra9.append(ra)
                dec9.append(dec)

            if Vhel <= 10000.:
                v10.append(Vhel)
                ra10.append(ra)
                dec10.append(dec)
            

            # update the counter
            percentComplete = round((float(count)/130759)*100,2)
            sys.stdout.write('Percent complete: {0}\r'.format(percentComplete))
            sys.stdout.flush()
            
            
    # scale EWs here:
    largestEW = max(all_Lya_Ws)
    smallestEW = min(all_Lya_Ws)
    maxSize = 500
    minSize = 80

    new_Ws = []
    for w in all_Lya_Ws:
        # cut off the super large outlier EWs so they don't mess up all the scaling
        if w > 1000:
            w = 1000
            
        new_W = ((float(w) - smallestEW)/(largestEW - smallestEW)) * (maxSize - minSize) + minSize
        new_Ws.append(new_W)
            
            
    for ra, dec, v, W in zip(all_RA_targets, all_Dec_targets, all_Lya_vs, new_Ws):
        v = float(v)
        W = float(W)
        ra = float(ra)
        dec = float(dec)

        
        if v <= 1000.:
            Lya_v1.append(v)
            Lya_ra1.append(ra)
            Lya_dec1.append(dec)
            Lya_W1.append(W)
        
        
        if v <= 2000.:
            Lya_v2.append(v)
            Lya_ra2.append(ra)
            Lya_dec2.append(dec)
            Lya_W2.append(W)


        if v <= 3000.:
            Lya_v3.append(v)
            Lya_ra3.append(ra)
            Lya_dec3.append(dec)
            Lya_W3.append(W)


        if v <= 4000.:
            Lya_v4.append(v)
            Lya_ra4.append(ra)
            Lya_dec4.append(dec)
            Lya_W4.append(W)

            
        if v <= 5000.:
            Lya_v5.append(v)
            Lya_ra5.append(ra)
            Lya_dec5.append(dec)
            Lya_W5.append(W)

            
        if v <= 6000.:
            Lya_v6.append(v)
            Lya_ra6.append(ra)
            Lya_dec6.append(dec)
            Lya_W6.append(W)

            
        if v <= 7000.:
            Lya_v7.append(v)
            Lya_ra7.append(ra)
            Lya_dec7.append(dec)
            Lya_W7.append(W)

            
        if v <= 8000.:
            Lya_v8.append(v)
            Lya_ra8.append(ra)
            Lya_dec8.append(dec)
            Lya_W8.append(W)

            
        if v <= 9000.:
            Lya_v9.append(v)
            Lya_ra9.append(ra)
            Lya_dec9.append(dec)
            Lya_W9.append(W)


        if v <= 10000.:
            Lya_v10.append(v)
            Lya_ra10.append(ra)
            Lya_dec10.append(dec)
            Lya_W10.append(W)

            
            
##########################################################################################
##########################################################################################

    saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/plotting_code/figs/all_sky_map_cumulative/'


    # plot it
#     colmap = cm.RdBu_r
#     colmap = cm.cool
#     colors = numpy.array(vhels)
#     
#     vmaxVal = max(colors)
#     vminVal = min(colors)
# 
#     norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
#     m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)

    alpha_galaxy = 0.5
    alpha_Lya = 0.8
    lw = 0.8
#     colmap = cm.cool
    colmap = cm.viridis


##########################################################################################
    # 0 - 1000 km/s

    colors = numpy.array(v1)
    
    vmaxVal = max(colors)
    vminVal = min(colors)

    norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)

    # galaxies first
    ras = numpy.array(ra1)
    decs = numpy.array(dec1)
#     colors = numpy.array(v1)
    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')

    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi

    # Now plot the data in Aitoff projection with a grid.
    fig = plt.figure(figsize=(12,6))
    ax = lab.subplot(111,projection="aitoff")

    plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
    vmax=vmaxVal,lw=0,cmap=colmap,s=4,alpha=alpha_galaxy)
    
    
    # now absorbers
    ras = numpy.array(Lya_ra1)
    decs = numpy.array(Lya_dec1)
    colors = numpy.array(Lya_v1)
    sizes = numpy.array(Lya_W1)
    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')

    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi

    plot2 = plt.scatter(ra_rad, dec_rad, marker='*', c=colors, vmin=vminVal,\
    vmax=vmaxVal, lw=lw, cmap=colmap, s=sizes, alpha=alpha_Lya, edgecolor='black')
    
    
    cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
    cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=530)
    ax.grid(color='k', linestyle='solid', linewidth=0.5)
    tight_layout()
    
#     plt.show()
    plt.savefig('{0}1.pdf'.format(saveDirectory),format='pdf')

##########################################################################################
    # 1000 - 2000 km/s
    
    colors = numpy.array(v2)
    
    vmaxVal = max(colors)
    vminVal = min(colors)

    norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)
    
    
    # first galaxies
    ras = numpy.array(ra2)
    decs = numpy.array(dec2)
#     colors = numpy.array(v2)
    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')

    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi

    # Now plot the data in Aitoff projection with a grid.
    fig = plt.figure(figsize=(12,6))
    ax = lab.subplot(111,projection="aitoff")
    # lab.subplot(111)

    plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
    vmax=vmaxVal,lw=0,cmap=colmap,s=4,alpha=alpha_galaxy)
    
    
    # now absorbers
    ras = numpy.array(Lya_ra2)
    decs = numpy.array(Lya_dec2)
    colors = numpy.array(Lya_v2)
    sizes = numpy.array(Lya_W2)
    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')

    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi

    plot2 = plt.scatter(ra_rad, dec_rad, marker='*', c=colors, vmin=vminVal,\
    vmax=vmaxVal, lw=lw, cmap=colmap, s=sizes, alpha=alpha_Lya, edgecolor='black')
    
    
    cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
    cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=530)
    ax.grid(color='k', linestyle='solid', linewidth=0.5)
    tight_layout()
    
#     plt.show()
    plt.savefig('{0}2.pdf'.format(saveDirectory),format='pdf')
    
##########################################################################################
    # 2000 - 3000 km/s
    
    colors = numpy.array(v3)
    
    vmaxVal = max(colors)
    vminVal = min(colors)

    norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)
    
    # first galaxies
    ras = numpy.array(ra3)
    decs = numpy.array(dec3)
#     colors = numpy.array(v3)
    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi

    # Now plot the data in Aitoff projection with a grid.
    fig = plt.figure(figsize=(12,6))
    ax = lab.subplot(111,projection="aitoff")
    # lab.subplot(111)

    plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
    vmax=vmaxVal,lw=0,cmap=colmap,s=4,alpha=alpha_galaxy)
    
    
    # now absorbers
    ras = numpy.array(Lya_ra3)
    decs = numpy.array(Lya_dec3)
    colors = numpy.array(Lya_v3)
    sizes = numpy.array(Lya_W3)
    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')

    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi

    plot2 = plt.scatter(ra_rad, dec_rad, marker='*', c=colors, vmin=vminVal,\
    vmax=vmaxVal, lw=lw, cmap=colmap, s=sizes, alpha=alpha_Lya, edgecolor='black')
    
    
    cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
    cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=530)
    ax.grid(color='k', linestyle='solid', linewidth=0.5)
    tight_layout()
    
#     plt.show()
    plt.savefig('{0}3.pdf'.format(saveDirectory),format='pdf')

##########################################################################################
    # 3000 - 4000 km/s
    
    colors = numpy.array(v4)
    
    vmaxVal = max(colors)
    vminVal = min(colors)

    norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)
    
    # first galaxies
    ras = numpy.array(ra4)
    decs = numpy.array(dec4)
#     colors = numpy.array(v4)
    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')

    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi

    # Now plot the data in Aitoff projection with a grid.
    fig = plt.figure(figsize=(12,6))
    ax = lab.subplot(111,projection="aitoff")
    # lab.subplot(111)

    plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
    vmax=vmaxVal,lw=0,cmap=colmap,s=4,alpha=alpha_galaxy)
    
    
    # now absorbers
    ras = numpy.array(Lya_ra4)
    decs = numpy.array(Lya_dec4)
    colors = numpy.array(Lya_v4)
    sizes = numpy.array(Lya_W4)
    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')

    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi

    plot2 = plt.scatter(ra_rad, dec_rad, marker='*', c=colors, vmin=vminVal,\
    vmax=vmaxVal, lw=lw, cmap=colmap, s=sizes, alpha=alpha_Lya, edgecolor='black')
    
    
    cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
    cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=530)
    ax.grid(color='k', linestyle='solid', linewidth=0.5)
    tight_layout()
    
#     plt.show()
    plt.savefig('{0}4.pdf'.format(saveDirectory),format='pdf')

##########################################################################################
    # 4000 - 5000 km/s
    
    colors = numpy.array(v5)
    
    vmaxVal = max(colors)
    vminVal = min(colors)

    norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)
    
    # first galaxies
    ras = numpy.array(ra5)
    decs = numpy.array(dec5)
#     colors = numpy.array(v5)
    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')

    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi

    # Now plot the data in Aitoff projection with a grid.
    fig = plt.figure(figsize=(12,6))
    ax = lab.subplot(111,projection="aitoff")
    # lab.subplot(111)

    plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
    vmax=vmaxVal,lw=0,cmap=colmap,s=4,alpha=alpha_galaxy)
    

    # now absorbers
    ras = numpy.array(Lya_ra5)
    decs = numpy.array(Lya_dec5)
    colors = numpy.array(Lya_v5)
    sizes = numpy.array(Lya_W5)
    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')

    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi

    plot2 = plt.scatter(ra_rad, dec_rad, marker='*', c=colors, vmin=vminVal,\
    vmax=vmaxVal, lw=lw, cmap=colmap, s=sizes, alpha=alpha_Lya, edgecolor='black')
    
    
    cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
    cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=530)
    ax.grid(color='k', linestyle='solid', linewidth=0.5)
    tight_layout()
    
#     plt.show()
    plt.savefig('{0}5.pdf'.format(saveDirectory),format='pdf')



##########################################################################################
    # 5000 - 6000 km/s
    
    colors = numpy.array(v6)
    
    vmaxVal = max(colors)
    vminVal = min(colors)

    norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)
    
    
    # first galaxies
    ras = numpy.array(ra6)
    decs = numpy.array(dec6)
#     colors = numpy.array(v6)
    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi

    # Now plot the data in Aitoff projection with a grid.
    fig = plt.figure(figsize=(12,6))
    ax = lab.subplot(111,projection="aitoff")
    # lab.subplot(111)

    plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
    vmax=vmaxVal,lw=0,cmap=colmap,s=4,alpha=alpha_galaxy)
    

    # now absorbers
    ras = numpy.array(Lya_ra6)
    decs = numpy.array(Lya_dec6)
    colors = numpy.array(Lya_v6)
    sizes = numpy.array(Lya_W6)
    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')

    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi

    plot2 = plt.scatter(ra_rad, dec_rad, marker='*', c=colors, vmin=vminVal,\
    vmax=vmaxVal, lw=lw, cmap=colmap, s=sizes, alpha=alpha_Lya, edgecolor='black')
    
    
    cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
    cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=530)
    ax.grid(color='k', linestyle='solid', linewidth=0.5)
    tight_layout()
    
#     plt.show()
    plt.savefig('{0}6.pdf'.format(saveDirectory),format='pdf')

##########################################################################################
    # 6000 - 7000 km/s
    
    colors = numpy.array(v7)
    
    vmaxVal = max(colors)
    vminVal = min(colors)

    norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)
    
    # first galaxies
    ras = numpy.array(ra7)
    decs = numpy.array(dec7)
#     colors = numpy.array(v7)
    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi

    # Now plot the data in Aitoff projection with a grid.
    fig = plt.figure(figsize=(12,6))
    ax = lab.subplot(111,projection="aitoff")
    # lab.subplot(111)

    plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
    vmax=vmaxVal,lw=0,cmap=colmap,s=4,alpha=alpha_galaxy)
    
    
    # now absorbers
    ras = numpy.array(Lya_ra7)
    decs = numpy.array(Lya_dec7)
    colors = numpy.array(Lya_v7)
    sizes = numpy.array(Lya_W7)
    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')

    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi

    plot2 = plt.scatter(ra_rad, dec_rad, marker='*', c=colors, vmin=vminVal,\
    vmax=vmaxVal, lw=lw, cmap=colmap, s=sizes, alpha=alpha_Lya, edgecolor='black')
    
    
    cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
    cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=530)
    ax.grid(color='k', linestyle='solid', linewidth=0.5)
    tight_layout()
    
#     plt.show()
    plt.savefig('{0}7.pdf'.format(saveDirectory),format='pdf')
    
    
    
##########################################################################################
    # 7000 - 8000 km/s
    
    colors = numpy.array(v8)
    
    vmaxVal = max(colors)
    vminVal = min(colors)

    norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)
    
    # first galaxies
    ras = numpy.array(ra8)
    decs = numpy.array(dec8)
#     colors = numpy.array(v8)
    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')

    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi

    # Now plot the data in Aitoff projection with a grid.
    fig = plt.figure(figsize=(12,6))
    ax = lab.subplot(111,projection="aitoff")
    # lab.subplot(111)

    plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
    vmax=vmaxVal,lw=0,cmap=colmap,s=4,alpha=alpha_galaxy)
    
    
    # now absorbers
    ras = numpy.array(Lya_ra8)
    decs = numpy.array(Lya_dec8)
    colors = numpy.array(Lya_v8)
    sizes = numpy.array(Lya_W8)
    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')

    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi

    plot2 = plt.scatter(ra_rad, dec_rad, marker='*', c=colors, vmin=vminVal,\
    vmax=vmaxVal, lw=lw, cmap=colmap, s=sizes, alpha=alpha_Lya, edgecolor='black')
    
    
    cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
    cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=530)
    ax.grid(color='k', linestyle='solid', linewidth=0.5)
    tight_layout()
    
#     plt.show()
    plt.savefig('{0}8.pdf'.format(saveDirectory),format='pdf')
    
    
##########################################################################################
    # 8000 - 9000 km/s
    
    colors = numpy.array(v9)
    
    vmaxVal = max(colors)
    vminVal = min(colors)

    norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)
    
    # first galaxies
    ras = numpy.array(ra9)
    decs = numpy.array(dec9)
#     colors = numpy.array(v9)
    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')


    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi

    # Now plot the data in Aitoff projection with a grid.
    fig = plt.figure(figsize=(12,6))
    ax = lab.subplot(111,projection="aitoff")
    # lab.subplot(111)

    plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
    vmax=vmaxVal,lw=0,cmap=colmap,s=4,alpha=alpha_galaxy)
    
    
    # now absorbers
    ras = numpy.array(Lya_ra9)
    decs = numpy.array(Lya_dec9)
    colors = numpy.array(Lya_v9)
    sizes = numpy.array(Lya_W9)
    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')

    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi

    plot2 = plt.scatter(ra_rad, dec_rad, marker='*', c=colors, vmin=vminVal,\
    vmax=vmaxVal, lw=lw, cmap=colmap, s=sizes, alpha=alpha_Lya, edgecolor='black')
    
    
    cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
    cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=530)
    ax.grid(color='k', linestyle='solid', linewidth=0.5)
    tight_layout()
    
#     plt.show()
    plt.savefig('{0}9.pdf'.format(saveDirectory),format='pdf')

##########################################################################################
    # 9000 - 10000 km/s
    
    colors = numpy.array(v10)
    
    vmaxVal = max(colors)
    vminVal = min(colors)

    norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)
    
    
    # first galaxies
    ras = numpy.array(ra10)
    decs = numpy.array(dec10)
#     colors = numpy.array(v10)
    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')

    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi

    # Now plot the data in Aitoff projection with a grid.
    fig = plt.figure(figsize=(12,6))
    ax = lab.subplot(111,projection="aitoff")
    # lab.subplot(111)

    plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
    vmax=vmaxVal,lw=0,cmap=colmap,s=4,alpha=alpha_galaxy)
    
    
    # now absorbers
    ras = numpy.array(Lya_ra10)
    decs = numpy.array(Lya_dec10)
    colors = numpy.array(Lya_v10)
    sizes = numpy.array(Lya_W10)
    # Transform the data into a SkyCoord object.
    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')

    # Matplotlib needs the coordinates in radians, so we have to convert them.
    # Furtermore matplotlib needs the RA coordinate in the range between -pi
    # and pi, not 0 and 2pi.
    ra_rad = c.ra.radian
    dec_rad = c.dec.radian
    ra_rad[ra_rad > np.pi] -= 2. * np.pi

    plot2 = plt.scatter(ra_rad, dec_rad, marker='*', c=colors, vmin=vminVal,\
    vmax=vmaxVal, lw=lw, cmap=colmap, s=sizes, alpha=alpha_Lya, edgecolor='black')
    
    
    cbar = plt.colorbar(plot1,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
    cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=530)
    ax.grid(color='k', linestyle='solid', linewidth=0.5)
    tight_layout()
    
#     plt.show()
    plt.savefig('{0}10.pdf'.format(saveDirectory),format='pdf')

##########################################################################################
##########################################################################################
    
    # close the files
    inFile.close()
    
    print "Done."
    print
    
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
