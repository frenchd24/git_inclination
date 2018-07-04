#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: all_sky_galaxies_absorbers2.py, v1 07/04/18

Plot an all sky map of galaxies and absorbers, split into 4 bins of velocity


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
        all_filename = '/Users/frenchd/Research/inclination/git_inclination/all8.p'

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

#         if flag !=1:
#             go = True
#         else:
#             go = False
            
        if flag == 0:
            go = True
        else:
            go = False
            
        if Vhel <= 450:
            go = False
            
        if go:
            if Vhel <= 2500.:
                v1.append(Vhel)
                ra1.append(ra)
                dec1.append(dec)
        
            if Vhel <= 5000. and Vhel >2500:
                v2.append(Vhel)
                ra2.append(ra)
                dec2.append(dec)

            if Vhel <= 7500. and Vhel >5000:
                v3.append(Vhel)
                ra3.append(ra)
                dec3.append(dec)

            if Vhel <= 10000. and Vhel >7500:
                v4.append(Vhel)
                ra4.append(ra)
                dec4.append(dec)
            

            # update the counter
            percentComplete = round((float(count)/130759)*100,2)
            sys.stdout.write('Percent complete: {0}\r'.format(percentComplete))
            sys.stdout.flush()
            
            
    # scale EWs here:
    largestEW = max(all_Lya_Ws)
    smallestEW = min(all_Lya_Ws)
    maxSize = 1500
    minSize = 20

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

        
        if v <= 2500.:
            Lya_v1.append(v)
            Lya_ra1.append(ra)
            Lya_dec1.append(dec)
            Lya_W1.append(W)
        
        
        if v <= 5000. and v > 2500:
            Lya_v2.append(v)
            Lya_ra2.append(ra)
            Lya_dec2.append(dec)
            Lya_W2.append(W)


        if v <= 7500. and v > 5000:
            Lya_v3.append(v)
            Lya_ra3.append(ra)
            Lya_dec3.append(dec)
            Lya_W3.append(W)


        if v <= 10000. and v > 7500:
            Lya_v4.append(v)
            Lya_ra4.append(ra)
            Lya_dec4.append(dec)
            Lya_W4.append(W)

            
            
##########################################################################################
##########################################################################################

    saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/plotting_code/figs/all_sky_map/'


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
    alpha_Lya = 0.9
    lw = 0.45
    size_galaxy = 2
#     colmap = cm.cool
#     colmap = cm.viridis
    colmap = cm.plasma

    include_galaxies = True

##########################################################################################
    # 0 - 2500 km/s
    # Now plot the data in Aitoff projection with a grid.
    fig = plt.figure(figsize=(12,6))
    ax = lab.subplot(111, projection="aitoff")
        
    colors = numpy.array(v1)
    vmaxVal = max(colors)
    vminVal = min(colors)

    norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)

    # galaxies first
    
    if include_galaxies:
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

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
    
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
    
    
    cbar = plt.colorbar(plot2,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
    cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=530)
    ax.grid(color='k', linestyle='solid', linewidth=0.5)
    tight_layout()
    
#     plt.show()
    plt.savefig('{0}2500kms_galaxies-{1}.pdf'.format(saveDirectory, include_galaxies),format='pdf')

##########################################################################################
    # 2500 - 5000 km/s

    # Now plot the data in Aitoff projection with a grid.
    fig = plt.figure(figsize=(12,6))
    ax = lab.subplot(111,projection="aitoff")
    
    colors = numpy.array(v2)
    vmaxVal = max(colors)
    vminVal = min(colors)

    norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)
    
    
    # first galaxies
    if include_galaxies:
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

        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
    
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
    
    
    cbar = plt.colorbar(plot2,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
    cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=530)
    ax.grid(color='k', linestyle='solid', linewidth=0.5)
    tight_layout()
    
#     plt.show()
    plt.savefig('{0}5000kms_galaxies-{1}.pdf'.format(saveDirectory, include_galaxies),format='pdf')
    
##########################################################################################
    # 5000 - 7500 km/s

    # Now plot the data in Aitoff projection with a grid.
    fig = plt.figure(figsize=(12,6))
    ax = lab.subplot(111,projection="aitoff")
    
    colors = numpy.array(v3)
    vmaxVal = max(colors)
    vminVal = min(colors)

    norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)
    
    # first galaxies
    if include_galaxies:
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

        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
    
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
    
    
    cbar = plt.colorbar(plot2,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
    cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=530)
    ax.grid(color='k', linestyle='solid', linewidth=0.5)
    tight_layout()
    
#     plt.show()
    plt.savefig('{0}7500kms_galaxies-{1}.pdf'.format(saveDirectory, include_galaxies),format='pdf')

##########################################################################################
    # 7500 - 10000 km/s
    
    # Now plot the data in Aitoff projection with a grid.
    fig = plt.figure(figsize=(12,6))
    ax = lab.subplot(111,projection="aitoff")
    
    colors = numpy.array(v4)
    vmaxVal = max(colors)
    vminVal = min(colors)

    norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)
    
    # first galaxies
    if include_galaxies:
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

        # lab.subplot(111)

        plot1 = plt.scatter(ra_rad, dec_rad,marker='o',c=colors,vmin=vminVal,\
        vmax=vmaxVal, lw=0, cmap=colmap, s=size_galaxy, alpha=alpha_galaxy)
    
    
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
    
    
    cbar = plt.colorbar(plot2,format=r'$\rm %d$',cmap=colmap,orientation='vertical',fraction=0.024, pad=0.03)
    cbar.set_label(r'$\rm Heliocentric Velocity ~[km ~s^{-1}]$')
    
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=530)
    ax.grid(color='k', linestyle='solid', linewidth=0.5)
    tight_layout()
    
#     plt.show()
    plt.savefig('{0}10000kms_galaxies-{1}.pdf'.format(saveDirectory, include_galaxies),format='pdf')

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
