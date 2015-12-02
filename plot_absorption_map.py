#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_absorption_map.py, v 1.0 09/11/2015

Plot absorption properties on a map showing their distribution around a central galaxy


'''

import sys
import os
import csv

from pylab import *
from scipy import stats
# import atpy
from math import *
from utilities import *
import getpass
import pickle

# from astropy.io.votable import parse,tree
# import numpy as np
# import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

from matplotlib import rc


def returnLinDiameters(major,minor,distance):
    # input major and minor in arcsec, distance in Mpc
    # outputs major and minor in kpc
    newMajor = math.tan(math.radians(float(major)))*(distance*1000)
    newMinor = math.tan(math.radians(float(minor)))*(distance*1000)
    return (newMajor,newMinor)
    

def returnAngDiameters(major,minor,distance):
    # input distances in mpc, major and minor is in kpc
    # outputs angular diameters in arcsec
    newMajor = math.atan((float(major)/1000)/float(distance))*(1/3600)
    newMinor = math.atan((float(minor)/1000)/float(distance))*(1/3600)
    return (newMajor,newMinor)
    
    

###########################################################################

    
def main():
        
    if getpass.getuser() == 'David':
        pickleFilename = '/Users/David/Research_Documents/inclination/pilotData2.p'
        saveDirectory = '/Users/David/Research_Documents/inclination/pilot_paper/figures'

    elif getpass.getuser() == 'frenchd':
        pickleFilename = '/usr/users/frenchd/inclination/pilotData2.p'
        saveDirectory = '/usr/users/frenchd/inclination/pilot_paper/figures'

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    pickleFile = open(pickleFilename,'rU')
    fullDict = pickle.load(pickleFile)
    
    pickleFile.close()
    
    
    # save each plot?
    save = False
    
    
    # overall structure: fullDict is a dictionary with all the lines and their data in it
    # separated into 'associated' and 'ambiguous' as the two keys. Associated contains
    # all the lists of data for lines associated with a galaxy. Ambiguous contains all
    # the lists of data for lines not unambiguously associated (could be many galaxies
    # or none)
    
##########################################################################################
##########################################################################################
    # all the lists to be used for associated lines
    
    lyaVList = fullDict['lyaVList']
    lyaWList = fullDict['lyaWList']
    lyaErrorList = fullDict['lyaErrorList']
    naList = fullDict['naList']
    bList = fullDict['bList']
    impactList = fullDict['impactList']
    azList = fullDict['azList']
    newAzList = fullDict['newAzList']
    incList = fullDict['incList']
    fancyIncList = fullDict['fancyIncList']
    cosIncList = fullDict['cosIncList']
    fancyCosIncList = fullDict['fancyCosIncList']
    paList = fullDict['paList']
    vcorrList = fullDict['vcorrList']
    majList = fullDict['majList']
    difList = fullDict['difList']
    envList = fullDict['envList']
    morphList = fullDict['morphList']
    galaxyNameList = fullDict['galaxyNameList']
    raList = fullDict['raList']
    decList = fullDict['decList']
    
    
    lyaV_blue = []
    lyaV_red = []
    lyaW_blue = []
    lyaW_red = []
    lyaErr_blue = []
    lyaErr_red = []
    na_blue = []
    na_red = []
    b_blue = []
    b_red = []
    impact_blue = []
    impact_red = []
    az_blue = []
    az_red = []
    newAz_blue = []
    newAz_red = []
    inc_blue = []
    inc_red = []
    fancyInc_blue = []
    fancyInc_red = []
    cosInc_blue = []
    cosInc_red = []
    fancyCosInc_blue = []
    fancyCosInc_red = []
    pa_blue = []
    pa_red = []
    vcorr_blue = []
    vcorr_red = []
    maj_blue = []
    maj_red = []
    dif_blue = []
    dif_red = []
    env_blue = []
    env_red = []
    morph_blue = []
    morph_red = []
    ra_blue = []
    ra_red = []
    dec_blue = []
    dec_red = []
    
    c = -1
    for d in difList:
        c +=1
        if d > 0:
            # blueshifted absorption
            lyaV_blue.append(lyaVList[c])
            lyaW_blue.append(lyaWList[c])
            lyaErr_blue.append(lyaErrorList[c])
            na_blue.append(naList[c])
            b_blue.append(bList[c])
            impact_blue.append(impactList[c])
            az_blue.append(azList[c])
            newAz_blue.append(newAzList[c])
            inc_blue.append(incList[c])
            fancyInc_blue.append(fancyIncList[c])
            cosInc_blue.append(cosIncList[c])
            fancyCosInc_blue.append(fancyCosIncList[c])
            pa_blue.append(paList[c])
            vcorr_blue.append(vcorrList[c])
            maj_blue.append(majList[c])
            dif_blue.append(difList[c])
            env_blue.append(envList[c])
            morph_blue.append(morphList[c])
            ra_blue.append(raList[c])
            dec_blue.append(decList[c])
            
        else:
            # redshifted absorption
            lyaV_red.append(lyaVList[c])
            lyaW_red.append(lyaWList[c])
            lyaErr_red.append(lyaErrorList[c])
            na_red.append(naList[c])
            b_red.append(bList[c])
            impact_red.append(impactList[c])
            az_red.append(azList[c])
            newAz_red.append(newAzList[c])
            inc_red.append(incList[c])
            fancyInc_red.append(fancyIncList[c])
            cosInc_red.append(cosIncList[c])
            fancyCosInc_red.append(fancyCosIncList[c])
            pa_red.append(paList[c])
            vcorr_red.append(vcorrList[c])
            maj_red.append(majList[c])
            dif_red.append(difList[c])
            env_red.append(envList[c])
            morph_red.append(morphList[c])
            ra_red.append(raList[c])
            dec_red.append(decList[c])
    
        
        
##########################################################################################
##########################################################################################
    # all the lists to be used for ambiguous lines
    
    lyaVListAmb = fullDict['lyaVListAmb']
    lyaWListAmb = fullDict['lyaWListAmb']
    lyaErrorListAmb = fullDict['lyaErrorListAmb']
    naListAmb = fullDict['naListAmb']
    bListAmb = fullDict['bListAmb']
    impactListAmb = fullDict['impactListAmb']
    azListAmb = fullDict['azListAmb']
    newAzListAmb = fullDict['newAzListAmb']
    incListAmb = fullDict['incListAmb']
    fancyIncListAmb = fullDict['fancyIncListAmb']
    cosIncListAmb = fullDict['cosIncListAmb']
    fancyCosIncListAmb = fullDict['fancyCosIncListAmb']
    paListAmb = fullDict['paListAmb']
    vcorrListAmb = fullDict['vcorrListAmb']
    majListAmb = fullDict['majListAmb']
    difListAmb = fullDict['difListAmb']
    envListAmb = fullDict['envListAmb']
    morphListAmb = fullDict['morphListAmb']
    galaxyNameListAmb = fullDict['galaxyNameListAmb']
    raListAmb = fullDict['raListAmb']
    decListAmb = fullDict['decListAmb']
    
    lyaV_blueAmb = []
    lyaV_redAmb = []
    lyaW_blueAmb = []
    lyaW_redAmb = []
    lyaErr_blueAmb = []
    lyaErr_redAmb = []
    na_blueAmb = []
    na_redAmb = []
    b_blueAmb = []
    b_redAmb = []
    impact_blueAmb = []
    impact_redAmb = []
    az_blueAmb = []
    az_redAmb = []
    newAz_blueAmb = []
    newAz_redAmb = []
    inc_blueAmb = []
    inc_redAmb = []
    fancyInc_blueAmb = []
    fancyInc_redAmb = []
    cosInc_blueAmb = []
    cosInc_redAmb = []
    fancyCosInc_blueAmb = []
    fancyCosInc_redAmb = []
    pa_blueAmb = []
    pa_redAmb = []
    vcorr_blueAmb = []
    vcorr_redAmb = []
    maj_blueAmb = []
    maj_redAmb = []
    dif_blueAmb = []
    dif_redAmb = []
    env_blueAmb = []
    env_redAmb = []
    morph_blueAmb = []
    morph_redAmb = []
    ra_blueAmb = []
    ra_redAmb = []
    dec_blueAmb = []
    dec_redAmb = []
    
    
    c = -1
    for d in difListAmb:
        c +=1
        if d > 0:
            # blueshifted absorption
            lyaV_blueAmb.append(lyaVListAmb[c])
            lyaW_blueAmb.append(lyaWListAmb[c])
            lyaErr_blueAmb.append(lyaErrorListAmb[c])
            na_blueAmb.append(naListAmb[c])
            b_blueAmb.append(bListAmb[c])
            impact_blueAmb.append(impactListAmb[c])
            az_blueAmb.append(azListAmb[c])
            newAz_blueAmb.append(newAzListAmb[c])
            inc_blueAmb.append(incListAmb[c])
            fancyInc_blueAmb.append(fancyIncListAmb[c])
            cosInc_blueAmb.append(cosIncListAmb[c])
            fancyCosInc_blueAmb.append(fancyCosIncListAmb[c])
            pa_blueAmb.append(paListAmb[c])
            vcorr_blueAmb.append(vcorrListAmb[c])
            maj_blueAmb.append(majListAmb[c])
            dif_blueAmb.append(difListAmb[c])
            env_blueAmb.append(envListAmb[c])
            morph_blueAmb.append(morphListAmb[c])
            ra_blueAmb.append(raListAmb[c])
            dec_blueAmb.append(decListAmb[c])
            
        else:
            # redshifted absorption
            lyaV_redAmb.append(lyaVListAmb[c])
            lyaW_redAmb.append(lyaWListAmb[c])
            lyaErr_redAmb.append(lyaErrorListAmb[c])
            na_redAmb.append(naListAmb[c])
            b_redAmb.append(bListAmb[c])
            impact_redAmb.append(impactListAmb[c])
            az_redAmb.append(azListAmb[c])
            newAz_redAmb.append(newAzListAmb[c])
            inc_redAmb.append(incListAmb[c])
            fancyInc_redAmb.append(fancyIncListAmb[c])
            cosInc_redAmb.append(cosIncListAmb[c])
            fancyCosInc_redAmb.append(fancyCosIncListAmb[c])
            pa_redAmb.append(paListAmb[c])
            vcorr_redAmb.append(vcorrListAmb[c])
            maj_redAmb.append(majListAmb[c])
            dif_redAmb.append(difListAmb[c])
            env_redAmb.append(envListAmb[c])
            morph_redAmb.append(morphListAmb[c])
            ra_redAmb.append(raListAmb[c])
            dec_redAmb.append(decListAmb[c])
    
##########################################################################################
##########################################################################################
    
    # lists for the full galaxy dataset
    
    allPA = fullDict['allPA']
    allInclinations = fullDict['allInclinations']
    allCosInclinations = fullDict['allCosInclinations']
    allFancyInclinations = fullDict['allFancyInclinations']
    allFancyCosInclinations = fullDict['allCosFancyInclinations']
    
    total = 0
    totalNo = 0
    totalYes = 0
    totalIsolated = 0
    totalGroup = 0
    

########################################################################################
########################################################################################

    # plot histograms of the inclinations for both associated galaxies and the 
    # full galaxy data set, combining both redshifted and blueshifted
    plot_edge_on = True
    
    if plot_edge_on:
        
        # new coordinates of absorption for plotting (w/ galaxy at (0,0))
        yList = []
        xList = []
        
        # new coordinates normalized for diameter
        y_dList = []
        x_dList = []
        
        # new coordinates normalized for virial radius
        y_vList = []
        x_vList = []
        
        # calculate the position on the sky for each absorption feature wrt to the galaxy
        for r,d,i,a,fInc,maj in zip(raList,decList,impactList,newAzList,fancyIncList,majList):
            if float(a) >= 0 and float(fInc)>=0 and float(fInc)<=90:
                
                # y coordinate
                y = float(i) * sin(a*pi/180.)
                yList.append(y)
                
                # x coordinate
                x = float(i) * cos(a*pi/180.)
                xList.append(x)
                
                # normalize by diameter
                # y coordinate
                y_d = (float(i)/float(maj)) * sin(a*pi/180.)
                y_dList.append(y_d)

                # x coordinate
                x_d = (float(i)/float(maj)) * cos(a*pi/180.)
                x_dList.append(x_d)
                
                # normalize by virial radius
                # 
                # this to be completed later
                                        
                                        
            else:
                print 'float(a) <0: ',r,d,i,a,fInc

        # calculate the average red vs blue azimuth line
        blueAvg = mean(newAz_blue)
        print 'blueAvg: ',blueAvg
        redAvg = mean(newAz_red)
        print 'redAvg: ',redAvg
        print
        
        xyBlueAvg = (500.* cos(blueAvg * pi/180.), 500.* sin(blueAvg * pi/180.))
        xyRedAvg = (500.* cos(redAvg * pi/180.), 500.* sin(redAvg * pi/180.))
        print 'xyBlueAvg: ',xyBlueAvg
        print 'xyRedAvg: ',xyRedAvg

        
        # plot the distributions 
        fig = figure(figsize=(8,8))
        subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)
        
        ax1 = fig.add_subplot(111)
        # plot the flat galaxy line
        plot1 = plot((-10,10),(0,0),c='black',linewidth=5)
        
        maxW = 120
        minW = 5
        newLyaWList = []
        for w in lyaWList:
            newW = ((float(w) - min(lyaWList)) / (max(lyaWList) - min(lyaWList)))*(maxW - minW)
            newLyaWList.append(newW)

#         colmap = cm.RdBu
#         norm = matplotlib.colors.Normalize(vmin = minW, vmax = maxW)
#         m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)
        
        # plot the absorption features
        for x,y,s,d in zip(x_dList,y_dList,newLyaWList,difList):
            if d>0:
                # blueshifted
                plot2 = scatter(x,y,marker='*',color='blue',s=s)
            else:
                # redshifted
                plot3 = scatter(x,y,marker='*',color='red',s=s)

#         plot4 = plot((0,xyBlueAvg[0]), (0,xyBlueAvg[1]), color = 'blue')
#         plot5 = plot((0,xyRedAvg[0]), (0,xyRedAvg[1]), color = 'red')
        show()



    
###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    