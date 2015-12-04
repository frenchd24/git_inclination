#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotAzHist.py, v 4.0 05/13/2015

This is the plotAzHist bit from histograms3.py. Now is separated, and loads in a pickle
file of the relevant data, as created by "buildDataLists.py"


Previous (from histograms3.py):
    Plot some stuff for the 100largest initial results

    Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

    Updated for the pilot paper (05/06/15)

v5: updated to work with the new, automatically updated LG_correlation_combined5.csv
    (12/04/15)

'''

import sys
import os
import csv

from pylab import *
# import atpy
from math import *
from utilities import *
import getpass
import pickle

# from astropy.io.votable import parse,tree

# from vo.table import parse
# import vo.tree


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# ## for Palatino and other serif fonts use:
# #rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)


    
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
    # assuming 'theFile' contains one name per line, read the file
    
    
    if getpass.getuser() == 'David':
        resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5.csv'
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots/'

    elif getpass.getuser() == 'frenchd':
        resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5.csv'
        saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots/'

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    
    # save each plot?
    save = False
    
    results = open(resultsFilename,'rU')
    reader = csv.DictReader(results)
    
    virInclude = False
    cusInclude = True
    finalInclude = False
    
    # if match, then the includes in the file have to MATCH the includes above. e.g., if 
    # virInclude = False, cusInclude = True, finalInclude = False, then only systems
    # matching those three would be included. Otherwise, all cusInclude = True would be included
    # regardless of the others
    match = False
    
    # all the lists to be used for associated lines
    lyaVList = []
    lyaWList = []
    naList = []
    bList = []
    impactList = []
    azList = []
    incList = []
    cosIncList = []
    paList = []
    vcorrList = []
    majList = []
    difList = []
    envList = []
    morphList = []
    m15List = []
    virList = []
    likeList = []
    likem15List = []
    
    
    for l in reader:
        include_vir = eval(l['include_vir'])
        include_cus = eval(l['include_custom'])
        include = l['include']
        
        go = False
        if match:
            if virInclude == include_vir and cusInclude == include_cus:
                go = True
            else:
                go = False
                
        else:
            if virInclude and include_vir:
                go = True
                
            if cusInclude and include_cus:
                go = True
            
            else:
                go = False
        
        if go:
            AGNra_dec = eval(l['degreesJ2000RA_DecAGN'])
            galaxyRA_Dec = eval(l['degreesJ2000RA_DecGalaxy'])
            lyaV = l['Lya_v']
            lyaW = l['Lya_W'].partition('pm')[0]
            lyaW_err = l['Lya_W'].partition('pm')[2]
            env = l['environment']
            galaxyName = l['galaxyName']
            impact = l['impactParameter (kpc)']
            galaxyDist = l['distGalaxy (Mpc)']
            pa = l['positionAngle (deg)']
            RC3pa = l['RC3pa (deg)']
            morph = l['morphology']
            vcorr = l['vcorrGalaxy (km/s)']
            maj = l['majorAxis (kpc)']
            min = l['minorAxis (kpc)']
            inc = l['inclination (deg)']
            az = l['azimuth (deg)']
            b = l['b'].partition('pm')[0]
            b_err = l['b'].partition('pm')[2]
            na = eval(l['Na'].partition(' pm ')[0])
            print "l['Na'].partition(' pm ')[2] : ",l['Na'].partition(' pm ')
            na_err = eval(l['Na'].partition(' pm ')[2])
            likelihood = l['likelihood']
            likelihoodm15 = l['likelihood_1.5']
            virialRadius = l['virialRadius']
            m15 = l['d^1.5']
            vel_diff = l['vel_diff']
            
            if isNumber(RC3pa) and not isNumber(pa):
                pa = RC3pa
            
            if isNumber(inc):
                cosInc = cos(float(inc) * pi/180.)
            else:
                cosInc = -99
                inc = -99
            
            # all the lists to be used for associated lines
            lyaVList.append(float(lyaV))
            lyaWList.append(float(lyaW))
            naList.append(na)
            bList.append(float(b))
            impactList.append(float(impact))
            azList.append(az)
            incList.append(float(inc))
            cosIncList.append(cosInc)
            paList.append(pa)
            vcorrList.append(vcorr)
            majList.append(maj)
            difList.append(float(vel_diff))
            envList.append(float(env))
            morphList.append(morph)
            m15List.append(m15)
            virList.append(virialRadius)
            likeList.append(likelihood)
            likem15List.append(likelihoodm15)

    results.close()
    
    total = 0
    totalNo = 0
    totalYes = 0
    totalIsolated = 0
    totalGroup = 0
    
    for a,n,g in zip(azList,newAzList,galaxyNameList):
        if a != -99 and n !=-99:
            if abs(a-n) >5:
                print "{0} : old = {1}, new = {2}".format(g,a,n)
                
    

########################################################################################

    # make a histogram of all the azimuth angles (old, hand measured ones)
    plotAzHist = False
    
    if plotAzHist:
        fig = figure()
        ax = fig.add_subplot(111)
        bins = [0,10,20,30,40,50,60,70,80,90]

        plot1 = hist(azList,bins=bins,histtype='bar')
        title('Distribution of old azimuths')
        xlabel('Azimuth (deg)')
        ylabel('Number')
        xlim(0,90)
        
        if save:
            savefig('{0}/hist(azimuth).pdf'.format(saveDirectory),format='pdf')
        else:
            show()
      
########################################################################################
 
    # make a histogram of all the azimuth angles (new, automatic ones)
    plotNewAzHist = False
    
    if plotNewAzHist:
        fig = figure()
        ax = fig.add_subplot(111)
        bins = [0,10,20,30,40,50,60,70,80,90]

        plot1 = hist(newAzList,bins=bins,histtype='bar')
        title("Distribution of new azimuths")
        xlabel('Azimuth (deg)')
        ylabel('Number')
        xlim(0,90)
        
        if save:
            savefig('{0}/hist(new_azimuth).pdf'.format(saveDirectory),format='pdf')
        else:
            show()

    
###############################################################################

    # plot both in the same window
    plotAzHist_both = True
    
    if plotAzHist_both:
        fig = figure()
        ax = fig.add_subplot(211)
        bins = [0,10,20,30,40,50,60,70,80,90]
#         bins = [0,5,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90]

#         bins = 20

        plot1 = hist(newAzList,bins=bins,histtype='bar')
        title("Distribution of new azimuths")
#         xlabel('Azimuth (deg)')
        ylabel('Number')
        xlim(0,90)
        ylim(0,10)
        
        ax = fig.add_subplot(212)
        plot1 = hist(azList,bins=bins,histtype='bar')
        title('Distribution of old azimuths')
        xlabel('Azimuth (deg)')
        ylabel('Number')
        xlim(0,90)
        ylim(0,10)
        
        
        
        if save:
            savefig('{0}/hist(both_azimuth).pdf'.format(saveDirectory),format='pdf')
        else:
            show()


###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    