#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotW_Az.py, v 5.1 04/21/2016

This is the plotW_Az_major bit from histograms3.py. Now is separated, and loads in a pickle
file of the relevant data, as created by "buildDataLists.py"

Plots equivalent width as a function of azimuth, also normalized by galaxy size, separated 
into red and blue shifted absorption samples

Previous (from histograms3.py):
    Plot some stuff for the 100largest initial results

    Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

    Updated for the pilot paper (05/06/15)

v5: updated to work with the new, automatically updated LG_correlation_combined5.csv
    (12/04/15)
    
v5.1: remake plots with v_hel instead of vcorr (4/21/16)
        
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

    

###########################################################################

    
def main():
    # assuming 'theFile' contains one name per line, read the file
    
    
    if getpass.getuser() == 'David':
        pickleFilename = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/pilotData2.p'
        resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_3.csv'
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots/'

    elif getpass.getuser() == 'frenchd':
        pickleFilename = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/pilotData2.p'
        resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_3.csv'
        saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots/'

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # use the old pickle file to get the full galaxy dataset info
    pickleFile = open(pickleFilename,'rU')
    fullDict = pickle.load(pickleFile)
    
    pickleFile.close()
    
    
    # save each plot?
    save = False
    
    results = open(resultsFilename,'rU')
    reader = csv.DictReader(results)
    
    virInclude = False
    cusInclude = False
    finalInclude = True
    
    # if match, then the includes in the file have to MATCH the includes above. e.g., if 
    # virInclude = False, cusInclude = True, finalInclude = False, then only systems
    # matching those three would be included. Otherwise, all cusInclude = True would be included
    # regardless of the others
    match = False
    
    # all the lists to be used for associated lines
    lyaVList = []
    lyaWList = []
    lyaErrList = []
    naList = []
    bList = []
    impactList = []
    azList = []
    incList = []
    fancyIncList = []
    cosIncList = []
    cosFancyIncList = []
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
    
    # for ambiguous lines
    lyaVAmbList = []
    lyaWAmbList = []
    envAmbList = []
    
    for l in reader:
        include_vir = eval(l['include_vir'])
        include_cus = eval(l['include_custom'])
        include = eval(l['include'])
        
        go = False
        if match:
            if virInclude == include_vir and cusInclude == include_cus:
                go = True
            else:
                go = False
                
        else:
            if virInclude and include_vir:
                go = True
                
            elif cusInclude and include_cus:
                go = True
                
            elif finalInclude and include:
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
            
            if isNumber(inc):
                cosInc = cos(float(inc) * pi/180.)
                
                if isNumber(maj) and isNumber(min):
                    q0 = 0.2
                    fancyInc = calculateFancyInclination(maj,min,q0)
                    cosFancyInc = cos(fancyInc * pi/180)
                else:
                    fancyInc = -99
                    cosFancyInc = -99
            else:
                cosInc = -99
                inc = -99
                fancyInc = -99
                cosFancyInc = -99
        
            if isNumber(pa):
                pa = float(pa)
            elif isNumber(RC3pa):
                pa = float(RC3pa)
            else:
                pa = -99
                
            if isNumber(az):
                az = float(az)
            else:
                az = -99
                
            if isNumber(maj):
                maj = float(maj)
                virialRadius = float(virialRadius)
            else:
                maj = -99
                virialRadius = -99
            
            # all the lists to be used for associated lines
            lyaVList.append(float(lyaV))
            lyaWList.append(float(lyaW))
            lyaErrList.append(float(lyaW_err))
            naList.append(na)
            bList.append(float(b))
            impactList.append(float(impact))
            azList.append(az)
            incList.append(float(inc))
            fancyIncList.append(fancyInc)
            cosIncList.append(cosInc)
            cosFancyIncList.append(cosFancyInc)
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
            
        else:
            lyaV = l['Lya_v']
            lyaW = l['Lya_W'].partition('pm')[0]
            lyaW_err = l['Lya_W'].partition('pm')[2]
            env = l['environment']
            
            lyaVAmbList.append(float(lyaV))
            lyaWAmbList.append(float(lyaW))
            envAmbList.append(float(env))

    results.close()
    
        
    # lists for the full galaxy dataset
    allPA = fullDict['allPA']
    allInclinations = fullDict['allInclinations']
    allCosInclinations = fullDict['allCosInclinations']
    allFancyInclinations = fullDict['allFancyInclinations']
    allCosFancyInclinations = fullDict['allCosFancyInclinations']
    
    total = 0
    totalNo = 0
    totalYes = 0
    totalIsolated = 0
    totalGroup = 0

########################################################################################
########################################################################################

    # plot equivalent width as a function of azimuth normalized by galaxy size, separated
    # into red and blue shifted absorption samples
    #
    
    plotW_Az_major = True
    save = True
    
    if plotW_Az_major:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        
        # give some stats:
        lessThan45 = 0
        for a in azList:
            if a <=45:
                lessThan45 +=1
        print '{0}/{1} have az <= 45 degrees'.format(lessThan45,len(azList))
        print 'average, median azimuth: ',average(azList),', ',median(azList)
        
        for d,a,w,m in zip(difList,azList,lyaWList,majList):
            # check if all the values are okay
            if isNumber(d) and isNumber(a) and isNumber(w) and isNumber(m):
                if d!=-99 and a!=-99 and w!=-99 and m!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(a/m,w,c='Blue',s=50,label= labelb)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(a/m,w,c='Red',s=50,label= labelr)
                
                    plot1 = scatter(a/m,w,c=color,s=50)
            
        title('W(azimuth/diameter) for red vs blue absorption')
        xlabel(r'Azimuth / Major Axis')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,max(lyaWList)+50)
        xlim(0,10)

        if save:
            savefig('{0}/W(azimuth_diameter)_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()


#########################################################################################
#########################################################################################

    # plot equivalent width as a function of azimuth normalized by galaxy size, separated
    # into red and blue shifted absorption samples
    #
    
    plotW_Az_vir = True
    save = True
    
    if plotW_Az_vir:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        
        # give some stats:
        lessThan45 = 0
        for a in azList:
            if a <=45:
                lessThan45 +=1
        print '{0}/{1} have az <= 45 degrees'.format(lessThan45,len(azList))
        print 'average, median azimuth: ',average(azList),', ',median(azList)
        
        for d,a,w,m in zip(difList,azList,lyaWList,virList):
            # check if all the values are okay
            if isNumber(d) and isNumber(a) and isNumber(w) and isNumber(m):
                if d!=-99 and a!=-99 and w!=-99 and m!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(a/m,w,c='Blue',s=50,label= labelb)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(a/m,w,c='Red',s=50,label= labelr)
                
                    plot1 = scatter(a/m,w,c=color,s=50)
            
        title('W(azimuth/R_vir) for red vs blue absorption')
        xlabel(r'$\rm Azimuth / R_{vir}$')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,max(lyaWList)+50)
        xlim(0,1)

        if save:
            savefig('{0}/W(azimuth_vir)_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()


#########################################################################################
#########################################################################################

    # plot equivalent width as a function of azimuth angle (old one) for red vs blue
    # shifted absorption
    #
    
    plotW_Az = True
    save = True
    
    if plotW_Az:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        print 'len(newAzList): ',len(azList)
        print 'len(azList): ',len(azList)
        
        for d,a,w,m in zip(difList,azList,lyaWList,majList):
            # check if all the values are good
            if isNumber(d) and isNumber(a) and isNumber(w) and isNumber(m):
                if d!=-99 and a!=-99 and w!=-99 and m!=-99:
                    count +=1
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(a,w,c='Blue',s=50,label= labelb)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(a,w,c='Red',s=50,label= labelr)
            
                    plot1 = scatter(a,w,c=color,s=50)
                    
        print 'countr: ',countr
        print 'countb: ',countb
        print 'count: ',count
        title('W(azimuth) for red vs blue shifted absorption')
        xlabel(r'Azimuth (deg)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,max(lyaWList)+50)
        xlim(0,90)

        if save:
            savefig('{0}/W(azimuth)_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()



###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    