#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotInc_vs_dif.py, v 4.0 05/13/2015

This is the plotInc_vs_dif bit from histograms3.py. Now is separated, and loads in a pickle
file of the relevant data, as created by "buildDataLists.py"


Previous (from histograms3.py):
    Plot some stuff for the 100largest initial results

    Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

    Updated for the pilot paper (05/06/15)


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
        pickleFilename = '/Users/David/Research_Documents/inclination/pilotData.p'
        resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5.csv'
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots/'

    elif getpass.getuser() == 'frenchd':
        pickleFilename = '/usr/users/inclination/pilotData.p'
        resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5.csv'
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

    # I have no idea what this is or was supposed to do.
    plotInc_vs_dif = True
    
    if plotInc_vs_dif:
        fig = figure()
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(211)
#         bins = [0,10,20,30,40,50,60,70,80,90]
#         bins = [5,15,25,35,45,55,65,75,85,95]

    #     bins = [5,15,25,35,45,55,65,75,85]
#         bins = [0,15,30,45,60,75,90]
#         bins = [0,30,60,90]
        bins = [0,30,60,90]
        incBlue = []
        incRed = []
        for i,d,l in zip(incList,difList,lyaWList):
            if d>0:
                incBlue.append(i)
            else:
                incRed.append(i)
                
        n, bins, patches = hist(incBlue, bins)
        setp(patches, 'facecolor', 'blue', 'alpha', 0.5)               
#         plot1 = hist(incBlue,bins=bins,histtype='bar',c='blue',alpha=0.5)
        title('Blue Shifted Absorbers')
        xlabel('Inclination (deg) / W (AA)')
        ylabel('Number')
#         xlim(0,90)
        
        ax = fig.add_subplot(212)
#         bins = [5,15,25,35,45,55,65,75,85,95]
    #     bins = [5,15,25,35,45,55,65,75,85]
#         bins = [0,15,30,45,60,75,90]
    
        n, bins, patches = hist(incRed, bins)
        setp(patches, 'facecolor', 'red', 'alpha', 0.5)
#         plot2 = hist(incRed,bins=bins,histtype='bar',c='red',alpha=0.5)
        title('Red Shifted Absorbers')
        xlabel('Inclination (deg) / W (AA)')
        ylabel('Number')
#         ylim(0,1)
#         xlim(0,90)
        tight_layout()


        # give me the stats:
        incList2 = []
        azList2 = []
        for i,a in zip(incList,azList):
            if i>=50:
                incList2.append(i)
            if a>50:
                azList2.append(a)
        print '{0} of {1} are inclined >=50%'.format(len(incList2),len(incList))
        print '{0} of {1} have azimuth >=50%'.format(len(azList2),len(azList))
        print 'a: ',a
        print 'len: ',len(incList)
        
        if save:
            savefig('{0}/inc_vs_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
    
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    