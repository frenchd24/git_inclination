#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotW_FancyCosInc_colorbar.py, v 4.0 05/18/2015

This is the plotW_FancyCosInc_colorbar bit from histograms3.py. Now is separated, and loads in a pickle
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

    # plot equivalent width as a function of cos(fancy-inclination) with red and blue shifted
    # absorption represented by a color bar
    plotW_FancyCosInc_colorbar= True
    
    if plotW_FancyCosInc_colorbar:
        # colormap the velocity difference of the absorber
        averaged =[]

        blueMap = cm.Blues
        redMap = cm.Reds
        
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
#         labelr = 'Red Shifted Absorber'
#         labelb = "Blue Shifted Absorber"
        rdif = []
        rcosInc = []
        rlyaW = []
        rMaj = []
        bdif = []
        bcosInc = []
        blyaW = []
        bMaj = []
        for d,i,w,m in zip(difList,fancyCosIncList,lyaWList,majList):
        
            # check if all the values are okay
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(m):
                if d !=-99 and i !=-99 and w!=-99 and m!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            bdif.append(d)
                            bcosInc.append(i)
                            blyaW.append(w)
                            bMaj.append(m)
        #                     plotb = ax.scatter(a, w, cmap=blueMap, c=d, s=50, vmin=0, vmax=400)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
        #                     plotr = ax.scatter(a, w, cmap=redMap, c=d, s=50, vmin=0, vmax=400)
                            rdif.append(d)
                            rcosInc.append(i)
                            rlyaW.append(w)
                            rMaj.append(m)


        print
        print 'average, median redshifts: ',average(rdif),', ',median(rdif)
        print 'average, median blueshifts: ',average(bdif),', ',median(bdif)
        print 'max, min red dif: ',max(rdif), ', ',min(rdif)
        print 'max, min blue dif: ',max(bdif), ', ', min(bdif)
        
        plotr = ax.scatter(rcosInc, rlyaW, cmap=redMap, c=rdif, s=50, vmin=0, vmax=-400)
        norm = matplotlib.colors.Normalize(vmin = 0, vmax = 400)
        cbarRed = plt.colorbar(plotr,cmap=redMap,orientation='vertical')
        cbarRed.set_label('vcorr - absorber velocity (km/s)')
        
        plotb = ax.scatter(bcosInc, blyaW, cmap=blueMap, c=bdif, s=50, vmin=0, vmax=400)
        cbarBlue = plt.colorbar(plotb,cmap=blueMap,orientation='vertical')
        cbarBlue.set_label('vcorr - absorber velocity (km/s)')
#       cbar.ax.set_yticklabels(ticks)
                    
        title('Equivalent width vs Cos(fancy_inclination)')
        xlabel(r'Cos(inclination) = b/a')
        ylabel(r'Equivalent Width ($\rm \AA$)')
#         legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(-0.02,1)

        if save:
            savefig('{0}/W(cos(fancy_inclination))_dif_colorbar.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
        
        
        
#         fig = figure()
#         ax = fig.add_subplot(211)
#         hist(rdif,color='red')
# 
#         ax = fig.add_subplot(212)
#         hist(bdif,color='blue')
#         show()
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    