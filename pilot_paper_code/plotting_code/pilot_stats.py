#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  pilot_stats.py, v 1.0 01/04/2016

Print out all the relevant stats on the dataset for the pilot paper

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
            minor = l['minorAxis (kpc)']
            inc = l['inclination (deg)']
            az = l['azimuth (deg)']
            b = l['b'].partition('pm')[0]
            b_err = l['b'].partition('pm')[2]
            na = eval(l['Na'].partition(' pm ')[0])
#             print "l['Na'].partition(' pm ')[2] : ",l['Na'].partition(' pm ')
            na_err = eval(l['Na'].partition(' pm ')[2])
            likelihood = l['likelihood']
            likelihoodm15 = l['likelihood_1.5']
            virialRadius = l['virialRadius']
            m15 = l['d^1.5']
            vel_diff = l['vel_diff']
            
            if isNumber(inc):
                cosInc = cos(float(inc) * pi/180.)
                
                if isNumber(maj) and isNumber(minor):
                    q0 = 0.2
                    fancyInc = calculateFancyInclination(maj,minor,q0)
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
#########################################################################################
    
    # print all the things
    #
    
    # absorber info lists
    blues = []
    reds = []
    blueW = []
    redW = []
    blueErr = []
    redErr = []
    blueV = []
    redV = []
    blueImpact = []
    redImpact = []
    
    # galaxy info lists
    blueInc = []
    redInc = []
    blueFancyInc = []
    redFancyInc = []
    blueAz = []
    redAz = []
    bluePA = []
    redPA = []
    blueVcorr = []
    redVcorr = []
    blueEnv = []
    redEnv = []
    blueVir = []
    redVir = []
    

    
    # for absorbers
    for d,w,e,v,i in zip(difList,lyaWList,lyaErrList,lyaVList,impactList):
        if d>=0:
            blues.append(float(d))
            blueW.append(float(w))
            blueErr.append(float(e))
            blueV.append(float(v))
            blueImpact.append(float(i))
        else:
            reds.append(float(d))
            redW.append(float(w))
            redErr.append(float(e))
            redV.append(float(v))
            redImpact.append(float(i))
            
    # for galaxies
    for d,inc,finc,az,pa,vcorr,e,vir in zip(difList,incList,fancyIncList,azList,paList,vcorrList,envList,virList):
        if d>=0:
            if inc !=-99:
                blueInc.append(float(inc))
            if finc !=-99:
                blueFancyInc.append(float(finc))
            if az !=-99:
                blueAz.append(float(az))
            if pa !=-99:
                bluePA.append(float(pa))
            if vcorr !=-99:
                blueVcorr.append(float(vcorr))
            blueEnv.append(float(e))
            if vir !=-99:
                blueVir.append(float(vir))
        else:
            if inc !=-99:
                redInc.append(float(inc))
            if finc !=-99:
                redFancyInc.append(float(finc))
            if az !=-99:
                redAz.append(float(az))
            if pa !=-99:
                redPA.append(float(pa))
            if vcorr !=-99:
                redVcorr.append(float(vcorr))
            redEnv.append(float(e))
            if vir !=-99:
                redVir.append(float(vir))
            
    print
    print '------------------------ Pilot Data ------------------------------'
    print
    print ' FOR THE FOLLOWING INCLUDE SET:'
    print ' Virial radius include = ',virInclude
    print ' Custom include =        ',cusInclude
    print ' Final include =         ',finalInclude
    print ' Match =                 ',match
    print
    print 'total number of lines: ', len(lyaWList) + len(lyaWAmbList)
    print 'total number of associated lines: ',len(difList)
    print '# of redshifted lines: ',len(reds)
    print reds
    print '# of blueshifted lines: ',len(blues)
    print
    print '----------------------- Absorber info ----------------------------'
    print
    print 'avg blueshifted EW: ',mean(blueW)
    print 'median blueshifted EW: ',median(blueW)

    print 'avg blueshifted vel_diff: ',mean(blues)
    print 'median blueshifted vel_diff: ',median(blues)

    print 'avg blue velocity: ',mean(blueV)
    print 'median blue velocity: ',median(blueV)

    print 'avg blue err: ',mean(blueErr)
    print 'median blue err: ',median(blueErr)
    
    print 'avg blue impact: ',mean(blueImpact)
    print 'median blue impact: ',median(blueImpact)

    print
    
    print 'avg redshifted EW: ',mean(redW)
    print 'median redshifted EW: ',median(redW)

    print 'avg redshifted vel_diff: ',mean(reds)
    print 'median redshifted vel_diff: ',median(reds)

    print 'avg red velocity: ',mean(redV)
    print 'median red velocity: ',median(redV)

    print 'avg red err: ',mean(redErr)
    print 'median red err: ',median(redErr)
    
    print 'avg red impact: ',mean(redImpact)
    print 'median red impact: ',median(redImpact)
    
    print
    print '----------------------- Galaxy info ----------------------------'
    print
    
    print 'avg blue inclination: ',mean(blueInc)
    print 'median blue inclination: ',median(blueInc)
    
    print 'avg blue fancy inclination: ',mean(blueFancyInc)
    print 'median blue fancy inclination: ',median(blueFancyInc)
    
    print 'avg blue azimuth: ',mean(blueAz)
    print 'median blue azimuth: ',median(blueAz)
    
    print 'avg blue PA: ',mean(bluePA)
    print 'median blue PA: ',median(bluePA)
    
    print 'avg blue vcorr: ',mean(blueVcorr)
    print 'median blue vcorr: ',median(blueVcorr)
    
    print 'avg blue environment: ',mean(blueEnv)
    print 'median blue environment: ',median(blueEnv)
    
    print 'avg blue R_vir: ',mean(blueVir)
    print 'median blue R_vir: ',median(blueVir)
    
    print

    print 'avg red inclination: ',mean(redInc)
    print 'median red inclination: ',median(redInc)
    
    print 'avg red fancy inclination: ',mean(redFancyInc)
    print 'median red fancy inclination: ',median(redFancyInc)
    
    print 'avg red azimuth: ',mean(redAz)
    print 'median red azimuth: ',median(redAz)
    
    print 'avg red PA: ',mean(redPA)
    print 'median red PA: ',median(redPA)
    
    print 'avg red vcorr: ',mean(redVcorr)
    print 'median red vcorr: ',median(redVcorr)
    
    print 'avg red environment: ',mean(redEnv)
    print 'median red environment: ',median(redEnv)
    
    print 'avg red R_vir: ',mean(redVir)
    print 'median red R_vir: ',median(redVir)
    
    print
    print 'Complete'
    
    

###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    