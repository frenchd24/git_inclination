#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  pilot_stats.py, v 1.4 09/23/16

Print out all the relevant stats on the dataset for the pilot paper (01/04/2016)

v1.1: updates for v_hel velocity and probably something else? (2/22/16)

v1.2: updates for the new large galaxy sample (07/14/16) -> /plots4/

v1.3: added ability to limit results by environment variable (7/14/16)
    also same for likelihood values, including some stats about them also
    
v1.4: updated to LG_correlation_combined5_11_25cut_edit4.csv and plots5 (9/23/16)

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
    
    if getpass.getuser() == 'frenchd':

#         pickleFilename = '/Users/frenchd/Research/inclination/git_inclination/picklePilot_plusSALT_14.p'
#         gtPickleFilename = '/Users/frenchd/Research/inclination/git_inclination/pickleGT.p'
#         saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/plotting_code/figs'

#         gtPickleFilename = '/Users/frenchd/Research/inclination/git_inclination/pickleGT_filteredAll.p'
        gtPickleFilename = '/Users/frenchd/Research/GT_update2/pickleGT_filteredAll.p'

        saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/plotting_code/figs/'
        
        isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/isolated4.p'
        L_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_isolated4.p'
        L_associated_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated4.p'
        L_associated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated4.p'
        L_nonassociated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated4.p'
        L_two_filename = '/Users/frenchd/Research/inclination/git_inclination/L_two4.p'
        L_two_plus_filename = '/Users/frenchd/Research/inclination/git_inclination/L_two_plus4.p'
        L_group_filename = '/Users/frenchd/Research/inclination/git_inclination/L_group4.p'


    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # pickle file for the whole galaxy table:
    gtPickleFile = open(gtPickleFilename,'rU')
    gtDict = pickle.load(gtPickleFile)
    gtPickleFile.close()
    
    
    # open all the pickle files
    isolated_file = open(isolated_filename,'r')
    L_isolated_file = open(L_isolated_filename,'r')
    L_associated_isolated_file = open(L_associated_isolated_filename,'r')
    L_associated_file = open(L_associated_filename,'r')
    L_nonassociated_file = open(L_nonassociated_filename,'r')
    L_two_file = open(L_two_filename,'r')
    L_two_plus_file = open(L_two_plus_filename,'r')
    L_group_file = open(L_group_filename,'r')

    # unload the data from them
    isolated = pickle.load(isolated_file)
    L_isolated = pickle.load(L_isolated_file)
    L_associated_isolated = pickle.load(L_associated_isolated_file)
    L_associated = pickle.load(L_associated_file)
    L_nonassociated = pickle.load(L_nonassociated_file)
    L_two = pickle.load(L_two_file)
    L_two_plus = pickle.load(L_two_plus_file)
    L_group = pickle.load(L_group_file)
    
    # close the files
    isolated_file.close()
    L_isolated_file.close()
    L_associated_isolated_file.close()
    L_associated_file.close()
    L_nonassociated_file.close()
    L_two_file.close()
    L_two_plus_file.close()
    L_group_file.close()
    
    
    # if match, then the includes in the file have to MATCH the includes above. e.g., if 
    # virInclude = False, cusInclude = True, finalInclude = False, then only systems
    # matching those three would be included. Otherwise, all cusInclude = True would be included
    # regardless of the others
    
    dataSet = L_associated_isolated
    
    Lya_vs = dataSet['Lya_vs']
    e_Lya_vs = dataSet['e_Lya_vs']
    Lya_Ws = dataSet['Lya_Ws']
    e_Lya_Ws = dataSet['e_Lya_Ws']
    Nas = dataSet['Nas']
    e_Nas = dataSet['e_Nas']
    bs = dataSet['bs']
    e_bs = dataSet['e_bs']
    Ws = dataSet['Ws']
    e_Ws = dataSet['e_Ws']
    targets = dataSet['targets']
    z_targets = dataSet['z_targets']
    RA_targets = dataSet['RA_targets']
    Dec_targets = dataSet['Dec_targets']
    Names = dataSet['Names']
    RA_galaxies = dataSet['RA_galaxies']
    Dec_galaxies = dataSet['Dec_galaxies']
    impacts = dataSet['impacts']
    azimuths = dataSet['azimuths']
    PAs = dataSet['PAs']
    incs = dataSet['incs']
    adjustedIncs = dataSet['adjustedIncs']
    ls = dataSet['ls']
    l_cuss = dataSet['l_cuss']
    R_virs = dataSet['R_virs']
    cuss = dataSet['cuss']
    MajDiams = dataSet['MajDiams']
    MTypes = dataSet['MTypes']
    Vhels = dataSet['Vhels']
    vcorrs = dataSet['vcorrs']
    bestDists = dataSet['bestDists']
    e_bestDists = dataSet['e_bestDists']
    group_nums = dataSet['group_nums']
    group_mems = dataSet['group_mems']
    group_dists = dataSet['group_dists']
    Lstar_meds = dataSet['Lstar_meds']
    e_Lstar_meds = dataSet['e_Lstar_meds']
    Bmags = dataSet['Bmags']



    majorAxisL = gtDict['majorAxis']
    incL = gtDict['inc']
    adjustedIncL = gtDict['adjustedInc']
    paL = gtDict['PA']
    BmagL = gtDict['Bmag']
#     Bmag_sdssL = gtDict['Bmag_sdss']
    RID_medianL = gtDict['RID_median']
    RID_meanL = gtDict['RID_mean']
    RID_stdL = gtDict['RID_std']
    VhelL = gtDict['Vhel']
    RAdegL = gtDict['RAdeg']
    DEdegL = gtDict['DEdeg']
    NameL= gtDict['Name']
    
    allPA = paL
    allInclinations = []
    allAdjustedIncs = []
    allCosInclinations = []

#     print 'type: ',type(incL)
    for i in incL:
        if i != -99:
            i = float(i)
            allInclinations.append(i)
            
            i2 = pi/180. * i
            cosi2 = cos(i)
            allCosInclinations.append(cosi2)
            
    allCosFancyCosInclinations = []
    for i in adjustedIncL:
        if str(i) != '-99':
            i = float(i)

            allAdjustedIncs.append(i)
            
            i2 = pi/180. * i
            cosi2 = cos(i)
            allCosFancyCosInclinations.append(cosi2)
            
    allDiameter = majorAxisL

    print 'finished with this shit'
    print 'len(allAdjustedIncs): ',len(allAdjustedIncs)
    print

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
    blueAbs = []
    redAbs = []
    blueW = []
    redW = []
    blueB = []
    redB = []
    e_blueB = []
    e_redB = []
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
    blueVir = []
    redVir = []
    blueLike = []
    redLike = []
    
    
    # for absorbers
    for Lya_v, w, e_w, Vhel, i, b, e_b in zip(Lya_vs, Lya_Ws, e_Lya_Ws, Vhels, impacts, bs, e_bs):
        vel_dif = Lya_v - Vhel
        if vel_dif >=0:
            reds.append(float(vel_dif))
            redW.append(float(w))
            redErr.append(float(e_w))
            redV.append(float(Vhel))
            redImpact.append(float(i))
            redAbs.append(abs(vel_dif))
            redB.append(float(b))
            e_redB.append(float(e_b))

        else:

            blues.append(float(vel_dif))
            blueW.append(float(w))
            blueErr.append(float(e_w))
            blueV.append(float(Vhel))
            blueImpact.append(float(i))
            blueAbs.append(abs(vel_dif))
            blueB.append(float(b))
            e_blueB.append(float(e_b))

            
##########################################################################################
##########################################################################################
                   
            
    nameDict = {}
    # for galaxies
    for Lya_v, Vhel, inc, adjustedInc, az, pa, vcorr, vir, l, name in zip(Lya_vs, Vhels, incs, adjustedIncs, azimuths, PAs, vcorrs, R_virs, ls, Names):
        vel_dif = Lya_v - Vhel

        if nameDict.has_key(name):
            i = nameDict[name]
            i+=1
            nameDict[name] = i
        else:
            nameDict[name] = 1
        
        if vel_dif>=0:
            if inc !=-99:
                redInc.append(float(inc))
            if adjustedInc !=-99:
                redFancyInc.append(float(adjustedInc))
            if az !=-99:
                redAz.append(float(az))
            if pa !=-99:
                redPA.append(float(pa))
            if vcorr !=-99:
                redVcorr.append(float(vcorr))
            if vir !=-99:
                redVir.append(float(vir))
            if l !=-99:
                redLike.append(float(l))
        else:

            if inc !=-99:
                blueInc.append(float(inc))
            if adjustedInc !=-99:
                blueFancyInc.append(float(adjustedInc))
            if az !=-99:
                blueAz.append(float(az))
            if pa !=-99:
                bluePA.append(float(pa))
            if vcorr !=-99:
                blueVcorr.append(float(vcorr))
            if vir !=-99:
                blueVir.append(float(vir))
            if l !=-99:
                blueLike.append(float(l))
                
                
    galaxyNames = nameDict.keys()
                
    # how many absorbers above vs below vel_cut?
    redVelCount200 = 0
    redVelCount100 = 0
    blueVelCount200 = 0
    blueVelCount100 = 0
    
    for b in blues:
        if b >=200:
            blueVelCount200 +=1
        if b >= 100:
            blueVelCount100 +=1
        
    for r in reds:
        if abs(r) >=200:
            redVelCount200 +=1
        if abs(r) >=100:
            redVelCount100 +=1
    

    assocFancyInc = adjustedIncs
    
    AGNnameDict = {}
    for i in targets:
        if AGNnameDict.has_key(i):
            c = AGNnameDict[i]
            c +=1
            AGNnameDict[i] = c
        else:
            AGNnameDict[i] = 1
        
    AGN_list = AGNnameDict.keys()


    # write out a file breaking down all this shit
#     out_directory = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/'
#     save_name = 'full_stats.txt'
#     stats_filename = '{0}/{1}'.format(out_directory, save_name)
#     stats_file = open(stats_filename,'wt')
    
    
    print
    print '------------------------ Pilot Data -----------------------------'
    print
#     print 'total number of lines: ', len(lyaWList) + len(lyaWAmbList)
    print 'total number of lines: ', len(Lya_vs)
    print 'total number of unique galaxies matched: ',len(galaxyNames)
    print 'total number of AGN: ',len(AGN_list)
    print '# of redshifted lines: ',len(reds)
    print '# of blueshifted lines: ',len(blues)
    print
    print
    print ' TARGETS '
    print
    print 'final target number: ',len(AGNnameDict.keys())
    for i in AGNnameDict.keys():
        print i
    print
    print
    print
    print
    print '----------------------- Absorber info ----------------------------'
    print
    print 'avg blueshifted EW: ',mean(blueW)
    print 'median blueshifted EW: ',median(blueW)
    print 'avg blue err: ',mean(blueErr)
    print 'median blue err: ',median(blueErr)
    print
    print 'std(blue EW): ',std(blueW)
    print 'stats.sem(blue EW): ',stats.sem(blueW)
    print 'stats.describe(blue EW): ',stats.describe(blueW)
    print
    print 'avg blueshifted vel_diff: ',mean(blues)
    print 'median blueshifted vel_diff: ',median(blues)
    print 'std(blueshifted vel_diff): ',std(blues)
    print 'stats.sem(blue vel_dif): ',stats.sem(blues)
    print 'stats.describe(blue vel_dif: ',stats.describe(blues)
    print
    print '% blueshifted which have vel_diff >= 200 km/s: {0}'.format(float(blueVelCount200)/len(blues))
    print 'total number with abs(vel_diff) >= 200 km/s: {0}'.format(blueVelCount200)
    print '% blueshifted which have vel_diff >= 100 km/s: {0}'.format(float(blueVelCount100)/len(blues))
    print 'total number with abs(vel_diff) >= 100 km/s: {0}'.format(blueVelCount100)
    print
    print
    print 'avg blue velocity: ',mean(blueV)
    print 'median blue velocity: ',median(blueV)
    print 'std(blue Velocity): ',std(blueV)
    print 'avg blue impact: ',mean(blueImpact)
    print 'median blue impact: ',median(blueImpact)
    print 'stats.sem(blue impact): ',stats.sem(blueImpact)
    print 'stats.describe(blue impact): ',stats.describe(blueImpact)

    print
    print
    
    print 'avg redshifted EW: ',mean(redW)
    print 'median redshifted EW: ',median(redW)
    print 'avg red err: ',mean(redErr)
    print 'median red err: ',median(redErr)
    print
    print 'std(red EW): ',std(redW)
    print 'stats.sem(red EW): ',stats.sem(redW)
    print 'stats.describe(red EW): ',stats.describe(redW)

    print
    print 'avg redshifted vel_diff: ',mean(reds)
    print 'median redshifted vel_diff: ',median(reds)
    print 'std(redshifted vel_dif): ',std(reds)
    print 'stats.sem(red vel_dif): ',stats.sem(reds)
    print 'stats.describe(red vel_dif): ',stats.describe(reds)
    print
    print '% redshifted which have abs(vel_diff) >= 200 km/s: {0}'.format(float(redVelCount200)/len(reds))
    print 'total number with abs(vel_diff) >= 200 km/s: {0}'.format(redVelCount200)
    print '% redshifted which have abs(vel_diff) >= 100 km/s: {0}'.format(float(redVelCount100)/len(reds))
    print 'total number with abs(vel_diff) >= 100 km/s: {0}'.format(redVelCount100)
    print

    print 'avg red velocity: ',mean(redV)
    print 'median red velocity: ',median(redV)
    print
    print 'avg red impact: ',mean(redImpact)
    print 'median red impact: ',median(redImpact)
    print 'stats.sem(red impact): ',stats.sem(redImpact)
    print 'stats.describe(red impact): ',stats.describe(redImpact)
    print 'std(red impact): ',std(redImpact)



    print
    print '----------------------- Galaxy info ----------------------------'
    print
    
    # regular inclinations
    incCut = 50
    totalBlueInc = len(blueInc)
    totalRedInc = len(redInc)
    
    print
    print
    print 'totalBlueInc: ',totalBlueInc
    print 'totalRedInc: ',totalRedInc
    print
    print "blueInc: ",blueInc
    print
    
    blueIncCount = 0
    for i in blueInc:
        if i >= incCut:
            blueIncCount +=1
            
    redIncCount = 0
    for i in redInc:
        if i >= incCut:
            redIncCount +=1
            
    totalInc = len(allInclinations)
    totalCount = 0
    for i in allInclinations:
        if i >= incCut:
            totalCount +=1
            
            
    # fancy inclinations
    totalBlueFancyInc = len(blueFancyInc)
    totalRedFancyInc = len(redFancyInc)
    
    blueFancyIncCount = 0
    for i in blueFancyInc:
        if i >= incCut:
            blueFancyIncCount +=1
            
    redFancyIncCount = 0
    for i in redFancyInc:
        if i >= incCut:
            redFancyIncCount +=1
            
    combinedCount = redFancyIncCount + blueFancyIncCount
    totalCombinedCount = totalRedFancyInc + totalBlueFancyInc
            
    totalFancyInc = len(allAdjustedIncs)
    totalFancyCount = 0
    for i in allAdjustedIncs:
        if i >= incCut:
            totalFancyCount +=1
    
    print
    print ' INCLINATIONS: '
    print 
    print 'Blue: {0} % of associated galaxies have >={1}% inclination'.format(float(blueIncCount)/float(totalBlueInc),incCut)
    print 'Red: {0} % of associated galaxies have >={1}% inclination'.format(float(redIncCount)/float(totalRedInc),incCut)
    print 'All: {0} % of ALL galaxies have >={1}% inclination'.format(float(totalCount)/float(totalInc),incCut)
    print
    print ' FANCY INCLINATIONS: '
    print
    print 'Blue: {0} % of associated galaxies have >={1}% fancy inclination'.format(float(blueFancyIncCount)/float(totalBlueFancyInc),incCut)
    print 'Red: {0} % of associated galaxies have >={1}% fancy inclination'.format(float(redFancyIncCount)/float(totalRedFancyInc),incCut)
    print 'All: {0} % of ALL galaxies have >={1}% fancy inclination'.format(float(totalFancyCount)/float(totalFancyInc),incCut)
    print 'Combined: {0} % of associated galaxies have >= {1} fancy inclination'.format(float(combinedCount)/float(totalCombinedCount),incCut)
    print
    print 'Average all fancy inclination: ',mean(allAdjustedIncs)
    print 'stats.sem(all): ',stats.sem(allAdjustedIncs)
    print    
    print 'avg blue inclination: ',mean(blueInc)
    print 'median blue inclination: ',median(blueInc)
    print 'avg blue fancy inclination: ',mean(blueFancyInc)
    print 'median blue fancy inclination: ',median(blueFancyInc)
    print
    print 'avg red inclination: ',mean(redInc)
    print 'median red inclination: ',median(redInc)
    print 'avg red fancy inclination: ',mean(redFancyInc)
    print 'median red fancy inclination: ',median(redFancyInc)
    
    print
    print 'mean associated: ',mean(assocFancyInc)
    print 'stats.sem(associated): ',stats.sem(assocFancyInc)
    print 'stats.describe(associated): ',stats.describe(assocFancyInc)
    print 'stats.sem(blue): ',stats.sem(blueFancyInc)
    print 'stats.describe(blue): ',stats.describe(blueFancyInc)
    print
    print 'stats.sem(red): ',stats.sem(redFancyInc)
    print 'stats.describe(red): ',stats.describe(redFancyInc)
    print
    print
    print "  AZIMUTHS and PA:  "
    print
    print 'avg blue azimuth: ',mean(blueAz)
    print 'median blue azimuth: ',median(blueAz)
    print 'stats.sem(blue az): ',stats.sem(blueAz)
    print 'stats.describe(blue az): ',stats.describe(blueAz)
    print
    print 'avg red azimuth: ',mean(redAz)
    print 'median red azimuth: ',median(redAz)
    print 'stats.sem(red az): ',stats.sem(redAz)
    print 'stats.describe(red az): ',stats.describe(redAz)
    print
    print 'avg blue PA: ',mean(bluePA)
    print 'median blue PA: ',median(bluePA)
    print
    print 'avg red PA: ',mean(redPA)
    print 'median red PA: ',median(redPA)
    
    print
    print ' VCORR : '
    print
    print 'avg blue vcorr: ',mean(blueVcorr)
    print 'median blue vcorr: ',median(blueVcorr)
    print
    print 'avg red vcorr: ',mean(redVcorr)
    print 'median red vcorr: ',median(redVcorr)
    
    
    print
    print ' R_vir: '
    print
    print 'avg blue R_vir: ',mean(blueVir)
    print 'median blue R_vir: ',median(blueVir)
    print 'stats.sem(blue R_vir): ',stats.sem(blueVir)
    print 'stats.describe(blue R_vir): ',stats.describe(blueVir)
    print
    print 'avg red R_vir: ',mean(redVir)
    print 'median red R_vir: ',median(redVir)
    print 'stats.sem(red R_vir): ',stats.sem(redVir)
    print 'stats.describe(red R_vir): ',stats.describe(redVir)

    print
    print ' LIKELIHOOD: '
    print
    print 'avg blue likelihood: ',mean(blueLike)
    print 'median blue likelihood: ',median(blueLike)
    print
    print 'avg red likelihood: ',mean(redLike)
    print 'median red likelihood: ',median(redLike)
    
    print
    print
    print '-------------------- Distribution analysis ----------------------'
    print
    print
    
    print ' FANCY INCLINATIONS: '
    
    # perform the K-S and AD tests for inclination
    ans1 = stats.ks_2samp(blueFancyInc, redFancyInc)
    ans1a = stats.anderson_ksamp([blueFancyInc,redFancyInc])

    print 'KS for blue vs red fancy inclinations: ',ans1
    print 'AD for blue vs red fancy inclinations: ',ans1a
    
    ans2 = stats.ks_2samp(blueFancyInc, allAdjustedIncs)
    print 'KS for blue vs all fancy inclinations: ',ans2
    
    ans3 = stats.ks_2samp(redFancyInc, allAdjustedIncs)
    print 'KS for red vs all fancy inclinations: ',ans3
    
    print
    z_statrb, p_valrb = stats.ranksums(blueFancyInc, redFancyInc)
    z_statall, p_valall = stats.ranksums(assocFancyInc, allAdjustedIncs)
    print 'ranksum red vs blue p-value: ',p_valrb
    print 'ranksum associated vs all: ',p_valall


    ans4 = stats.ks_2samp(assocFancyInc, allAdjustedIncs)
    ans4a = stats.anderson_ksamp([assocFancyInc,allAdjustedIncs])

    print 'KS for all associated vs all fancy inclinations: ',ans4
    print 'AD for all associated vs all fancy inclinations: ',ans4a
    
    print

#     ans5 = stats.ks_2samp(spiralIncList, allSpiralIncList)
#     ans5a = stats.anderson_ksamp([spiralIncList,allSpiralIncList])
# 
#     print 'KS for all spiral associated vs all spiral fancy inclinations: ',ans5
#     print 'AD for all spiral associated vs all spiral fancy inclinations: ',ans5a
    
    print
    print ' INCLINATIONS: '
    print
    
    # perform the K-S and AD tests for inclination
    ans1 = stats.ks_2samp(blueInc, redInc)
    ans1a = stats.anderson_ksamp([blueInc,redInc])

    print 'KS for blue vs red inclinations: ',ans1
    print 'AD for blue vs red inclinations: ',ans1a
    
    ans2 = stats.ks_2samp(blueInc, allInclinations)
    print 'KS for blue vs all inclinations: ',ans2
    
    ans3 = stats.ks_2samp(redInc, allInclinations)
    print 'KS for red vs all inclinations: ',ans3
    
    assocInc = incs
    ans4 = stats.ks_2samp(assocInc, allInclinations)
    print 'KS for associated vs all inclinations: ',ans4
    
    print
    print ' EW Distributions: '
    print
    
    # perform the K-S and AD tests for EW
    ans1 = stats.ks_2samp(blueW, redW)
    ans1a = stats.anderson_ksamp([blueW,redW])
    print 'KS for blue vs red EW: ',ans1
    print 'AD for blue vs red EW: ',ans1a
    

    print
    print ' Impact parameter Distributions: '
    print
    
    # perform the K-S and AD tests for impact parameter
    ans1 = stats.ks_2samp(blueImpact, redImpact)
    ans1a = stats.anderson_ksamp([blueImpact,redImpact])
    print 'KS for blue vs red impact parameters: ',ans1
    print 'AD for blue vs red impact parameters: ',ans1a
    
    print
    print ' \Delta v Distributions: '
    print
    
    # perform the K-S and AD tests for \delta v
    ans1 = stats.ks_2samp(blueAbs, redAbs)
    ans1a = stats.anderson_ksamp([blueAbs,redAbs])
    print 'KS for blue vs red \Delta v: ',ans1
    print 'AD for blue vs red \Delta v: ',ans1a
    
    print
    print ' Azimuth Distributions: '
    print
    
    # perform the K-S and AD tests for azimuth
    ans1 = stats.ks_2samp(blueAz, redAz)
    ans1a = stats.anderson_ksamp([blueAz,redAz])
    print 'KS for blue vs red azimuth: ',ans1
    print 'AD for blue vs red azimuth: ',ans1a
    print
    
    # now against a flat distribution
    flatRed = arange(0,90,1)
    flatBlue = arange(0,90,1)

    ans1 = stats.ks_2samp(blueAz, flatBlue)
    ans1a = stats.anderson_ksamp([blueAz,flatBlue])
    print 'KS for blue vs flat azimuth: ',ans1
    print 'AD for blue vs flat azimuth: ',ans1a
    print
    ans1 = stats.ks_2samp(redAz, flatRed)
    ans1a = stats.anderson_ksamp([redAz,flatRed])
    print 'KS for red vs flat azimuth: ',ans1
    print 'AD for erd vs flat azimuth: ',ans1a
    print
    
    
    print
    print ' R_vir Distributions: '
    print
    
    # perform the K-S and AD tests for r_vir
    ans1 = stats.ks_2samp(blueVir, redVir)
    ans1a = stats.anderson_ksamp([blueVir,redVir])
    print 'KS for blue vs red R_vir: ',ans1
    print 'AD for blue vs red R_vir: ',ans1a
    
    print
    print ' Doppler parameter Distributions: '
    print
    
    # perform the K-S and AD tests for doppler parameter
    ans1 = stats.ks_2samp(blueB, redB)
    ans1a = stats.anderson_ksamp([blueB,redB])
    print 'KS for blue vs red doppler parameter: ',ans1
    print 'AD for blue vs red doppler parameter: ',ans1a
    
    print
    print ' Likelihood Distributions: '
    print
    
    # perform the K-S and AD tests for doppler parameter
    ans1 = stats.ks_2samp(blueLike, redLike)
    ans1a = stats.anderson_ksamp([blueLike,redLike])
    print 'KS for blue vs red likelihood: ',ans1
    print 'AD for blue vs red likelihood: ',ans1a
    
    print
    print ' COMPLETED. '

    

###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    