#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotW_impact_nonassociated.py, v 1 01/05/18

Plot W vs impact to nearest galaxy




Based on plotW_impact_final.py:
Plot EW as a function of impact parameter, and impact parameter/diameter and /R_vir
    (01/04/2016)


This is the plotW_b_diam bit from histograms3.py. Now is separated, and loads in a pickle
file of the relevant data, as created by "buildDataLists.py"

Previous (from histograms3.py):
    Plot some stuff for the 100largest initial results

    Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

    Updated for the pilot paper (05/06/15)


v5.1: updated for LG_correlation_combined5_8_edit2.csv for l_min = 0.001 (02/24/2016)

v5.2: remake plots with v_hel instead of vcorr (4/21/16)

v5.3: remake plots with new large galaxy sample (7/13/16) -> /plots4/

v5.4: add the ability to limit results based on 'environment' number (7/14/16)
        also add a likelihood limit 


v5.5: major edits to structure and functions included. Same ideas, but better formatting
    and removed some duplicate functions. Made plots4/ for new pilot paper (8/05/16)
    
v5.6: update with LG_correlation_combined5_11_25cut_edit4.csv and /plots5/
    (9/26/16)

v5.7: final version used after first referee report (12/09/16)
    - make marker points larger (60), make markers for imp/R_vir <=1 systems open

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
from scipy import stats


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

fontScale = 18
rc('text', usetex=True)
rc('font', size=18, family='serif', weight='normal')
rc('xtick.major',size=8,width=0.6)
rc('xtick.minor',size=5,width=0.6)
rc('ytick.major',size=8,width=0.6)
rc('ytick.minor',size=5,width=0.6)
rc('xtick',labelsize = fontScale)
rc('ytick',labelsize = fontScale)
rc('axes',labelsize = fontScale)
rc('xtick', labelsize = fontScale)
rc('ytick',labelsize = fontScale)
# rc('font', weight = 450)
# rc('axes',labelweight = 'bold')
rc('axes',linewidth = 1)
    

###########################################################################

def perc90(a):
    if len(a)>0:
        return percentile(a,90)
    else:
        return 0
        
def perc10(a):
    if len(a)>0:
        return percentile(a,10)
    else:
        return 0
        
def perc70(a):
    if len(a)>0:
        return percentile(a,70)
    else:
        return 0

def errors(a):
    # return the standard error in the mean for the input array
    return stats.sem(a)

    
def main():
    # assuming 'theFile' contains one name per line, read the file
    
    if getpass.getuser() == 'frenchd':
#         pickleFilename = '/Users/frenchd/Research/inclination/git_inclination/pilot_paper_code/pilotData2.p'
#         resultsFilename = '/Users/frenchd/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit4.csv'
#         saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/pilot_paper_code/plots6/'
#         WS09data = '/Users/frenchd/Research/inclination/git_inclination/WS2009_lya_data.tsv'
        
#         pickleFilename = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/pickleSALT.p'
#         saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/figures/'

#         pickleFilename = '/Users/frenchd/Research/inclination/git_inclination/picklePilot_plusSALTcut.p'
        pickleFilename = '/Users/frenchd/Research/inclination/git_inclination/picklePilot_plusSALT_14.p'
        gtPickleFilename = '/Users/frenchd/Research/inclination/git_inclination/pickleGT.p'
        saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/plotting_code/figs/'


    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # use the old pickle file to get the full galaxy dataset info
    pickleFile = open(pickleFilename,'rU')
    fullDict = pickle.load(pickleFile)
    pickleFile.close()
    
    # for the whole galaxy table:
    gtPickleFile = open(gtPickleFilename,'rU')
    gtDict = pickle.load(gtPickleFile)
    gtPickleFile.close()
    
    
    # save each plot?
    save = False
    
#     results = open(resultsFilename,'rU')
#     reader = csv.DictReader(results)
    
#     WS = open(WS09data,'rU')
#     WSreader = csv.DictReader(WS,delimiter=';')
    
    virInclude = False
    cusInclude = False
    finalInclude = 1
    
    maxEnv = 100
    minL = 0.001
    
    # if match, then the includes in the file have to MATCH the includes above. e.g., if 
    # virInclude = False, cusInclude = True, finalInclude = False, then only systems
    # matching those three would be included. Otherwise, all cusInclude = True would be included
    # regardless of the others
    match = False
    
    # all the lists to be used for associated lines
    raList = []
    decList = []
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
    
    AGNnameList = []
    nameList = []
    
    # for ambiguous lines (include = 0)
    lyaVAmbList = []
    lyaWAmbList = []
    envAmbList = []
    ambAGNnameList = []
    
    # for include = 2 lines
    lyaV_2List = []
    lyaW_2List = []
    env_2List = []
    vir_2List = []
    impact_2List = []
    like_2List = []
    
    # for include = 3 lines
    lyaV_3List = []
    lyaW_3List = []
    env_3List = []
    vir_3List = []
    impact_3List = []
    like_3List = []
    
    
    
    # for all lines with a galaxy within 500 kpc
    lyaV_nearestList = []
    lyaW_nearestList = []
    env_nearestList = []
    impact_nearestList = []
    diam_nearestList = []
    vir_nearestList = []
    cus_nearestList = []
    
    # WS lists
#     WSvcorr = []
#     WSdiam = []
#     WSimpact =[]
#     WSew = []
#     WSvel = []
#     WSlya = []
#     WSvel_dif = []
#     WSvir = []
#     WSlike = []
#     
#     l_min = 0.001
# 
#     for w in WSreader:
#         vcorr = w['HV']
#         diam = w['Diam']
#         rho = w['rho']
#         ew = w['EWLya']
#         vel = w['LyaVel']
#         lya = w['Lya']
#         
#         if lya == 'Lya  ' and isNumber(diam) and isNumber(ew) and isNumber(rho):
#             if float(rho) <=500.0:
#                 # this is a single galaxy association
#                 vir = calculateVirialRadius(float(diam))
#                 
#                 vel_dif = float(vcorr) - float(vel)
#     
#                 # try this "sphere of influence" value instead
#                 m15 = float(diam)**1.5
# 
#                 # first for the virial radius
#                 likelihood = math.exp(-(float(rho)/vir)**2) * math.exp(-(vel_dif/200.)**2)
#                 
#                 if vir>= float(rho):
#                     likelihood = likelihood*2
#                     
#                 # then for the second 'virial like' m15 radius
#                 likelihoodm15 = math.exp(-(float(rho)/m15)**2) * math.exp(-(vel_dif/200.)**2)
#                 
#                 if m15>= float(rho):
#                     likelihoodm15 = likelihoodm15*2
#                     
#                 if likelihood <= likelihoodm15:
#                     likelihood = likelihoodm15
#                     
#                 WSlike.append(likelihood)
#                 
# #                 l_min=0
#                 
#                 if likelihood >= l_min:
#                 
#                     WSvcorr.append(float(vcorr))
#                     WSdiam.append(float(diam))
#                     WSvir.append(vir)
#                     WSimpact.append(float(rho))
#                     WSew.append(float(ew))
#                     WSvel.append(float(vel))
#                     WSlya.append(lya)
#                     WSvel_dif.append(vel_dif)
    
    
    
    targetNameL= fullDict['targetName']
    galaxyNameL = fullDict['galaxyName']
    environmentL = fullDict['environment']
    RA_agnL = fullDict['RA_agn']
    Dec_agnL = fullDict['Dec_agn']
    RA_galL = fullDict['RA_gal']
    Dec_galL = fullDict['Dec_gal']
    likelihoodL = fullDict['likelihood']
    likelihood_cusL = fullDict['likelihood_cus']
    virialRadiusL = fullDict['virialRadius']
    cusL = fullDict['cus']
    impactParameterL = fullDict['impact']
    vcorrL = fullDict['vcorr']
    radialVelocityL = fullDict['radialVelocity']
    vel_diffL = fullDict['vel_diff']
    distGalaxyL = fullDict['distGalaxy']
    majorAxisL = fullDict['majorAxis']
    minorAxisL = fullDict['minorAxis']
    inclinationL = fullDict['inclination']
    positionAngleL = fullDict['PA']
    azimuthL = fullDict['azimuth']
    RC3flagL = fullDict['RC3flag']
    RC3typeL = fullDict['RC3type']
    RC3incL = fullDict['RC3inc']
    RC3paL = fullDict['RC3pa']
    final_morphologyL = fullDict['final_morphology']
    includeL = fullDict['include']
    include_virL = fullDict['include_vir']
    include_customL = fullDict['include_custom']
    Lya_vL = fullDict['Lya_v']
    vlimitsL = fullDict['vlimits']
    Lya_WL = fullDict['Lya_W']
    NaL = fullDict['Na']
    bL = fullDict['b']
    identifiedL = fullDict['identified']
    sourceL = fullDict['source']
    
    print 'initial len(Lya_vL): ',len(Lya_vL)
#     print
#     print 'type(includeL): ',type(includeL)
#     print 'type(includeL[0]): ',type(includeL[0])
    includeL = [int(i) for i in includeL]

    i = -1
    for include,include_vir,include_cus in zip(includeL,include_virL,include_customL):
        i+=1
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

        galaxyName = galaxyNameL[i]
        targetName = targetNameL[i]
        RA_agn = RA_agnL[i]
        Dec_agn = Dec_agnL[i]
        RA_gal = RA_galL[i]
        Dec_gal = Dec_galL[i]
        lyaV = Lya_vL[i]
        lyaW = Lya_WL[i]
        lyaW_err = lyaW*0.1
        env = environmentL[i]
        impact = impactParameterL[i]
        galaxyDist = distGalaxyL[i]
        pa = positionAngleL[i]
        RC3pa = RC3paL[i]
        morph = final_morphologyL[i]
        vcorr = vcorrL[i]
        maj = majorAxisL[i]
        minor = minorAxisL[i]
        inc = inclinationL[i]
        az = azimuthL[i]
        b = bL[i]
        b_err = b*0.1
        na = NaL[i]
        na_err = na*0.1
        likelihood = likelihoodL[i]
        likelihoodm15 = likelihood_cusL[i]
        virialRadius = virialRadiusL[i]
        m15 = cusL[i]
        vel_diff = vel_diffL[i]
        source = sourceL[i]
        
        AGNnameList.append(targetName)
        
        # for ambiguous lines
        if include == 0:
            lyaVAmbList.append(float(lyaV))
            lyaWAmbList.append(float(lyaW))
            envAmbList.append(float(env))
            ambAGNnameList.append(targetName)
            
        print 'include = ', include
        if include == 2:
            print 'include2 = ',include
            # for include = 2 lines
            lyaV_2List.append(float(lyaV))
            lyaW_2List.append(float(lyaW))
            env_2List.append(float(env))
            vir_2List.append(float(virialRadius))
            impact_2List.append(float(impact))
            like_2List.append(float(likelihood))
    
        if include == 3:
            # for include = 3 lines
            lyaV_3List.append(float(lyaV))
            lyaW_3List.append(float(lyaW))
            env_3List.append(float(env))
            vir_3List.append(float(virialRadius))
            impact_3List.append(float(impact))
            like_3List.append(float(likelihood))
            
        # for all absorbers with a galaxy within 500kpc
        if isNumber(impact):
            lyaV_nearestList.append(float(lyaV))
            lyaW_nearestList.append(float(lyaW))
            env_nearestList.append(float(env))
            impact_nearestList.append(float(impact))
            diam_nearestList.append(float(maj))
            nameList.append(galaxyName)
            vir_nearestList.append(float(virialRadius))
            cus_nearestList.append(float(m15))

            
#         if go and source == 'salt':
#         if go and source == 'pilot':
        if go and env <=maxEnv:
#         if go:
            if isNumber(RC3pa) and not isNumber(pa):
                pa = RC3pa
            
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
            
            # all the lists to be used for associated lines
            if float(env) <= maxEnv and float(likelihood) >= minL:
                raList.append(RA_gal)
                decList.append(Dec_gal)
                lyaVList.append(float(lyaV))
                lyaWList.append(float(lyaW))
                lyaErrList.append(float(lyaW_err))
                naList.append(na)
                bList.append(float(b))
                impactList.append(float(impact))
                azList.append(float(az))
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
                nameList.append(galaxyName)


    # lists for the full galaxy dataset
    majorAxisL = gtDict['majorAxis']
    incL = gtDict['inc']
    adjustedIncL = gtDict['adjustedInc']
    paL = gtDict['PA']
    BmagL = gtDict['Bmag']
    Bmag_sdssL = gtDict['Bmag_sdss']
    RID_medianL = gtDict['RID_median']
    RID_meanL = gtDict['RID_mean']
    RID_stdL = gtDict['RID_std']
    VhelL = gtDict['Vhel']
    RAdegL = gtDict['RAdeg']
    DEdegL = gtDict['DEdeg']
    NameL= gtDict['Name']
    
    allPA = paL
    allInclinations = []
    allCosInclinations = []

#     print 'type: ',type(incL)
    for i in incL:
        if i != -99:
            i = float(i)
            allInclinations.append(i)
            
            i2 = pi/180. * i
            cosi2 = cos(i)
            allCosInclinations.append(cosi2)
            
    allFancyInclinations = []
    allCosFancyCosInclinations = []
    for i in adjustedIncL:
        if i != -99:
            i = float(i)

            allFancyInclinations.append(i)
            
            i2 = pi/180. * i
            cosi2 = cos(i)
            allCosFancyCosInclinations.append(cosi2)
            
    allDiameter = majorAxisL

    print 'finished with this shit'
    print 'len(lyaV_2List): ',len(lyaV_2List)

    total = 0
    totalNo = 0
    totalYes = 0
    totalIsolated = 0
    totalGroup = 0




##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter, split between
    # absorbers at impact < 1 R_vir, overplot median histograms for each
    #
    # FOR ANY ABSORBER WITH A GALAXY WITHIN 500 KPC
    
    plotW_impact_virseparate_med = True
    save = True
    
    if plotW_impact_virseparate_med:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        alpha = 0.7
        alphaInside = 0.7
        binSize = 75
        bins = arange(0,625,binSize)
        markerSize = 60
        color = 'blue'
        symbol = 'o'
        maxW = 1500.


        labelr = r'$\rm Redshifted ~Absorber$'
        labelb = r'$\rm Blueshifted ~Absorber$'
#         bSymbol = 'D'
        
        xVals = []
        yVals = []
        redX = []
        redY = []
        blueX = []
        blueY = []
        
        for i,w,v in zip(impact_nearestList,lyaW_nearestList,vir_nearestList):
            # check if all the values are good
            if isNumber(i) and isNumber(w) and isNumber(v):
                if i!=-99 and w!=-99 and v!=-99 and w <=maxW:
                    xVal = float(i)
                    yVal = float(w)
                    
                    xVals.append(xVal)
                    yVals.append(yVal)
                        
                    if float(i) > float(v):
                        # impact parameter > virial radius
                         a = alpha
                         fc = color
                         ec = 'black'
                    
                    if float(i) <= float(v):
                        # impact parameter <= virial radius
                        a = alphaInside
                        fc = 'none'
                        ec = color

                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=markerSize,\
                    facecolor=fc,edgecolor=ec,alpha=a)
     
    
        # histograms
        bin_means,edges,binNumber = stats.binned_statistic(array(xVals), array(yVals), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dashed',color='black',lw=1.7,alpha=alpha+0.1,label=r'$\rm Mean ~EW$')
        
        
        # x-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        xlabel(r'$\rm \rho ~[kpc]$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ax.legend(scatterpoints=1,prop={'size':14},loc=1,fancybox=True)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0,500)

        if save:
            savefig('{0}/W(impact)_mean_{1}_virsep_nonassociated.pdf'.format(saveDirectory,binSize,maxEnv),format='pdf',bbox_inches='tight')
        else:
            show()
            
            
            
##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter/R_vir, split between
    # red and blue shifted absorption, overplot median histograms for red and blue
    #
    # plot shaded error regions around histogram
    # 
    # FOR ANY ABSORBER WITH A GALAXY WITHIN 500 KPC

    
    plotW_impact_vir_med = True
    save = True
    
    if plotW_impact_vir_med:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1

        alpha = 0.7
        alphaInside = 0.7
        markerSize = 60
        
        binSize = 0.25
        bins = arange(0,2.5,binSize)
        
        labelr = r'$\rm Redshifted ~Absorber$'
        labelb = r'$\rm Blueshifted ~Absorber$'
        color = 'blue'
        symbol = 'o'
        maxW = 1500.


                    
        xVals = []
        yVals = []
        redX = []
        redY = []
        blueX = []
        blueY = []
        
        for i,w,v in zip(impact_nearestList,lyaW_nearestList,vir_nearestList):
            # check if all the values are good
            if isNumber(i) and isNumber(w) and isNumber(v):
                if i!=-99 and w!=-99 and v!=-99 and w <=maxW:
                    xVal = float(i)/float(v)
                    yVal = float(w)
                    
                    xVals.append(xVal)
                    yVals.append(yVal)
                        
                    if float(i) > float(v):
                        # impact parameter > virial radius
                         a = alpha
                         fc = color
                         ec = 'black'
                    
                    if float(i) <= float(v):
                        # impact parameter <= virial radius
                        a = alphaInside
                        fc = 'none'
                        ec = color

                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=markerSize,\
                    facecolor=fc,edgecolor=ec,alpha=a)
     
    
        # histograms
        bin_means,edges,binNumber = stats.binned_statistic(array(xVals), array(yVals), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dashed',color='black',lw=1.7,alpha=alpha+0.1,label=r'$\rm Mean ~EW$')
    
        
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        xlabel(r'$\rm \rho / R_{vir}$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        leg = ax.legend(scatterpoints=1,prop={'size':14},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0,2.0)

        if save:
            savefig('{0}/W(impact_vir)_mean_{1}_virsep_nonassociated.pdf'.format(saveDirectory,binSize,maxEnv),format='pdf',bbox_inches='tight')
        else:
            show()
            
            
#########################################################################################
#########################################################################################



##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter, split between
    # absorbers at impact < 1 R_vir, overplot median histograms for each
    #
    # FOR INCLUDE = 2 ABSORBERS

    
    plotW_impact_virseparate_include2 = True
    save = True
    
    if plotW_impact_virseparate_include2:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        alpha = 0.7
        alphaInside = 0.7
        binSize = 100
        bins = arange(0,625,binSize)
        markerSize = 60
        color = 'blue'
        symbol = 'o'
        maxW = 1500.


        labelr = r'$\rm Redshifted ~Absorber$'
        labelb = r'$\rm Blueshifted ~Absorber$'
#         bSymbol = 'D'
        
        xVals = []
        yVals = []
        redX = []
        redY = []
        blueX = []
        blueY = []
        
        for i,w,v in zip(impact_2List,lyaW_2List,vir_2List):
            # check if all the values are good
            if isNumber(i) and isNumber(w) and isNumber(v):
                if i!=-99 and w!=-99 and v!=-99 and w <=maxW:
                    xVal = float(i)
                    yVal = float(w)
                    
                    xVals.append(xVal)
                    yVals.append(yVal)
                        
                    if float(i) > float(v):
                        # impact parameter > virial radius
                         a = alpha
                         fc = color
                         ec = 'black'
                    
                    if float(i) <= float(v):
                        # impact parameter <= virial radius
                        a = alphaInside
                        fc = 'none'
                        ec = color

                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=markerSize,\
                    facecolor=fc,edgecolor=ec,alpha=a)
     
    
        # histograms
        bin_means,edges,binNumber = stats.binned_statistic(array(xVals), array(yVals), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dashed',color='black',lw=1.7,alpha=alpha+0.1,label=r'$\rm Mean ~EW$')
        
        
        # x-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        xlabel(r'$\rm \rho ~[kpc]$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ax.legend(scatterpoints=1,prop={'size':14},loc=1,fancybox=True)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0,500)

        if save:
            savefig('{0}/W(impact)_mean_{1}_virsep_include2.pdf'.format(saveDirectory,binSize,maxEnv),format='pdf',bbox_inches='tight')
        else:
            show()
            
            
            
##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter/R_vir, split between
    # red and blue shifted absorption, overplot median histograms for red and blue
    #
    # plot shaded error regions around histogram
    # 
    # FOR INCLUDE = 2 ABSORBERS


    plotW_impact_vir_median_include2 = True
    save = True
    
    if plotW_impact_vir_median_include2:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1

        alpha = 0.7
        alphaInside = 0.7
        markerSize = 60
        
        binSize = 0.5
        bins = arange(0,2.5,binSize)
        
        labelr = r'$\rm Redshifted ~Absorber$'
        labelb = r'$\rm Blueshifted ~Absorber$'
        color = 'blue'
        symbol = 'o'
        maxW = 1500.


                    
        xVals = []
        yVals = []
        redX = []
        redY = []
        blueX = []
        blueY = []
        
        for i,w,v in zip(impact_2List,lyaW_2List,vir_2List):
            # check if all the values are good
            if isNumber(i) and isNumber(w) and isNumber(v):
                if i!=-99 and w!=-99 and v!=-99 and w <=maxW:
                    xVal = float(i)/float(v)
                    yVal = float(w)
                    
                    xVals.append(xVal)
                    yVals.append(yVal)
                        
                    if float(i) > float(v):
                        # impact parameter > virial radius
                         a = alpha
                         fc = color
                         ec = 'black'
                    
                    if float(i) <= float(v):
                        # impact parameter <= virial radius
                        a = alphaInside
                        fc = 'none'
                        ec = color

                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=markerSize,\
                    facecolor=fc,edgecolor=ec,alpha=a)
     
    
        # histograms
        bin_means,edges,binNumber = stats.binned_statistic(array(xVals), array(yVals), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dashed',color='black',lw=1.7,alpha=alpha+0.1,label=r'$\rm Mean ~EW$')
    
        
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        xlabel(r'$\rm \rho / R_{vir}$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        leg = ax.legend(scatterpoints=1,prop={'size':14},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0,2.0)

        if save:
            savefig('{0}/W(impact_vir)_mean_{1}_virsep_include2.pdf'.format(saveDirectory,binSize,maxEnv),format='pdf',bbox_inches='tight')
        else:
            show()



#########################################################################################
#########################################################################################



##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter, split between
    # absorbers at impact < 1 R_vir, overplot median histograms for each
    #
    # FOR INCLUDE = 3 ABSORBERS

    
    plotW_impact_virseparate_include3 = True
    save = True
    
    if plotW_impact_virseparate_include3:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        alpha = 0.7
        alphaInside = 0.7
        binSize = 100
        bins = arange(0,625,binSize)
        markerSize = 60
        color = 'blue'
        symbol = 'o'
        maxW = 1500.


        labelr = r'$\rm Redshifted ~Absorber$'
        labelb = r'$\rm Blueshifted ~Absorber$'
#         bSymbol = 'D'
        
        xVals = []
        yVals = []
        redX = []
        redY = []
        blueX = []
        blueY = []
        
        for i,w,v in zip(impact_3List,lyaW_3List,vir_3List):
            # check if all the values are good
            if isNumber(i) and isNumber(w) and isNumber(v):
                if i!=-99 and w!=-99 and v!=-99 and w <=maxW:
                    xVal = float(i)
                    yVal = float(w)
                    
                    xVals.append(xVal)
                    yVals.append(yVal)
                        
                    if float(i) > float(v):
                        # impact parameter > virial radius
                         a = alpha
                         fc = color
                         ec = 'black'
                    
                    if float(i) <= float(v):
                        # impact parameter <= virial radius
                        a = alphaInside
                        fc = 'none'
                        ec = color

                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=markerSize,\
                    facecolor=fc,edgecolor=ec,alpha=a)
     
    
        # histograms
        bin_means,edges,binNumber = stats.binned_statistic(array(xVals), array(yVals), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dashed',color='black',lw=1.7,alpha=alpha+0.1,label=r'$\rm Mean ~EW$')
        
        
        # x-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        xlabel(r'$\rm \rho ~[kpc]$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ax.legend(scatterpoints=1,prop={'size':14},loc=1,fancybox=True)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0,500)

        if save:
            savefig('{0}/W(impact)_mean_{1}_virsep_include3.pdf'.format(saveDirectory,binSize,maxEnv),format='pdf',bbox_inches='tight')
        else:
            show()
            
            
            
##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter/R_vir, split between
    # red and blue shifted absorption, overplot median histograms for red and blue
    #
    # plot shaded error regions around histogram
    # 
    # FOR INCLUDE = 3 ABSORBERS


    plotW_impact_vir_median_include3 = True
    save = True
    
    if plotW_impact_vir_median_include3:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1

        alpha = 0.7
        alphaInside = 0.7
        markerSize = 60
        
        binSize = 0.5
        bins = arange(0,2.5,binSize)
        
        labelr = r'$\rm Redshifted ~Absorber$'
        labelb = r'$\rm Blueshifted ~Absorber$'
        color = 'blue'
        symbol = 'o'
        maxW = 1500.


                    
        xVals = []
        yVals = []
        redX = []
        redY = []
        blueX = []
        blueY = []
        
        for i,w,v in zip(impact_3List,lyaW_3List,vir_3List):
            # check if all the values are good
            if isNumber(i) and isNumber(w) and isNumber(v):
                if i!=-99 and w!=-99 and v!=-99 and w <=maxW:
                    xVal = float(i)/float(v)
                    yVal = float(w)
                    
                    xVals.append(xVal)
                    yVals.append(yVal)
                        
                    if float(i) > float(v):
                        # impact parameter > virial radius
                         a = alpha
                         fc = color
                         ec = 'black'
                    
                    if float(i) <= float(v):
                        # impact parameter <= virial radius
                        a = alphaInside
                        fc = 'none'
                        ec = color

                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=markerSize,\
                    facecolor=fc,edgecolor=ec,alpha=a)
     
    
        # histograms
        bin_means,edges,binNumber = stats.binned_statistic(array(xVals), array(yVals), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dashed',color='black',lw=1.7,alpha=alpha+0.1,label=r'$\rm Mean ~EW$')
    
        
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        xlabel(r'$\rm \rho / R_{vir}$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        leg = ax.legend(scatterpoints=1,prop={'size':14},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0,2.0)

        if save:
            savefig('{0}/W(impact_vir)_mean_{1}_virsep_include3.pdf'.format(saveDirectory,binSize,maxEnv),format='pdf',bbox_inches='tight')
        else:
            show()



#########################################################################################
#########################################################################################




#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


if __name__=="__main__":
    # do the work
    main()
    