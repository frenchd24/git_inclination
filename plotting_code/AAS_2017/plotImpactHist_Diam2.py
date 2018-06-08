#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotImpactHist_Diam2.py, v 5.6 7/18/16

Plots histograms of impact parameter, and various versions of normalized impact parameters
also plots impact parameter vs R_vir



This is the plotImpactHist_Diam bit from histograms3.py. Now is separated, and loads in a pickle
file of the relevant data, as created by "buildDataLists.py"

Previous (from histograms3.py):
    Plot some stuff for the 100largest initial results

    Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

    Updated for the pilot paper (05/06/15)

v5: updated to work with the new, automatically updated LG_correlation_combined5.csv
    (12/04/15) - original updates to the individual files

v5.1: included a second option, to normalize by virial radius or d^1.5 instead (12/29/15)

v5.2: updated for LG_correlation_combined5_8_edit2.csv with l_min = 0.001 (02/22/2016)

v5.3: included a function to plot JUST impact parameter alone (no normalizing)
    - (03/02/2016)
    
v5.4: remake plots with v_hel instead of vcorr (4/21/16)

v5.5: remake plots with new large galaxy sample (7/13/16) -> /plots4/

v5.6: add ability to limit results by environment or likelihood (7/14/16)

    - some edits. Changed the old 'plotImpactHist_Vir_dif' function to plot a single
    marginal histogram style plot with both red and blue impact/R_vir overplotted 
    (7/18/16)

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
# rc('text', usetex=True)
# rc('font',size=16,weight='bold')
    
fontScale = 15
rc('text', usetex=True)
rc('font', size=15, family='serif', weight='normal')
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
    
    maxEnv = 1000
    minL = 0.001
    maxLyaW = 1500
    
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
    WSvcorr = []
    WSdiam = []
    WSimpact =[]
    WSew = []
    WSvel = []
    WSlya = []
    WSvel_dif = []
    WSvir = []
    WSlike = []
    
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
            if float(env) <= maxEnv and float(likelihood) >= minL and float(lyaW) <= maxLyaW:
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
    # make a histogram of the distribution of impact parameters for associated galaxies 
    # normalized by major diameter
    #
    # I'm not seeing this as useful
    
    plotImpactHist_Diam = False
    save = False
    
    if plotImpactHist_Diam:
        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)

        bins = arange(0,20.0,1.0)
        
        impactArray = np.array([])
        majArray = np.array([])
        
        # make sure all the values are okay
        for i,m in zip(impactList,majList):
            if isNumber(i) and isNumber(m):
                if i !=-99 and m !=-99:
                    impactArray = append(impactArray,float(i))
                    majArray = append(majArray,float(m))
        
        normalizedImpactArray = impactArray/majArray
        
        # x-axis
        majorLocator   = MultipleLocator(1)
        majorFormatter = FormatStrFormatter(r'$\rm %s$')
        minorLocator   = MultipleLocator(0.5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        plot1 = hist(normalizedImpactArray,bins=bins,histtype='bar')
        
        xlabel(r'$\rm Impact ~Parameter / Diameter$')
        ylabel(r'$\rm Number$')
#         ax.tick_params(axis='x', labelsize=0)
#         ax.tick_params(axis='y',labelsize=8)

        if save:
            savefig('{0}/hist(Impact_Diameter).pdf'.format(saveDirectory),format='pdf')
        else:
            show()


##########################################################################################
##########################################################################################
    # make a histogram of the distribution of impact parameters for associated galaxies 
    # normalized by virial radius
    #
    
    plotImpactHist_Vir = True
    save = True
    
    if plotImpactHist_Vir:
        fig = figure(figsize=(10,4))
        ax = fig.add_subplot(111)
        alpha = 0.75
        
#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
#         bins = [0,5,10,15,20,25,30,35,40]
        bins = arange(0,2.2,0.2)
        
        impactArray = np.array([])
        virArray = np.array([])
        
        # make sure all the values are okay
        for i,v in zip(impactList,virList):
            if isNumber(i) and isNumber(v):
                if i !=-99 and v !=-99:
                    impactArray = append(impactArray,float(i))
                    virArray = append(virArray,float(v))
        
        normalizedImpactArray = impactArray/virArray
        
        print 'median: ',median(normalizedImpactArray)
        print 'mean: ',mean(normalizedImpactArray)
        
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %s$')
        minorLocator   = MultipleLocator(0.25)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        
        plot1 = hist(normalizedImpactArray,bins=bins,histtype='bar',alpha=alpha)
        
#         title('Distribution of impact parameters')
        xlabel(r'$\rm \rho / R_{vir}$')
        ylabel(r'$\rm Number$')
#         grid(True)
    
#         ax.tick_params(axis='x', labelsize=10)
#         ax.tick_params(axis='y',labelsize=10)
        tight_layout()
#         ylim(0,5)

        if save:
            savefig('{0}/hist(Impact_vir).pdf'.format(saveDirectory),format='pdf')
        else:
            show()


##########################################################################################
##########################################################################################
    # make a histogram of the distribution of impact parameters for associated galaxies 
    # normalized by virial radius, overplotting red and blue samples
    #
    # This is specifically to make a marginal histogram for the Berlin poster
    
    plotImpactHist_Vir_dif_marginal = False
    save = False
    
    if plotImpactHist_Vir_dif_marginal:
        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        subplots_adjust(hspace=0.200)

        binsize = 0.2
        bins = arange(0,2.0,binsize)
        alpha = 0.7
        
        reds = np.array([])
        blues = np.array([])
        valList = np.array([])
        
        # make sure all the values are okay
        for i,v,d in zip(impactList,virList,difList):
            if isNumber(i) and isNumber(v):
                if i !=-99 and v !=-99:
                    val = float(i)/float(v)
                    valList = append(valList,val)
                    
                    # if blueshifted
                    if d>=0:
                        blues = append(blues,val)
                    
                    # for redshifted
                    else:
                        reds = append(reds,val)
        
        

        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %s$')
        minorLocator   = MultipleLocator(0.25)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # don't want any axis here though!
        ax.xaxis.set_major_formatter(plt.NullFormatter())
        
        # y axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)


        plot1 = hist(valList,bins=bins,histtype='step',lw=2.0,alpha=0.9,color='black',label=r'$\rm All$')

        hist(blues,bins=bins,histtype='bar',color='Blue',alpha = alpha,lw=1.5,ls='dashed',label=r'$\rm Blueshifted$')
        hist(reds,bins=bins,histtype='bar',color='red',alpha = alpha,lw=1.5,label=r'$\rm Redshifted$')

#         xlabel(r'$\rm \rho / R_{vir}$')
        ylabel(r'$\rm Number$')
        tight_layout()
        xlim(0,2.0)
        
#         legend(scatterpoints=1,prop={'size':14},fancybox=True,loc=2)
        

        if save:
            savefig('{0}/hist(Impact_vir_dif)_bin_{1}.pdf'.format(saveDirectory,binsize),format='pdf')
        else:
            show()


##########################################################################################
##########################################################################################
    # make a histogram of the distribution of impact parameters for associated galaxies 
    # normalized by (major diameter)^1.5
    #
    
    plotImpactHist_m15 = False
    save = False
    
    if plotImpactHist_m15:
        fig = figure(figsize=(10,2))
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(111)
        
#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
#         bins = [0,10,20,30,40,50,60,70,80,90]
        bins = arange(0,35,0.5)

        impactArray = np.array([])
        m15Array = np.array([])
        
        # make sure all the values are okay
        for i,m in zip(impactList,m15List):
            if isNumber(i) and isNumber(m):
                if i !=-99 and m !=-99:
                    impactArray = append(impactArray,float(i))
                    m15Array = append(m15Array,float(m))
        
        normalizedImpactArray = impactArray/m15Array
        
        plot1 = hist(normalizedImpactArray,bins=bins,histtype='bar')
        
        title('Distribution of impact parameters')
        xlabel('Impact Parameter/m15')
        ylabel('Number')
        ax.tick_params(axis='x', labelsize=0)
        ax.tick_params(axis='y',labelsize=8)
#         ylim(0,5)

        if save:
            savefig('{0}/hist(Impact_m15).pdf'.format(saveDirectory),format='pdf')
        else:
            show()



##########################################################################################
##########################################################################################
    # make a histogram of the distribution of impact parameters for associated galaxies, 
    # one each for red and blue samples
    #
    
    plotImpactHist_dif = False
    save = False
    
    if plotImpactHist_dif:
        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        subplots_adjust(hspace=0.200)
        alpha = 0.6

        binsize = 50
        bins = arange(0,500+binsize,binsize)
        
        reds = np.array([])
        blues = np.array([])
        
        # make sure all the values are okay
        for i,v,d in zip(impactList,virList,difList):
            if isNumber(i) and isNumber(v):
                if i !=-99 and v !=-99:
                    val = float(i)
                    
                    # if blueshifted
                    if d>=0:
                        blues = append(blues,val)
                    
                    # for redshifted
                    else:
                        reds = append(reds,val)
        
#         title('Distribution of impact parameters')

        plot1 = hist(blues,bins=bins,histtype='bar',color='blue',label='Blueshifted',alpha=alpha)
#         xlabel('Impact Parameter/R_vir (blueshifted)')
        ylabel('Number')
        xlabel(r'$\rho$ (kpc)')
        legend(scatterpoints=1,prop={'size':12})
        ax.tick_params(axis='y',labelsize=11)
        ax.tick_params(axis='x', labelsize=0)
        ylim(0,10)
#         tight_layout()
        
        if save:
            savefig('{0}/hist(Impact_blue)_bin_{1}.pdf'.format(saveDirectory,binsize),format='pdf')
        else:
            show()

        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        plot2 = hist(reds,bins=bins,histtype='bar',color="red",label='Redshifted',alpha=alpha)
        xlabel(r'$\rho$ (kpc)')
        ylabel('Number')
        legend(scatterpoints=1,prop={'size':12})
        
        ax.tick_params(axis='x', labelsize=11)
        ax.tick_params(axis='y',labelsize=11)
        ylim(0,10)
#         tight_layout()
        
        if save:
            savefig('{0}/hist(Impact_red)_bin_{1}.pdf'.format(saveDirectory,binsize),format='pdf')
        else:
            show()



##########################################################################################
##########################################################################################
    # make a histogram of the distribution of impact parameters for associated galaxies 
    # normalized by virial radius, overplotting red and blue samples
    #
    
    plotImpactHist_Vir_dif = False
    save = False
    
    if plotImpactHist_Vir_dif:
        fig = figure(figsize=(10,4))
        ax = fig.add_subplot(111)
        subplots_adjust(hspace=0.200)

        binSize = 50
        bins = arange(0,550,binSize)
        alpha = 0.6
        
        reds = np.array([])
        blues = np.array([])
        valList = np.array([])
        
        # make sure all the values are okay
        for i,v,d in zip(impactList,virList,difList):
            if isNumber(i) and isNumber(v):
                if i !=-99 and v !=-99:
                    val = float(i)
                    valList = append(valList,val)
                    
                    # if blueshifted
                    if d>=0:
                        blues = append(blues,val)
                    
                    # for redshifted
                    else:
                        reds = append(reds,val)
        

        # x-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)


        plot1 = hist(valList,bins=bins,histtype='step',lw=2.0,alpha=0.9,color='black',label=r'$\rm All$')

        hist(blues,bins=bins,histtype='bar',color='Blue',alpha = alpha,lw=1.5,ls='dashed',label=r'$\rm Blueshifted$')
        hist(reds,bins=bins,histtype='bar',color='red',alpha = alpha,lw=1.5,label=r'$\rm Redshifted$')

        xlabel(r'$\rm \rho ~[kpc]$')
        ylabel(r'$\rm Number$')
#         tight_layout()
#         xlim(0,2.0)
        
#         legend(scatterpoints=1,prop={'size':14},fancybox=True,loc=2)
        
        if save:
            savefig('{0}/hist(Impact_dif)_{1}.pdf'.format(saveDirectory,binSize),format='pdf')
        else:
            show()


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


if __name__=="__main__":
    # do the work
    main()
    