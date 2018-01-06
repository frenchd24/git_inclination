#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotAzHist2.py, v 5.6 01/05/18

This is the plotAzHist bit from histograms3.py. Now is separated, and loads in a pickle
file of the relevant data, as created by "buildDataLists.py"

Plots a histogram of azimuth angles for red and blue shifted absorbers

Previous (from histograms3.py):
    Plot some stuff for the 100largest initial results

    Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

    Updated for the pilot paper (05/06/15)

v5: updated to work with the new, automatically updated LG_correlation_combined5.csv
    (12/04/15) more on (12/28/2015)
    
v5.1: updated to work with LG_correlation_combined5_8_edit2.csv and l_min = 0.001
    (02/22/2016)

v5.2: remake plots with v_hel instead of vcorr (4/21/16)

v5.3: remake plots with newest large galaxy sample (07/13/16) -> /plots4/

v5.4: include ability to limit results based on environment number (7/14/16)

v5.5: update for LG_correlation_combined5_11_25cut_edit4.csv -> /plots5/ (9/26/16)

v5.6: updated for AAS_2017 (01/05/18)
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

# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# ## for Palatino and other serif fonts use:
# #rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)

from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})

# plt.clf()
# plt.rcdefaults()
# fontScale = 18
# params = {'axes.labelsize': fontScale,
# 'axes.titlesize': fontScale,
# 'font.size': fontScale,
# 'xtick.labelsize': fontScale,
# 'ytick.labelsize': fontScale,
# 'font.weight': 450,
# 'axes.labelweight': 450,
# 'xtick.major.size': 4,
# 'xtick.major.width': 1.2,
# 'xtick.minor.size': 3,
# 'xtick.minor.width': 1.2,
# 'ytick.major.size': 4,
# 'ytick.major.width': 1.2,
# 'ytick.minor.size': 3,
# 'ytick.minor.width': 1.2
# }


fontScale = 24
rc('text', usetex=True)
rc('font', size=24, family='serif', weight='normal')
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
    
    maxEnv = 100
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
    # make a histogram of all the azimuth angles
    #
    
    plotAzHist = False
    save = False
    
    #font = {'family':'serif','size':16}
#     font = {'family':'serif','size':18, 'serif': ['computer modern roman']}
#     rc('font',**font)
#     rc('legend',**{'fontsize':18})
#     rcParams['text.latex.preamble']=[r'\usepackage{amsmath}']
    
    
    if plotAzHist:
        fig = figure(figsize=(10,8))
        ax = fig.add_subplot(111)
        alpha = 0.85
        
#         bins = [0,10,20,30,40,50,60,70,80,90]
#         bins = array([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90])
#         bins +=4
#         bins = 45

        bins = arange(0,100,10)
        print 'bins:' ,bins
        print 'stdev: ',std(azList)
        
        bot = 0
        top = 45
        count = 0
        for a in azList:
            if a <=top and a>=bot:
                count+=1
                
        print 'between {0} and {1} = {2}'.format(top,bot,count)
        print 'total number: ',len(azList)
        print 'ratio = ',float(count)/len(azList)
        
        print 'mean: ',mean(azList)
        print 'median: ',median(azList)
        
#         rc('text', usetex=True)
#         rc('font', family='serif',weight='bold')
        
        plot1 = hist(azList,bins=bins,histtype='bar',alpha=alpha)
#         ax.tick_params(axis='both', which='major', labelsize=20)
        
#         print 'azList: ',azList
#         hist(azList)
#         title('Distribution of azimuths')

#         xlabel(r'Azimuth (deg)',fontsize=20,weight='bold')
#         ylabel(r'Number',fontsize=20,weight='bold')
        xlabel(r'Azimuth (deg)')
        ylabel(r'Number')
        xlim(0,90)
        ylim(0,10)
        tight_layout()
        
        if save:
            savefig('{0}/hist(azimuth).pdf'.format(saveDirectory),format='pdf')
        else:
            show()
      
#########################################################################################
#########################################################################################
    # make histograms for red and blue shifted azimuth angles
    #
    
    plotAzHist_dif = False
    save = False
    
    if plotAzHist_dif:
    
        fig = figure(figsize=(10,6))
        subplots_adjust(hspace=0.200)

#         bins = [0,15,30,45,60,75,90]
#         bins = arange(0,90,10)
        binsize = 10
        bins = arange(0,100,binsize)
        blue = []
        red = []
        alpha = 0.80
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        for d,a,l,e in zip(difList,azList,lyaWList,lyaErrList):
            if a != -99:
                if d >=0:
                    # blue shifted absorber, but galaxy is REDSHIFTED
                    print 'd: ',d
                    blue.append(a)
                    blueLya.append(l)
                    blueLyaErr.append(e)
                else:
                    # red shifted absorber, but galaxy is BLUESHIFTED
                    red.append(a)
                    redLya.append(l)
                    redLyaErr.append(e)
        
        print 'stats-----'
        print 'max red: ',max(red)
        print 'min red: ',min(red)
        print 'max blue:' ,max(blue)
        print 'min blue: ',min(blue)
        
        ax = fig.add_subplot(211)        
        hist(red,bins=bins,histtype='bar',color='red',alpha = alpha,label='Redshifted absorbers')
#         title('Red shifted absorption: Galaxies')    
        ylabel("Number")
        ylim(0,7)
        xlim(0,90)
        legend()

        ax = fig.add_subplot(212)
        hist(blue,bins=bins,histtype='bar',color='Blue',alpha = alpha,label='Blueshifted absorbers')
#         title('Blue shifted absorption: Galaxies')
        ylabel('Number')
        xlabel("Azimuth (deg)")
        xlim(0,90)
        ylim(0,7)
        legend()

#         tight_layout()

        if save:
            savefig('{0}/hist(azimuth)_dif_{1}.pdf'.format(saveDirectory,binsize),format='pdf')
        else:
            show()

#########################################################################################
#########################################################################################
    # make histograms for red and blue shifted azimuth angles overlaid on each other
    #
    
    plotAzHist_over_dif = False
    save = False
    
    if plotAzHist_over_dif:
    
        fig = figure(figsize=(10,6))
        subplots_adjust(hspace=0.200)
        
        alpha = 0.75

        bins = arange(0,100,10)
        blue = []
        red = []
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        for d,a,l,e in zip(difList,azList,lyaWList,lyaErrList):
            if a != -99:
                if d >=0:
                    # blue shifted absorber, but galaxy is REDSHIFTED
                    print 'd: ',d
                    blue.append(a)
                    blueLya.append(l)
                    blueLyaErr.append(e)
                else:
                    # red shifted absorber, but galaxy is BLUESHIFTED
                    red.append(a)
                    redLya.append(l)
                    redLyaErr.append(e)
        
        ax = fig.add_subplot(111)        
        hist(red,bins=bins,histtype='bar',color='red',alpha = alpha,hatch = '/',label='Redshifted absorbers')
        hist(blue,bins=bins,histtype='bar',color='Blue',alpha = alpha,label='Blueshifted absorbers')

        ylabel('Number')
        xlabel("Azimuth (deg)")
        xlim(0,90)
        ylim(0,7)
        legend()

#         tight_layout()

        if save:
            savefig('{0}/hist(azimuth)_overlaid_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
            
            
#########################################################################################
#########################################################################################
    # make histograms for all azimuth angles, with red and blue shifted overlaid
    #
    
    plotAzHist_all_over = True
    save = True
    
    if plotAzHist_all_over:
    
        fig = figure(figsize=(10,5))
#         subplots_adjust(hspace=0.200)
        
        alpha = 0.6

        bins = arange(0,100,10)
        blue = []
        red = []
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        for d,a,l,e in zip(difList,azList,lyaWList,lyaErrList):
            if a != -99:
                if d >=0:
                    # blue shifted absorber, but galaxy is REDSHIFTED
                    print 'd: ',d
                    blue.append(a)
                    blueLya.append(l)
                    blueLyaErr.append(e)
                else:
                    # red shifted absorber, but galaxy is BLUESHIFTED
                    red.append(a)
                    redLya.append(l)
                    redLyaErr.append(e)
        
        
        ax = fig.add_subplot(111)
        
        # all first
        hist(blue,bins=bins,histtype='bar',color='Blue',alpha=alpha+0.25,lw=1.7,hatch='//',label=r'$\rm Blueshifted$')
        hist(red,bins=bins,histtype='bar',color='red',alpha=alpha,lw=1.7,label=r'$\rm Redshifted$')
        plot1 = hist(azList,bins=bins,histtype='step',lw=2.1,alpha=0.9,color='black',label=r'$\rm All$')


        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
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

        ylabel(r'$\rm Number$')
        xlabel(r'$\rm Azimuth ~ [deg]$')
        xlim(0,90)
        ylim(0,20)
        legend(scatterpoints=1,prop={'size':16},loc=2,fancybox=True)


#         tight_layout()

        if save:
            savefig('{0}/hist(azimuth)_overlaid_all.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()

#########################################################################################
#########################################################################################
    # make histograms for red and blue shifted azimuth angles vs flat distributions
    #
    
    plotAzHist_dif_flat = False
    save = False
    
    if plotAzHist_dif_flat:
    
        fig = figure(figsize=(10,6))
        subplots_adjust(hspace=0.200)

#         bins = [0,15,30,45,60,75,90]
#         bins = arange(0,90,10)
        bins = arange(0,90,12)
        blue = []
        red = []
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        for d,a,l,e in zip(difList,azList,lyaWList,lyaErrList):
            if a != -99:
                if d >=0:
                    # blue shifted absorber, but galaxy is REDSHIFTED
                    blue.append(a)
                    blueLya.append(l)
                    blueLyaErr.append(e)
                else:
                    # red shifted absorber, but galaxy is BLUESHIFTED
                    red.append(a)
                    redLya.append(l)
                    redLyaErr.append(e)
        
        
#         flatRedNum = len(red)/9.0
#         print 'flatRedNum: ',flatRedNum
#         flatRed = ones(9)*flatRedNum
#     
#         flatBlueNum = len(blue)/9.0
#         flatBlue = ones(9)*flatBlueNum

        flatRed = arange(0,90,5)
        flatBlue = arange(0,90,5)
        
        print 'flatRed: ',flatRed
        print
        print 'flatBlue: ',flatBlue
        print
        
        ax = fig.add_subplot(411)        
        hist(red,bins=bins,histtype='bar',color='red',alpha = 0.80,label='Redshifted absorbers')
#         title('Red shifted absorption: Galaxies')    
        ylabel("Number")
        ylim(0,6)
        xlim(0,90)
        legend()
        
        ax = fig.add_subplot(412)        
        hist(flatRed,bins=bins,histtype='bar',color='red',alpha = 0.80,label='Flat Redshifted absorbers')
#         title('Red shifted absorption: Galaxies')    
        ylabel("Number")
        ylim(0,6)
        xlim(0,90)
        legend()

        ax = fig.add_subplot(413)
        hist(blue,bins=bins,histtype='bar',color='Blue',alpha = 0.80,label='Blueshifted absorbers')
#         title('Blue shifted absorption: Galaxies')
        ylabel('Number')
        xlabel("Azimuth (deg)")
        xlim(0,90)
        ylim(0,6)
        legend()

        ax = fig.add_subplot(414)
        hist(flatBlue,bins=bins,histtype='bar',color='Blue',alpha = 0.80,label='Flat Blueshifted absorbers')
#         title('Blue shifted absorption: Galaxies')
        ylabel('Number')
        xlabel("Azimuth (deg)")
        xlim(0,90)
        ylim(0,6)
        legend()
        
#         tight_layout()

        if save:
            savefig('{0}/hist(azimuth)_dif_flat.pdf'.format(saveDirectory),format='pdf')
        else:
            show()

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


if __name__=="__main__":
    # do the work
    main()
    