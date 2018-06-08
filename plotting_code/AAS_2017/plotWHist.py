#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotWHist.py, v 6.1 01/07/2018

Make histograms and CDFs of Lya equivalent widths. This is an evolution of the original
code -> plotLyaWHist_both2.py, v 5.2 04/21/2016

    - which is itself an evolution of the original histograms3.py, detailed below:


    This is the plotLyaWHist_both bit from histograms3.py. Now is separated, and loads in
    a pickle file of the relevant data, as created by "buildDataLists.py"

    Previous (from histograms3.py):
        Plot some stuff for the 100largest initial results

        Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

        Updated for the pilot paper (05/06/15)


    v5: updated to work with the new, automatically updated LG_correlation_combined5.csv
        (12/04/15) - original updates to the individual files
    
    v5.1: update for LG_correlation_combined5_8_edit2.csv and l_min = 0.001 (02/23/2016)

    v5.2: remake plots with v_hel instead of vcorr (4/21/16)

v6.0: Name change to better match the rest of the plotting codes. Include CDF plotting 
    functions as well. Output into /plots4 now
    
v6.1: update for AAS_2017 (01/07/18)

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

# fontScale = 16
# rc('text', usetex=True)
# rc('font',size=18)
# rc('xtick.major',size=5,width=1.2)
# rc('xtick.minor',size=3,width=1.2)
# rc('ytick.major',size=5,width=1.2)
# rc('ytick.minor',size=3,width=1.2)
# rc('xtick',labelsize=16)
# rc('ytick',labelsize=16)
# rc('axes',labelsize=16)
# rc('xtick', labelsize = fontScale)
# rc('ytick',labelsize = fontScale)
# # rc('font', weight = 450)
# # rc('axes',labelweight = 'bold')
# rc('axes',linewidth = 1)

fontScale = 15
rc('text', usetex=True)
# rc('font', size=15, family='serif', weight=450)
rc('font', size=15, family='serif',weight='normal')
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
    
    maxEnv = 1
    minL = 0.001
    maxLyaW = 10000
    
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

#########################################################################################
#########################################################################################
    # make a histogram of the distribution of Ly-alpha equivalent widths for both the 
    # associated and ambiguous samples
    #
    #
    # normByEnv doesn't work for the ambigous lines, because most of them have env = 0
    #
    
    plotLyaWHist_both = False
    save = False
    
    if plotLyaWHist_both:
#         fig = figure(figsize=(2,8))
        fig = figure()
        ax = fig.add_subplot(211)
#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
#         bins = arange(0,max(max(lyaWList),max(lyaWAmbList)),20)
        bins = 15
        
        lyaWArray = array(lyaWList)
        lyaWAmbArray = array(lyaWAmbList)
        
        envArray = array(envList)
        envAmbArray = array(envAmbList)
    
        print 'lyaWAmbArray: ',lyaWAmbArray
        print 'envAmbArray: ',envAmbArray
    
        # see above, does not work the for the ambigous ones
        normByEnv = False
        
#         envArray2 = []
#         for i in envArray:
#             if int(i) == 0:
#                 envArray2.append(1)
#             else:
#                 envArray2.append(int(i))
#         
#         envAmbArray2 = []
#         for i in envAmbArray:
#             if int(i) == 0:
#                 envAmbArray2.append(1)
#             else:
#                 envAmbArray2.append(int(i))
        
        if normByEnv:
            pass
#             plot1 = hist(lyaWArray/envArray2,bins=bins,histtype='bar',orientation = 'vertical')
#             title('Distribution of Lya W - Associated')
# 
#             ax = fig.add_subplot(212)
#             plot1 = hist(lyaWAmbArray/envAmbArray2,bins=bins,histtype='bar',orientation = 'vertical')

        else:
            plot1 = hist(lyaWList,bins=bins,histtype='bar',orientation = 'vertical',label='Associated')
            legend(scatterpoints=1,prop={'size':12},loc=1)
            ylabel('Number')
            xlim(0,1000)

            ax = fig.add_subplot(212)
            plot1 = hist(lyaWAmbList,bins=bins,histtype='bar',orientation = 'vertical',label='Ambiguous')
            legend(scatterpoints=1,prop={'size':12},loc=1)
            xlim(0,1000)


        
        xlabel(r'Equivalent Width ($\rm m\AA$)')
        ylabel('Number')
        ax.tick_params(axis='x', labelsize=11)
        ax.tick_params(axis='y',labelsize=11)
#         tight_layout()
        
        if save:
            savefig('{0}/hist(lyaW_assoc_vs_ambig)_2.pdf'.format(saveDirectory),format='pdf')
        else:
            show()


#########################################################################################
#########################################################################################
    # make a histogram of the distribution of Ly-alpha equivalent widths for the 
    # associated sample, splitting on red vs blue shifted absorbers
    #
    #
    
    plotWHist_dif = False
    save = False
    
    if plotWHist_dif:
#         fig = figure(figsize=(2,8))
        fig = figure()

        bins = arange(0,1300,100)
        
        lyaWArray = array(lyaWList)
        lyaWAmbArray = array(lyaWAmbList)
        
        envArray = array(envList)
        envAmbArray = array(envAmbList)
        
        blues = []
        reds = []
        for d,l in zip(difList,lyaWList):
            if float(d) >0:
                # blueshifted ABSORBER: dif = v_gal - v_absorber
                blues.append(float(l))
            else:
                reds.append(float(l))
    
        ax = fig.add_subplot(211)
        plot1 = hist(blues,bins=bins,histtype='bar',orientation = 'vertical',label=r'$\rm Blueshifted$',color="Blue",alpha = 0.85)
        legend(scatterpoints=1,prop={'size':12},loc=1)
        ylabel(r'$\rm Number$')
        ylim(0,10)
        
        # X-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)

        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # Y-axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        

        ax = fig.add_subplot(212)
        plot1 = hist(reds,bins=bins,histtype='bar',orientation = 'vertical',label=r'$\rm Redshifted$',color="Red",alpha = 0.85)
        legend(scatterpoints=1,prop={'size':12},loc=1)
        ylabel(r'$\rm Number$')
        ylim(0,10)
        

        # X-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)

        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # Y-axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        
        xlabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
#         ax.tick_params(axis='x', labelsize=11)
#         ax.tick_params(axis='y',labelsize=11)
#         xlim(0,1200)
        
        if save:
            savefig('{0}/hist(lyaW_blue_vs_red).pdf'.format(saveDirectory),format='pdf')
        else:
            show()
            

#########################################################################################
#########################################################################################
    # make a histogram of the distribution of Ly-alpha equivalent widths for the 
    # associated sample, overplotting red and blue shifted absorbers as well as the
    # combined set.
    #
    # This makes a marginal histogram for the Berlin poster
    
    plotWHist_dif_all = False
    save = False
    
    if plotWHist_dif_all:
        fig = figure(figsize=(10,2))
        
        alpha = 0.7
        
        
        bins = arange(0,1200,100)
        
        lyaWArray = array(lyaWList)
        
        blues = []
        reds = []
        
        for d,l in zip(difList,lyaWList):
            if float(d) >0:
                # blueshifted ABSORBER: dif = v_gal - v_absorber
                blues.append(float(l))
            else:
                reds.append(float(l))
    
        ax = fig.add_subplot(111)


        # x-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # don't want any axis here though!
        ax.xaxis.set_major_formatter(plt.NullFormatter())
        
        # y axis
        majorLocator   = MultipleLocator(4)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)


        plot1 = hist(lyaWList,bins=bins,histtype='step',lw=2.0,alpha=0.9,color='black',label=r'$\rm All$')

        hist(blues,bins=bins,histtype='bar',color='Blue',alpha = alpha,lw=1.5,ls='dashed',label=r'$\rm Blueshifted$')
        hist(reds,bins=bins,histtype='bar',color='red',alpha = alpha,lw=1.5,label=r'$\rm Redshifted$')

        ylabel(r'$\rm Number$')
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
        ax.yaxis.set_ticks_position('both')
        plt.gca().invert_xaxis()
        labels = ax.get_yticklabels()
        plt.setp(labels, rotation=90)

#         ylim(0,10)

#         xlabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')

#         xlim(0,1200)
        
        if save:
            savefig('{0}/hist(lyaW_dif_all).pdf'.format(saveDirectory),format='pdf')
        else:
            show()


#########################################################################################
#########################################################################################
    # make a CDF of the distribution of Ly-alpha equivalent widths for the 
    # associated sample, splitting on red vs blue shifted absorbers
    #
    #
    
    plotWCDF_dif = False
    save = False
    
    if plotWCDF_dif:
#         fig = figure(figsize=(2,8))
        fig = figure()
        ax = fig.add_subplot(111)

        bins = arange(0,1300,100)
        
        lyaWArray = array(lyaWList)
        lyaWAmbArray = array(lyaWAmbList)
        
        envArray = array(envList)
        envAmbArray = array(envAmbList)
        
        blues = []
        reds = []
        for d,l in zip(difList,lyaWList):
            if float(d) >0:
                # blueshifted ABSORBER: dif = v_gal - v_absorber
                blues.append(float(l))
            else:
                reds.append(float(l))

        sorted_blues = sort(blues)
        sorted_reds = sort(reds)

    
        # compute the CDF y-values
        yvals_blues = np.arange(len(sorted_blues)) / float(len(sorted_blues))
        yvals_reds = np.arange(len(sorted_reds)) / float(len(sorted_reds))

        
#         plot1 = hist(blues,bins=bins,histtype='bar',orientation = 'vertical',label='Blueshifted',color="Blue",alpha = 0.85)
#         legend(scatterpoints=1,prop={'size':12},loc=1)
#         ylabel(r'Number')
#         ylim(0,10)
# 
#         ax = fig.add_subplot(212)
#         plot1 = hist(reds,bins=bins,histtype='bar',orientation = 'vertical',label='Redshifted',color="Red",alpha = 0.85)
#         legend(scatterpoints=1,prop={'size':12},loc=1)
#         ylabel(r'Number')
#         ylim(0,10)
    
        # plot the y-values against the sorted data
        plt.plot(sorted_reds,yvals_reds,color='red',lw=3,label=r'$\rm Redshifted ~CDF$')
        plt.plot(sorted_blues,yvals_blues,color='blue',lw=3,label=r'$\rm Blueshifted ~CDF$')
        
        xlabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ylabel(r'$\rm Number$')
        
        
        # format all the axis and stuff
        # X-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)

        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # Y-axis
        majorLocator   = MultipleLocator(0.1)
        majorFormatter = FormatStrFormatter(r'$\rm %s$')
        minorLocator   = MultipleLocator(0.05)
        
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        tight_layout()
        legend(scatterpoints=1,prop={'size':12},loc=2)
#         ylim(0,5)

        if save:
            savefig('{0}/CDF(lyaW_blue_vs_red).pdf'.format(saveDirectory),format='pdf')
        else:
            show()
    
#########################################################################################
#########################################################################################
    # make a CDF of the distribution of Ly-alpha equivalent widths for the 
    # associated sample, splitting on red vs blue shifted absorbers
    #
    # USE THE BUILT IN CDF FUNCTION INSTEAD
    
    plotWCDF_alt_dif = True
    save = True
    
    if plotWCDF_alt_dif:
#         fig = figure(figsize=(2,8))
        fig = figure()
        ax = fig.add_subplot(111)
        
        lyaWArray = array(lyaWList)
        lyaWAmbArray = array(lyaWAmbList)
        
        envArray = array(envList)
        envAmbArray = array(envAmbList)
        
        bins = arange(0,max(lyaWArray),5)

        
        blues = []
        reds = []
        for d,l in zip(difList,lyaWList):
            if float(d) >0:
                # blueshifted ABSORBER: dif = v_gal - v_absorber
                blues.append(float(l))
            else:
                reds.append(float(l))
        
#         n_reds, bins_reds, patches_reds = hist(reds, bins, normed=1, histtype="step",\
#         cumulative=True,color='red',lw=1,label=r'$\rm Redshifted ~CDF$')
#         
#         n_blues, bins_bluess, patches_bluess = hist(blues, bins, normed=1, histtype="step",\
#         cumulative=True,lw=1,label=r'$\rm Blueshifted ~CDF$',color='blue')

        n_1, bins_1, patches_1 = hist(lyaWList, bins, normed=1, histtype="step",\
        cumulative=True,color='black',lw=1,label=r'$\rm Associated ~Ly\alpha$')
        
        n_amb, bins_amb, patches_amb = hist(lyaWAmbList, bins, normed=1, histtype="step",\
        cumulative=True,lw=1,label=r'$\rm Ambiguous ~Ly\alpha$',color='grey')

#         n_2, bins_2, patches_2 = hist(lyaW_2List, bins, normed=1, histtype="step",\
#         cumulative=True,lw=1,label=r'$\rm Include=2 ~CDF$',color='purple')
#         
#         n_3, bins_3, patches_3 = hist(lyaW_3List, bins, normed=1, histtype="step",\
#         cumulative=True,lw=1,label=r'$\rm Include=3 ~CDF$',color='red')
    
#         # plot the y-values against the sorted data
#         plt.plot(sorted_reds,yvals_reds,color='red',lw=3,label=r'$\rm Redshifted ~CDF$')
#         plt.plot(sorted_blues,yvals_blues,color='blue',lw=3,label=r'$\rm Blueshifted ~CDF$')
        
        xlabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ylabel(r'$\rm CDF$')
        
        # format all the axis and stuff
        # X-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)

        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # Y-axis
        majorLocator   = MultipleLocator(0.1)
        majorFormatter = FormatStrFormatter(r'$\rm %s$ ')
        minorLocator   = MultipleLocator(0.05)
        
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        tight_layout()
        legend(scatterpoints=1,prop={'size':12},loc='lower right')
#         xlim(0,max(lyaWArray)+1)
        xlim(0,1000)
    
        if save:
            savefig('{0}/CDF(lyaW_assoc_vs_amb)_maxEnv{1}.pdf'.format(saveDirectory,maxEnv),format='pdf')
        else:
            show()
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    