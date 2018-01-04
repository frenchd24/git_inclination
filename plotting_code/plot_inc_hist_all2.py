#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_inc_hist_all2.py, v6 01/03/18

Plot histograms of the inclinations of the associated galaxies (and all galaxies too)

    - this now combines plotCosIncDifHist_full2.py, plotFancyIncHist_full2.py, 
    plotFancyCosIncDifHist_full2.py, plotCosIncHist_full2.py, plotIncHist_full2.py,
    and plotFancyCosIncHist_full2.py all into this one program (split into sections below)

v5: updated to work with the new, automatically updated LG_correlation_combined5.csv
    (12/04/15) - original updates to the individual files
    
v5.1: updated for LG_correlation_combined5_8_edit2.csv
    (2/17/2016)
    
v5.2: remake plots with v_hel instead of vcorr (4/22/16)

v5.3: remake plots with new large galaxy sample (7/13/16) -> /plots4/

v5.4: include ability to limit results based on environment number (7/14/16)
    also likelihood limits

v6: Updated for AAS_2017 results (01/03/18)

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

        pickleFilename = '/Users/frenchd/Research/inclination/git_inclination/picklePilot_plusSALTcut.p'
        gtPickleFilename = '/Users/frenchd/Research/inclination/git_inclination/pickleGT.p'
        saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/plotting_code/figs'


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
    
    maxEnv = 3000
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
    
    print 'initial len(Lya_vL): ',len(Lya_vL)
    print

    i = -1
    for inc,inc_vir,inc_cus in zip(includeL,include_virL,include_customL):
        i+=1
        go = False
        if match:
            if virInclude == inc_vir and cusInclude == inc_cus:
                go = True
            else:
                go = False
                
        else:
            if virInclude and inc_vir:
                go = True
                
            elif cusInclude and inc_cus:
                go = True
                
            elif finalInclude and inc:
                go = True
            
            else:
                go = False
        
        if go:
            RA_agn = RA_agnL[i]
            Dec_agn = Dec_agnL[i]
            RA_gal = RA_galL[i]
            Dec_gal = Dec_galL[i]
            lyaV = Lya_vL[i]
            lyaW = Lya_WL[i]
            lyaW_err = lyaW*0.1
            env = environmentL[i]
            galaxyName = galaxyNameL[i]
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
                print az
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

        
    # lists for the full galaxy dataset    
#     allPA = gtDict['allPA']
#     allInclinations = gtDict['allInclinations']
#     allCosInclinations = gtDict['allCosInclinations']
#     allFancyInclinations = gtDict['allFancyInclinations']
#     allCosFancyInclinations = gtDict['allCosFancyInclinations']
#     allDiameter = gtDict['allDiameters']


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

    total = 0
    totalNo = 0
    totalYes = 0
    totalIsolated = 0
    totalGroup = 0
    

########################################################################################
########################################################################################
    # cos(inclination) histograms for redshifted and blueshifted distributions of absorbers
    # as well as the full data table
    # 
    # from original file: plotCosIncDifHist_full2 (12/28/15)
    
    plotCosIncDifHist_full = False
    save = False
    
    if plotCosIncDifHist_full:
    
#         fig = figure(figsize=(12,10))
#         fig = figure()
#         subplots_adjust(hspace=0.200)

#         subplots_adjust(hspace=0.7)
#         ax = fig.add_subplot(311)
#         subplot(311)

        bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        blue = []
        red = []
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []

        
        for d,i,l,e in zip(difList,cosIncList,lyaWList,lyaErrList):
            if d >=0:
                # blue shifted absorber, but galaxy is REDSHIFTED
                print 'd: ',d
                blue.append(i)
                blueLya.append(l)
                blueLyaErr.append(e)
            else:
                # red shifted absorber, but galaxy is BLUESHIFTED
                red.append(i)
                redLya.append(l)
                redLyaErr.append(e)
        
#         plot1 = hist([red,blue],bins=bins,histtype='bar',color=['Red','blue'],alpha=0.5)
 #        title('Red shifted aborption: Galaxies')
#         xlabel('Inclination (deg)')
#         ylim(0,6)
#         ylabel('Number')

        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        subplots_adjust(top = 0.85,bottom = 0.25)
        
        hist(red,bins=bins,histtype='bar',color='red',alpha = 0.9)
        print 'average red: ',average(redLya)
        print 'median red: ', median(redLya)
        print 'average(error): ',average(redLyaErr)
        print 'median(error): ',median(redLyaErr)
        print
        title('Red shifted absorption: Galaxies')
        xlabel('Cos(inclination) = b/a')
#         ylim(0,6)
        ylabel('Number')
        
        if save:
            savefig('{0}/hist(cos(inclination))_dif_red.pdf'.format(saveDirectory),format='pdf')
        else:
            show()


#         fig.add_subplot(212)
#         subplot(212)
        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        subplots_adjust(top = 0.85,bottom=0.25)

        hist(blue,bins=bins,histtype='bar',color='Blue',alpha = 0.9)
        print 'average blue: ',average(blueLya)
        print 'median blue: ',median(blueLya)
        print 'avereage(error): ',average(blueLyaErr)
        print 'median(error): ',median(blueLyaErr)
        
        title('Blue shifted absorption: Galaxies')
        xlabel('Cos(inclination) = b/a')
#         ylim(0,5)
        ylabel('Number')
        
        if save:
            savefig('{0}/hist(cos(inclination))_dif_blue.pdf'.format(saveDirectory),format='pdf')
        else:
            show()

#         ax = fig.add_subplot(312)
#         subplot(312)    
        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        subplots_adjust(top = 0.85,bottom=0.25)

        hist(allCosInclinations,bins=bins,histtype='bar',color = 'green',alpha=0.9)
        title('Full galaxy sample inclinations')
        xlabel('Cos(inclination) = b/a')
        ylabel('Number')
#         tight_layout()

        if save:
            savefig('{0}/hist(cos(inclination))_fulldataset.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
    

#########################################################################################
#########################################################################################
    # plot histograms of the fancy inclination for both associated galaxies and the 
    # full galaxy data set
    #
    # originally from: plotFancyIncHist_full2.py
    #
    # All this shows is that the associated galaxies sample the full distribution pretty 
    # well
    
    plotFancyIncHist_full = False
    save = False
    
    if plotFancyIncHist_full:
        fig = figure(figsize=(10,3))
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(111)
        alpha=0.75
        
        bins = [0,10,20,30,40,50,60,70,80,90]
        subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)

    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
    
#         x-axis
#         majorLocator   = MultipleLocator(10)
#         majorFormatter = FormatStrFormatter('%d')
#         minorLocator   = MultipleLocator(2)
#         ax.xaxis.set_major_locator(majorLocator)
#         ax.xaxis.set_major_formatter(majorFormatter)
#         ax.xaxis.set_minor_locator(minorLocator)
#         
#         y-axis
#         majorLocator   = MultipleLocator(.08)
#         majorFormatter = FormatStrFormatter('%d')
#         minorLocator   = MultipleLocator(.01)
#         ax.yaxis.set_major_locator(majorLocator)
#         ax.yaxis.set_major_formatter(majorFormatter)
#         ax.yaxis.set_minor_locator(minorLocator)
#     
#         plot1 = hist(fancyIncList,bins=bins,histtype='bar',label='Associated Galaxies',\
#         normed=True,alpha=alpha)
#         xlabel(r'Galaxy Inclination (deg)')
#         ylabel(r'Number')
#         legend(loc=2,fontsize=16)
# 
#         tight_layout()
# 
#         if save:
#             savefig('{0}/hist(fancy_inclination)_associated2.pdf'.format(saveDirectory),format='pdf')
#         else:
#             show()

        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator   = MultipleLocator(2)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(5000)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator   = MultipleLocator(1000)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        plot1 = hist(allFancyInclinations,bins=bins,histtype='bar',label="All Galaxies",alpha=alpha,color='green')
        
        xlabel(r'Galaxy Inclination (deg)')
        ylabel(r'Number')
        legend(loc=2,fontsize=16)
        
        tight_layout()

        if save:
            savefig('{0}/hist(fancy_inclination)_all2.pdf'.format(saveDirectory),format='pdf')
        else:
            show()


#########################################################################################
#########################################################################################
    # cos(inclination) histograms for redshifted vs blueshifted distributions of absorbers
    #
    # originally from plotFancyCosIncDifHist_full2.py
    #
    # looks about the same as the above (1st) function
    #
    
    plotCosFancyIncDifHist_full = False
    save = False
    
    if plotCosFancyIncDifHist_full:
    
#         fig = figure(figsize=(12,10))
#         fig = figure()
#         subplots_adjust(hspace=0.200)

#         subplots_adjust(hspace=0.7)
#         ax = fig.add_subplot(311)
#         subplot(311)
        bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        blue = []
        red = []
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []

        
        for d,i,l,e in zip(difList,cosFancyIncList,lyaWList,lyaErrList):
            if d >=0:
                # blue shifted absorber, but galaxy is REDSHIFTED
                print 'd: ',d
                blue.append(i)
                blueLya.append(l)
                blueLyaErr.append(e)
            else:
                # red shifted absorber, but galaxy is BLUESHIFTED
                red.append(i)
                redLya.append(l)
                redLyaErr.append(e)
        
#         plot1 = hist([red,blue],bins=bins,histtype='bar',color=['Red','blue'],alpha=0.5)
 #        title('Red shifted aborption: Galaxies')
#         xlabel('Inclination (deg)')
#         ylim(0,6)
#         ylabel('Number')
        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        subplots_adjust(top = 0.85,bottom = 0.25)
        
        hist(red,bins=bins,histtype='bar',color='red',alpha = 0.9)
        print 'average red: ',average(redLya)
        print 'median red: ', median(redLya)
        print 'average(error): ',average(redLyaErr)
        print 'median(error): ',median(redLyaErr)
        print
        title('Red shifted absorption: Galaxies')
        xlabel('Cos(inclination) = b/a')
#         ylim(0,6)
        ylabel('Number')
        
        if save:
            savefig('{0}/hist(cos(fancy_inclination))_dif_red.pdf'.format(saveDirectory),format='pdf')
        else:
            show()


#         fig.add_subplot(212)
#         subplot(212)
        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        subplots_adjust(top = 0.85,bottom=0.25)

        hist(blue,bins=bins,histtype='bar',color='Blue',alpha = 0.9)
        print 'average blue: ',average(blueLya)
        print 'median blue: ',median(blueLya)
        print 'avereage(error): ',average(blueLyaErr)
        print 'median(error): ',median(blueLyaErr)
        
        title('Blue shifted absorption: Galaxies')
        xlabel('Cos(inclination) = b/a')
#         ylim(0,5)
        ylabel('Number')
        
        if save:
            savefig('{0}/hist(cos(fancy_inclination))_dif_blue.pdf'.format(saveDirectory),format='pdf')
        else:
            show()

#         ax = fig.add_subplot(312)
#         subplot(312)    
        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        subplots_adjust(top = 0.85,bottom=0.25)

        hist(allCosFancyInclinations,bins=bins,histtype='bar',color = 'green',alpha=0.9)
        title('Full galaxy sample inclinations')
        xlabel('Cos(inclination) = b/a')
        ylabel('Number')
#         tight_layout()

        if save:
            savefig('{0}/hist(cos(fancy_inclination))_fulldataset.pdf'.format(saveDirectory),format='pdf')
        else:
            show()

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
    # fancyInclination histograms for redshifted vs blueshifted distributions of absorbers
    # as well as the whole table
    #
    
    plotFancyIncDifHist_full = False
    save = False
    
    if plotFancyIncDifHist_full:
    
        fig = figure(figsize=(10,6))
        subplots_adjust(hspace=0.300)

#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        bins = arange(0,100,10)
        blue = []
        red = []
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        for d,i,l,e in zip(difList,fancyIncList,lyaWList,lyaErrList):
            if d >=0:
                # blue shifted absorber, but galaxy is REDSHIFTED
                print 'd: ',d
                blue.append(i)
                blueLya.append(l)
                blueLyaErr.append(e)
            else:
                # red shifted absorber, but galaxy is BLUESHIFTED
                red.append(i)
                redLya.append(l)
                redLyaErr.append(e)
                
        
        ax = fig.add_subplot(311)        
        hist(red,bins=bins,histtype='bar',color='red',alpha = 0.9,label='Redshifted')
        legend(scatterpoints=1,prop={'size':12},loc=2)

        ax = fig.add_subplot(312)
        hist(blue,bins=bins,histtype='bar',color='Blue',alpha = 0.9,label='Blueshifted')
        legend(scatterpoints=1,prop={'size':12},loc=2)
        ylabel('Number')

        ax = fig.add_subplot(313)
        hist(allFancyInclinations,bins=bins,histtype='bar',color = 'green',alpha=0.9,label='All')
        legend(scatterpoints=1,prop={'size':12},loc=2)
        xlabel('Inclination (deg)')
#         tight_layout()

        if save:
            savefig('{0}/hist(fancy_inclination)_red_blue_all.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
            
            
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
    # fancyInclination histograms for redshifted vs blueshifted distributions of absorbers
    # plotted together and separately all in one, then the whole table as a subplot
    #
    
    plotFancyIncDifHist_full_all = True
    save = True
    
    if plotFancyIncDifHist_full_all:
    
        fig = figure(figsize=(10,6))
        subplots_adjust(hspace=0.200)

#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        bins = arange(0,100,10)
        blue = []
        red = []
        
        alpha = 0.6
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        for d,i,l,e in zip(difList,fancyIncList,lyaWList,lyaErrList):
            if d >=0:
                # blue shifted absorber, but galaxy is REDSHIFTED
#                 print 'd: ',d
                blue.append(i)
                blueLya.append(l)
                blueLyaErr.append(e)
            else:
                # red shifted absorber, but galaxy is BLUESHIFTED
                red.append(i)
                redLya.append(l)
                redLyaErr.append(e)
                
                
        # just associated
        ax = fig.add_subplot(211)
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(5)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(2.5)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

#         hist(red,bins=bins,histtype='bar',color='red',hatch='\\',lw=1.5,alpha = alpha,label=r'$\rm Redshifted$')        
#         hist(blue,bins=bins,histtype='bar',color='Blue',lw=1.5,alpha = alpha,label=r'$\rm Blueshifted$')
        hist(blue,bins=bins,histtype='bar',color='Blue',hatch='\\',lw=1.7,alpha = alpha+0.25,label=r'$\rm Blueshifted$')
        hist(red,bins=bins,histtype='bar',color='red',lw=1.7,alpha = alpha,label=r'$\rm Redshifted$')        
        hist(fancyIncList,bins=bins,histtype='step',color='Black',lw=2.1,alpha = 0.9,label=r'$\rm All ~ Associated$')
        
        legend(scatterpoints=1,prop={'size':15},loc=2,fancybox=True)
        ylabel(r'$\rm Number$')
        

        # full table
        ax = fig.add_subplot(212)
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(5000)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(2500)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        hist(allFancyInclinations,bins=bins,histtype='bar',lw=1.5,color = 'green',alpha=0.9,label=r'$\rm All$')
        legend(scatterpoints=1,prop={'size':15},loc=2,fancybox=True)
        xlabel(r'$\rm Galaxy ~ Inclination ~ [deg]$')
        ylabel(r'$\rm Number$')
#         tight_layout()

        if save:
            savefig('{0}/hist(fancy_inclination)_red_blue_full_all.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()
            

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
    # cos(fancyInclination) histograms for redshifted vs blueshifted distributions of absorbers
    # plotted together and separately all in one, then the whole table as a subplot
    #
    
    plotCosFancyIncDifHist_full_all = False
    save = False
    
    if plotCosFancyIncDifHist_full_all:
    
        fig = figure(figsize=(10,6))
        subplots_adjust(hspace=0.200)

#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        bins = arange(0,1.0,0.1)
        blue = []
        red = []
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        for d,i,l,e in zip(difList,cosFancyIncList,lyaWList,lyaErrList):
            if d >=0:
                # blue shifted absorber, but galaxy is REDSHIFTED
                print 'd: ',d
                blue.append(i)
                blueLya.append(l)
                blueLyaErr.append(e)
            else:
                # red shifted absorber, but galaxy is BLUESHIFTED
                red.append(i)
                redLya.append(l)
                redLyaErr.append(e)
                
                
        # just associated
        ax = fig.add_subplot(211)
        
        # x-axis
        majorLocator   = MultipleLocator(0.1)
        majorFormatter = FormatStrFormatter(r'$\rm %s$')
        minorLocator   = MultipleLocator(0.05)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(5)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        hist(blue,bins=bins,histtype='bar',color='Blue',lw=1.5,alpha = 0.7,label=r'$\rm Blueshifted$')
        hist(red,bins=bins,histtype='bar',color='red',lw=1.5,alpha = 0.7,label=r'$\rm Redshifted$')        
        hist(cosFancyIncList,bins=bins,histtype='step',color='Black',lw=2.5,alpha = 0.9,label=r'$\rm All ~ Associated$')
        
        legend(scatterpoints=1,prop={'size':14},loc=2,fancybox=True)
        ylabel(r'$\rm Number$')
        

        # full table
        ax = fig.add_subplot(212)
        
        # x-axis
        majorLocator   = MultipleLocator(0.1)
        majorFormatter = FormatStrFormatter(r'$\rm %s$')
        minorLocator   = MultipleLocator(0.05)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(5000)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1000)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        hist(allCosFancyInclinations,bins=bins,histtype='bar',lw=1.5,color = 'green',alpha=0.9,label=r'$\rm All$')
        legend(scatterpoints=1,prop={'size':14},loc=2,fancybox=True)
        xlabel(r'$\rm Cos[inc]$')
        ylabel(r'$\rm Number$')
#         tight_layout()

        if save:
            savefig('{0}/hist(cos(fancy_inclination))_red_blue_full_all.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()
            
            
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
    # fancyInclination CDF for redshifted vs blueshifted distributions of absorbers
    # plotted together and separately all in one, then the whole table as a subplot
    #
    
    plotFancyIncDifCDF_full_all = False
    save = False
    
    if plotFancyIncDifCDF_full_all:
    
        fig = figure(figsize=(10,6))
        subplots_adjust(hspace=0.200)

        bins = arange(0,90,0.5)
        blues = []
        reds = []
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        for d,i,l,e in zip(difList,fancyIncList,lyaWList,lyaErrList):
            if d >=0:
                # blue shifted absorber, but galaxy is REDSHIFTED
                print 'd: ',d
                blues.append(i)
                blueLya.append(l)
                blueLyaErr.append(e)
            else:
                # red shifted absorber, but galaxy is BLUESHIFTED
                reds.append(i)
                redLya.append(l)
                redLyaErr.append(e)
                
                
        # just associated
        ax = fig.add_subplot(111)
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(0.1)
        majorFormatter = FormatStrFormatter(r'$\rm %s$')
        minorLocator   = MultipleLocator(0.05)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        n_reds, bins_reds, patches_reds = hist(reds, bins, normed=1, histtype="step",\
        cumulative=True,color='red',lw=1,label=r'$\rm Redshifted ~CDF$')
        
        n_blues, bins_blues, patches_blues = hist(blues, bins, normed=1, histtype="step",\
        cumulative=True,lw=1,label=r'$\rm Blueshifted ~CDF$',color='blue')

#         n_blues, bins_blues, patches_blues = hist(fancyIncList, bins=bins, normed=1, histtype="step",\
#         cumulative=True,lw=2,label=r'$\rm Associated$',color='green')
        
        n_all, bins_all, patches_blues = hist(allFancyInclinations, bins=bins, normed=1, \
        histtype='step',cumulative=True, color='Black',lw=1.5,alpha = 0.9,label=r'$\rm All$')

    
        legend(scatterpoints=1,prop={'size':14},loc=2,fancybox=True)
        xlabel(r'$\rm Galaxy ~ Inclination ~ [deg]$')
        ylabel(r'$\rm Number$')
#         tight_layout()

        if save:
            savefig('{0}/CDF(fancy_inclination)_red_blue_full_all.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()


#########################################################################################
#########################################################################################
    # inclination histograms for redshifted vs blueshifted distributions of absorbers
    # as well as the whole table
    #
    
    plotIncDifHist_all = False
    save = False
    
    if plotIncDifHist_all:
    
        fig = figure(figsize=(10,6))
        subplots_adjust(hspace=0.300)

#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        bins = arange(0,100,10)
        blue = []
        red = []
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        print 'incList: ',incList
        
        for d,i,l,e in zip(difList,incList,lyaWList,lyaErrList):
            if d >=0:
                # blue shifted absorber, but galaxy is REDSHIFTED
                print 'd: ',d
                blue.append(i)
                blueLya.append(l)
                blueLyaErr.append(e)
            else:
                # red shifted absorber, but galaxy is BLUESHIFTED
                red.append(i)
                redLya.append(l)
                redLyaErr.append(e)
                
        
        ax = fig.add_subplot(311)        
        hist(red,bins=bins,histtype='bar',color='red',alpha = 0.9,label='Redshifted')
        legend(scatterpoints=1,prop={'size':12},loc=2)

        ax = fig.add_subplot(312)
        hist(blue,bins=bins,histtype='bar',color='Blue',alpha = 0.9,label='Blueshifted')
        legend(scatterpoints=1,prop={'size':12},loc=2)
        ylabel('Number')

        ax = fig.add_subplot(313)
        hist(allFancyInclinations,bins=bins,histtype='bar',color = 'green',alpha=0.9,label='All')
        legend(scatterpoints=1,prop={'size':12},loc=2)
        xlabel('Inclination (deg)')
#         tight_layout()

        if save:
            savefig('{0}/hist(inclination)_red_blue_all.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
            



#########################################################################################
#########################################################################################
    # plot histograms of the fancy inclination for associated galaxies, normalizing by
    # the full galaxy table
    #
    # originally from: plotFancyIncHist_full2.py
    #
    # All this shows is that the associated galaxies sample the full distribution pretty 
    # well
    
    plotFancyIncHist_full_norm = False
    save = False
    
    if plotFancyIncHist_full_norm:
        fig = figure()
        ax = fig.add_subplot(211)
        bins = [0,10,20,30,40,50,60,70,80,90]
        subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)
        
        normCounts, normBins= histogram(allFancyInclinations, bins=bins,normed= True)
        print "normCounts: ",normCounts
        
        dataCounts, dataBins = histogram(fancyIncList,bins=bins)
        
        normedCounts = dataCounts/normCounts
        print 'normedCounts: ',normedCounts
        print len(normedCounts)
        print len(bins)
        
#         plot1 = hist(normedCounts,bins=bins,histtype='bar')
        plot1 = bar([0,10,20,30,40,50,60,70,80],normedCounts)
        ylabel('Number')

        ax = fig.add_subplot(212)
        plot2 = bar([0,10,20,30,40,50,60,70,80],dataCounts)
        xlabel('Galaxy Inclination (deg)')
        ylabel('Number')


#         plot1 = hist(fancyIncList,bins=bins,histtype='bar',weights=normCounts)
#         title('Absorber-associated galaxies')
# #         xlabel('Inclination (deg)')
#         ylabel('Number')

#         ax = fig.add_subplot(212)
#         plot1 = hist(allFancyInclinations,bins=bins,histtype='bar')
#         title('Total galaxy inclination distribution')
#         xlabel('Galaxy Inclination (deg)')
#         ylabel('Number')


        if save:
            savefig('{0}/hist(fancy_inclination_norm).pdf'.format(saveDirectory),format='pdf')
        else:
            show()
            
            
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
    # fancyInclination histograms for redshifted vs blueshifted distributions of absorbers
    # plotted together and separately all in one, then the whole table as a subplot, but
    # only those galaxies with D>25kpc
    #
    
    plotFancyIncDifHist_full_largeOnly = True
    save = False
    
    if plotFancyIncDifHist_full_largeOnly:
    
        fig = figure(figsize=(10,6))
        subplots_adjust(hspace=0.200)

#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        bins = arange(0,100,10)
        blue = []
        red = []
        
        alpha = 0.6
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        for d,i,l,e in zip(difList,fancyIncList,lyaWList,lyaErrList):
            if d >=0:
                # blue shifted absorber, but galaxy is REDSHIFTED
                print 'd: ',d
                blue.append(i)
                blueLya.append(l)
                blueLyaErr.append(e)
            else:
                # red shifted absorber, but galaxy is BLUESHIFTED
                red.append(i)
                redLya.append(l)
                redLyaErr.append(e)
                
                
        # just associated
        ax = fig.add_subplot(211)
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(5)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(2.5)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

#         hist(red,bins=bins,histtype='bar',color='red',hatch='\\',lw=1.5,alpha = alpha,label=r'$\rm Redshifted$')        
#         hist(blue,bins=bins,histtype='bar',color='Blue',lw=1.5,alpha = alpha,label=r'$\rm Blueshifted$')
        hist(blue,bins=bins,histtype='bar',color='Blue',hatch='\\',lw=1.7,alpha = alpha+0.25,label=r'$\rm Blueshifted$')
        hist(red,bins=bins,histtype='bar',color='red',lw=1.7,alpha = alpha,label=r'$\rm Redshifted$')        
        hist(fancyIncList,bins=bins,histtype='step',color='Black',lw=2.1,alpha = 0.9,label=r'$\rm All ~ Associated$')
        
        legend(scatterpoints=1,prop={'size':15},loc=2,fancybox=True)
        ylabel(r'$\rm Number$')
        

        # full table
        ax = fig.add_subplot(212)
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(5000)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(2500)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        hist(allFancyInclinations,bins=bins,histtype='bar',lw=1.5,color = 'green',alpha=0.9,label=r'$\rm All ~(D>25~kpc)$')
        legend(scatterpoints=1,prop={'size':15},loc=2,fancybox=True)
        xlabel(r'$\rm Galaxy ~ Inclination ~ [deg]$')
        ylabel(r'$\rm Number$')
#         tight_layout()

        if save:
            savefig('{0}/hist(fancy_inclination)_red_blue_full_largeOnly.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()
            
            
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


if __name__=="__main__":
    # do the work
    main()
    