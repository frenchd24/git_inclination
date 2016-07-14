#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotImpactHist_Diam2.py, v 5.6 7/14/16

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
    
    
    if getpass.getuser() == 'David':
        pickleFilename = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/pilotData2.p'
#         resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_8_edit2.csv'
#         saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots2/'
        resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit.csv'
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots4/'

    elif getpass.getuser() == 'frenchd':
        pickleFilename = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/pilotData2.p'
#         resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_8_edit2.csv'
#         saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots2/'
        resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit.csv'
        saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots4/'

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
    
    maxEnv = 300
    minL = 0.3
    
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
    

##########################################################################################
##########################################################################################
    # make a histogram of the distribution of impact parameters for associated galaxies 
    # normalized by major diameter
    #
    
    plotImpactHist_Diam = False
    save = False
    
    if plotImpactHist_Diam:
        fig = figure(figsize=(10,2))
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(111)
        
#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        bins = [0,10,20,30,40,50,60,70,80,90]
        
        impactArray = np.array([])
        majArray = np.array([])
        
        # make sure all the values are okay
        for i,m in zip(impactList,majList):
            if isNumber(i) and isNumber(m):
                if i !=-99 and m !=-99:
                    impactArray = append(impactArray,float(i))
                    majArray = append(majArray,float(m))
        
        normalizedImpactArray = impactArray/majArray
        
        plot1 = hist(normalizedImpactArray,bins=bins,histtype='bar')
        
        title('Distribution of impact parameters')
        xlabel('Impact Parameter (kpc)')
        ylabel('Number')
        ax.tick_params(axis='x', labelsize=0)
        ax.tick_params(axis='y',labelsize=8)
#         ylim(0,5)

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
    save = False
    
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
    # make a histogram of the distribution of impact parameters for associated galaxies 
    # normalized by virial radius, and split into red and blue shifted bins
    #
    
    plotImpactHist_Vir_dif = False
    save = False
    
    if plotImpactHist_Vir_dif:
        fig = figure(figsize=(10,5))
        ax = fig.add_subplot(211)
        subplots_adjust(hspace=0.200)

#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
#         bins = [0,5,10,15,20,25,30,35,40]

        binsize = 0.2
        bins = arange(0,3.0,binsize)
        
        reds = np.array([])
        blues = np.array([])
        
        # make sure all the values are okay
        for i,v,d in zip(impactList,virList,difList):
            if isNumber(i) and isNumber(v):
                if i !=-99 and v !=-99:
                    val = float(i)/float(v)
                    
                    # if blueshifted
                    if d>=0:
                        blues = append(blues,val)
                    
                    # for redshifted
                    else:
                        reds = append(reds,val)
        
#         title('Distribution of impact parameters')

        plot1 = hist(blues,bins=bins,histtype='bar',color='blue',label='Blueshifted',alpha=0.75)
#         xlabel('Impact Parameter/R_vir (blueshifted)')
        ylabel('Number')
        legend(scatterpoints=1,prop={'size':12})
        ax.tick_params(axis='y',labelsize=11)
        ax.tick_params(axis='x', labelsize=0)

        ax = fig.add_subplot(212)
        plot2 = hist(reds,bins=bins,histtype='bar',color="red",label='Redshifted',alpha=0.75)
        xlabel(r'$\rho / R_{vir}$')
        ylabel('Number')
        
        legend(scatterpoints=1,prop={'size':12})
        
        ax.tick_params(axis='x', labelsize=11)
        ax.tick_params(axis='y',labelsize=11)
#         ylim(0,5)

        if save:
            savefig('{0}/hist(Impact_vir_dif)_bin_{1}.pdf'.format(saveDirectory,binsize),format='pdf')
        else:
            show()


##########################################################################################
##########################################################################################
    # make a histogram of the distribution of impact parameters for associated galaxies, 
    # and split into red and blue shifted bins
    #
    
    plotImpactHist_dif = False
    save = False
    
    if plotImpactHist_dif:
        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        subplots_adjust(hspace=0.200)
        alpha = 0.75

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
        tight_layout()
        
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
        tight_layout()
        
        if save:
            savefig('{0}/hist(Impact_red)_bin_{1}.pdf'.format(saveDirectory,binsize),format='pdf')
        else:
            show()


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


if __name__=="__main__":
    # do the work
    main()
    