#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_inc_hist_all2.py, v 5.1 02/17/2016

    - this now combines plotCosIncDifHist_full2.py, plotFancyIncHist_full2.py, 
    plotFancyCosIncDifHist_full2.py, plotCosIncHist_full2.py, plotIncHist_full2.py,
    and plotFancyCosIncHist_full2.py all into this one program (split into sections below)

v5: updated to work with the new, automatically updated LG_correlation_combined5.csv
    (12/04/15) - original updates to the individual files
    
v5.1: updated for LG_correlation_combined5_8_edit2.csv
    (2/17/2016)
    
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
        resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_8_edit2.csv'
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots2/'

    elif getpass.getuser() == 'frenchd':
        pickleFilename = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/pilotData2.p'
        resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_8_edit2.csv'
        saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots2/'

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
                
                if isNumber(maj) and isNumber(min):
                    q0 = 0.2
                    fancyInc = calculateFancyInclination(maj,min,q0)
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
        fig = figure()
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(211)
        bins = [0,10,20,30,40,50,60,70,80,90]
        subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)

    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        plot1 = hist(fancyIncList,bins=bins,histtype='bar')
        title('Absorber-associated galaxies fancy inclination')
        xlabel('Inclination (deg)')
        ylabel('Number')
#         ylim(0,5)

        ax = fig.add_subplot(212)
        plot1 = hist(allFancyInclinations,bins=bins,histtype='bar')
        title('Full galaxy sample fancy inclination')
        xlabel('Inclination (deg)')
        ylabel('Number')
#         ylim(0,8000)
#         tight_layout()

        if save:
            savefig('{0}/hist(fancy_inclination).pdf'.format(saveDirectory),format='pdf')
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
    # inclination histograms for redshifted vs blueshifted distributions of absorbers
    # as well as the whole table
    #
    
    plotIncDifHist_all = True
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
    # plot histograms of the cos(inclinations) for both associated galaxies and the 
    # full galaxy data set
    #
    # originally from plotCosIncHist_full2.py
    #
    
    plotCosIncHist_full = False
    save = False
    
    if plotCosIncHist_full:
        fig = figure()
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(211)
        bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)

    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        plot1 = hist(cosIncList,bins=bins,histtype='bar')
        title('Absorber-associated galaxies cos(inclination)')
        xlabel('Inclination (deg)')
        ylabel('Number')
#         ylim(0,5)

        ax = fig.add_subplot(212)
        bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        plot1 = hist(allCosInclinations,bins=bins,histtype='bar')
        title('Full galaxy sample cos(inclination)')
        xlabel('Inclination (deg)')
        ylabel('Number')
#         ylim(0,8000)
#         tight_layout()

        if save:
            savefig('{0}/hist(cos(inclination)).pdf'.format(saveDirectory),format='pdf')
        else:
            show()


#########################################################################################
#########################################################################################
    # plot histograms of the cos(inclinations) for both associated galaxies and the 
    # full galaxy data set
    #
    # originally from plotFancyCosIncHist_full2.py
    #
    
    plotFancyCosIncHist_full = False
    save = False
    
    if plotFancyCosIncHist_full:
        fig = figure()
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(211)
        bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)

    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        plot1 = hist(cosFancyIncList,bins=bins,histtype='bar')
        title('Absorber-associated galaxies cos(fancy_inclination)')
        xlabel('Inclination (deg)')
        ylabel('Number')
#         ylim(0,5)

        ax = fig.add_subplot(212)
        bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        plot1 = hist(allCosFancyInclinations,bins=bins,histtype='bar')
        title('Full galaxy sample cos(fancy_inclination)')
        xlabel('Inclination (deg)')
        ylabel('Number')
#         ylim(0,8000)
#         tight_layout()

        if save:
            savefig('{0}/hist(cos(fancy_inclination)).pdf'.format(saveDirectory),format='pdf')
        else:
            show()


########################################################################################
########################################################################################
    # plot histograms of the associated galaxies' inclinations along with that of the full
    # galaxy set
    #
    
    plotIncHist_full = False
    save = False
    
    if plotIncHist_full:
        fig = figure()
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(211)
        bins = [0,10,20,30,40,50,60,70,80,90]
        subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)

        plot1 = hist(incList,bins=bins,histtype='bar')
        title('Absorber-associated galaxies inclinations')
        xlabel('Inclination (deg)')
        ylabel('Number')
#         ylim(0,5)
        
        ax = fig.add_subplot(212)
        bins = [0,10,20,30,40,50,60,70,80,90]

        plot1 = hist(allInclinations,bins=bins,histtype='bar')
        title('Full Galaxy Sample inclinations')
        xlabel('Inclination (deg)')
        ylabel('Number')
#         ylim(0,8000)
#         tight_layout()

        if save:
            savefig('{0}/hist(inclination).pdf'.format(saveDirectory),format='pdf')
        else:
            show()


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


if __name__=="__main__":
    # do the work
    main()
    