#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotW_impact.py, v 5.2 04/21/2016

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

fontScale = 15
rc('text', usetex=True)
rc('font', size=15, family='serif', weight=450)
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
        resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_9_edit2.csv'
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots3/'

    elif getpass.getuser() == 'frenchd':
        pickleFilename = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/pilotData2.p'
#         resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_8_edit2.csv'
#         saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots2/'
        resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_9_edit2.csv'
        saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots3/'

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
########################################################################################
    
    # plot equivalent width as a function of impact parameter, splitting up red and 
    # blue shifted absorption
    #
    
    plotW_b = False
    save = False
    
    if plotW_b:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        alpha = 0.85
        
        for d,i,w in zip(difList,impactList,lyaWList):
            count +=1
            if d>0:
                # galaxy is behind absorber, so gas is blue shifted
                color = 'Blue'
                if countb == 0:
                    countb +=1
                    plotb = ax.scatter(i,w,c='Blue',s=50,label= labelb,alpha=alpha)
            if d<0:
                # gas is red shifted compared to galaxy
                color = 'Red'
                if countr == 0:
                    countr +=1
                    plotr = ax.scatter(i,w,c='Red',s=50,label= labelr,alpha=alpha)
                    
            plot1 = scatter(i,w,c=color,s = 50,alpha = alpha)
        
#         title('W(impact parameter) for red and blue shifted absorption')
        xlabel(r'$\rm \rho$ (kpc)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(-1,501)
        ax.legend(scatterpoints=1,prop={'size':12},loc=1)
        
        if save:
            savefig('{0}/W(impact)_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()

##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter, split between
    # red and blue shifted absorption, overplot average histograms
    #
    
    plotW_impact_hist = True
    save = True
    
    if plotW_impact_hist:
        fig = figure()
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        binSize = 100
        alpha = 0.85
        
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        
        placeArrayr = zeros(5)
        placeCountr = zeros(5)
        placeArrayb = zeros(5)
        placeCountb = zeros(5)
        
        xVals = []
        
        for d,i,w,v in zip(difList,impactList,lyaWList,virList):
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(v):
                if d!=-99 and i!=-99 and w!=-99 and v!=-99:
                    xVal = float(i)
                    yVal = float(w)
                    
                    xVals.append(xVal)
                    
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        
                        # which bin does it belong too?
                        place = xVal/binSize
                        print 'place: ',place
                        placeArrayb[place] += yVal
                        print 'placeArrayb: ',placeArrayb
                        placeCountb[place] +=1.
                        print 'placecountb: ',placeCountb
                        
                        if countb == 0:
                            countb +=1
#                             plotb = ax.scatter(v,w,c='Blue',s=50,label= labelb)
                            plotb = ax.scatter(xVal,yVal,c='Blue',s=50,alpha=alpha)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        
                        # which bin does it belong too?
                        place = xVal/binSize
                        placeArrayr[place] += yVal
                        placeCountr[place] +=1.
                        
                        if countr == 0:
                            countr +=1
#                             plotr = ax.scatter(v,w,c='Red',s=50,label= labelr)
                            plotr = ax.scatter(xVal,yVal,c='Red',s=50,alpha=alpha)

                    plot1 = scatter(xVal,yVal,c=color,s=50,alpha=alpha)
                    
        rHist = placeArrayr/placeCountr
        print 'rHist: ',rHist
        bHist = placeArrayb/placeCountb
        print 'bHist: ',bHist
        
        totalrHist = []
        totalrVir = []
        totalbHist = []
        totalbVir = []
        
        print 'xvals: ',xVals
        
        for r,v in zip(rHist,arange(0,max(xVals),binSize)):
            if not isNumber(r):
                r = 0
            
            totalrHist.append(r)
            totalrHist.append(r)

            totalrVir.append(v)
            totalrVir.append(v+binSize)
            
        for b,v in zip(bHist,arange(0,max(xVals),binSize)):
            if not isNumber(b):
                b = 0
                
            totalbHist.append(b)
            totalbHist.append(b)

            totalbVir.append(v)
            totalbVir.append(v+binSize)
        
        print 'totalrVir: ',totalrVir
        print 'totalrHist: ',totalrHist
        print
        print 'totalbVir: ',totalbVir
        print 'totalbHist: ',totalbHist
        print
        
        # x-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator   = MultipleLocator(50)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        
        plot2 = ax.plot(totalrVir,totalrHist,c='Red',lw=2.5,ls='dotted',\
        label=r'$\rm Average ~ Redshifted ~ EW$',alpha=alpha)
        
        plot3 = ax.plot(totalbVir,totalbHist,c='Blue',lw=1.5,ls='dashed',\
        label=r'$\rm Average ~ Blueshifted ~ EW$',alpha=alpha)
        
        xlabel(r'$\rm \rho ~ [kpc]$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ax.legend(scatterpoints=1,prop={'size':14},loc=1,fancybox=True)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1200)
        xlim(0,500)

        if save:
            savefig('{0}/W(impact)_avgHistograms.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()


########################################################################################
########################################################################################

    # plot equivalent width as a function of impact parameter/diameter, split between
    # red and blue shifted absorption
    #
    
    plotW_b_diam= False
    save = False
    
    if plotW_b_diam:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        alpha = 0.85
        for d,i,w,m in zip(difList,impactList,lyaWList,majList):
        
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(m):
                if d !=-99 and i !=-99 and w!=-99 and m!=-99:
                    count +=1
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(i/m,w,c='Blue',s=50,label= labelb,alpha=alpha)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(i/m,w,c='Red',s=50,label= labelr,alpha=alpha)
                    
                    plot1 = scatter(i/m,w,c=color,s = 50,alpha=alpha)
            
        # make the legend work properly
#         labelr = 'Red Shifted Absorber'
#         labelb = "Blue Shifted Absorber"
#         plotb = scatter(i[countb]/m[countb],w[countb],c='Blue',s=50,label= labelb)
#         plotr = scatter(i[countr]/m[countr],w[countr],c='Red',s=50,label= labelr)
        
#         title('W(impact/diameter) for red and blue shifted absorption')
        xlabel('Impact Parameter / Diameter')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,max(lyaWList)+100)
#         xlim(-1,150)

        ax.legend(scatterpoints=1)
        
        if save:
            savefig('{0}/W(impact_diam)_dif_cut.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
    
    
##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter/R_vir, split between
    # red and blue shifted absorption
    #
    
    plotW_b_vir= False
    save = False
    alpha = 0.75
    
    if plotW_b_vir:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,i,w,m in zip(difList,impactList,lyaWList,virList):
        
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(m):
                if d !=-99 and i !=-99 and w!=-99 and m!=-99:
                    count +=1
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(i/m,w,c='Blue',s=50,label= labelb,alpha=alpha)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(i/m,w,c='Red',s=50,label= labelr,alpha=alpha)
                    
                    plot1 = scatter(i/m,w,c=color,s = 50,alpha=alpha)
            
        # make the legend work properly
#         labelr = 'Red Shifted Absorber'
#         labelb = "Blue Shifted Absorber"
#         plotb = scatter(i[countb]/m[countb],w[countb],c='Blue',s=50,label= labelb)
#         plotr = scatter(i[countr]/m[countr],w[countr],c='Red',s=50,label= labelr)
        
#         title('W(impact/R_vir) for red and blue shifted absorption')
        xlabel(r'$\rm Impact Parameter / R_{vir}$')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,max(lyaWList)+100)
        
        # cut out the one outlier?
#         xlim(0,13)

        # or plot the whole range?
#         xlim(-0.2,30)
        
        ax.legend(scatterpoints=1)
        
        if save:
            # cut off the one outlier?
#             savefig('{0}/W(impact_vir)_dif_cut.pdf'.format(saveDirectory),format='pdf')
            
            # or plot the whole range?
            savefig('{0}/W(impact_vir)_dif_cut_lighter.pdf'.format(saveDirectory),format='pdf')

        else:
            show()


##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter/R_vir, split between
    # red and blue shifted absorption, overplot average histograms
    #
    
    plotW_impact_vir_hist = True
    save = True
    
    if plotW_impact_vir_hist:
        fig = figure()
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        binSize = 0.5
        alpha = 0.85
        
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        
        placeArrayr = zeros(6)
        placeCountr = zeros(6)
        placeArrayb = zeros(6)
        placeCountb = zeros(6)
        
        xVals = []
        
        for d,i,w,v in zip(difList,impactList,lyaWList,virList):
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(v):
                if d!=-99 and i!=-99 and w!=-99 and v!=-99:
                    xVal = float(i)/float(v)
                    yVal = float(w)
                    
                    xVals.append(xVal)
                    
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        
                        # which bin does it belong too?
                        place = xVal/binSize
                        print 'place: ',place
                        placeArrayb[place] += yVal
                        print 'placeArrayb: ',placeArrayb
                        placeCountb[place] +=1.
                        print 'placecountb: ',placeCountb
                        
                        if countb == 0:
                            countb +=1
#                             plotb = ax.scatter(v,w,c='Blue',s=50,label= labelb)
                            plotb = ax.scatter(xVal,yVal,c='Blue',s=50,alpha=alpha)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        
                        # which bin does it belong too?
                        place = xVal/binSize
                        placeArrayr[place] += yVal
                        placeCountr[place] +=1.
                        
                        if countr == 0:
                            countr +=1
#                             plotr = ax.scatter(v,w,c='Red',s=50,label= labelr)
                            plotr = ax.scatter(xVal,yVal,c='Red',s=50,alpha=alpha)

                    plot1 = scatter(xVal,yVal,c=color,s=50,alpha=alpha)
                    
        rHist = placeArrayr/placeCountr
        print 'rHist: ',rHist
        bHist = placeArrayb/placeCountb
        print 'bHist: ',bHist
        
        totalrHist = []
        totalrVir = []
        totalbHist = []
        totalbVir = []
        
        print 'xvals: ',xVals
        
        for r,v in zip(rHist,arange(0,max(xVals),binSize)):
            if not isNumber(r):
                r = 0
            
            totalrHist.append(r)
            totalrHist.append(r)

            totalrVir.append(v)
            totalrVir.append(v+binSize)
            
        for b,v in zip(bHist,arange(0,max(xVals),binSize)):
            if not isNumber(b):
                b = 0
                
            totalbHist.append(b)
            totalbHist.append(b)

            totalbVir.append(v)
            totalbVir.append(v+binSize)
        
        print 'totalrVir: ',totalrVir
        print 'totalrHist: ',totalrHist
        print
        print 'totalbVir: ',totalbVir
        print 'totalbHist: ',totalbHist
        print
        
        
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter('%0.1f')
        minorLocator   = MultipleLocator(0.25)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        
        plot2 = ax.plot(totalrVir,totalrHist,c='Red',lw=2.5,ls='dotted',\
        label=r'$\rm Average ~ Redshifted ~ EW$',alpha=alpha)
        
        plot3 = ax.plot(totalbVir,totalbHist,c='Blue',lw=1.5,ls='dashed',\
        label=r'$\rm Average ~ Blueshifted ~ EW$',alpha=alpha)
        
        xlabel(r'$\rm \rho / R_{vir}$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ax.legend(scatterpoints=1,prop={'size':14},loc=1,fancybox=True)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1200)
        xlim(0,3.0)

        if save:
            savefig('{0}/W(impact_vir)_avgHistograms.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


if __name__=="__main__":
    # do the work
    main()
    