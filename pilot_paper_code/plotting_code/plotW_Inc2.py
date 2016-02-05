#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotW_Inc2.py, v 5.0 12/04/2015

This is the plotW_Inc bit from histograms3.py. Now is separated, and loads in a pickle
file of the relevant data, as created by "buildDataLists.py"

Plot EW as a function of inclination, fancy(inc), cos(inc) and cos(fancy(inc)). Also plot
with a colormap.


Previous (from histograms3.py):
    Plot some stuff for the 100largest initial results

    Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

    Updated for the pilot paper (05/06/15)


v5: Updated for the final pilot paper results (12/04/15)
    No longer uses pilotData.p or whatever, just loads the files here
    
    - all the versions of W_vs_inc are here now: plotW_Inc2.py, plotW_CosInc_colorbar2.py,
    plotW_fancyInc2.py, plotW_FancyCosInc2.py, plotW_FancyCosInc_colorbar2.py,
    plotW_CosInc2.py
    
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


    
# def returnLinDiameters(major,minor,distance):
#     # input major and minor in arcsec, distance in Mpc
#     # outputs major and minor in kpc
#     newMajor = math.tan(math.radians(float(major)))*(distance*1000)
#     newMinor = math.tan(math.radians(float(minor)))*(distance*1000)
#     return (newMajor,newMinor)
#     
#     
# 
# def returnAngDiameters(major,minor,distance):
#     # input distances in mpc, major and minor is in kpc
#     # outputs angular diameters in arcsec
#     newMajor = math.atan((float(major)/1000)/float(distance))*(1/3600)
#     newMinor = math.atan((float(minor)/1000)/float(distance))*(1/3600)
#     return (newMajor,newMinor)
    
    

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
            print "l['Na'].partition(' pm ')[2] : ",l['Na'].partition(' pm ')
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
    
    
##########################################################################################
##########################################################################################

    # plot equivalent width as a function of inclination for red and blue shifted
    # absorption
    #
    
    plotW_Inc = False
    save = False
    
    if plotW_Inc:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        
        for d,i,w,m in zip(difList,incList,lyaWList,majList):
            # check if all the values are okay
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(m):
                if d!=-99 and i!=-99 and w!=-99 and m!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(i,w,c='Blue',s=50,label= labelb)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(i,w,c='Red',s=50,label= labelr)
                
                    plot1 = scatter(i,w,c=color,s=50)
            
#         title('W(inclination) for red vs blue shifted absorption')
        xlabel(r'Inclination (deg)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(0,90)
        
        if save:
            savefig('{0}/W(inclination)_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
    
    
##########################################################################################
##########################################################################################

    # plot equivalent width as a function of cos(inclination) with red and blue shifted
    # absorption represented by a color bar
    #
    
    plotW_CosInc_colorbar= False
    save = False
    
    if plotW_CosInc_colorbar:
        # colormap the velocity difference of the absorber
        averaged =[]

        blueMap = cm.Blues
        redMap = cm.Reds
        
        fig = figure()
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1

        rdif = []
        rcosInc = []
        rlyaW = []
        rMaj = []
        bdif = []
        bcosInc = []
        blyaW = []
        bMaj = []
        
        for d,i,w,m in zip(difList,cosIncList,lyaWList,majList):
            # check if all the values are okay
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(m):
                if d!=-99 and i!=-99 and w!=-99 and m!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            bdif.append(d)
                            bcosInc.append(i)
                            blyaW.append(w)
                            bMaj.append(m)
        #                     plotb = ax.scatter(a, w, cmap=blueMap, c=d, s=50, vmin=0, vmax=400)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
        #                     plotr = ax.scatter(a, w, cmap=redMap, c=d, s=50, vmin=0, vmax=400)
                            rdif.append(d)
                            rcosInc.append(i)
                            rlyaW.append(w)
                            rMaj.append(m)


        print
        print 'rdif: ',rdif
        print 'average, median redshifts: ',average(rdif),', ',median(rdif)
        print 'average, median blueshifts: ',average(bdif),', ',median(bdif)

        print 'max, min red dif: ',max(rdif), ', ',min(rdif)
        print 'max, min blue dif: ',max(bdif), ', ', min(bdif)
        
        plotr = ax.scatter(rcosInc, rlyaW, cmap=redMap, c=rdif, s=50, vmin=-400, vmax=0)
        norm = matplotlib.colors.Normalize(vmin = 0, vmax = 400)
        cbarRed = plt.colorbar(plotr,cmap=redMap,orientation='vertical')
        cbarRed.set_label('galaxy - absorber velocity (km/s)')
        
        plotb = ax.scatter(bcosInc, blyaW, cmap=blueMap, c=bdif, s=50, vmin=0, vmax=400)
        cbarBlue = plt.colorbar(plotb,cmap=blueMap,orientation='vertical')
        cbarBlue.set_label('galaxy - absorber velocity (km/s)')
#       cbar.ax.set_yticklabels(ticks)
                    
        title("Equivalent Width vs Cos(inclination)")
        xlabel(r'Cos(inclination) = b/a')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
#         legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(-0.02,1)

        if save:
            savefig('{0}/W(cos(inclination))_colorbar.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
        
        
        
#         fig = figure()
#         ax = fig.add_subplot(211)
#         hist(rdif,color='red')
#        
#         ax = fig.add_subplot(212)
#         hist(bdif,color='blue')
#         show()


##########################################################################################
##########################################################################################

    # plot equivalent width as a function of fancy_inclination for red and blue shifted
    # absorption
    #
    
    plotW_fancyInc = True
    save = True
    alpha = 0.75
    
    if plotW_fancyInc:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,i,w,m in zip(difList,fancyIncList,lyaWList,majList):
            # check if all the values are okay
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(m):
                if d!=-99 and i!=-99 and w!=-99 and m!=-99:
                    if i ==90:
                        print d, i, w, m
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
                
                    plot1 = scatter(i,w,c=color,s=50,alpha=alpha)
            
#         title('W(fancy_inclination) for red vs blue shifted absorption')
        xlabel(r'Inclination (deg)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(0,90)
        
        if save:
            savefig('{0}/W(fancy_inclination)_dif2.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
            
            
##########################################################################################
##########################################################################################          
            
    # plot equivalent width as a function of cos(fancy_inclination) for red and blue shifted
    # absorption
    #
    
    plotW_CosFancyInc = False
    save = False
    
    if plotW_CosFancyInc:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,i,w,m in zip(difList,cosFancyIncList,lyaWList,majList):
        
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(m):
                if d!=-99 and i!=-99 and w!=-99 and m!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(i,w,c='Blue',s=50,label= labelb)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(i,w,c='Red',s=50,label= labelr)
                
                    plot1 = scatter(i,w,c=color,s=50)
            
        title('W(cos(fancy_inclination)) for red vs blue shifted absorption')
        xlabel(r'Cos(fancy_inclination)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(-0.01,1)
        
        if save:
            savefig('{0}/W(cos(fancy_inclination))_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
            
            
##########################################################################################
##########################################################################################

    # plot equivalent width as a function of cos(fancy_inclination) with velocity dif
    # represented by a color bar
    #
    
    plotW_FancyCosInc_colorbar= False
    save = False
    
    if plotW_FancyCosInc_colorbar:
        # colormap the velocity difference of the absorber
        averaged =[]

        blueMap = cm.Blues
        redMap = cm.Reds
        
        fig = figure()
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1

        rdif = []
        rcosInc = []
        rlyaW = []
        rMaj = []
        bdif = []
        bcosInc = []
        blyaW = []
        bMaj = []
        
        for d,i,w,m in zip(difList,cosFancyIncList,lyaWList,majList):
            # check if all the values are okay
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(m):
                if d !=-99 and i !=-99 and w!=-99 and m!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            bdif.append(d)
                            bcosInc.append(i)
                            blyaW.append(w)
                            bMaj.append(m)
        #                     plotb = ax.scatter(a, w, cmap=blueMap, c=d, s=50, vmin=0, vmax=400)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
        #                     plotr = ax.scatter(a, w, cmap=redMap, c=d, s=50, vmin=0, vmax=400)
                            rdif.append(d)
                            rcosInc.append(i)
                            rlyaW.append(w)
                            rMaj.append(m)


        print
        print 'average, median redshifts: ',average(rdif),', ',median(rdif)
        print 'average, median blueshifts: ',average(bdif),', ',median(bdif)
        print 'max, min red dif: ',max(rdif), ', ',min(rdif)
        print 'max, min blue dif: ',max(bdif), ', ', min(bdif)
        
        plotr = ax.scatter(rcosInc, rlyaW, cmap=redMap, c=rdif, s=50, vmin=-400, vmax=0)
        norm = matplotlib.colors.Normalize(vmin = 0, vmax = 400)
        cbarRed = plt.colorbar(plotr,cmap=redMap,orientation='vertical')
        cbarRed.set_label('galaxy - absorber velocity (km/s)')
        
        plotb = ax.scatter(bcosInc, blyaW, cmap=blueMap, c=bdif, s=50, vmin=0, vmax=400)
        cbarBlue = plt.colorbar(plotb,cmap=blueMap,orientation='vertical')
        cbarBlue.set_label('galaxy - absorber velocity (km/s)')
#       cbar.ax.set_yticklabels(ticks)
                    
        title('Equivalent width vs Cos(fancy_inclination) colormap')
        xlabel(r'Cos(fancy_inclination) = b/a')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
#         legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(-0.02,1)

        if save:
            savefig('{0}/W(cos(fancy_inclination))_colorbar.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
        
        
#         fig = figure()
#         ax = fig.add_subplot(211)
#         hist(rdif,color='red')
# 
#         ax = fig.add_subplot(212)
#         hist(bdif,color='blue')
#         show()


##########################################################################################
##########################################################################################

    # plot equivalent width as a function of cos(inclination) for red and blue shifted
    # absorption
    #
    
    plotW_CosInc = False
    save = False
    
    if plotW_CosInc:
    
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,i,w,m in zip(difList,cosIncList,lyaWList,majList):
        
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(m):
                if d!=-99 and i!=-99 and w!=-99 and m!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(i,w,c='Blue',s=50,label= labelb)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(i,w,c='Red',s=50,label= labelr)
                
                    plot1 = scatter(i,w,c=color,s=50)
            
        title('W(cos(inclination)) for red vs blue shifted absorption')
        xlabel(r'Cos(inclination) = b/a')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(0,1)
        
        if save:
            savefig('{0}/W(cos(inclination))_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()



##########################################################################################
##########################################################################################

    # plot equivalent width as a function of cos(inclination) for red and blue shifted
    # absorption
    #
    
    plotW_inc_az = False
    save = False
    
    if plotW_inc_az:
    
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        labelrless = 'Red Shifted Absorber <45 az'
        labelbless = "Blue Shifted Absorber <45 az"
        
        for d,i,w,a in zip(difList,fancyIncList,lyaWList,azList):
        
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(a):
                if d!=-99 and i!=-99 and w!=-99 and a!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'

                        if a >=45:
                            plotb = ax.scatter(i,w,c='Blue',s=50,label= labelb)
                        else:
                            plotb = ax.scatter(i,w,c='Blue',s=50,label= labelbless,marker='*',lw=0)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'

                        if a >=45:
                            plotr = ax.scatter(i,w,c='Red',s=50,label= labelr)
                        else:
                            plotr = ax.scatter(i,w,c='Red',s=50,label= labelrless,marker='*',lw=0)

                
#                     plot1 = scatter(i,w,c=color,s=10)
            
#         title('W(cos(inclination)) for red vs blue shifted absorption')
        xlabel(r'Inclination')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
#         legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
#         ylim(-1,1200)
#         xlim(0,1)
        
        if save:
            savefig('{0}/W(inclination)_dif_az.pdf'.format(saveDirectory),format='pdf')
        else:
            show()

##########################################################################################
##########################################################################################



##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


if __name__=="__main__":
    # do the work
    main()
    