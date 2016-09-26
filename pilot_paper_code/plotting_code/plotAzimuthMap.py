#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotAzimuthMap.py, v 1.0 09/19/2016

Plot absorption properties on a map showing their distribution around a central galaxy


'''

import sys
import os
import csv
from scipy import stats

from pylab import *
# import atpy
from math import *
from utilities import *
import getpass
import pickle
from matplotlib.patches import Ellipse


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

    
def main():
    # assuming 'theFile' contains one name per line, read the file
    
    
    if getpass.getuser() == 'David':
        pickleFilename = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/pilotData2.p'
#         resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_8_edit2.csv'
#         saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots2/'
        resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit4.csv'
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots4/'
        WS09data = '/Users/David/Research_Documents/inclination/git_inclination/WS2009_lya_data.tsv'

    elif getpass.getuser() == 'frenchd':
        pickleFilename = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/pilotData2.p'
#         resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_8_edit2.csv'
#         saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots2/'
        resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit4.csv'
        saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots4/'
        WS09data = '/usr/users/frenchd/inclination/git_inclination/WS2009_lya_data.tsv'

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
    
    WS = open(WS09data,'rU')
    WSreader = csv.DictReader(WS,delimiter=';')
    
    virInclude = False
    cusInclude = False
    finalInclude = True
    
    maxEnv = 300
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
    WSvcorr = []
    WSdiam = []
    WSimpact =[]
    WSew = []
    WSvel = []
    WSlya = []
    WSvel_dif = []
    WSvir = []
    WSlike = []
    
    l_min = 0.001

    for w in WSreader:
        vcorr = w['HV']
        diam = w['Diam']
        rho = w['rho']
        ew = w['EWLya']
        vel = w['LyaVel']
        lya = w['Lya']
        
        if lya == 'Lya  ' and isNumber(diam) and isNumber(ew) and isNumber(rho):
            if float(rho) <=500.0:
                # this is a single galaxy association
                vir = calculateVirialRadius(float(diam))
                
                vel_dif = float(vcorr) - float(vel)
    
                # try this "sphere of influence" value instead
                m15 = float(diam)**1.5

                # first for the virial radius
                likelihood = math.exp(-(float(rho)/vir)**2) * math.exp(-(vel_dif/200.)**2)
                
                if vir>= float(rho):
                    likelihood = likelihood*2
                    
                # then for the second 'virial like' m15 radius
                likelihoodm15 = math.exp(-(float(rho)/m15)**2) * math.exp(-(vel_dif/200.)**2)
                
                if m15>= float(rho):
                    likelihoodm15 = likelihoodm15*2
                    
                if likelihood <= likelihoodm15:
                    likelihood = likelihoodm15
                    
                WSlike.append(likelihood)
                
#                 l_min=0
                
                if likelihood >= l_min:
                
                    WSvcorr.append(float(vcorr))
                    WSdiam.append(float(diam))
                    WSvir.append(vir)
                    WSimpact.append(float(rho))
                    WSew.append(float(ew))
                    WSvel.append(float(vel))
                    WSlya.append(lya)
                    WSvel_dif.append(vel_dif)
    
    
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
                raList.append(galaxyRA_Dec[0])
                decList.append(galaxyRA_Dec[1])
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

    results.close()
    WS.close()
            
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

    # plot histograms of the inclinations for both associated galaxies and the 
    # full galaxy data set, combining both redshifted and blueshifted
    plot_edge_on = True
    
    if plot_edge_on:
        
        azBlue_1 = []
        azRed_1 = []
        azBlue_2 = []
        azRed_2 = []
        azBlue_3 = []
        azRed_3 = []
        
        # new coordinates of absorption for plotting (w/ galaxy at (0,0))
        yList_1 = []
        xList_1 = []
        yList_2 = []
        xList_2 = []
        yList_3 = []
        xList_3 = []
        
        # new coordinates normalized for diameter
        y_dList = []
        x_dList = []
        
        # new coordinates normalized for virial radius
        y_vList = []
        x_vList = []
        
        # calculate the position on the sky for each absorption feature wrt to the galaxy
        for r,d,i,a,fInc,maj,vir,dif,inc in zip(raList,decList,impactList,azList,fancyIncList,majList,virList,difList,incList):
            if float(a) >= 0 and float(fInc)>=0 and float(fInc)<=90:
                
                if inc <=40:
                    if dif >0:
                        #blue absorber
                        azBlue_1.append(a)
                    else:
                        azRed_1.append(a)
                    
                    # y coordinate
                    y = (float(i)/float(vir)) * sin((a*pi)/180.)
                    yList_1.append(y)
                
                    # x coordinate
                    x = (float(i)/float(vir)) * cos(a*pi/180.)
                    xList_1.append(x)

                if inc > 40 and inc <=65:
                    if dif >0:
                        #blue absorber
                        azBlue_2.append(a)
                    else:
                        azRed_2.append(a)
                    
                    # y coordinate
                    y = (float(i)/float(vir)) * sin((a*pi)/180.)
                    yList_2.append(y)
                
                    # x coordinate
                    x = (float(i)/float(vir)) * cos(a*pi/180.)
                    xList_2.append(x)
                    
                if inc > 65:
                    if dif >0:
                        #blue absorber
                        azBlue_3.append(a)
                    else:
                        azRed_3.append(a)
                    
                    # y coordinate
                    y = (float(i)/float(vir)) * sin((a*pi)/180.)
                    yList_3.append(y)
                
                    # x coordinate
                    x = (float(i)/float(vir)) * cos(a*pi/180.)
                    xList_3.append(x)
                                        
            else:
                print 'float(a) <0: ',r,d,i,a,fInc

        # calculate the average red vs blue azimuth line
        blueAvg1 = mean(azBlue_1)
        print 'blueAvg1: ',blueAvg1
        redAvg1 = mean(azRed_1)
        print 'redAvg1: ',redAvg1
        print
        blueAvg2 = mean(azBlue_2)
        print 'blueAvg2: ',blueAvg2
        redAvg2 = mean(azRed_2)
        print 'redAvg2: ',redAvg2
        print
        blueAvg3 = mean(azBlue_3)
        print 'blueAvg3: ',blueAvg3
        redAvg3 = mean(azRed_3)
        print 'redAvg3: ',redAvg3
        print
        
        xyBlueAvg1 = (500.* cos(blueAvg1 * pi/180.), 500.* sin(blueAvg1 * pi/180.))
        xyRedAvg1 = (500.* cos(redAvg1 * pi/180.), 500.* sin(redAvg1 * pi/180.))
        print 'xyBlueAvg1: ',xyBlueAvg1
        print 'xyRedAvg1: ',xyRedAvg1
        xyBlueAvg2 = (500.* cos(blueAvg2 * pi/180.), 500.* sin(blueAvg2 * pi/180.))
        xyRedAvg2 = (500.* cos(redAvg2 * pi/180.), 500.* sin(redAvg2 * pi/180.))
        print 'xyBlueAvg2: ',xyBlueAvg2
        print 'xyRedAvg2: ',xyRedAvg2
        xyBlueAvg3 = (500.* cos(blueAvg3 * pi/180.), 500.* sin(blueAvg3 * pi/180.))
        xyRedAvg3 = (500.* cos(redAvg3 * pi/180.), 500.* sin(redAvg3 * pi/180.))
        print 'xyBlueAvg3: ',xyBlueAvg3
        print 'xyRedAvg3: ',xyRedAvg3

        
        # plot the distributions 
        fig = figure(figsize=(4,10))
        subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)
        
        ax1 = fig.add_subplot(311)
        # plot the flat galaxy line
        plot1 = plot((-0.5,0.5),(0,0),c='black',linewidth=5)
        
        maxW = 120
        minW = 5
        minLya = min(lyaWList)
        maxLya = max(lyaWList)
        
        plot2 = scatter(xList_1,yList_1)
        xlim(0,2.5)
        ylim(0,2.5)

        
        ax2 = fig.add_subplot(312)
        # plot the flat galaxy line
        e = Ellipse(xy=(0,0), width=0.5, height=0.2, angle=0)
                        
        # no transparency
        e.set_alpha(0.3)
        ax2.add_artist(e)
        e.set_facecolor('black')
        e.set_edgecolor('black')
        
        plot3 = scatter(xList_2,yList_2)
        xlim(0,2.5)
        ylim(0,2.5)

        
        ax3 = fig.add_subplot(313)
        # plot the flat galaxy line
        e = Ellipse(xy=(0,0), width=0.5, height=0.4, angle=0)
                        
        # no transparency
        e.set_alpha(0.3)
        ax3.add_artist(e)
        e.set_facecolor('black')
        e.set_edgecolor('black')
        
        plot4 = scatter(xList_3,yList_3)
        xlim(0,2.5)
        ylim(0,2.5)


#             newW = ((float(w) - min(lyaWList)) / (max(lyaWList) - min(lyaWList)))*(maxW - minW)
        
        # plot the absorption features
#         for x,y,s,d in zip(xList_1,yList_1):
#             newW = ((float(w) - min(lyaWList)) / (max(lyaWList) - min(lyaWList)))*(maxW - minW)
#             if d>0:
#                 # blueshifted
#                 plot2 = scatter(x,y,marker='*',color='blue',s=s)
#             else:
#                 # redshifted
#                 plot3 = scatter(x,y,marker='*',color='red',s=s)


        show()



    
###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    