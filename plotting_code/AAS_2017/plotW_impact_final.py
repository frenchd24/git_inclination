#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotW_impact_final.py, v 5.7 12/09/16

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
    
    maxEnv = 1000
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
    
    # for ambiguous lines
    lyaVAmbList = []
    lyaWAmbList = []
    envAmbList = []
    ambAGNnameList = []
    
    # for all lines with a galaxy within 500 kpc
    lyaV_nearestList = []
    lyaW_nearestList = []
    env_nearestList = []
    impact_nearestList = []
    diam_nearestList = []
    
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
    print

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
            
        # for all absorbers with a galaxy within 500kpc
        if isNumber(impact):
            lyaV_nearestList.append(float(lyaV))
            lyaW_nearestList.append(float(lyaW))
            env_nearestList.append(float(env))
            impact_nearestList.append(float(impact))
            diam_nearestList.append(float(maj))
            nameList.append(galaxyName)

            
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

    total = 0
    totalNo = 0
    totalYes = 0
    totalIsolated = 0
    totalGroup = 0
    
    

##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter, split between
    # red and blue shifted absorption, overplot median histograms for red and blue
    #
    # include standard error in the mean shaded regions for each bar
    #
    
    plotW_impact_hist_errors = False
    save = False
    
    if plotW_impact_hist_errors:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        alpha = 0.7
        errorAlpha=0.25
        
        binSize = 125
        bins = arange(0,625,binSize)
        
        labelr = r'$\rm Redshifted ~Absorber$'
        labelb = r'$\rm Blueshifted ~Absorber$'
        bSymbol = 'D'
        rSymbol = 'o'
        
        xVals = []
        redX = []
        redY = []
        blueX = []
        blueY = []
        
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
                        symbol = bSymbol
                        
                        blueX.append(xVal)
                        blueY.append(yVal)
                        
                        if countb == 0:
                            countb +=1
#                             plotb = ax.scatter(v,w,c='Blue',s=50,label= labelb)
                            plotb = ax.scatter(xVal,yVal,marker=symbol,c='Blue',s=50,\
                            alpha=alpha,label=labelb)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        symbol = rSymbol
                        
                        redX.append(xVal)
                        redY.append(yVal)
                        
                        if countr == 0:
                            countr +=1
#                             plotr = ax.scatter(v,w,c='Red',s=50,label= labelr)
                            plotr = ax.scatter(xVal,yVal,marker=symbol,c='Red',s=50,\
                            alpha=alpha,label=labelr)

                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=50,alpha=alpha)

        
        # avg red
        bin_means,edges,binNumber = stats.binned_statistic(array(redX), array(redY), \
        statistic='mean', bins=bins)
        
        bin_errors,edges_e,binNumber_e = stats.binned_statistic(array(redX), array(redY), \
        statistic=lambda y: errors(y), bins=bins)
        
        bin_std,edges_std,binNumber_std = stats.binned_statistic(array(redX), array(redY), \
        statistic=lambda y: std(y), bins=bins)
        
        print 'bin_means,edges,binNumber : ',bin_means,edges,binNumber
        print
        print 'bin_errors, edges_e,binNumber_e : ',bin_errors,edges_e,binNumber_e
        print
        print 'bin_std,edges_std,binNumber_std : ',bin_std,edges_std,binNumber_std
        
        # the mean
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        
        # the errors
        left_e,right_e = edges_e[:-1],edges_e[1:]        
        X_e = array([left_e,right_e]).T.flatten()
        Y_e = array([nan_to_num(bin_errors),nan_to_num(bin_errors)]).T.flatten()
        
        yErrorsTop = Y + Y_e
        yErrorsBot = Y - Y_e

        plot(X_e,yErrorsBot, ls='solid',color='red',lw=1,alpha=errorAlpha)
        plot(X_e,yErrorsTop, ls='solid',color='red',lw=1,alpha=errorAlpha)
        fill_between(X_e, yErrorsBot, yErrorsTop, facecolor='red', interpolate=True,alpha=errorAlpha)

        plot(X,Y, ls='dotted',color='red',lw=2.1,alpha=alpha+0.2,label=r'$\rm Mean~ Redshifted ~EW$')
    
    
    
    
        # avg blue
        bin_means,edges,binNumber = stats.binned_statistic(array(blueX), array(blueY), \
        statistic='mean', bins=bins)
        
        bin_errors,edges_e,binNumber_e = stats.binned_statistic(array(blueX), array(blueY), \
        statistic=lambda y: errors(y), bins=bins)
        
        bin_std,edges_std,binNumber_std = stats.binned_statistic(array(blueX), array(blueY), \
        statistic=lambda y: std(y), bins=bins)
        
        print
        print 'bin_means,edges,binNumber : ',bin_means,edges,binNumber
        print
        print 'bin_errors, edges_e,binNumber_e : ',bin_errors,edges_e,binNumber_e
        print
        print 'bin_std,edges_std,binNumber_std : ',bin_std,edges_std,binNumber_std
        
        # the mean
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        
        # the errors
        left_e,right_e = edges_e[:-1],edges_e[1:]        
        X_e = array([left_e,right_e]).T.flatten()
        Y_e = array([nan_to_num(bin_errors),nan_to_num(bin_errors)]).T.flatten()
        
        yErrorsTop = Y + Y_e
        yErrorsBot = Y - Y_e

        plot(X_e,yErrorsBot, ls='solid',color='blue',lw=1,alpha=errorAlpha)
        plot(X_e,yErrorsTop, ls='solid',color='blue',lw=1,alpha=errorAlpha)
        fill_between(X_e, yErrorsBot, yErrorsTop, facecolor='blue', interpolate=True,alpha=errorAlpha)
        
        plot(X,Y, ls='dashed',color='blue',lw=1.7,alpha=alpha+0.1,label=r'$\rm Mean~ Blueshifted ~EW$')    
        
        
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
        leg = ax.legend(scatterpoints=1,prop={'size':13},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0,500)

        if save:
            savefig('{0}/W(impact)_mean_{1}_difHistograms2.pdf'.format(saveDirectory,binSize),format='pdf',bbox_inches='tight')
        else:
            show()


##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter/R_vir, split between
    # red and blue shifted absorption, overplot median histograms for red and blue
    #
    # plot shaded error regions around histogram
    # 
    
    plotW_impact_vir_hist_errors = False
    save = False
    
    if plotW_impact_vir_hist_errors:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        alpha = 0.7
        
        binSize = 0.5
        bins = arange(0,2.5,binSize)
        
        labelr = r'$\rm Redshifted ~Absorber$'
        labelb = r'$\rm Blueshifted ~Absorber$'
        bSymbol = 'D'
        rSymbol = 'o'
        
        
        xVals = []
        redX = []
        redY = []
        blueX = []
        blueY = []
        
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
                        symbol = bSymbol
                        
                        blueX.append(xVal)
                        blueY.append(yVal)
                        
                        if countb == 0:
                            countb +=1
#                             plotb = ax.scatter(v,w,c='Blue',s=50,label= labelb)
                            plotb = ax.scatter(xVal,yVal,marker=symbol,c='Blue',s=50,\
                            alpha=alpha,label=labelb)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        symbol = rSymbol
                        
                        redX.append(xVal)
                        redY.append(yVal)
                        
                        if countr == 0:
                            countr +=1
#                             plotr = ax.scatter(v,w,c='Red',s=50,label= labelr)
                            plotr = ax.scatter(xVal,yVal,marker=symbol,c='Red',s=50,\
                            alpha=alpha,label=labelr)

                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=50,alpha=alpha)
        
        
        # avg red
        bin_means,edges,binNumber = stats.binned_statistic(array(redX), array(redY), \
        statistic='mean', bins=bins)
        
        bin_errors,edges_e,binNumber_e = stats.binned_statistic(array(redX), array(redY), \
        statistic=lambda y: errors(y), bins=bins)
        
        bin_std,edges_std,binNumber_std = stats.binned_statistic(array(redX), array(redY), \
        statistic=lambda y: std(y), bins=bins)
        
        print 'bin_means,edges,binNumber : ',bin_means,edges,binNumber
        print
        print 'bin_errors, edges_e,binNumber_e : ',bin_errors,edges_e,binNumber_e
        print
        print 'bin_std,edges_std,binNumber_std : ',bin_std,edges_std,binNumber_std
        
        # the mean
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        
        # the errors
        left_e,right_e = edges_e[:-1],edges_e[1:]        
        X_e = array([left_e,right_e]).T.flatten()
        Y_e = array([nan_to_num(bin_errors),nan_to_num(bin_errors)]).T.flatten()
        
        yErrorsTop = Y + Y_e
        yErrorsBot = Y - Y_e

        plot(X_e,yErrorsBot, ls='solid',color='red',lw=1,alpha=errorAlpha)
        plot(X_e,yErrorsTop, ls='solid',color='red',lw=1,alpha=errorAlpha)
        fill_between(X_e, yErrorsBot, yErrorsTop, facecolor='red', interpolate=True,alpha=errorAlpha)

        plot(X,Y, ls='dotted',color='red',lw=2.1,alpha=alpha+0.2,label=r'$\rm Mean~ Redshifted ~EW$')
    
    
    
        # avg blue
        bin_means,edges,binNumber = stats.binned_statistic(array(blueX), array(blueY), \
        statistic='mean', bins=bins)
        
        bin_errors,edges_e,binNumber_e = stats.binned_statistic(array(blueX), array(blueY), \
        statistic=lambda y: errors(y), bins=bins)
        
        bin_std,edges_std,binNumber_std = stats.binned_statistic(array(blueX), array(blueY), \
        statistic=lambda y: std(y), bins=bins)
        
        print
        print 'bin_means,edges,binNumber : ',bin_means,edges,binNumber
        print
        print 'bin_errors, edges_e,binNumber_e : ',bin_errors,edges_e,binNumber_e
        print
        print 'bin_std,edges_std,binNumber_std : ',bin_std,edges_std,binNumber_std
        
        # the mean
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        
        # the errors
        left_e,right_e = edges_e[:-1],edges_e[1:]        
        X_e = array([left_e,right_e]).T.flatten()
        Y_e = array([nan_to_num(bin_errors),nan_to_num(bin_errors)]).T.flatten()
        
        yErrorsTop = Y + Y_e
        yErrorsBot = Y - Y_e

        plot(X_e,yErrorsBot, ls='solid',color='blue',lw=1,alpha=errorAlpha)
        plot(X_e,yErrorsTop, ls='solid',color='blue',lw=1,alpha=errorAlpha)
        fill_between(X_e, yErrorsBot, yErrorsTop, facecolor='blue', interpolate=True,alpha=errorAlpha)
        
        plot(X,Y, ls='dashed',color='blue',lw=1.7,alpha=alpha+0.1,label=r'$\rm Mean~ Blueshifted ~EW$')
        
        
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
        leg = ax.legend(scatterpoints=1,prop={'size':13},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0,2.0)

        if save:
            savefig('{0}/W(impact_vir)_mean_{1}_difHistograms2.pdf'.format(saveDirectory,binSize),format='pdf',bbox_inches='tight')
        else:
            show()


##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter/R_vir, split between
    # red and blue shifted absorption, overplot single median histogram (total EW)
    #
    
    plotW_impact_vir_medHist = False
    save = False
    
    if plotW_impact_vir_medHist:
        fig = figure()
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        alpha = 0.7
        
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        bSymbol = 'D'
        rSymbol = 'o'
        
        xVals = []
        yVals = []
        
        for d,i,w,v in zip(difList,impactList,lyaWList,virList):
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(v):
                if d!=-99 and i!=-99 and w!=-99 and v!=-99:
                    xVal = float(i)/float(v)
                    yVal = float(w)
                    
                    xVals.append(xVal)
                    yVals.append(yVal)
                    
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        symbol = bSymbol

                        
                        if countb == 0:
                            countb +=1
#                             plotb = ax.scatter(v,w,c='Blue',s=50,label= labelb)
                            plotb = ax.scatter(xVal,yVal,marker=bSymbol,c='Blue',s=50,alpha=alpha)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        symbol = rSymbol
                        
                        if countr == 0:
                            countr +=1
#                             plotr = ax.scatter(v,w,c='Red',s=50,label= labelr)
                            plotr = ax.scatter(xVal,yVal,marker=rSymbol,c='Red',s=50,alpha=alpha)

                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=50,alpha=alpha)
        
        
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
    
        binSize = 0.5
        bins = arange(0,2.5,binSize)

        # 50% percentile
        bin_means,edges,binNumber = stats.binned_statistic(array(xVals), array(yVals), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='solid',color='black',lw=1.7,alpha=alpha+0.1,label=r'$\rm Mean ~EW$')
        
#         # 90% percentile
#         bin_means,edges,binNumber = stats.binned_statistic(array(xVals), array(yVals), \
#         statistic=lambda y: perc90(y), bins=bins)
#         left,right = edges[:-1],edges[1:]
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plot(X,Y, ls='dashed',color='dimgray',lw=1.7,alpha=alpha+0.1,label=r'$\rm 90th\% ~EW$')
        
#         # 10th percentile
#         bin_means,edges,binNumber = stats.binned_statistic(array(xVals), array(yVals), \
#         statistic=lambda y: perc10(y), bins=bins)
#         left,right = edges[:-1],edges[1:]        
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plot(X,Y, ls='dotted',color='green',lw=1.7,alpha=alpha+0.1,label=r'$\rm 10th\% ~EW$')
        
        
        xlabel(r'$\rm \rho / R_{vir}$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ax.legend(scatterpoints=1,prop={'size':15},loc=1,fancybox=True)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1200)
        xlim(0,2.0)

        if save:
            savefig('{0}/W(impact_vir)_median_{1}_Histograms_maxLyaW{2}.pdf'.format(saveDirectory,binSize,maxLyaW),format='pdf',bbox_inches='tight')
        else:
            show()
            

##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter, split between
    # red and blue shifted absorption, overplot single median histogram (total EW)
    #
    
    plotW_impact_medHist = False
    save = False
    
    if plotW_impact_medHist:
        fig = figure()
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        alpha = 0.7
        
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        bSymbol = 'D'
        rSymbol = 'o'
        
        xVals = []
        yVals = []
        
        for d,i,w,v in zip(difList,impactList,lyaWList,virList):
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(v):
                if d!=-99 and i!=-99 and w!=-99 and v!=-99:
                    xVal = float(i)
                    yVal = float(w)
                    
                    xVals.append(xVal)
                    yVals.append(yVal)
                    
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        symbol = bSymbol
                        
                        if countb == 0:
                            countb +=1
#                             plotb = ax.scatter(v,w,c='Blue',s=50,label= labelb)
                            plotb = ax.scatter(xVal,yVal,marker=bSymbol,c='Blue',s=50,alpha=alpha)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        symbol = rSymbol
                        
                        if countr == 0:
                            countr +=1
#                             plotr = ax.scatter(v,w,c='Red',s=50,label= labelr)
                            plotr = ax.scatter(xVal,yVal,marker=rSymbol,c='Red',s=50,alpha=alpha)

                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=50,alpha=alpha)
        
        
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
    
        binSize = 100
        bins = arange(0,600,binSize)
        
        # 50% percentile
        bin_means,edges,binNumber = stats.binned_statistic(array(xVals), array(yVals), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='solid',color='black',lw=1.5,alpha=alpha,label=r'$\rm Mean ~EW$')
        
#         # 90% percentile
#         bin_means,edges,binNumber = stats.binned_statistic(array(xVals), array(yVals), \
#         statistic=lambda y: perc90(y), bins=bins)
#         left,right = edges[:-1],edges[1:]        
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plot(X,Y, ls='dashed',color='dimgray',lw=1.5,alpha=alpha,label=r'$\rm 90th\% ~EW$')
        
#         # 10th percentile
#         bin_means,edges,binNumber = stats.binned_statistic(array(xVals), array(yVals), \
#         statistic=lambda y: perc10(y), bins=bins)
#         left,right = edges[:-1],edges[1:]        
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plot(X,Y, ls='dotted',color='green',lw=1.5,alpha=alpha,label=r'$\rm 10th\% ~EW$')
        
        xlabel(r'$\rm \rho ~[kpc]$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ax.legend(scatterpoints=1,prop={'size':15},loc=1,fancybox=True)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1200)
        xlim(0,500)

        if save:
            savefig('{0}/W(impact)_mean_{1}_Histograms_maxLyaW{2}.pdf'.format(saveDirectory,binSize,maxLyaW),format='pdf',bbox_inches='tight')
        else:
            show()




##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter, split between
    # absorbers at impact < 1 R_vir, overplot median histograms for each
    #
    
    plotW_impact_virseparate = False
    save = False
    
    if plotW_impact_virseparate:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        alpha = 0.7
        alphaInside = 0.7
        binSize = 125
        bins = arange(0,625,binSize)
        markerSize = 60
        
        labelr = r'$\rm Redshifted ~Absorber$'
        labelb = r'$\rm Blueshifted ~Absorber$'
        bSymbol = 'D'
        rSymbol = 'o'
        
        xVals = []
        redX = []
        redY = []
        blueX = []
        blueY = []
        
        for d,i,w,v,inc in zip(difList,impactList,lyaWList,virList,incList):
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(v) and isNumber(inc):
                if d!=-99 and i!=-99 and w!=-99 and v!=-99 and inc!=-99:
                    xVal = float(i)
                    yVal = float(w)
                    
                    xVals.append(xVal)
                        
                    if d>0:
                        # blueshifted
                        color = 'Blue'
                        symbol = bSymbol

                        blueX.append(xVal)
                        blueY.append(yVal)
                        
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
                    
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(xVal,yVal,marker=symbol,c='Blue',\
                            facecolor=fc,edgecolor=ec,s=markerSize,alpha=a,label=labelb)

                    if d<0:
                        # redshifted
                        color = 'Red'
                        symbol = rSymbol
                        
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
                    
                        redX.append(xVal)
                        redY.append(yVal)
                    
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(xVal,yVal,marker=symbol,c='Red',\
                            s=markerSize,facecolor=fc,edgecolor=ec,alpha=a,label=labelr)

                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=markerSize,\
                    facecolor=fc,edgecolor=ec,alpha=a)
     
        
        # redshifted
        bin_means,edges,binNumber = stats.binned_statistic(array(redX), array(redY), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dotted',color='red',lw=1.7,alpha=alpha+0.1,label=r'$\rm Mean~ Redshifted ~EW$')
    
        # blueshifted
        bin_means,edges,binNumber = stats.binned_statistic(array(blueX), array(blueY), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dashed',color='blue',lw=1.7,alpha=alpha+0.1,label=r'$\rm Mean~ Blueshifted ~EW$')
        
        
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
            savefig('{0}/W(impact)_mean_{1}_virsep_maxEnv{2}_maxLyaW{3}.pdf'.format(saveDirectory,binSize,maxEnv,maxLyaW),format='pdf',bbox_inches='tight')
        else:
            show()
            
            
            
##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter/R_vir, split between
    # red and blue shifted absorption, overplot median histograms for red and blue
    #
    # plot shaded error regions around histogram
    # 
    
    plotW_impact_vir_hist_errors = False
    save = False
    
    if plotW_impact_vir_hist_errors:
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
        bSymbol = 'D'
        rSymbol = 'o'
        
        
        xVals = []
        redX = []
        redY = []
        blueX = []
        blueY = []
        
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
                        symbol = bSymbol
                        
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
                        
                        
                        blueX.append(xVal)
                        blueY.append(yVal)
                        
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(xVal,yVal,marker=symbol,c='Blue',\
                            facecolor=fc,edgecolor=ec,s=markerSize,alpha=a,label=labelb)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        symbol = rSymbol
                        
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
                        
                        redX.append(xVal)
                        redY.append(yVal)
                        
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(xVal,yVal,marker=symbol,c='Red',\
                            s=markerSize,facecolor=fc,edgecolor=ec,alpha=a,label=labelr)

                    plot1 = scatter(xVal,yVal,marker=symbol,c=color,s=markerSize,\
                    facecolor=fc,edgecolor=ec,alpha=a)
        
        
        # avg red
        bin_means,edges,binNumber = stats.binned_statistic(array(redX), array(redY), \
        statistic='mean', bins=bins)
        
#         bin_errors,edges_e,binNumber_e = stats.binned_statistic(array(redX), array(redY), \
#         statistic=lambda y: errors(y), bins=bins)
#         
#         bin_std,edges_std,binNumber_std = stats.binned_statistic(array(redX), array(redY), \
#         statistic=lambda y: std(y), bins=bins)
#         
#         print 'bin_means,edges,binNumber : ',bin_means,edges,binNumber
#         print
#         print 'bin_errors, edges_e,binNumber_e : ',bin_errors,edges_e,binNumber_e
#         print
#         print 'bin_std,edges_std,binNumber_std : ',bin_std,edges_std,binNumber_std
        
        # the mean
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        
        # the errors
#         left_e,right_e = edges_e[:-1],edges_e[1:]        
#         X_e = array([left_e,right_e]).T.flatten()
#         Y_e = array([nan_to_num(bin_errors),nan_to_num(bin_errors)]).T.flatten()
#         
#         yErrorsTop = Y + Y_e
#         yErrorsBot = Y - Y_e
# 
#         plot(X_e,yErrorsBot, ls='solid',color='red',lw=1,alpha=errorAlpha)
#         plot(X_e,yErrorsTop, ls='solid',color='red',lw=1,alpha=errorAlpha)
#         fill_between(X_e, yErrorsBot, yErrorsTop, facecolor='red', interpolate=True,alpha=errorAlpha)

        plot(X,Y, ls='dotted',color='red',lw=2.1,alpha=alpha+0.2,label=r'$\rm Mean~ Redshifted ~EW$')
    
    
    
        # avg blue
        bin_means,edges,binNumber = stats.binned_statistic(array(blueX), array(blueY), \
        statistic='mean', bins=bins)
        
#         bin_errors,edges_e,binNumber_e = stats.binned_statistic(array(blueX), array(blueY), \
#         statistic=lambda y: errors(y), bins=bins)
#         
#         bin_std,edges_std,binNumber_std = stats.binned_statistic(array(blueX), array(blueY), \
#         statistic=lambda y: std(y), bins=bins)
#         
#         print
#         print 'bin_means,edges,binNumber : ',bin_means,edges,binNumber
#         print
#         print 'bin_errors, edges_e,binNumber_e : ',bin_errors,edges_e,binNumber_e
#         print
#         print 'bin_std,edges_std,binNumber_std : ',bin_std,edges_std,binNumber_std
        
        # the mean
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        
        # the errors
#         left_e,right_e = edges_e[:-1],edges_e[1:]        
#         X_e = array([left_e,right_e]).T.flatten()
#         Y_e = array([nan_to_num(bin_errors),nan_to_num(bin_errors)]).T.flatten()
#         
#         yErrorsTop = Y + Y_e
#         yErrorsBot = Y - Y_e

#         plot(X_e,yErrorsBot, ls='solid',color='blue',lw=1,alpha=errorAlpha)
#         plot(X_e,yErrorsTop, ls='solid',color='blue',lw=1,alpha=errorAlpha)
#         fill_between(X_e, yErrorsBot, yErrorsTop, facecolor='blue', interpolate=True,alpha=errorAlpha)
        
        plot(X,Y, ls='dashed',color='blue',lw=1.7,alpha=alpha+0.1,label=r'$\rm Mean~ Blueshifted ~EW$')
        
        
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
            savefig('{0}/W(impact_vir)_mean_{1}_virsep_maxEnv{2}_maxLyaW{3}.pdf'.format(saveDirectory,binSize,maxEnv,maxLyaW),format='pdf',bbox_inches='tight')
        else:
            show()


##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter, split between
    # absorbers at impact < 1 R_vir, overplot median histograms for each
    #
    # FOR INCLUDE = 1 ABSORBERS

    
    plotW_impact_virseparate_include1 = True
    save = True
    
    if plotW_impact_virseparate_include1:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        alpha = 0.6
        alphaInside = 0.6
        binSize = 125
        bins = arange(0,625,binSize)
        markerSize = 60
        color = 'Blue'
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
        
        for i,w,v in zip(impactList,lyaWList,virList):
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
        statistic='median', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dashed',color='black',lw=1.7,alpha=alpha+0.2,label=r'$\rm Median ~EW$')
        
        
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
            savefig('{0}/W(impact)_median_binSize{1}_virsep_include1_maxEnv{2}.pdf'.format(saveDirectory,binSize,maxEnv),format='pdf',bbox_inches='tight')
        else:
            show()
            
            
            
##########################################################################################
##########################################################################################
    # plot equivalent width as a function of impact parameter/R_vir
    # 
    # FOR INCLUDE = 1 ABSORBERS


    plotW_impact_vir_median_include1 = True
    save = True
    
    if plotW_impact_vir_median_include1:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1

        alpha = 0.6
        alphaInside = 0.6
        markerSize = 60
        
        binSize = 0.5
        bins = arange(0,2.5,binSize)
        
        labelr = r'$\rm Redshifted ~Absorber$'
        labelb = r'$\rm Blueshifted ~Absorber$'
        color = 'Blue'
        symbol = 'o'
        maxW = 1500.


                    
        xVals = []
        yVals = []
        redX = []
        redY = []
        blueX = []
        blueY = []
        
        for i,w,v in zip(impactList,lyaWList,virList):
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
        statistic='median', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,Y, ls='dashed',color='black',lw=1.7,alpha=alpha+0.2,label=r'$\rm Median ~EW$')
    
        
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
            savefig('{0}/W(impact_vir)_median_binSize{1}_virsep_include1_maxEnv{2}.pdf'.format(saveDirectory,binSize,maxEnv),format='pdf',bbox_inches='tight')
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
    