#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotW_Inc2_final.py, v 5.6 01/05/18

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
    
v5.1: updated for LG_correlation_combined5_8_edit2.csv and l_min = 0.001
    (02/18/2016)
    
v5.2: remake plots with v_hel instead of vcorr (4/21/16)

v5.3: remake plots with new large galaxy sample (7/13/16) -> /plots4/

v5.4: remake plots after adjusting some targets (9/20/16)
        -> LG_correlation_combined5_11_25cut_edit2.csv now ->/plots5/
        
v5.5: remake after first referee report (12/12/16)
    - change facecolor for impact parameter < /> R_vir systems
    
v5.6: updated for AAS_2017 (01/05/18)
'''

import sys
import os
import csv

from pylab import *
from scipy import stats
from math import *
from utilities import *
import getpass
import pickle

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

def perc50(a):
    if len(a)>0:
        return percentile(a,50)
    else:
        return 0


    
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
##########################################################################################
##########################################################################################
    # plot equivalent width as a function of fancy_inclination for red and blue shifted
    # absorption, include an overlying 50th and 90th percentile histograms
    #
    # Open and closed symbols for impact > /< R_vir systems
    
    plotW_inc_percentile = True
    save = True
    
    if plotW_inc_percentile:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        binSize = 10
        numBins = 9
        alpha = 0.7
        markerSize = 60

        binSize = 12
        bins = arange(0,90+binSize,binSize)

        labelr = r'$\rm Redshifted$'
        labelb = r'$\rm Blueshifted$'
        bSymbol = 'D'
        rSymbol = 'o'
        
        allInc = []
        allW = []
        allVir = []
        
        rInc = []
        bInc = []
        rW = []
        bW = []
        rVir = []
        bVir = []
        
        for d,i,w,v,imp in zip(difList,fancyIncList,lyaWList,virList,impactList):
            # check if all the values are good
            if isNumber(d) and isNumber(i) and isNumber(w) and isNumber(v):
                if d!=-99 and i!=-99 and w!=-99 and v!=-99:
                    
                    allInc.append(float(i))
                    allW.append(float(w))
                    allVir.append(float(v))
                        
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        symbol = bSymbol
                        
                        if float(imp) > float(v):
                            # impact parameter > R_vir
                            fc = color
                            ec = 'black'
                            
                        if float(imp) <= float(v):
                            # impact parameter <= R_vir
                            fc = 'none'
                            ec = color
                        
                        bInc.append(float(i))
                        bW.append(float(w))
                        bVir.append(float(v))
                                            
                        if countb == 0:
                            countb +=1
#                             plotb = ax.scatter(v,w,c='Blue',s=50,label= labelb)
                            plotb = ax.scatter(i,w,marker=bSymbol,c='Blue',s=markerSize,\
                            facecolor=fc,edgecolor=ec,alpha=alpha)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        symbol = rSymbol
                        
                        if float(imp) > float(v):
                            # impact parameter > R_vir
                            fc = color
                            ec = 'black'
                            
                        if float(imp) <= float(v):
                            # impact parameter <= R_vir
                            fc = 'none'
                            ec = color
                        
                        rInc.append(float(i))
                        rW.append(float(w))
                        rVir.append(float(v))
                        
                        if countr == 0:
                            countr +=1
#                             plotr = ax.scatter(v,w,c='Red',s=50,label= labelr)
                            plotr = ax.scatter(i,w,marker=rSymbol,c='Red',s=markerSize,\
                            facecolor=fc,edgecolor=ec,alpha=alpha)

                    plot1 = scatter(i,w,marker=symbol,c=color,s=markerSize,\
                    facecolor=fc,edgecolor=ec,alpha=alpha)
                            

#         totals
#         bin_means,edges,binNumber = stats.binned_statistic(array(allInc), array(allW), statistic='mean', bins=bins)
#         left,right = edges[:-1],edges[1:]
#         X = np.array([left,right]).T.flatten()
#         Y = np.array([bin_means,bin_means]).T.flatten()
#         plt.plot(X,Y, c='Black',ls='dashed',lw=2,alpha=alpha,label=r'$\rm Average ~ EW$')

        # mean
        bin_means,edges,binNumber = stats.binned_statistic(array(allInc), array(allW), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='solid',color='black',lw=2.0,alpha=alpha+0.1,label=r'$\rm Mean ~EW$')
        
        # 90% percentile
        bin_means,edges,binNumber = stats.binned_statistic(array(allInc), array(allW), \
        statistic=lambda y: perc90(y), bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='dashed',color='dimgrey',lw=2.0,alpha=alpha+0.1,label=r'$\rm 90th\% ~EW$')

#         # 50% percentile
#         bin_means,edges,binNumber = stats.binned_statistic(array(allInc), array(allW), \
#         statistic=lambda y: perc50(y), bins=bins)
#         left,right = edges[:-1],edges[1:]        
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plt.plot(X,Y, ls='solid',color='black',lw=2.0,alpha=alpha,label=r'$\rm Median ~EW$')
        
        
#         10th percentile
#         bin_means,edges,binNumber = stats.binned_statistic(array(allInc), array(allW), \
#         statistic=lambda y: perc10(y), bins=bins)
#         left,right = edges[:-1],edges[1:]        
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plt.plot(X,Y, ls='dashed',color='green',lw=1.5,alpha=alpha,label=r'$\rm 10th\% ~EW$')


        # x axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        
        xlabel(r'$\rm Galaxy ~ Inclination ~ [deg]$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ax.legend(scatterpoints=1,prop={'size':14},loc=2,fancybox=True)
        ax.grid(b=None,which='major',axis='both')
#         ylim(-5,max(lyaWList)+100)
        ylim(-2,1000)
        xlim(0,90.5)

        if save:
            savefig('{0}/W(fancy_inc)_mean_90_hist_sep.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


if __name__=="__main__":
    # do the work
    main()
    