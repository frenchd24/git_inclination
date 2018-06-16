#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_EW_impact.py, v 6.0 06/11/18


Based on:
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
    # plot equivalent width as a function of impact parameter/R_vir, split red and blue
    # shifted absorbers
    plot_EW_impact_vir_median = False
    plot_EW_impact_vir_median_save = False
    
    # plot EW as a function of impact/R_vir, also include L_associated and median
    # EW histograms for both
    plot_EW_impact_vir_median_plus_associated = False
    plot_EW_impact_vir_median_plus_associated_save = False

    # plot EW as a function of impact, also include L_associated and median
    # EW histograms for both
    plot_EW_impact_median_plus_associated = False
    plot_EW_impact_median_plus_associated_save = False

    # plot EW as a function of impact, also include L_two and median
    # EW histograms for both
    plot_EW_impact_median_plus_two = False
    plot_EW_impact_median_plus_two_save = False
    
    # plot EW as a function of impact, also include L_two_plus and median
    # EW histograms for both
    plot_EW_impact_median_plus_three = False
    plot_EW_impact_median_plus_three_save = False
    
    # plot EW as a function of impact, also include L_group and median
    # EW histograms for both
    plot_EW_impact_median_plus_group = True
    plot_EW_impact_median_plus_group_save = True


    # some colors
    color_blue = '#436bad'      # french blue
    color_red = '#ec2d01'     # tomato red

    if getpass.getuser() == 'frenchd':

#         pickleFilename = '/Users/frenchd/Research/inclination/git_inclination/picklePilot_plusSALT_14.p'
#         gtPickleFilename = '/Users/frenchd/Research/inclination/git_inclination/pickleGT.p'
#         saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/plotting_code/figs'

#         gtPickleFilename = '/Users/frenchd/Research/inclination/git_inclination/pickleGT_filteredAll.p'
        gtPickleFilename = '/Users/frenchd/Research/GT_update2/pickleGT_filteredAll.p'

        saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/plotting_code/figs/'
        
        isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/isolated6.p'
        L_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_isolated6.p'
        L_associated_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated6.p'
        L_associated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated6.p'
        L_nonassociated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated6.p'
        L_two_filename = '/Users/frenchd/Research/inclination/git_inclination/L_two6.p'
        L_two_plus_filename = '/Users/frenchd/Research/inclination/git_inclination/L_two_plus6.p'
        L_group_filename = '/Users/frenchd/Research/inclination/git_inclination/L_group6.p'


    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # pickle file for the whole galaxy table:
    gtPickleFile = open(gtPickleFilename,'rU')
    gtDict = pickle.load(gtPickleFile)
    gtPickleFile.close()
    
    
    # open all the pickle files
    isolated_file = open(isolated_filename,'r')
    L_isolated_file = open(L_isolated_filename,'r')
    L_associated_isolated_file = open(L_associated_isolated_filename,'r')
    L_associated_file = open(L_associated_filename,'r')
    L_nonassociated_file = open(L_nonassociated_filename,'r')
    L_two_file = open(L_two_filename,'r')
    L_two_plus_file = open(L_two_plus_filename,'r')
    L_group_file = open(L_group_filename,'r')

    # unload the data from them
    isolated = pickle.load(isolated_file)
    L_isolated = pickle.load(L_isolated_file)
    L_associated_isolated = pickle.load(L_associated_isolated_file)
    L_associated = pickle.load(L_associated_file)
    L_nonassociated = pickle.load(L_nonassociated_file)
    L_two = pickle.load(L_two_file)
    L_two_plus = pickle.load(L_two_plus_file)
    L_group = pickle.load(L_group_file)
    
    # close the files
    isolated_file.close()
    L_isolated_file.close()
    L_associated_isolated_file.close()
    L_associated_file.close()
    L_nonassociated_file.close()
    L_two_file.close()
    L_two_plus_file.close()
    L_group_file.close()
    
    # which dataset to use for plotting?
    dataSet = L_associated_isolated
    
    Lya_vs = dataSet['Lya_vs']
    e_Lya_vs = dataSet['e_Lya_vs']
    Lya_Ws = dataSet['Lya_Ws']
    e_Lya_Ws = dataSet['e_Lya_Ws']
    Nas = dataSet['Nas']
    e_Nas = dataSet['e_Nas']
    bs = dataSet['bs']
    e_bs = dataSet['e_bs']
    Ws = dataSet['Ws']
    e_Ws = dataSet['e_Ws']
    targets = dataSet['targets']
    z_targets = dataSet['z_targets']
    RA_targets = dataSet['RA_targets']
    Dec_targets = dataSet['Dec_targets']
    Names = dataSet['Names']
    RA_galaxies = dataSet['RA_galaxies']
    Dec_galaxies = dataSet['Dec_galaxies']
    impacts = dataSet['impacts']
    azimuths = dataSet['azimuths']
    PAs = dataSet['PAs']
    incs = dataSet['incs']
    adjustedIncs = dataSet['adjustedIncs']
    ls = dataSet['ls']
    l_cuss = dataSet['l_cuss']
    R_virs = dataSet['R_virs']
    cuss = dataSet['cuss']
    MajDiams = dataSet['MajDiams']
    MTypes = dataSet['MTypes']
    Vhels = dataSet['Vhels']
    vcorrs = dataSet['vcorrs']
    bestDists = dataSet['bestDists']
    e_bestDists = dataSet['e_bestDists']
    group_nums = dataSet['group_nums']
    group_mems = dataSet['group_mems']
    group_dists = dataSet['group_dists']
    Lstar_meds = dataSet['Lstar_meds']
    e_Lstar_meds = dataSet['e_Lstar_meds']
    Bmags = dataSet['Bmags']



    majorAxisL = gtDict['majorAxis']
    incL = gtDict['inc']
    adjustedIncL = gtDict['adjustedInc']
    paL = gtDict['PA']
    BmagL = gtDict['Bmag']
#     Bmag_sdssL = gtDict['Bmag_sdss']
    RID_medianL = gtDict['RID_median']
    RID_meanL = gtDict['RID_mean']
    RID_stdL = gtDict['RID_std']
    VhelL = gtDict['Vhel']
    RAdegL = gtDict['RAdeg']
    DEdegL = gtDict['DEdeg']
    NameL= gtDict['Name']
    
    allPA = paL
    allInclinations = []
    allAdjustedIncs = []
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

            allAdjustedIncs.append(i)
            
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
    # plot equivalent width as a function of impact parameter, split between
    # absorbers at impact < 1 R_vir, overplot median histograms for each
    #
    # FOR INCLUDE = 1 ABSORBERS

    
    plotW_impact_virseparate_include1 = False
    save = False
    
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
    
    if plot_EW_impact_vir_median:
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
        color = 'Blue'
        bSymbol = 'D'
        rSymbol = 'o'
        maxW = 1500.

        xVals = []
        yVals = []
        redX = []
        redY = []
        blueX = []
        blueY = []
        
        for i, w, rvir, Lya_v, Vhel in zip(impacts, Lya_Ws, R_virs, Lya_vs, Vhels):
            vel_dif = float(Lya_v) - float(Vhel)
        
            # check if all the values are good
            xVal = float(i)/float(rvir)
            yVal = float(w)
            
            xVals.append(xVal)
            yVals.append(yVal)
                
                
            if vel_dif >= 0:
                # redshifted
                color = color_red
                symbol = rSymbol
                
                if float(i) > float(rvir):
                    # impact parameter > virial radius
                    a = alpha
                    fc = color
                    ec = 'black'
            
                if float(i) <= float(rvir):
                    # impact parameter <= virial radius
                    a = alphaInside
                    fc = 'none'
                    ec = color
            
                redX.append(xVal)
                redY.append(yVal)
                
                if countr == 0:
                    countr +=1
                    plotr = ax.scatter(xVal,
                                        yVal,
                                        marker=symbol,
                                        c=color_red,
                                        s=markerSize,
                                        facecolor=fc,
                                        edgecolor=ec,
                                        alpha=a,
                                        label=labelr)
                
            else:
                # blueshifted
                color = color_blue
                symbol = bSymbol

                blueX.append(xVal)
                blueY.append(yVal)
                
                if float(i) > float(rvir):
                    # impact parameter > virial radius
                     a = alpha
                     fc = color
                     ec = 'black'
                
                if float(i) <= float(rvir):
                    # impact parameter <= virial radius
                    a = alphaInside
                    fc = 'none'
                    ec = color
            
                if countb == 0:
                    countb +=1
                    plotb = ax.scatter(xVal,
                                        yVal,
                                        marker=symbol,
                                        c=color_blue,
                                        facecolor=fc,
                                        edgecolor=ec,
                                        s=markerSize,
                                        alpha=a,
                                        label=labelb)
            
            plot1 = scatter(xVal,
                            yVal,
                            marker=symbol,
                            c=color,
                            s=markerSize,
                            facecolor=fc,
                            edgecolor=ec,
                            alpha=a)
     
    
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
#         ylim(0,1300)
        xlim(0, 2.0)

        if plot_EW_impact_vir_median_save:
            savefig('{0}/W(impact_vir)_median_binSize{1}_virsep.pdf'.format(saveDirectory, binSize),format='pdf',bbox_inches='tight')
        else:
            show()



#########################################################################################
#########################################################################################

##########################################################################################
##########################################################################################
    
    if plot_EW_impact_vir_median_plus_associated:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1

        alpha_isolated = 0.7
        alpha_associated = 0.7
        alpha_bins = 0.8
        markerSize = 60
        
        binSize = 0.5
        bins = arange(0,2.5,binSize)
        
        label_isolated = r'$\rm Isolated$'
        label_associated = r'$\rm Associated$'
        symbol_isolated = 'D'
        symbol_associated = 'o'
        color_isolated = color_blue
        color_associated = color_red
        
        maxW = 1500.


        associated_Lya_Ws = L_associated['Lya_Ws']
        associated_R_virs = L_associated['R_virs']
        associated_impacts = L_associated['impacts']
        
        isolated_xs = np.array(impacts)/np.array(R_virs)
        isolated_ys = Lya_Ws
        associated_xs = np.array(associated_impacts)/np.array(associated_R_virs)
        associated_ys = associated_Lya_Ws
        
        
        # isolated
        plot1 = scatter(isolated_xs,
                        isolated_ys,
                        marker=symbol_isolated,
                        c=color_isolated,
                        s=markerSize,
                        edgecolor='black',
                        alpha=alpha_isolated)
                        
        # associated
        plot1 = scatter(associated_xs,
                        associated_ys,
                        marker=symbol_associated,
                        c=color_associated,
                        s=markerSize,
                        edgecolor='black',
                        alpha=alpha_associated)

    
        # histogram isolated
        bin_means, edges, binNumber = stats.binned_statistic(isolated_xs,
                                                            isolated_ys,
                                                            statistic='median',
                                                            bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,
            Y,
            ls='dashed',
            color=color_isolated,
            lw=1.7,
            alpha=alpha_bins,
            label=r'$\rm Isolated~ Median ~EW$')
        
        
        # histogram associated
        bin_means, edges, binNumber = stats.binned_statistic(associated_xs,
                                                            associated_ys,
                                                            statistic='median',
                                                            bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,
            Y,
            ls='dashed',
            color=color_associated,
            lw=1.7,
            alpha=alpha_bins,
            label=r'$\rm Assoc. ~Median ~EW$')
        
    
        
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
#         xlim(0, 2.0)

        if plot_EW_impact_vir_median_plus_associated_save:
            savefig('{0}/W(impact_vir)_median_binSize{1}_plus_associated.pdf'.format(saveDirectory, binSize),format='pdf',bbox_inches='tight')
        else:
            show()



#########################################################################################
#########################################################################################
##########################################################################################
##########################################################################################
    
    if plot_EW_impact_median_plus_associated:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1

        alpha_isolated = 0.35
        alpha_associated = 0.35
        alpha_bins = 0.93
        markerSize = 60
        
        binSize = 50
        bins = arange(0, 550, binSize)
        
        label_isolated = r'$\rm Isolated$'
        label_associated = r'$\rm Associated$'
        symbol_isolated = 'D'
        symbol_associated = 'o'
        color_isolated = color_blue
        color_associated = color_red
        
        maxW = 1500.

        associated_Lya_Ws = L_associated['Lya_Ws']
        associated_R_virs = L_associated['R_virs']
        associated_impacts = L_associated['impacts']
        
        isolated_xs = np.array(impacts)
        isolated_ys = Lya_Ws
        associated_xs = np.array(associated_impacts)
        associated_ys = associated_Lya_Ws
        
        
        # isolated
        plot1 = scatter(isolated_xs,
                        isolated_ys,
                        marker=symbol_isolated,
                        c=color_isolated,
                        s=markerSize,
                        edgecolor='black',
                        alpha=alpha_isolated,
                        label=label_isolated)
                        
        # associated
        plot1 = scatter(associated_xs,
                        associated_ys,
                        marker=symbol_associated,
                        c=color_associated,
                        s=markerSize,
                        edgecolor='black',
                        alpha=alpha_associated,
                        label=label_associated)

    
        # histogram isolated
        bin_means, edges, binNumber = stats.binned_statistic(isolated_xs,
                                                            isolated_ys,
                                                            statistic='median',
                                                            bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,
            Y,
            ls='solid',
            color=color_isolated,
            lw=2.0,
            alpha=alpha_bins,
            label=r'$\rm Isolated~ Median ~EW$')
        
        
        # histogram associated
        bin_means, edges, binNumber = stats.binned_statistic(associated_xs,
                                                            associated_ys,
                                                            statistic='median',
                                                            bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,
            Y,
            ls='dashed',
            color=color_associated,
            lw=2.0,
            alpha=alpha_bins,
            label=r'$\rm Assoc. ~Median ~EW$')
        
    
        
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
        
        
        xlabel(r'$\rm \rho$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        leg = ax.legend(scatterpoints=1,prop={'size':14},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0, 500)

        if plot_EW_impact_median_plus_associated_save:
            savefig('{0}/W(impact)_median_binSize{1}_plus_associated.pdf'.format(saveDirectory, binSize),format='pdf',bbox_inches='tight')
        else:
            show()
            
            
##########################################################################################
##########################################################################################
    
    if plot_EW_impact_median_plus_two:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1

        alpha_isolated = 0.35
        alpha_associated = 0.35
        alpha_two = 0.35
        alpha_bins = 0.93
        markerSize = 60
        
        binSize = 50
        bins = arange(0, 550, binSize)
        
        label_isolated = r'$\rm Isolated$'
        label_associated = r'$\rm Associated$'
        label_two = r'$\rm Two$'
        symbol_isolated = 'D'
        symbol_associated = 'o'
        symbol_two = 'o'

        color_isolated = color_blue
        color_associated = color_red
        color_two = color_red
        
        maxW = 1500.

        associated_Lya_Ws = L_associated['Lya_Ws']
        associated_R_virs = L_associated['R_virs']
        associated_impacts = L_associated['impacts']
        
        two_Lya_Ws = L_two['Lya_Ws']
        two_R_virs = L_two['R_virs']
        two_impacts = L_two['impacts']
                
        isolated_xs = np.array(impacts)
        isolated_ys = Lya_Ws
        associated_xs = np.array(associated_impacts)
        associated_ys = associated_Lya_Ws
        
        two_xs = np.array(two_impacts)
        two_ys = np.array(two_Lya_Ws)
        
        # isolated
        plot1 = scatter(isolated_xs,
                        isolated_ys,
                        marker=symbol_isolated,
                        c=color_isolated,
                        s=markerSize,
                        edgecolor='black',
                        alpha=alpha_isolated,
                        label=label_isolated)
                        
        # two
        plot1 = scatter(two_xs,
                        two_ys,
                        marker=symbol_two,
                        c=color_two,
                        s=markerSize,
                        edgecolor='black',
                        alpha=alpha_two,
                        label=label_two)

    
        # histogram isolated
        bin_means, edges, binNumber = stats.binned_statistic(isolated_xs,
                                                            isolated_ys,
                                                            statistic='median',
                                                            bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,
            Y,
            ls='solid',
            color=color_isolated,
            lw=2.0,
            alpha=alpha_bins,
            label=r'$\rm Isolated~ Median ~EW$')
        
        
        # histogram two
        bin_means, edges, binNumber = stats.binned_statistic(two_xs,
                                                            two_ys,
                                                            statistic='median',
                                                            bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,
            Y,
            ls='dashed',
            color=color_two,
            lw=2.0,
            alpha=alpha_bins,
            label=r'$\rm Two ~Median ~EW$')
        
    
        
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
        
        
        xlabel(r'$\rm \rho$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        leg = ax.legend(scatterpoints=1,prop={'size':14},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0, 500)

        if plot_EW_impact_median_plus_two_save:
            savefig('{0}/W(impact)_median_binSize{1}_plus_two.pdf'.format(saveDirectory, binSize),format='pdf',bbox_inches='tight')
        else:
            show()



##########################################################################################
##########################################################################################
    
    if plot_EW_impact_median_plus_three:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        
        color_green = '#1b9e77'
        color_purple = '#7570b3'
        color_orange = '#d95f02'

        alpha_isolated = 0.4
        alpha_associated = 0.35
        alpha_two = 0.35
        alpha_three = 0.35
        alpha_bins = 0.93
        markerSize = 60
        
        binSize = 50
        bins = arange(0, 550, binSize)
        
        label_isolated = r'$\rm Isolated$'
        label_associated = r'$\rm Associated$'
        label_two = r'$\rm Two$'
        label_three = r'$\rm Three+$'
        symbol_isolated = 'D'
        symbol_associated = 'o'
        symbol_two = 'o'
        symbol_three = 'o'

        color_isolated = color_orange
        color_associated = color_red
        color_two = color_red
        color_three = color_purple
        
        maxW = 1500.

        associated_Lya_Ws = L_associated['Lya_Ws']
        associated_R_virs = L_associated['R_virs']
        associated_impacts = L_associated['impacts']
        
        two_Lya_Ws = L_two['Lya_Ws']
        two_R_virs = L_two['R_virs']
        two_impacts = L_two['impacts']
        
        three_Lya_Ws = L_two_plus['Lya_Ws']
        three_R_virs = L_two_plus['R_virs']
        three_impacts = L_two_plus['impacts']
                
        isolated_xs = np.array(impacts)
        isolated_ys = Lya_Ws
        
        associated_xs = np.array(associated_impacts)
        associated_ys = associated_Lya_Ws
        
        two_xs = np.array(two_impacts)
        two_ys = np.array(two_Lya_Ws)

        three_xs = np.array(three_impacts)
        three_ys = np.array(three_Lya_Ws)
        
        
        # isolated
        plot1 = scatter(isolated_xs,
                        isolated_ys,
                        marker=symbol_isolated,
                        c=color_isolated,
                        s=markerSize,
                        edgecolor='black',
                        alpha=alpha_isolated,
                        label=label_isolated)
                        
        # three
        plot1 = scatter(three_xs,
                        three_ys,
                        marker=symbol_three,
                        c=color_three,
                        s=markerSize,
                        edgecolor='black',
                        alpha=alpha_three,
                        label=label_three)

    
        # histogram isolated
        bin_means, edges, binNumber = stats.binned_statistic(isolated_xs,
                                                            isolated_ys,
                                                            statistic='median',
                                                            bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,
            Y,
            ls='solid',
            color=color_isolated,
            lw=2.0,
            alpha=alpha_bins,
            label=r'$\rm Isolated~ Median ~EW$')
        
        
        # histogram two
        bin_means, edges, binNumber = stats.binned_statistic(three_xs,
                                                            three_ys,
                                                            statistic='median',
                                                            bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,
            Y,
            ls='dashed',
            color=color_three,
            lw=2.0,
            alpha=alpha_bins,
            label=r'$\rm Three+ ~Median ~EW$')
        
    
        
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
        
        
        xlabel(r'$\rm \rho$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        leg = ax.legend(scatterpoints=1,prop={'size':14},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0, 500)

        if plot_EW_impact_median_plus_three_save:
            savefig('{0}/W(impact)_median_binSize{1}_plus_three.pdf'.format(saveDirectory, binSize),format='pdf',bbox_inches='tight')
        else:
            show()



##########################################################################################
##########################################################################################
    
    if plot_EW_impact_median_plus_group:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        
        color_green = '#1b9e77'
        color_purple = '#7570b3'
        color_orange = '#d95f02'
        color_purple2 = '#984ea3'

        alpha_isolated = 0.4
        alpha_associated = 0.35
        alpha_two = 0.35
        alpha_three = 0.35
        alpha_group = 0.35
        alpha_bins = 0.93
        markerSize = 60
        
        binSize = 50
        bins = arange(0, 550, binSize)
        
        label_isolated = r'$\rm Isolated$'
        label_associated = r'$\rm Associated$'
        label_two = r'$\rm Two$'
        label_three = r'$\rm Three+$'
        label_group = r'$\rm Group$'

        symbol_isolated = 'D'
        symbol_associated = 'o'
        symbol_two = 'o'
        symbol_three = 'o'
        symbol_group = 'o'

        color_isolated = 'black'
        color_associated = color_red
        color_two = color_red
        color_three = color_purple
        color_group = color_orange
        
        maxW = 1500.

        # define the x and y data for the isolated set
        isolated_xs = np.array(impacts)
        isolated_ys = Lya_Ws

        # grab the associated data and define the x and y data
        associated_Lya_Ws = L_associated['Lya_Ws']
        associated_R_virs = L_associated['R_virs']
        associated_impacts = L_associated['impacts']
        associated_xs = np.array(associated_impacts)
        associated_ys = associated_Lya_Ws
        
        
        # grab the two data and define the x and y data
        two_Lya_Ws = L_two['Lya_Ws']
        two_R_virs = L_two['R_virs']
        two_impacts = L_two['impacts']
        two_xs = np.array(two_impacts)
        two_ys = np.array(two_Lya_Ws)
        
        
        # grab the two_plus data and define the x and y data
        three_Lya_Ws = L_two_plus['Lya_Ws']
        three_R_virs = L_two_plus['R_virs']
        three_impacts = L_two_plus['impacts']
        three_xs = np.array(three_impacts)
        three_ys = np.array(three_Lya_Ws)
        
        
        # grab the group data and define the x and y data
        group_Lya_Ws = L_group['Lya_Ws']
        group_R_virs = L_group['R_virs']
        group_impacts = L_group['impacts']
        group_xs = np.array(group_impacts)
        group_ys = np.array(group_Lya_Ws)
        
        
        # isolated
        plot1 = scatter(isolated_xs,
                        isolated_ys,
                        marker=symbol_isolated,
                        c=color_isolated,
                        s=markerSize,
                        edgecolor='black',
                        alpha=alpha_isolated,
                        label=label_isolated)
                        
        # three
        plot1 = scatter(group_xs,
                        group_ys,
                        marker=symbol_group,
                        c=color_group,
                        s=markerSize,
                        edgecolor='black',
                        alpha=alpha_group,
                        label=label_group)

    
        # histogram isolated
        bin_means, edges, binNumber = stats.binned_statistic(isolated_xs,
                                                            isolated_ys,
                                                            statistic='median',
                                                            bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,
            Y,
            ls='solid',
            color=color_isolated,
            lw=2.0,
            alpha=alpha_bins,
            label=r'$\rm Isolated~ Median ~EW$')
        
        
        # histogram two
        bin_means, edges, binNumber = stats.binned_statistic(group_xs,
                                                            group_ys,
                                                            statistic='median',
                                                            bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,
            Y,
            ls='dashed',
            color=color_group,
            lw=2.0,
            alpha=alpha_bins,
            label=r'$\rm Group ~Median ~EW$')
        
    
        
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
        
        
        xlabel(r'$\rm \rho$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        leg = ax.legend(scatterpoints=1,prop={'size':14},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0, 500)

        if plot_EW_impact_median_plus_group_save:
            savefig('{0}/W(impact)_median_binSize{1}_plus_group.pdf'.format(saveDirectory, binSize),format='pdf',bbox_inches='tight')
        else:
            show()



##########################################################################################
##########################################################################################
    
    if plot_EW_impact_median_plus_group:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        
        color_green = '#1b9e77'
        color_purple = '#7570b3'
        color_orange = '#d95f02'
        color_purple2 = '#984ea3'

        alpha_isolated = 0.4
        alpha_associated = 0.35
        alpha_two = 0.35
        alpha_three = 0.35
        alpha_group = 0.35
        alpha_bins = 0.93
        markerSize = 60
        
        binSize = 50
        bins = arange(0, 550, binSize)
        
        label_isolated = r'$\rm Isolated$'
        label_associated = r'$\rm Associated$'
        label_two = r'$\rm Two$'
        label_three = r'$\rm Three+$'
        label_group = r'$\rm Group$'

        symbol_isolated = 'D'
        symbol_associated = 'o'
        symbol_two = 'o'
        symbol_three = 'o'
        symbol_group = 'o'

        color_isolated = 'black'
        color_associated = color_red
        color_two = color_red
        color_three = color_purple
        color_group = color_orange
        
        maxW = 1500.

        # define the x and y data for the isolated set
        isolated_xs = np.array(impacts)
        isolated_ys = Lya_Ws

        # grab the associated data and define the x and y data
        associated_Lya_Ws = L_associated['Lya_Ws']
        associated_R_virs = L_associated['R_virs']
        associated_impacts = L_associated['impacts']
        associated_xs = np.array(associated_impacts)
        associated_ys = associated_Lya_Ws
        
        
        # grab the two data and define the x and y data
        two_Lya_Ws = L_two['Lya_Ws']
        two_R_virs = L_two['R_virs']
        two_impacts = L_two['impacts']
        two_xs = np.array(two_impacts)
        two_ys = np.array(two_Lya_Ws)
        
        
        # grab the two_plus data and define the x and y data
        three_Lya_Ws = L_two_plus['Lya_Ws']
        three_R_virs = L_two_plus['R_virs']
        three_impacts = L_two_plus['impacts']
        three_xs = np.array(three_impacts)
        three_ys = np.array(three_Lya_Ws)
        
        
        # grab the group data and define the x and y data
        group_Lya_Ws = L_group['Lya_Ws']
        group_R_virs = L_group['R_virs']
        group_impacts = L_group['impacts']
        group_xs = np.array(group_impacts)
        group_ys = np.array(group_Lya_Ws)
        
        
        # isolated
        plot1 = scatter(isolated_xs,
                        isolated_ys,
                        marker=symbol_isolated,
                        c=color_isolated,
                        s=markerSize,
                        edgecolor='black',
                        alpha=alpha_isolated,
                        label=label_isolated)
                        
        # three
        plot1 = scatter(group_xs,
                        group_ys,
                        marker=symbol_group,
                        c=color_group,
                        s=markerSize,
                        edgecolor='black',
                        alpha=alpha_group,
                        label=label_group)

    
        # histogram isolated
        bin_means, edges, binNumber = stats.binned_statistic(isolated_xs,
                                                            isolated_ys,
                                                            statistic='median',
                                                            bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,
            Y,
            ls='solid',
            color=color_isolated,
            lw=2.0,
            alpha=alpha_bins,
            label=r'$\rm Isolated~ Median ~EW$')
        
        
        # histogram two
        bin_means, edges, binNumber = stats.binned_statistic(group_xs,
                                                            group_ys,
                                                            statistic='median',
                                                            bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,
            Y,
            ls='dashed',
            color=color_group,
            lw=2.0,
            alpha=alpha_bins,
            label=r'$\rm Group ~Median ~EW$')
        
    
        
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
        
        
        xlabel(r'$\rm \rho$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        leg = ax.legend(scatterpoints=1,prop={'size':14},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0, 500)

        if plot_EW_impact_median_plus_group_save:
            savefig('{0}/W(impact)_median_binSize{1}_plus_group.pdf'.format(saveDirectory, binSize),format='pdf',bbox_inches='tight')
        else:
            show()


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


if __name__=="__main__":
    # do the work
    main()
    