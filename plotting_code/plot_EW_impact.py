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
    
    # plot EW as a function of impact, also include L_three_plus and median
    # EW histograms for both
    plot_EW_impact_median_plus_three = False
    plot_EW_impact_median_plus_three_save = False
    
    # plot EW as a function of impact, also include L_group and median
    # EW histograms for both
    plot_EW_impact_median_plus_group = False
    plot_EW_impact_median_plus_group_save = False
    
    # plot EW as a function of impact/r_vir w/o splitting, add other sets if you want
    plot_EW_impact_vir_mean_plus = False
    plot_EW_impact_vir_mean_plus_save = False
    
    # plot EW as a function of impact/r_vir w/o splitting, add other sets if you want
    plot_EW_impact_vir_MType = True
    plot_EW_impact_vir_MType_save = True
    
    # plot EW as a function of impact w/o splitting, add other sets if you want
    plot_EW_impact_MType = True
    plot_EW_impact_MType_save = True
    
    min_EW = 0
    max_EW = 10000
    
    # which data set to use
    data_set = '_double'

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
        
        isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/isolated8{0}.p'.format(data_set)
        L_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_isolated8{0}.p'.format(data_set)
        L_associated_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated8{0}.p'.format(data_set)
        L_associated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated8{0}.p'.format(data_set)
        L_nonassociated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated8{0}.p'.format(data_set)
        L_two_filename = '/Users/frenchd/Research/inclination/git_inclination/L_two8{0}.p'.format(data_set)
        L_three_plus_filename = '/Users/frenchd/Research/inclination/git_inclination/L_three_plus8{0}.p'.format(data_set)
        L_group_filename = '/Users/frenchd/Research/inclination/git_inclination/L_group8{0}.p'.format(data_set)
        L_summed_filename = '/Users/frenchd/Research/inclination/git_inclination/L_summed8{0}.p'.format(data_set)
        all_filename = '/Users/frenchd/Research/inclination/git_inclination/all8{0}.p'.format(data_set)

        morph_dict_filename = '/Users/frenchd/Research/inclination/git_inclination/morph_dict4.p'
    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # pickle file for the whole galaxy table:
    gtPickleFile = open(gtPickleFilename,'rU')
    gtDict = pickle.load(gtPickleFile)
    gtPickleFile.close()
    
    
    # pickle file for the morph dictionary
    morph_dict_file = open(morph_dict_filename,'rU')
    morph_dict = pickle.load(morph_dict_file)
    morph_dict_file.close()
    
    morph_counts = {}
    morph_values = morph_dict.values()
    for morph in morph_values:
        if morph_counts.has_key(morph):
            c = morph_counts[morph]
            c+=1
            morph_counts[morph] = c
        else:
            morph_counts[morph] = 1
            
    k = morph_counts.keys()
    v = morph_counts.values()
    for k, v in zip(k, v):
        print '{0} = {1}'.format(k,v)
        
    print 'done'
    
    
    # open all the pickle files
    isolated_file = open(isolated_filename,'r')
    L_isolated_file = open(L_isolated_filename,'r')
    L_associated_isolated_file = open(L_associated_isolated_filename,'r')
    L_associated_file = open(L_associated_filename,'r')
    L_nonassociated_file = open(L_nonassociated_filename,'r')
    L_two_file = open(L_two_filename,'r')
    L_three_plus_file = open(L_three_plus_filename,'r')
    L_group_file = open(L_group_filename,'r')

    # unload the data from them
    isolated = pickle.load(isolated_file)
    L_isolated = pickle.load(L_isolated_file)
    L_associated_isolated = pickle.load(L_associated_isolated_file)
    L_associated = pickle.load(L_associated_file)
    L_nonassociated = pickle.load(L_nonassociated_file)
    L_two = pickle.load(L_two_file)
    L_three_plus = pickle.load(L_three_plus_file)
    L_group = pickle.load(L_group_file)
    
    # close the files
    isolated_file.close()
    L_isolated_file.close()
    L_associated_isolated_file.close()
    L_associated_file.close()
    L_nonassociated_file.close()
    L_two_file.close()
    L_three_plus_file.close()
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
        
        three_Lya_Ws = L_three_plus['Lya_Ws']
        three_R_virs = L_three_plus['R_virs']
        three_impacts = L_three_plus['impacts']
                
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
        
        
        # grab the three_plus data and define the x and y data
        three_Lya_Ws = L_three_plus['Lya_Ws']
        three_R_virs = L_three_plus['R_virs']
        three_impacts = L_three_plus['impacts']
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
        
        
        # grab the three_plus data and define the x and y data
        three_Lya_Ws = L_three_plus['Lya_Ws']
        three_R_virs = L_three_plus['R_virs']
        three_impacts = L_three_plus['impacts']
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
    
    if plot_EW_impact_vir_mean_plus:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        
        color_purple = '#7570b3'
        color_purple2 = '#984ea3'
        
        color_green = '#1b9e77'
        color_orange = '#d95f02'
        color_purple3 = '#7570b3'
        color_pink = '#e7298a'
        color_lime = '#66a61e'
        color_yellow = '#e6ab02'
        color_brown = '#a6761d'
        color_coal = '#666666'
        

        alpha_isolated = 0.5
        alpha_assoc = 0.5
        alpha_two = 0.5
        alpha_three = 0.5
        alpha_group = 0.5
        alpha_bins = 0.99
        markerSize = 30
        
#         binSize = 50
#         bins = arange(0, 550, binSize)
        binSize = 0.5
        bins = arange(0, 3.0, binSize)

        
        label_isolated = r'$\rm Isolated$'
        label_assoc = r'$\rm Associated$'
        label_two = r'$\rm Two$'
        label_three = r'$\rm Three+$'
        label_group = r'$\rm Group$'

        symbol_isolated = 'D'
        symbol_assoc = 'o'
        symbol_two = 'o'
        symbol_three = 'o'
        symbol_group = 'o'

#         color_isolated = color_green
#         color_assoc = color_orange
#         color_two = color_purple3
#         color_group = color_yellow

        color_isolated = 'black'
        color_assoc = color_green
        color_two = color_purple2
        color_group = color_orange
        
        # define the x and y data for the isolated set
        
        Lya_Ws2 = []
        R_virs2 = []
        impacts2 = []
        for w, r, i in zip(Lya_Ws, R_virs, impacts):
            if float(w) <= max_EW:
                Lya_Ws2.append(w)
                R_virs2.append(r)
                impacts2.append(i)
        
        isolated_xs = np.array(impacts2)/np.array(R_virs2)
        isolated_ys = np.array(Lya_Ws2)
        
        

        # grab the associated data and define the x and y data
        associated_Lya_Ws = L_associated['Lya_Ws']
        associated_R_virs = L_associated['R_virs']
        associated_impacts = L_associated['impacts']
        
        associated_Lya_Ws2 = []
        associated_R_virs2 = []
        associated_impacts2 = []
        for w, r, i in zip(associated_Lya_Ws, associated_R_virs, associated_impacts):
            if float(w) <= max_EW:
                associated_Lya_Ws2.append(w)
                associated_R_virs2.append(r)
                associated_impacts2.append(i)
        
        associated_xs = np.array(associated_impacts2)/np.array(associated_R_virs2)
        associated_ys = np.array(associated_Lya_Ws2)
        
        
        # grab the two data and define the x and y data
        two_Lya_Ws = L_two['Lya_Ws']
        two_R_virs = L_two['R_virs']
        two_impacts = L_two['impacts']

        two_Lya_Ws2 = []
        two_R_virs2 = []
        two_impacts2 = []
        for w, r, i in zip(two_Lya_Ws, two_R_virs, two_impacts):
            if float(w) <= max_EW:
                two_Lya_Ws2.append(w)
                two_R_virs2.append(r)
                two_impacts2.append(i)
        
        two_xs = np.array(two_impacts2)/np.array(two_R_virs2)
        two_ys = np.array(two_Lya_Ws2)
        
        
        # grab the three_plus data and define the x and y data
        three_Lya_Ws = L_three_plus['Lya_Ws']
        three_R_virs = L_three_plus['R_virs']
        three_impacts = L_three_plus['impacts']
        
        three_Lya_Ws2 = []
        three_R_virs2 = []
        three_impacts2 = []
        for w, r, i in zip(three_Lya_Ws, three_R_virs, three_impacts):
            if float(w) <= max_EW:
                three_Lya_Ws2.append(w)
                three_R_virs2.append(r)
                three_impacts2.append(i)
        
        three_xs = np.array(three_impacts2)/np.array(three_R_virs2)
        three_ys = np.array(three_Lya_Ws2)
        
        
        
        # grab the group data and define the x and y data
        group_Lya_Ws = L_group['Lya_Ws']
        group_R_virs = L_group['R_virs']
        group_impacts = L_group['impacts']
        group_mems = L_group['group_mems']

        group_Lya_Ws2 = []
        group_R_virs2 = []
        group_impacts2 = []
        for w, r, i, group in zip(group_Lya_Ws, group_R_virs, group_impacts, group_mems):
            if float(group) >= 2 and float(w) <= max_EW:
                group_Lya_Ws2.append(w)
                group_R_virs2.append(r)
                group_impacts2.append(i)
        
        group_xs = np.array(group_impacts2)/np.array(group_R_virs2)
        group_ys = np.array(group_Lya_Ws2)
        
        
        
##########################################################################################
        # do the plotting 
        
        # isolated
        plot1 = scatter(isolated_xs,
                        isolated_ys,
                        marker=symbol_isolated,
                        c=color_isolated,
                        s=markerSize,
                        edgecolor='black',
                        alpha=alpha_isolated,
                        label=label_isolated)
                        
                        
                        
        # histogram isolated
        bin_means, edges, binNumber = stats.binned_statistic(isolated_xs,
                                                            isolated_ys,
                                                            statistic='mean',
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
            
            
        # associated
        plot1 = scatter(associated_xs,
                        associated_ys,
                        marker=symbol_assoc,
                        c=color_assoc,
                        s=markerSize,
                        edgecolor='black',
                        alpha=alpha_assoc,
                        label=label_assoc)
        
        # histogram associated
        bin_means, edges, binNumber = stats.binned_statistic(associated_xs,
                                                            associated_ys,
                                                            statistic='mean',
                                                            bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,
            Y,
            ls='solid',
            color=color_assoc,
            lw=2.0,
            alpha=alpha_bins,
            label=r'$\rm Assoc. ~Median ~EW$')

           
        # group
        plot1 = scatter(two_xs,
                        two_ys,
                        marker=symbol_two,
                        c=color_two,
                        s=markerSize,
                        edgecolor='black',
                        alpha=alpha_two,
                        label=label_two)
        
        # histogram group
        bin_means, edges, binNumber = stats.binned_statistic(two_xs,
                                                            two_ys,
                                                            statistic='mean',
                                                            bins=bins)
        left,right = edges[:-1],edges[1:]
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,
            Y,
            ls='solid',
            color=color_two,
            lw=2.0,
            alpha=alpha_bins,
            label=r'$\rm Two ~Mean ~EW$')
           
                        
        # group
        plot1 = scatter(group_xs,
                        group_ys,
                        marker=symbol_group,
                        c=color_group,
                        s=markerSize,
                        edgecolor='black',
                        alpha=alpha_group,
                        label=label_group)
        
        # histogram group
        bin_means, edges, binNumber = stats.binned_statistic(group_xs,
                                                            group_ys,
                                                            statistic='mean',
                                                            bins=bins)
        left,right = edges[:-1],edges[1:]
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plot(X,
            Y,
            ls='solid',
            color=color_group,
            lw=2.0,
            alpha=alpha_bins,
            label=r'$\rm Group ~Mean ~EW$')
        
    
        
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %.2f$')
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
        
        
        xlabel(r'$\rm \rho / R_{{vir}}$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        leg = ax.legend(scatterpoints=1,prop={'size':12},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,1500)
        xlim(0, 2.5)

        if plot_EW_impact_vir_mean_plus_save:
            savefig('{0}/W(impact_vir)_mean_binSize{1}_plus4_EWcut{2}-{3}_dataset{4}.pdf'.format(saveDirectory, binSize, min_EW, max_EW, data_set),format='pdf',bbox_inches='tight')
        else:
            show()
            
            
##########################################################################################
##########################################################################################
    
    if plot_EW_impact_vir_MType:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        
        color_purple = '#7570b3'
        color_purple2 = '#984ea3'
        
        color_green = '#1b9e77'
        color_orange = '#d95f02'
        color_purple3 = '#7570b3'
        color_pink = '#e7298a'
        color_lime = '#66a61e'
        color_yellow = '#e6ab02'
        color_brown = '#a6761d'
        color_coal = '#666666'
        
        color_purple4 = '#810f7c'

        
        include_histograms = False
        include_fit = True

        alpha_S = 0.4
        alpha_S0 = 0.99
        alpha_E = 0.99
        alpha_I = 0.99
        alpha_other = 0.99
        alpha_bins = 0.99
        markerSize_S = 30
        markerSize_E = 50
        markerSize_I = 50
        markerSize_S0 = 50
        markerSize_other = 50
                
#         binSize = 50
#         bins = arange(0, 550, binSize)
        binSize = 0.5
        bins = arange(0, 3.0, binSize)

        lw = 0.6
        lw_other = 1.2
        
        label_isolated = r'$\rm Isolated$'
        label_assoc = r'$\rm Associated$'
        label_two = r'$\rm Two$'
        label_three = r'$\rm Three+$'
        label_group = r'$\rm Group$'

        label_S = r'$\rm S$'
        label_S0 = r'$\rm S0$'
        label_E = r'$\rm E$'
        label_I = r'$\rm I$'
        label_other = r'$\rm ?$'

        symbol_E = 'X'
        symbol_S = 'o'
        symbol_S0 = 'D'
        symbol_I = 's'
        symbol_other = 'o'

        color_E = color_red
        color_S = color_blue
        color_S0 = color_purple4
        color_I = color_green
        color_other = 'black'

        color_isolated = 'black'
        color_assoc = color_green
        color_two = color_purple4
        color_group = color_orange
        
        # define the x and y data for the isolated set
        Lya_Ws2_S = []
        R_virs2_S = []
        impacts2_S = []

        Lya_Ws2_S0 = []
        R_virs2_S0 = []
        impacts2_S0 = []

        Lya_Ws2_E = []
        R_virs2_E = []
        impacts2_E = []

        Lya_Ws2_I = []
        R_virs2_I = []
        impacts2_I = []
        
        Lya_Ws2_other = []
        R_virs2_other = []
        impacts2_other = []
        for w, r, i, name in zip(Lya_Ws, R_virs, impacts, Names):
            mtype = morph_dict[name]
            
            if float(w) <= max_EW and float(w) >= min_EW:
                if mtype == 'e':
                    Lya_Ws2_E.append(w)
                    R_virs2_E.append(r)
                    impacts2_E.append(i)
                    print 'e: ',name
                
                elif mtype == 's0':
                    Lya_Ws2_S0.append(w)
                    R_virs2_S0.append(r)
                    impacts2_S0.append(i)
                    print 's0: ',name
                    
                elif mtype == 'i':
                    Lya_Ws2_I.append(w)
                    R_virs2_I.append(r)
                    impacts2_I.append(i)
                    
                elif mtype == 'sa' or mtype == 'sb':
                    # combine sa, sb here
                    Lya_Ws2_S.append(w)
                    R_virs2_S.append(r)
                    impacts2_S.append(i)
                    
                else:
                    # unknown type
                    Lya_Ws2_other.append(w)
                    R_virs2_other.append(r)
                    impacts2_other.append(i)
        
        
        isolated_xs_E = np.array(impacts2_E)/np.array(R_virs2_E)
        isolated_ys_E = np.array(Lya_Ws2_E)
        
        isolated_xs_S = np.array(impacts2_S)/np.array(R_virs2_S)
        isolated_ys_S = np.array(Lya_Ws2_S)
        
        isolated_xs_S0 = np.array(impacts2_S0)/np.array(R_virs2_S0)
        isolated_ys_S0 = np.array(Lya_Ws2_S0)
        
        isolated_xs_I = np.array(impacts2_I)/np.array(R_virs2_I)
        isolated_ys_I = np.array(Lya_Ws2_I)

        isolated_xs_other = np.array(impacts2_other)/np.array(R_virs2_other)
        isolated_ys_other = np.array(Lya_Ws2_other)

        # grab the associated data and define the x and y data
        associated_Lya_Ws = L_associated['Lya_Ws']
        associated_R_virs = L_associated['R_virs']
        associated_impacts = L_associated['impacts']
        associated_MTypes = L_associated['MTypes']
        associated_Names = L_associated['Names']

        associated_Lya_Ws2_S = []
        associated_R_virs2_S = []
        associated_impacts2_S = []

        associated_Lya_Ws2_S0 = []
        associated_R_virs2_S0 = []
        associated_impacts2_S0 = []
        
        associated_Lya_Ws2_E = []
        associated_R_virs2_E = []
        associated_impacts2_E = []
        
        associated_Lya_Ws2_I = []
        associated_R_virs2_I = []
        associated_impacts2_I = []
        
        associated_Lya_Ws2_other = []
        associated_R_virs2_other = []
        associated_impacts2_other = []
        for w, r, i, name in zip(associated_Lya_Ws, associated_R_virs, associated_impacts, associated_Names):
            mtype = morph_dict[name]
            
            if float(w) <= max_EW and float(w) >= min_EW:            
                if mtype == 'e':
                    associated_Lya_Ws2_E.append(w)
                    associated_R_virs2_E.append(r)
                    associated_impacts2_E.append(i)
                    print 'e: ', name
        
                elif mtype == 's0':
                    associated_Lya_Ws2_S0.append(w)
                    associated_R_virs2_S0.append(r)
                    associated_impacts2_S0.append(i)
                    print 's0: ', name
        
                elif mtype == 'i':
                    associated_Lya_Ws2_I.append(w)
                    associated_R_virs2_I.append(r)
                    associated_impacts2_I.append(i)
                    
                elif mtype == 'sa' or mtype == 'sb':
                    # combine sa, sb here
                    associated_Lya_Ws2_S.append(w)
                    associated_R_virs2_S.append(r)
                    associated_impacts2_S.append(i)
                    
                else:
                    # unknown type
                    associated_Lya_Ws2_other.append(w)
                    associated_R_virs2_other.append(r)
                    associated_impacts2_other.append(i)
        
        associated_xs_E = np.array(associated_impacts2_E)/np.array(associated_R_virs2_E)
        associated_ys_E = np.array(associated_Lya_Ws2_E)
        
        associated_xs_S = np.array(associated_impacts2_S)/np.array(associated_R_virs2_S)
        associated_ys_S = np.array(associated_Lya_Ws2_S)
        
        associated_xs_S0 = np.array(associated_impacts2_S0)/np.array(associated_R_virs2_S0)
        associated_ys_S0 = np.array(associated_Lya_Ws2_S0)
        
        associated_xs_I = np.array(associated_impacts2_I)/np.array(associated_R_virs2_I)
        associated_ys_I = np.array(associated_Lya_Ws2_I)
        
        associated_xs_other = np.array(associated_impacts2_other)/np.array(associated_R_virs2_other)
        associated_ys_other = np.array(associated_Lya_Ws2_other)
        
        
        # grab the two data and define the x and y data
#         two_Lya_Ws = L_two['Lya_Ws']
#         two_R_virs = L_two['R_virs']
#         two_impacts = L_two['impacts']
#         two_MTypes = L_two['MTypes']
#         two_Names = L_two['Names']
# 
#         two_Lya_Ws2_S = []
#         two_R_virs2_S = []
#         two_impacts2_S = []
#         
#         two_Lya_Ws2_S0 = []
#         two_R_virs2_S0 = []
#         two_impacts2_S0 = []
#         
#         two_Lya_Ws2_E = []
#         two_R_virs2_E = []
#         two_impacts2_E = []
#         
#         two_Lya_Ws2_I = []
#         two_R_virs2_I = []
#         two_impacts2_I = []
#         for w, r, i, name in zip(two_Lya_Ws, two_R_virs, two_impacts, two_Names):
#             mtype = morph_dict[name]
#             
#             if float(w) <= max_EW and float(w) >= min_EW:
#                 if mtype == 'e':
#                     two_Lya_Ws2_E.append(w)
#                     two_R_virs2_E.append(r)
#                     two_impacts2_E.append(i)
#                     
#                 if mtype == 's0':
#                     two_Lya_Ws2_S0.append(w)
#                     two_R_virs2_S0.append(r)
#                     two_impacts2_S0.append(i)
#         
#                 elif mtype == 'i':
#                     two_Lya_Ws2_I.append(w)
#                     two_R_virs2_I.append(r)
#                     two_impacts2_I.append(i)
#                     
#                 else:
#                     # combine sa, sb here
#                     two_Lya_Ws2_S.append(w)
#                     two_R_virs2_S.append(r)
#                     two_impacts2_S.append(i)
#         
#         two_xs_E = np.array(two_impacts2_E)/np.array(two_R_virs2_E)
#         two_ys_E = np.array(two_Lya_Ws2_E)
# 
#         two_xs_S = np.array(two_impacts2_S)/np.array(two_R_virs2_S)
#         two_ys_S = np.array(two_Lya_Ws2_S)
#         
#         two_xs_S0 = np.array(two_impacts2_S0)/np.array(two_R_virs2_S0)
#         two_ys_S0 = np.array(two_Lya_Ws2_S0)
#         
#         two_xs_I = np.array(two_impacts2_I)/np.array(two_R_virs2_I)
#         two_ys_I = np.array(two_Lya_Ws2_I)
        
        
        # grab the three_plus data and define the x and y data
#         three_Lya_Ws = L_three_plus['Lya_Ws']
#         three_R_virs = L_three_plus['R_virs']
#         three_impacts = L_three_plus['impacts']
#         three_Names = L_three_plus['Names']
# 
#         three_Lya_Ws2_S = []
#         three_R_virs2_S = []
#         three_impacts2_S = []
#         
#         three_Lya_Ws2_S0 = []
#         three_R_virs2_S0 = []
#         three_impacts2_S0 = []
#         
#         three_Lya_Ws2_E = []
#         three_R_virs2_E = []
#         three_impacts2_E = []
#         
#         three_Lya_Ws2_I = []
#         three_R_virs2_I = []
#         three_impacts2_I = []
#         for w, r, i, name in zip(three_Lya_Ws, three_R_virs, three_impacts, three_Names):
#             mtype = morph_dict[name]
#             
#             if float(w) <= max_EW and float(w) >= min_EW:
#                 if mtype == 'e':
#                     three_Lya_Ws2_E.append(w)
#                     three_R_virs2_E.append(r)
#                     three_impacts2_E.append(i)
#                     
#                 if mtype == 's0':
#                     three_Lya_Ws2_S0.append(w)
#                     three_R_virs2_S0.append(r)
#                     three_impacts2_S0.append(i)
#         
#                 elif mtype == 'i':
#                     three_Lya_Ws2_I.append(w)
#                     three_R_virs2_I.append(r)
#                     three_impacts2_I.append(i)
#                     
#                 else:
#                     # combine sa, sb here
#                     three_Lya_Ws2_S.append(w)
#                     three_R_virs2_S.append(r)
#                     three_impacts2_S.append(i)
#         
#         three_xs_E = np.array(three_impacts2_E)/np.array(three_R_virs2_E)
#         three_ys_E = np.array(three_Lya_Ws2_E)
# 
#         three_xs_S = np.array(three_impacts2_S)/np.array(three_R_virs2_S)
#         three_ys_S = np.array(three_Lya_Ws2_S)
#         
#         three_xs_S0 = np.array(three_impacts2_S0)/np.array(three_R_virs2_S0)
#         three_ys_S0 = np.array(three_Lya_Ws2_S0)
#         
#         three_xs_I = np.array(three_impacts2_I)/np.array(three_R_virs2_I)
#         three_ys_I = np.array(three_Lya_Ws2_I)
        
        
        
##########################################################################################
        # do the plotting 
        
        all_Rvir_E = np.array(list(associated_R_virs2_E) + list(R_virs2_E))
        all_Rvir_I = np.array(list(associated_R_virs2_I) + list(R_virs2_I))
        all_Rvir_S = np.array(list(associated_R_virs2_S) + list(R_virs2_S))
        all_Rvir_S0 = np.array(list(associated_R_virs2_S0) + list(R_virs2_S0))
        all_Rvir_other = np.array(list(associated_R_virs2_other) + list(R_virs2_other))
        
        from scipy import stats

        print 'all_Rvir_E: ',stats.describe(all_Rvir_E)
        print 'all_Rvir_S: ',stats.describe(all_Rvir_S)
        print 'all_Rvir_I: ',stats.describe(all_Rvir_I)
        print 'all_Rvir_S0: ',stats.describe(all_Rvir_S0)
        print 'all_Rvir_other: ',stats.describe(all_Rvir_other)

        
        all_xs_E = np.array(list(associated_xs_E) + list(isolated_xs_E))
        all_ys_E = np.array(list(associated_ys_E) + list(isolated_ys_E))

        all_xs_S = np.array(list(associated_xs_S) + list(isolated_xs_S))
        all_ys_S = np.array(list(associated_ys_S) + list(isolated_ys_S))
        
        all_xs_S0 = np.array(list(associated_xs_S0) + list(isolated_xs_S0))
        all_ys_S0 = np.array(list(associated_ys_S0) + list(isolated_ys_S0))

        all_xs_I = np.array(list(associated_xs_I) + list(isolated_xs_I))
        all_ys_I = np.array(list(associated_ys_I) + list(isolated_ys_I))

        all_xs_other = np.array(list(associated_xs_other) + list(isolated_xs_other))
        all_ys_other = np.array(list(associated_ys_other) + list(isolated_ys_other))            

###########
        # Spiral
        plot1 = scatter(all_xs_S,
                        all_ys_S,
                        marker=symbol_S,
                        c=color_S,
                        s=markerSize_S,
                        edgecolor='black',
                        lw = lw,
                        alpha=alpha_S,
                        label=label_S)
                        
        # histogram spiral
        if include_histograms:
            bin_means, edges, binNumber = stats.binned_statistic(all_xs_S,
                                                                all_ys_S,
                                                                statistic='mean',
                                                                bins=bins)
            left,right = edges[:-1],edges[1:]        
            X = array([left,right]).T.flatten()
            Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
            plot(X,
                Y,
                ls='solid',
                color=color_S,
                lw=2.0,
                alpha=alpha_bins,
                label=r'$\rm S-type ~EW$')
                
###########
        # Elliptical
        plot1 = scatter(all_xs_E,
                        all_ys_E,
                        marker=symbol_E,
                        c=color_E,
                        s=markerSize_E,
                        edgecolor='black',
                        lw = lw,
                        alpha=alpha_E,
                        label=label_E)
                        
        # histogram elliptical
        if include_histograms:
            bin_means, edges, binNumber = stats.binned_statistic(all_xs_E,
                                                                all_ys_E,
                                                                statistic='mean',
                                                                bins=bins)
            left,right = edges[:-1],edges[1:]        
            X = array([left,right]).T.flatten()
            Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
            plot(X,
                Y,
                ls='solid',
                color=color_E,
                lw=2.0,
                alpha=alpha_bins,
                label=r'$\rm E-type ~EW$')
            
###########
        # S0 type
        plot1 = scatter(all_xs_S0,
                        all_ys_S0,
                        marker=symbol_S0,
                        c=color_S0,
                        s=markerSize_S0,
                        edgecolor='black',
                        lw = lw,
                        alpha=alpha_S0,
                        label=label_S0)
                        
        # histogram S0
        if include_histograms:
            bin_means, edges, binNumber = stats.binned_statistic(all_xs_S0,
                                                                all_ys_S0,
                                                                statistic='mean',
                                                                bins=bins)
            left,right = edges[:-1],edges[1:]        
            X = array([left,right]).T.flatten()
            Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
            plot(X,
                Y,
                ls='solid',
                color=color_S0,
                lw=2.0,
                alpha=alpha_bins,
                label=r'$\rm S0-type ~EW$')
                
                
###########
        # Irregular
        plot1 = scatter(all_xs_I,
                        all_ys_I,
                        marker=symbol_I,
                        c=color_I,
                        s=markerSize_I,
                        edgecolor='black',
                        lw = lw,
                        alpha=alpha_I,
                        label=label_I)

        # histogram irregular
        if include_histograms:
            bin_means, edges, binNumber = stats.binned_statistic(all_xs_I,
                                                                all_ys_I,
                                                                statistic='mean',
                                                                bins=bins)
            left,right = edges[:-1],edges[1:]        
            X = array([left,right]).T.flatten()
            Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
            plot(X,
                Y,
                ls='solid',
                color=color_I,
                lw=2.0,
                alpha=alpha_bins,
                label=r'$\rm I-type ~EW$')
                
###########
        # unknown
        plot1 = scatter(all_xs_other,
                        all_ys_other,
                        marker=symbol_other,
                        edgecolors=color_other,
                        s=markerSize_other,
                        lw=lw_other,
                        facecolors='none',
                        alpha=alpha_other,
                        label=label_other)
                        
        # histogram unknown
        if include_histograms:
            bin_means, edges, binNumber = stats.binned_statistic(all_xs_other,
                                                                all_ys_other,
                                                                statistic='mean',
                                                                bins=bins)
            left,right = edges[:-1],edges[1:]        
            X = array([left,right]).T.flatten()
            Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
            plot(X,
                Y,
                ls='solid',
                color=color_other,
                lw=2.0,
                alpha=alpha_bins,
                label=r'$\rm ?-type ~EW$')
        
        
###########
        if include_fit:
            from scipy.interpolate import interp1d
            from scipy import interpolate, optimize
            
            def line_func(x, m, b):
                y = m * x + b
                return y
    
            m = -1
            b = 500.
    
            x = np.linspace(0, 2.5, num=500)
            
            x_all =  np.array(list(all_xs_E) + list(all_xs_S) + list(all_xs_S0) + list(all_xs_I) + list(all_xs_other))
            y_all = np.array(list(all_ys_E) + list(all_ys_S) + list(all_ys_S0) + list(all_ys_I) + list(all_ys_other))
    
            popt, pcov = optimize.curve_fit(line_func, x_all, y_all)
            perr = np.sqrt(np.diag(pcov))
            print
            print 'popt: ',popt
            print
            print 'pcov: ',pcov
            print
            print 'perr: ',perr
            print
    
    
            m = int(round(popt[0], 0))
            b = int(round(popt[1], 0))
            m_err = int(perr[0])
            b_err = int(perr[1])
        
#             plt.plot(x, line_func(x, *popt), ls='dashed', color='black', alpha = 0.9, label=r'$\rm Fit = {0}\pm {1}~(\rho/R_{{vir}}) + {2}\pm{3}$'.format(m, m_err, b, b_err))
            plt.plot(x, line_func(x, *popt), ls='dashed', color='black', alpha = 0.9, label=r'$\rm Fit = {0}(\rho/R_{{vir}}) + {1}$'.format(m, b))

            from scipy.stats import linregress
            
            EW_min300 = []
            impact_vir_min300 = []
            for w, imp in zip(y_all, x_all):
                if w >= 300:
                    EW_min300.append(w)
                    impact_vir_min300.append(imp)
            
            print "rho/rvir linregress(a, b) = ",linregress(x_all, y_all)
            print
            print "min300 rho/rvir linregress(a, b) = ",linregress(impact_vir_min300, EW_min300)
            print
        
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %.2f$')
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
        
        
        xlabel(r'$\rm \rho / R_{{vir}}$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        leg = ax.legend(scatterpoints=1,prop={'size':12},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,1400)
        xlim(0, 2.5)

        if plot_EW_impact_vir_MType_save:
            savefig('{0}/W(impact_vir)_MType_binSize{1}_EWcut{2}-{3}_dataset{4}.pdf'.format(saveDirectory, binSize, min_EW, max_EW, data_set),format='pdf',bbox_inches='tight')
        else:
            show()
            
            
##########################################################################################
##########################################################################################
    
    if plot_EW_impact_MType:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        
        color_purple = '#7570b3'
        color_purple2 = '#984ea3'
        
        color_green = '#1b9e77'
        color_orange = '#d95f02'
        color_purple3 = '#7570b3'
        color_pink = '#e7298a'
        color_lime = '#66a61e'
        color_yellow = '#e6ab02'
        color_brown = '#a6761d'
        color_coal = '#666666'
        
        color_purple4 = '#810f7c'
        
        include_histograms = False
        include_fit = True

        alpha_S = 0.4
        alpha_S0 = 0.99
        alpha_E = 0.99
        alpha_I = 0.99
        alpha_other = 0.99
        alpha_bins = 0.99
        markerSize_S = 25
        markerSize_E = 55
        markerSize_I = 55
        markerSize_S0 = 55
        markerSize_other = 55


        lw = 0.6
        lw_other = 1.2
        
        binSize = 50
        bins = arange(0, 550, binSize)
#         binSize = 0.5
#         bins = arange(0, 3.0, binSize)

        
        label_isolated = r'$\rm Isolated$'
        label_assoc = r'$\rm Associated$'
        label_two = r'$\rm Two$'
        label_three = r'$\rm Three+$'
        label_group = r'$\rm Group$'

        label_S = r'$\rm S$'
        label_S0 = r'$\rm S0$'
        label_E = r'$\rm E$'
        label_I = r'$\rm I$'
        label_other = r'$\rm ?$'

        symbol_E = 'X'
        symbol_S = 'o'
        symbol_S0 = 'D'
        symbol_I = 's'
        symbol_other = 'o'

        color_E = color_red
        color_S = color_blue
        color_S0 = color_purple4
        color_I = color_green
        color_other = 'black'

        color_isolated = 'black'
        color_assoc = color_green
        color_two = color_purple4
        color_group = color_orange
        
        # define the x and y data for the isolated set
        Lya_Ws2_S = []
        R_virs2_S = []
        impacts2_S = []

        Lya_Ws2_S0 = []
        R_virs2_S0 = []
        impacts2_S0 = []

        Lya_Ws2_E = []
        R_virs2_E = []
        impacts2_E = []

        Lya_Ws2_I = []
        R_virs2_I = []
        impacts2_I = []
        
        Lya_Ws2_other = []
        R_virs2_other = []
        impacts2_other = []
        for w, r, i, name in zip(Lya_Ws, R_virs, impacts, Names):
            mtype = morph_dict[name]
            
            if float(w) <= max_EW and float(w) >= min_EW:
                
                if mtype == 'e':
                    Lya_Ws2_E.append(w)
                    R_virs2_E.append(r)
                    impacts2_E.append(i)
                
                elif mtype == 's0':
                    Lya_Ws2_S0.append(w)
                    R_virs2_S0.append(r)
                    impacts2_S0.append(i)
                    
                elif mtype == 'i':
                    Lya_Ws2_I.append(w)
                    R_virs2_I.append(r)
                    impacts2_I.append(i)
                    
                elif mtype == 'sa' or mtype == 'sb':
                    # combine sa, sb here
                    Lya_Ws2_S.append(w)
                    R_virs2_S.append(r)
                    impacts2_S.append(i)
                    
                else:
                    # unknown type
                    Lya_Ws2_other.append(w)
                    R_virs2_other.append(r)
                    impacts2_other.append(i)
                    
        
        isolated_xs_E = np.array(impacts2_E)
        isolated_ys_E = np.array(Lya_Ws2_E)
        
        isolated_xs_S = np.array(impacts2_S)
        isolated_ys_S = np.array(Lya_Ws2_S)
        
        isolated_xs_S0 = np.array(impacts2_S0)
        isolated_ys_S0 = np.array(Lya_Ws2_S0)
        
        isolated_xs_I = np.array(impacts2_I)
        isolated_ys_I = np.array(Lya_Ws2_I)
        
        isolated_xs_other = np.array(impacts2_other)
        isolated_ys_other = np.array(Lya_Ws2_other)

        # grab the associated data and define the x and y data
        associated_Lya_Ws = L_associated['Lya_Ws']
        associated_R_virs = L_associated['R_virs']
        associated_impacts = L_associated['impacts']
        associated_MTypes = L_associated['MTypes']
        associated_Names = L_associated['Names']

        associated_Lya_Ws2_S = []
        associated_R_virs2_S = []
        associated_impacts2_S = []

        associated_Lya_Ws2_S0 = []
        associated_R_virs2_S0 = []
        associated_impacts2_S0 = []
        
        associated_Lya_Ws2_E = []
        associated_R_virs2_E = []
        associated_impacts2_E = []
        
        associated_Lya_Ws2_I = []
        associated_R_virs2_I = []
        associated_impacts2_I = []
        
        associated_Lya_Ws2_other = []
        associated_R_virs2_other = []
        associated_impacts2_other = []
        for w, r, i, name in zip(associated_Lya_Ws, associated_R_virs, associated_impacts, associated_Names):
            mtype = morph_dict[name]
            
            if float(w) <= max_EW and float(w) >= min_EW:            
                if mtype == 'e':
                    associated_Lya_Ws2_E.append(w)
                    associated_R_virs2_E.append(r)
                    associated_impacts2_E.append(i)
        
                elif mtype == 's0':
                    associated_Lya_Ws2_S0.append(w)
                    associated_R_virs2_S0.append(r)
                    associated_impacts2_S0.append(i)
        
                elif mtype == 'i':
                    associated_Lya_Ws2_I.append(w)
                    associated_R_virs2_I.append(r)
                    associated_impacts2_I.append(i)
                    
                elif mtype == 'sa' or mtype == 'sb':
                    # combine sa, sb here
                    associated_Lya_Ws2_S.append(w)
                    associated_R_virs2_S.append(r)
                    associated_impacts2_S.append(i)
                    
                else:
                    # unknown types
                    associated_Lya_Ws2_other.append(w)
                    associated_R_virs2_other.append(r)
                    associated_impacts2_other.append(i)
                    
        
        associated_xs_E = np.array(associated_impacts2_E)
        associated_ys_E = np.array(associated_Lya_Ws2_E)
        
        associated_xs_S = np.array(associated_impacts2_S)
        associated_ys_S = np.array(associated_Lya_Ws2_S)
        
        associated_xs_S0 = np.array(associated_impacts2_S0)
        associated_ys_S0 = np.array(associated_Lya_Ws2_S0)
        
        associated_xs_I = np.array(associated_impacts2_I)
        associated_ys_I = np.array(associated_Lya_Ws2_I)
        
        associated_xs_other = np.array(associated_impacts2_other)
        associated_ys_other = np.array(associated_Lya_Ws2_other)
        
        # grab the two data and define the x and y data
#         two_Lya_Ws = L_two['Lya_Ws']
#         two_R_virs = L_two['R_virs']
#         two_impacts = L_two['impacts']
#         two_MTypes = L_two['MTypes']
#         two_Names = L_two['Names']
# 
#         two_Lya_Ws2_S = []
#         two_R_virs2_S = []
#         two_impacts2_S = []
#         
#         two_Lya_Ws2_S0 = []
#         two_R_virs2_S0 = []
#         two_impacts2_S0 = []
#         
#         two_Lya_Ws2_E = []
#         two_R_virs2_E = []
#         two_impacts2_E = []
#         
#         two_Lya_Ws2_I = []
#         two_R_virs2_I = []
#         two_impacts2_I = []
#         for w, r, i, name in zip(two_Lya_Ws, two_R_virs, two_impacts, two_Names):
#             mtype = morph_dict[name]
#             
#             if float(w) <= max_EW and float(w) >= min_EW:
#                 if mtype == 'e':
#                     two_Lya_Ws2_E.append(w)
#                     two_R_virs2_E.append(r)
#                     two_impacts2_E.append(i)
#                     
#                 if mtype == 's0':
#                     two_Lya_Ws2_S0.append(w)
#                     two_R_virs2_S0.append(r)
#                     two_impacts2_S0.append(i)
#         
#                 elif mtype == 'i':
#                     two_Lya_Ws2_I.append(w)
#                     two_R_virs2_I.append(r)
#                     two_impacts2_I.append(i)
#                     
#                 else:
#                     # combine sa, sb here
#                     two_Lya_Ws2_S.append(w)
#                     two_R_virs2_S.append(r)
#                     two_impacts2_S.append(i)
#         
#         two_xs_E = np.array(two_impacts2_E)
#         two_ys_E = np.array(two_Lya_Ws2_E)
# 
#         two_xs_S = np.array(two_impacts2_S)
#         two_ys_S = np.array(two_Lya_Ws2_S)
#         
#         two_xs_S0 = np.array(two_impacts2_S0)
#         two_ys_S0 = np.array(two_Lya_Ws2_S0)
#         
#         two_xs_I = np.array(two_impacts2_I)
#         two_ys_I = np.array(two_Lya_Ws2_I)
        
        
        # grab the three_plus data and define the x and y data
#         three_Lya_Ws = L_three_plus['Lya_Ws']
#         three_R_virs = L_three_plus['R_virs']
#         three_impacts = L_three_plus['impacts']
#         three_Names = L_three_plus['Names']
# 
#         three_Lya_Ws2_S = []
#         three_R_virs2_S = []
#         three_impacts2_S = []
#         
#         three_Lya_Ws2_S0 = []
#         three_R_virs2_S0 = []
#         three_impacts2_S0 = []
#         
#         three_Lya_Ws2_E = []
#         three_R_virs2_E = []
#         three_impacts2_E = []
#         
#         three_Lya_Ws2_I = []
#         three_R_virs2_I = []
#         three_impacts2_I = []
#         for w, r, i, name in zip(three_Lya_Ws, three_R_virs, three_impacts, three_Names):
#             mtype = morph_dict[name]
#             
#             if float(w) <= max_EW and float(w) >= min_EW:
#                 if mtype == 'e':
#                     three_Lya_Ws2_E.append(w)
#                     three_R_virs2_E.append(r)
#                     three_impacts2_E.append(i)
#                     
#                 if mtype == 's0':
#                     three_Lya_Ws2_S0.append(w)
#                     three_R_virs2_S0.append(r)
#                     three_impacts2_S0.append(i)
#         
#                 elif mtype == 'i':
#                     three_Lya_Ws2_I.append(w)
#                     three_R_virs2_I.append(r)
#                     three_impacts2_I.append(i)
#                     
#                 else:
#                     # combine sa, sb here
#                     three_Lya_Ws2_S.append(w)
#                     three_R_virs2_S.append(r)
#                     three_impacts2_S.append(i)
#         
#         three_xs_E = np.array(three_impacts2_E)
#         three_ys_E = np.array(three_Lya_Ws2_E)
# 
#         three_xs_S = np.array(three_impacts2_S)
#         three_ys_S = np.array(three_Lya_Ws2_S)
#         
#         three_xs_S0 = np.array(three_impacts2_S0)
#         three_ys_S0 = np.array(three_Lya_Ws2_S0)
#         
#         three_xs_I = np.array(three_impacts2_I)
#         three_ys_I = np.array(three_Lya_Ws2_I)
        
        
        
##########################################################################################
        # do the plotting 
        
        all_xs_E = np.array(list(associated_xs_E) + list(isolated_xs_E))
        all_ys_E = np.array(list(associated_ys_E) + list(isolated_ys_E))

        all_xs_S = np.array(list(associated_xs_S) + list(isolated_xs_S))
        all_ys_S = np.array(list(associated_ys_S) + list(isolated_ys_S))
        
        all_xs_S0 = np.array(list(associated_xs_S0) + list(isolated_xs_S0))
        all_ys_S0 = np.array(list(associated_ys_S0) + list(isolated_ys_S0))

        all_xs_I = np.array(list(associated_xs_I) + list(isolated_xs_I))
        all_ys_I = np.array(list(associated_ys_I) + list(isolated_ys_I))

        all_xs_other = np.array(list(associated_xs_other) + list(isolated_xs_other))
        all_ys_other = np.array(list(associated_ys_other) + list(isolated_ys_other))

            
###########
        # Spiral
        plot1 = scatter(all_xs_S,
                        all_ys_S,
                        marker=symbol_S,
                        c=color_S,
                        s=markerSize_S,
                        edgecolor='black',
                        lw = lw,
                        alpha=alpha_S,
                        label=label_S)
                        
        # histogram spiral
        if include_histograms:
            bin_means, edges, binNumber = stats.binned_statistic(all_xs_S,
                                                                all_ys_S,
                                                                statistic='mean',
                                                                bins=bins)
            left,right = edges[:-1],edges[1:]        
            X = array([left,right]).T.flatten()
            Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
            plot(X,
                Y,
                ls='solid',
                color=color_S,
                lw=2.0,
                alpha=alpha_bins,
                label=r'$\rm S-type ~EW$')
            
###########
        # Elliptical
        plot1 = scatter(all_xs_E,
                        all_ys_E,
                        marker=symbol_E,
                        c=color_E,
                        s=markerSize_E,
                        edgecolor='black',
                        lw = lw,
                        alpha=alpha_E,
                        label=label_E)
                        
        # histogram elliptical
        if include_histograms:
            bin_means, edges, binNumber = stats.binned_statistic(all_xs_E,
                                                                all_ys_E,
                                                                statistic='mean',
                                                                bins=bins)
            left,right = edges[:-1],edges[1:]        
            X = array([left,right]).T.flatten()
            Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
            plot(X,
                Y,
                ls='solid',
                color=color_E,
                lw=2.0,
                alpha=alpha_bins,
                label=r'$\rm E-type ~EW$')
            
###########
        # S0 type
        plot1 = scatter(all_xs_S0,
                        all_ys_S0,
                        marker=symbol_S0,
                        c=color_S0,
                        s=markerSize_S0,
                        edgecolor='black',
                        lw = lw,
                        alpha=alpha_S0,
                        label=label_S0)
                        
        # histogram S0
        if include_histograms:
            bin_means, edges, binNumber = stats.binned_statistic(all_xs_S0,
                                                                all_ys_S0,
                                                                statistic='mean',
                                                                bins=bins)
            left,right = edges[:-1],edges[1:]        
            X = array([left,right]).T.flatten()
            Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
            plot(X,
                Y,
                ls='solid',
                color=color_S0,
                lw=2.0,
                alpha=alpha_bins,
                label=r'$\rm S0-type ~EW$')


###########
        # Irregular
        plot1 = scatter(all_xs_I,
                        all_ys_I,
                        marker=symbol_I,
                        c=color_I,
                        s=markerSize_I,
                        edgecolor='black',
                        lw = lw,
                        alpha=alpha_I,
                        label=label_I)
                          
        # histogram irregular
        if include_histograms:
            bin_means, edges, binNumber = stats.binned_statistic(all_xs_I,
                                                                all_ys_I,
                                                                statistic='mean',
                                                                bins=bins)
            left,right = edges[:-1],edges[1:]        
            X = array([left,right]).T.flatten()
            Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
            plot(X,
                Y,
                ls='solid',
                color=color_I,
                lw=2.0,
                alpha=alpha_bins,
                label=r'$\rm I-type ~EW$')
        
###########
        # Unknown
        plot1 = scatter(all_xs_other,
                        all_ys_other,
                        marker=symbol_other,
                        edgecolors=color_other,
                        s=markerSize_other,
                        lw=lw_other,
                        facecolors='none',
                        alpha=alpha_other,
                        label=label_other)
                        
        # histogram unknown
        if include_histograms:
            bin_means, edges, binNumber = stats.binned_statistic(all_xs_other,
                                                                all_ys_other,
                                                                statistic='mean',
                                                                bins=bins)
            left,right = edges[:-1],edges[1:]        
            X = array([left,right]).T.flatten()
            Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
            plot(X,
                Y,
                ls='solid',
                color=color_other,
                lw=2.0,
                alpha=alpha_bins,
                label=r'$\rm ?-type ~EW$')




        if include_fit:
            from scipy.interpolate import interp1d
            from scipy import interpolate, optimize
            
            def line_func(x, m, b):
                y = m * x + b
                return y
    
            m = -1
            b = 500.
    
            x = np.linspace(0, 500, num=500)
            
            x_all =  np.array(list(all_xs_E) + list(all_xs_S) + list(all_xs_S0) + list(all_xs_I) + list(all_xs_other))
            y_all = np.array(list(all_ys_E) + list(all_ys_S) + list(all_ys_S0) + list(all_ys_I) + list(all_ys_other))
    
            popt, pcov = optimize.curve_fit(line_func, x_all, y_all)
            perr = np.sqrt(np.diag(pcov))
            print
            print 'popt: ',popt
            print
            print 'pcov: ',pcov
            print
            print 'perr: ',perr
            print
    
    
            m = round(popt[0],2)
            b = int(round(popt[1],0))
        
            plt.plot(x, line_func(x, *popt), ls='dashed', color='black', alpha = 0.9, label=r'$\rm Fit = {0}(\rho) + {1}$'.format(m, b))

            from scipy.stats import linregress
            print "impact only linregress(a, b) = ",linregress(x_all,y_all)
            print
            
            EW_min300 = []
            impact_min300 = []
            for w, imp in zip(y_all, x_all):
                if w >= 300:
                    EW_min300.append(w)
                    impact_min300.append(imp)
            
            print "impact linregress(a, b) = ",linregress(x_all, y_all)
            print
            print "min300 impact linregress(a, b) = ",linregress(impact_min300, EW_min300)
            print




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
        leg = ax.legend(scatterpoints=1,prop={'size':12},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,1400)
        xlim(0, 500)

        if plot_EW_impact_MType_save:
            savefig('{0}/W(impact)_MType_binSize{1}_EWcut{2}-{3}_dataset{4}.pdf'.format(saveDirectory, binSize, min_EW, max_EW, data_set),format='pdf',bbox_inches='tight')
        else:
            show()
            

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


if __name__=="__main__":
    # do the work
    main()
    