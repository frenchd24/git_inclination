#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_EW_az_cus.py, v 5.4 06/28/18

Plot equivalent width as a function of azimuth. Use the likelihood_cus sorted pickle files


Based on:
$Id:  plotW_Az.py, v 5.3 01/05/18

This is the plotW_Az_major bit from histograms3.py. Now is separated, and loads in a pickle
file of the relevant data, as created by "buildDataLists.py"

Plots equivalent width as a function of azimuth, also normalized by galaxy size, separated 
into red and blue shifted absorption samples

Previous (from histograms3.py):
    Plot some stuff for the 100largest initial results

    Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

    Updated for the pilot paper (05/06/15)

v5: updated to work with the new, automatically updated LG_correlation_combined5.csv
    (12/04/15)
    
v5.1: remake plots with v_hel instead of vcorr (4/21/16)
        
v5.2: remake plots with new large galaxy sample (7/14/16) -> /plots4/
        also included ability to limit results by environment number
        
v5.3: remake with LG_correlation_combined5_11_25cut_edit4.csv (9/23/16) -> /plots5

v5.4: update for AAS_2017 (01/05/18)

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

fontScale = 16
rc('text', usetex=True)
rc('font', size=16, family='serif', weight='normal')
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
    # plot equivalent width as a function of azimuth angle for red vs blue
    # shifted absorption
    plot_EW_az = False
    plot_EW_az_save = False
    
    plot_EW_az_allsubsets = False
    plot_EW_az_allsubsets_save = False
    
    # plot combined L_isolated_associated + L_associated vs combined L_two + L_three_plus
    # median histograms ONLY, no data points
    plot_EW_az_assoc_vs_not = True
    plot_EW_az_assoc_vs_not_save = True
    
    # plot EW vs azimuth for L_isolated_associated sample. Median 
    # histograms included
    plot_EW_az_isolated = False
    plot_EW_az_isolated_save = False
    
    # plot EW vs azimuth for L_associated sample. Median 
    # histograms included
    plot_EW_az_assoc = False
    plot_EW_az_assoc_save = False
    
    # plot EW vs azimuth for L_two sample. Median 
    # histograms included
    plot_EW_az_two = False
    plot_EW_az_two_save = False
    
    # plot EW vs azimuth for L_isolated_associated and L_associated samples. Median 
    # histograms included
    plot_EW_az_isolated_assoc = False
    plot_EW_az_isolated_assoc_save = False
    
    

    # some colors
    color_blue = '#436bad'      # french blue
    color_red = '#ec2d01'     # tomato red

    if getpass.getuser() == 'frenchd':

#         pickleFilename = '/Users/frenchd/Research/inclination/git_inclination/picklePilot_plusSALT_14.p'
#         gtPickleFilename = '/Users/frenchd/Research/inclination/git_inclination/pickleGT.p'
#         saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/plotting_code/figs'

#         gtPickleFilename = '/Users/frenchd/Research/inclination/git_inclination/pickleGT_filteredAll.p'
        gtPickleFilename = '/Users/frenchd/Research/GT_update2/pickleGT_filteredAll.p'

        saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/plotting_code/likelihood_cus/figs/'
        
        isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/isolated7_cus.p'
        L_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_isolated7_cus.p'
        L_associated_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated7_cus.p'
        L_associated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated7_cus.p'
        L_nonassociated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated7_cus.p'
        L_two_filename = '/Users/frenchd/Research/inclination/git_inclination/L_two7_cus.p'
        L_three_plus_filename = '/Users/frenchd/Research/inclination/git_inclination/L_three_plus7_cus.p'
        L_group_filename = '/Users/frenchd/Research/inclination/git_inclination/L_group7_cus.p'


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

########################################################################################
########################################################################################

    # plot equivalent width as a function of azimuth normalized by galaxy size, separated
    # into red and blue shifted absorption samples
    #
    
    plotW_Az_major = False
    save = False
    
    if plotW_Az_major:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        
        # give some stats:
        lessThan45 = 0
        for a in azList:
            if a <=45:
                lessThan45 +=1
        print '{0}/{1} have az <= 45 degrees'.format(lessThan45,len(azList))
        print 'average, median azimuth: ',average(azList),', ',median(azList)
        
        for d,a,w,m in zip(difList,azList,lyaWList,majList):
            # check if all the values are okay
            if isNumber(d) and isNumber(a) and isNumber(w) and isNumber(m):
                if d!=-99 and a!=-99 and w!=-99 and m!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(a/m,w,c='Blue',s=50,label= labelb)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(a/m,w,c='Red',s=50,label= labelr)
                
                    plot1 = scatter(a/m,w,c=color,s=50)
            
        title('W(azimuth/diameter) for red vs blue absorption')
        xlabel(r'Azimuth / Major Axis')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,max(lyaWList)+50)
        xlim(0,10)

        if save:
            savefig('{0}/W(azimuth_diameter)_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()


#########################################################################################
#########################################################################################

    # plot equivalent width as a function of azimuth normalized by galaxy size, separated
    # into red and blue shifted absorption samples
    #
    
    plotW_Az_vir = False
    save = False
    
    if plotW_Az_vir:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        
        # give some stats:
        lessThan45 = 0
        for a in azList:
            if a <=45:
                lessThan45 +=1
        print '{0}/{1} have az <= 45 degrees'.format(lessThan45,len(azList))
        print 'average, median azimuth: ',average(azList),', ',median(azList)
        
        for d,a,w,m in zip(difList,azList,lyaWList,virList):
            # check if all the values are okay
            if isNumber(d) and isNumber(a) and isNumber(w) and isNumber(m):
                if d!=-99 and a!=-99 and w!=-99 and m!=-99:
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(a/m,w,c='Blue',s=50,label= labelb)
                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(a/m,w,c='Red',s=50,label= labelr)
                
                    plot1 = scatter(a/m,w,c=color,s=50)
            
        title('W(azimuth/R_vir) for red vs blue absorption')
        xlabel(r'$\rm Azimuth / R_{vir}$')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,max(lyaWList)+50)
        xlim(0,1)

        if save:
            savefig('{0}/W(azimuth_vir)_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()


#########################################################################################
#########################################################################################
    
    if plot_EW_az:
        fig = figure(figsize=(7.7, 5.7))
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        bSymbol = 'D'
        rSymbol = 'o'
        alpha = 0.7
#         bins = 9.
#         bins = 10.
        bins = arange(0,100, 10)

        
        allAz  = []
        allW = []
        
        for Lya_v, Vhel, i, r, w, e, a in zip(Lya_vs, Vhels, impacts, R_virs, Lya_Ws, e_Lya_Ws, azimuths):
            vel_dif = float(Lya_v) - float(Vhel)
            
            count +=1
            allAz.append(float(a))
            allW.append(float(w))
                
            if vel_dif >= 0:
                # gas is red shifted compared to galaxy
                color = color_red
                symbol = rSymbol

                if countr == 0:
                    countr +=1
                    plotr = ax.scatter(a,
                                        w,
                                        c=color_red,
                                        s=50,
                                        label=labelr,
                                        marker=symbol,
                                        alpha=alpha)

            else:
                # galaxy is behind absorber, so gas is blue shifted
                color = color_blue
                symbol = bSymbol

                if countb == 0:
                    countb +=1
                    plotb = ax.scatter(a,
                                        w,
                                        c=color_blue,
                                        s=50,
                                        label=labelb,
                                        marker=symbol,
                                        alpha=alpha)

        
            plot1 = scatter(a,
                            w,
                            c=color,
                            s=50,
                            marker=symbol,
                            alpha=alpha)

                    
        print 'countr: ',countr
        print 'countb: ',countb
        print 'count: ',count
        
        # mean
        bin_means,edges,binNumber = stats.binned_statistic(array(allAz), array(allW), \
        statistic='mean', bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='solid',color='black',lw=2.0,alpha=alpha+0.1,label=r'$\rm Mean ~EW$')
        
        # 90% percentile
        bin_means,edges,binNumber = stats.binned_statistic(array(allAz), array(allW), \
        statistic=lambda y: perc90(y), bins=bins)
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='dashed',color='dimgrey',lw=2.0,alpha=alpha+0.1,label=r'$\rm 90th\% ~EW$')
        
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
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
        
#         title('W(azimuth) for red vs blue shifted absorption')
        xlabel(r'$\rm Azimuth ~[deg]$')
        ylabel(r'$\rm Equivalent ~Width ~[m\AA]$')
        
#         legend(scatterpoints=1,fancybox=True,prop={'size':15},loc=2)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1000)
        xlim(0,90)
        tight_layout()
    
    
        print
        print 'len(azimuths): ',len(azimuths)
        print
        print 'max(Lya_Ws): ',max(Lya_Ws)

        if plot_EW_az_save:
            savefig('{0}/W(azimuth).pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()



#########################################################################################
#########################################################################################
    
    if plot_EW_az_allsubsets:
        fig = figure(figsize=(7.7, 5.7))
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        bSymbol = 'D'
        rSymbol = 'o'
        alpha = 0.7
#         bins = 9.
#         bins = 10.
        bins = arange(0,100, 10)

        L_associated_azimuths = L_associated['azimuths']
        L_nonassociated_azimuths = L_nonassociated['azimuths']
        L_two_azimuths = L_two['azimuths']
        L_three_plus_azimuths = L_three_plus['azimuths']
        L_group_azimuths = L_group['azimuths']
        
        L_associated_Lya_Ws = L_associated['Lya_Ws']
        L_nonassociated_Lya_Ws = L_nonassociated['Lya_Ws']
        L_two_Lya_Ws = L_two['Lya_Ws']
        L_three_plus_Lya_Ws = L_three_plus['Lya_Ws']
        L_group_Lya_Ws = L_group['Lya_Ws']
        
        
        alpha_L_associated_isolated = 1.0
        alpha_L_associated = 0.8
        alpha_L_two = 1.0
        alpha_L_three_plus = 1.0
        alpha_L_group = 1.0
        
        lw = 2.0
        
        color_L_associated_isolated = '#377eb8' # blue
        color_L_associated = 'black'
        color_L_two = '#ff7f00' # orange
        color_L_three_plus = '#984ea3' # purple
        color_L_group = 'grey'
    
        
        # mean L_associated_isolated
        bin_means,edges,binNumber = stats.binned_statistic(array(azimuths), 
                                                            array(Lya_Ws),
                                                            statistic='mean',
                                                            bins=bins)
                                                            
        left,right = edges[:-1],edges[1:]
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='solid',color='black',lw=2.0,alpha=alpha+0.1,label=r'$\rm Isolated$')
        
        
        # mean L_associated
        bin_means,edges,binNumber = stats.binned_statistic(array(L_associated_azimuths), 
                                                            array(L_associated_Lya_Ws),
                                                            statistic='mean',
                                                            bins=bins)
                                                            
        left,right = edges[:-1],edges[1:]
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='solid',color='black',lw=2.0,alpha=alpha+0.1,label=r'$\rm Associated$')


        # mean L_two
        bin_means,edges,binNumber = stats.binned_statistic(array(L_two_azimuths), 
                                                            array(L_two_Lya_Ws),
                                                            statistic='mean',
                                                            bins=bins)
                                                            
        left,right = edges[:-1],edges[1:]
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='solid',color='black',lw=2.0,alpha=alpha+0.1,label=r'$\rm Two$')
        
        
        # mean L_three_plus
        bin_means,edges,binNumber = stats.binned_statistic(array(L_three_plus_azimuths), 
                                                            array(L_three_plus_Lya_Ws),
                                                            statistic='mean',
                                                            bins=bins)
                                                            
        left,right = edges[:-1],edges[1:]
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='solid',color='black',lw=2.0,alpha=alpha+0.1,label=r'$\rm Three +$')
        
        
        # mean L_three_plus
        bin_means,edges,binNumber = stats.binned_statistic(array(L_group_azimuths), 
                                                            array(L_group_Lya_Ws),
                                                            statistic='mean',
                                                            bins=bins)
                                                            
        left,right = edges[:-1],edges[1:]
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='solid',color='black',lw=2.0,alpha=alpha+0.1,label=r'$\rm Group$')
        
        
        # 90% percentile
#         bin_means,edges,binNumber = stats.binned_statistic(array(allAz), array(allW), \
#         statistic=lambda y: perc90(y), bins=bins)
#         left,right = edges[:-1],edges[1:]        
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plt.plot(X,Y, ls='dashed',color='dimgrey',lw=2.0,alpha=alpha+0.1,label=r'$\rm 90th\% ~EW$')
        
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
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
        
#         title('W(azimuth) for red vs blue shifted absorption')
        xlabel(r'$\rm Azimuth ~[deg]$')
        ylabel(r'$\rm Equivalent ~Width ~[m\AA]$')
        
#         legend(scatterpoints=1,fancybox=True,prop={'size':15},loc=2)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1000)
        xlim(0,90)
        tight_layout()
    
    
        print
        print 'len(azimuths): ',len(azimuths)
        print
        print 'max(Lya_Ws): ',max(Lya_Ws)

        if plot_EW_az_allsubsets_save:
            savefig('{0}/W(azimuth)_allsubsets.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()       
            
            
#########################################################################################
#########################################################################################
    
    if plot_EW_az_isolated:
        fig = figure(figsize=(7.7, 5.7))
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        symbol_diamond = 'D'
        symbol_circle = 'o'
        alpha = 0.7
#         bins = 9.
#         bins = 10.
        bins = arange(0,100, 10)
        
        
        color_green = '#1b9e77'
        color_orange = '#d95f02'
        color_purple3 = '#7570b3'
        color_pink = '#e7298a'
        color_lime = '#66a61e'
        color_yellow = '#e6ab02'
        color_brown = '#a6761d'
        color_coal = '#666666'
        
        alpha_isolated = 0.4
        alpha_assoc = 0.4
        
        size_isolated = 30
        size_assoc = 30
        

        L_associated_azimuths = L_associated['azimuths']
        L_nonassociated_azimuths = L_nonassociated['azimuths']
        L_two_azimuths = L_two['azimuths']
        L_three_plus_azimuths = L_three_plus['azimuths']
        L_group_azimuths = L_group['azimuths']
        
        L_associated_Lya_Ws = L_associated['Lya_Ws']
        L_nonassociated_Lya_Ws = L_nonassociated['Lya_Ws']
        L_two_Lya_Ws = L_two['Lya_Ws']
        L_three_plus_Lya_Ws = L_three_plus['Lya_Ws']
        L_group_Lya_Ws = L_group['Lya_Ws']
        
        
        alpha_L_associated_isolated = 1.0
        alpha_L_associated = 0.8
        alpha_L_two = 1.0
        alpha_L_three_plus = 1.0
        alpha_L_group = 1.0
        
        lw = 2.0
        
        color_L_associated_isolated = '#377eb8' # blue
        color_L_associated = 'black'
        color_L_two = '#ff7f00' # orange
        color_L_three_plus = '#984ea3' # purple
        color_L_group = 'grey'
    
    
        all_assoc_azimuths = azimuths + L_associated_azimuths
        all_assoc_EWs = Lya_Ws + L_associated_Lya_Ws
        
        all_not_assoc_azimuths = L_two_azimuths + L_three_plus_azimuths
        all_not_assoc_EWs = L_two_Lya_Ws + L_three_plus_Lya_Ws
        
        
        # mean L_isolated_associated
        
        plot1 = ax.scatter(azimuths,
                            Lya_Ws,
                            c=color_green,
                            s=size_isolated,
                            label=r'$\rm Isolated$',
                            marker=symbol_diamond,
                            alpha=alpha_isolated)
        
        bin_means,edges,binNumber = stats.binned_statistic(array(azimuths), 
                                                            array(Lya_Ws),
                                                            statistic='mean',
                                                            bins=bins)
                                                            
        left,right = edges[:-1],edges[1:]
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='solid',color='black',lw=2.0,alpha=1.0,label=r'$\rm Isolated$')
        
        # 90% percentile
        bin_means,edges,binNumber = stats.binned_statistic(array(azimuths), array(Lya_Ws), \
        statistic=lambda y: perc90(y), bins=bins)
        left,right = edges[:-1],edges[1:]
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='dotted',color='dimgrey',lw=2.0,alpha=1.0,label=r'$\rm Isolated ~90th\% ~EW$')
        
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
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
        
#         title('W(azimuth) for red vs blue shifted absorption')
        xlabel(r'$\rm Azimuth ~[deg]$')
        ylabel(r'$\rm Equivalent ~Width ~[m\AA]$')
        
#         legend(scatterpoints=1,fancybox=True,prop={'size':15},loc=2)
        ax.grid(b=None,which='major',axis='both')
        legend(scatterpoints=1, prop={'size':14}, loc='upper right', fancybox=True)
        ylim(0,1300)
        xlim(0,90)
        tight_layout()
    
    
        print
        print 'len(azimuths): ',len(azimuths)
        print
        print 'max(Lya_Ws): ',max(Lya_Ws)

        if plot_EW_az_isolated_save:
            savefig('{0}/W(azimuth)_isolated.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()
            
            
            
#########################################################################################
#########################################################################################
    
    if plot_EW_az_assoc:
        fig = figure(figsize=(7.7, 5.7))
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        symbol_diamond = 'D'
        symbol_circle = 'o'
        alpha = 0.7
#         bins = 9.
#         bins = 10.
        bins = arange(0,100, 10)
        
        
        color_green = '#1b9e77'
        color_orange = '#d95f02'
        color_purple3 = '#7570b3'
        color_pink = '#e7298a'
        color_lime = '#66a61e'
        color_yellow = '#e6ab02'
        color_brown = '#a6761d'
        color_coal = '#666666'
        
        alpha_isolated = 0.4
        alpha_assoc = 0.4
        
        size_isolated = 30
        size_assoc = 30
        

        L_associated_azimuths = L_associated['azimuths']
        L_nonassociated_azimuths = L_nonassociated['azimuths']
        L_two_azimuths = L_two['azimuths']
        L_three_plus_azimuths = L_three_plus['azimuths']
        L_group_azimuths = L_group['azimuths']
        
        L_associated_Lya_Ws = L_associated['Lya_Ws']
        L_nonassociated_Lya_Ws = L_nonassociated['Lya_Ws']
        L_two_Lya_Ws = L_two['Lya_Ws']
        L_three_plus_Lya_Ws = L_three_plus['Lya_Ws']
        L_group_Lya_Ws = L_group['Lya_Ws']
        
        
        alpha_L_associated_isolated = 1.0
        alpha_L_associated = 0.8
        alpha_L_two = 1.0
        alpha_L_three_plus = 1.0
        alpha_L_group = 1.0
        
        lw = 2.0
        
        color_L_associated_isolated = '#377eb8' # blue
        color_L_associated = 'black'
        color_L_two = '#ff7f00' # orange
        color_L_three_plus = '#984ea3' # purple
        color_L_group = 'grey'
    
    
        all_assoc_azimuths = azimuths + L_associated_azimuths
        all_assoc_EWs = Lya_Ws + L_associated_Lya_Ws
        
        all_not_assoc_azimuths = L_two_azimuths + L_three_plus_azimuths
        all_not_assoc_EWs = L_two_Lya_Ws + L_three_plus_Lya_Ws
        
        
        # mean L_isolated_associated
        
        plot1 = ax.scatter(L_associated_azimuths,
                            L_associated_Lya_Ws,
                            c=color_orange,
                            s=size_assoc,
                            label=r'$\rm Assoc.$',
                            marker=symbol_diamond,
                            alpha=alpha_assoc)
        
        bin_means,edges,binNumber = stats.binned_statistic(array(L_associated_azimuths), 
                                                            array(L_associated_Lya_Ws),
                                                            statistic='mean',
                                                            bins=bins)
                                                            
        left,right = edges[:-1],edges[1:]
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='solid',color='black',lw=2.0,alpha=alpha+0.1,label=r'$\rm Assoc.$')
        
        # 90% percentile
        bin_means,edges,binNumber = stats.binned_statistic(array(L_associated_azimuths), 
        array(L_associated_Lya_Ws), statistic=lambda y: perc90(y), bins=bins)
        left,right = edges[:-1],edges[1:]
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='dotted',color='dimgrey',lw=2.0,alpha=1.0,label=r'$\rm Assoc. ~90th\% ~EW$')
        
        
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
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
        
#         title('W(azimuth) for red vs blue shifted absorption')
        xlabel(r'$\rm Azimuth ~[deg]$')
        ylabel(r'$\rm Equivalent ~Width ~[m\AA]$')
        
#         legend(scatterpoints=1,fancybox=True,prop={'size':15},loc=2)
        ax.grid(b=None,which='major',axis='both')
        legend(scatterpoints=1, prop={'size':14}, loc='upper right', fancybox=True)
        ylim(0,1300)
        xlim(0,90)
        tight_layout()
    
    
        print
        print 'len(azimuths): ',len(azimuths)
        print
        print 'max(Lya_Ws): ',max(Lya_Ws)

        if plot_EW_az_assoc_save:
            savefig('{0}/W(azimuth)_assoc.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()
            


#########################################################################################
#########################################################################################
    
    if plot_EW_az_two:
        fig = figure(figsize=(7.7, 5.7))
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        symbol_diamond = 'D'
        symbol_circle = 'o'
        alpha = 0.7
#         bins = 9.
#         bins = 10.
        bins = arange(0,100, 10)
        
        
        color_green = '#1b9e77'
        color_orange = '#d95f02'
        color_purple3 = '#7570b3'
        color_pink = '#e7298a'
        color_lime = '#66a61e'
        color_yellow = '#e6ab02'
        color_brown = '#a6761d'
        color_coal = '#666666'
        
        alpha_isolated = 0.4
        alpha_assoc = 0.4
        alpha = 0.4
        
        size_isolated = 30
        size_assoc = 30
        size = 30
        

        L_associated_azimuths = L_associated['azimuths']
        L_nonassociated_azimuths = L_nonassociated['azimuths']
        L_two_azimuths = L_two['azimuths']
        L_three_plus_azimuths = L_three_plus['azimuths']
        L_group_azimuths = L_group['azimuths']
        
        L_associated_Lya_Ws = L_associated['Lya_Ws']
        L_nonassociated_Lya_Ws = L_nonassociated['Lya_Ws']
        L_two_Lya_Ws = L_two['Lya_Ws']
        L_three_plus_Lya_Ws = L_three_plus['Lya_Ws']
        L_group_Lya_Ws = L_group['Lya_Ws']
        
        
        alpha_L_associated_isolated = 1.0
        alpha_L_associated = 0.8
        alpha_L_two = 1.0
        alpha_L_three_plus = 1.0
        alpha_L_group = 1.0
        
        lw = 2.0
        
        color_L_associated_isolated = '#377eb8' # blue
        color_L_associated = 'black'
        color_L_two = '#ff7f00' # orange
        color_L_three_plus = '#984ea3' # purple
        color_L_group = 'grey'
    
    
        all_assoc_azimuths = azimuths + L_associated_azimuths
        all_assoc_EWs = Lya_Ws + L_associated_Lya_Ws
        
        all_not_assoc_azimuths = L_two_azimuths + L_three_plus_azimuths
        all_not_assoc_EWs = L_two_Lya_Ws + L_three_plus_Lya_Ws
        
        
        # mean L_isolated_associated
        
        plot1 = ax.scatter(L_two_azimuths,
                            L_two_Lya_Ws,
                            c=color_purple3,
                            s=size,
                            label=r'$\rm Two$',
                            marker=symbol_circle,
                            alpha=alpha)
        
        bin_means,edges,binNumber = stats.binned_statistic(array(L_two_azimuths), 
                                                            array(L_two_Lya_Ws),
                                                            statistic='mean',
                                                            bins=bins)
                                                            
        left,right = edges[:-1],edges[1:]
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='solid',color='black',lw=2.0,alpha=alpha+0.1,label=r'$\rm Two$')
        
        # 90% percentile
        bin_means,edges,binNumber = stats.binned_statistic(array(L_two_azimuths), 
        array(L_two_Lya_Ws), statistic=lambda y: perc90(y), bins=bins)
        left,right = edges[:-1],edges[1:]
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='dotted',color='dimgrey',lw=2.0,alpha=1.0,label=r'$\rm Two ~90th\% ~EW$')
        
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
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
        
#         title('W(azimuth) for red vs blue shifted absorption')
        xlabel(r'$\rm Azimuth ~[deg]$')
        ylabel(r'$\rm Equivalent ~Width ~[m\AA]$')
        
#         legend(scatterpoints=1,fancybox=True,prop={'size':15},loc=2)
        ax.grid(b=None,which='major',axis='both')
        legend(scatterpoints=1, prop={'size':14}, loc='upper right', fancybox=True)
        ylim(0,1300)
        xlim(0,90)
        tight_layout()
    
    
        print
        print 'len(azimuths): ',len(azimuths)
        print
        print 'max(Lya_Ws): ',max(Lya_Ws)

        if plot_EW_az_two_save:
            savefig('{0}/W(azimuth)_two.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()

#########################################################################################
#########################################################################################
    
    if plot_EW_az_isolated_assoc:
        fig = figure(figsize=(7.7, 5.7))
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        symbol_diamond = 'D'
        symbol_circle = 'o'
        alpha = 0.7
#         bins = 9.
#         bins = 10.
        bins = arange(0,100, 10)
        
        
        color_green = '#1b9e77'
        color_orange = '#d95f02'
        color_purple3 = '#7570b3'
        color_pink = '#e7298a'
        color_lime = '#66a61e'
        color_yellow = '#e6ab02'
        color_brown = '#a6761d'
        color_coal = '#666666'
        
        alpha_isolated = 0.4
        alpha_assoc = 0.4
        
        size_isolated = 30
        size_assoc = 30
        

        L_associated_azimuths = L_associated['azimuths']
        L_nonassociated_azimuths = L_nonassociated['azimuths']
        L_two_azimuths = L_two['azimuths']
        L_three_plus_azimuths = L_three_plus['azimuths']
        L_group_azimuths = L_group['azimuths']
        
        L_associated_Lya_Ws = L_associated['Lya_Ws']
        L_nonassociated_Lya_Ws = L_nonassociated['Lya_Ws']
        L_two_Lya_Ws = L_two['Lya_Ws']
        L_three_plus_Lya_Ws = L_three_plus['Lya_Ws']
        L_group_Lya_Ws = L_group['Lya_Ws']
        
        
        alpha_L_associated_isolated = 1.0
        alpha_L_associated = 0.8
        alpha_L_two = 1.0
        alpha_L_three_plus = 1.0
        alpha_L_group = 1.0
        
        lw = 2.0
        
        color_L_associated_isolated = '#377eb8' # blue
        color_L_associated = 'black'
        color_L_two = '#ff7f00' # orange
        color_L_three_plus = '#984ea3' # purple
        color_L_group = 'grey'
    
    
        all_assoc_azimuths = azimuths + L_associated_azimuths
        all_assoc_EWs = Lya_Ws + L_associated_Lya_Ws
        
        all_not_assoc_azimuths = L_two_azimuths + L_three_plus_azimuths
        all_not_assoc_EWs = L_two_Lya_Ws + L_three_plus_Lya_Ws
        
        
        # mean L_isolated_associated
        
        plot1 = ax.scatter(azimuths,
                            Lya_Ws,
                            c=color_green,
                            s=size_isolated,
                            label=r'$\rm Isolated$',
                            marker=symbol_diamond,
                            alpha=alpha_isolated)
        
        bin_means,edges,binNumber = stats.binned_statistic(array(azimuths), 
                                                            array(Lya_Ws),
                                                            statistic='mean',
                                                            bins=bins)
                                                            
        left,right = edges[:-1],edges[1:]
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='dashed',color='black',lw=2.0,alpha=alpha+0.1,label=r'$\rm Isolated$')
        
        
        # mean L_associated
        plot1 = ax.scatter(L_associated_azimuths,
                            L_associated_Lya_Ws,
                            c=color_orange,
                            s=size_assoc,
                            label=r'$\rm Isolated$',
                            marker=symbol_circle,
                            alpha=alpha_assoc)
        
        bin_means,edges,binNumber = stats.binned_statistic(array(L_associated_azimuths), 
                                                            array(L_associated_Lya_Ws),
                                                            statistic='mean',
                                                            bins=bins)
                                                            
        left,right = edges[:-1],edges[1:]
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='dotted',color='black',lw=2.0,alpha=alpha+0.1,label=r'$\rm Assoc$')
        
        
        
        
#         bin_means,edges,binNumber = stats.binned_statistic(array(all_assoc_azimuths), 
#                                                             array(all_assoc_EWs),
#                                                             statistic='mean',
#                                                             bins=bins)
#                                                             
#         left,right = edges[:-1],edges[1:]
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plt.plot(X,Y, ls='solid',color='black',lw=2.0,alpha=alpha+0.1,label=r'$\rm All ~Assoc.$')
        
        
        


        # mean L_two
#         bin_means,edges,binNumber = stats.binned_statistic(array(all_not_assoc_azimuths), 
#                                                             array(all_not_assoc_EWs),
#                                                             statistic='mean',
#                                                             bins=bins)
#                                                             
#         left,right = edges[:-1],edges[1:]
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plt.plot(X,Y, ls='solid',color=color_red,lw=2.0,alpha=alpha+0.1,label=r'$\rm Two+$')
        
        
        
        # 90% percentile
#         bin_means,edges,binNumber = stats.binned_statistic(array(allAz), array(allW), \
#         statistic=lambda y: perc90(y), bins=bins)
#         left,right = edges[:-1],edges[1:]        
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plt.plot(X,Y, ls='dashed',color='dimgrey',lw=2.0,alpha=alpha+0.1,label=r'$\rm 90th\% ~EW$')
        
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
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
        
#         title('W(azimuth) for red vs blue shifted absorption')
        xlabel(r'$\rm Azimuth ~[deg]$')
        ylabel(r'$\rm Equivalent ~Width ~[m\AA]$')
        
#         legend(scatterpoints=1,fancybox=True,prop={'size':15},loc=2)
        ax.grid(b=None,which='major',axis='both')
        legend(scatterpoints=1, prop={'size':14}, loc='upper right', fancybox=True)
        ylim(0,1300)
        xlim(0,90)
        tight_layout()
    
    
        print
        print 'len(azimuths): ',len(azimuths)
        print
        print 'max(Lya_Ws): ',max(Lya_Ws)

        if plot_EW_az_isolated_assoc_save:
            savefig('{0}/W(azimuth)_isolated_vs_assoc.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()
            
            
            
#########################################################################################
#########################################################################################
    
    if plot_EW_az_assoc_vs_not:
        fig = figure(figsize=(7.7, 5.7))
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        symbol_diamond = 'D'
        symbol_circle = 'o'
        alpha = 0.7
#         bins = 9.
#         bins = 10.
        bins = arange(0,100, 10)
        
        
        color_green = '#1b9e77'
        color_orange = '#d95f02'
        color_purple3 = '#7570b3'
        color_pink = '#e7298a'
        color_lime = '#66a61e'
        color_yellow = '#e6ab02'
        color_brown = '#a6761d'
        color_coal = '#666666'
        
        alpha_isolated = 0.8
        alpha_assoc = 0.8
        

        L_associated_azimuths = L_associated['azimuths']
        L_nonassociated_azimuths = L_nonassociated['azimuths']
        L_two_azimuths = L_two['azimuths']
        L_three_plus_azimuths = L_three_plus['azimuths']
        L_group_azimuths = L_group['azimuths']
        
        L_associated_Lya_Ws = L_associated['Lya_Ws']
        L_nonassociated_Lya_Ws = L_nonassociated['Lya_Ws']
        L_two_Lya_Ws = L_two['Lya_Ws']
        L_three_plus_Lya_Ws = L_three_plus['Lya_Ws']
        L_group_Lya_Ws = L_group['Lya_Ws']
        
        
        alpha_L_associated_isolated = 1.0
        alpha_L_associated = 0.8
        alpha_L_two = 1.0
        alpha_L_three_plus = 1.0
        alpha_L_group = 1.0
        
        lw = 2.0
        
        color_L_associated_isolated = '#377eb8' # blue
        color_L_associated = 'black'
        color_L_two = '#ff7f00' # orange
        color_L_three_plus = '#984ea3' # purple
        color_L_group = 'grey'
    
    
        all_assoc_azimuths = azimuths + L_associated_azimuths
        all_assoc_EWs = Lya_Ws + L_associated_Lya_Ws
        
        all_not_assoc_azimuths = L_two_azimuths + L_three_plus_azimuths
        all_not_assoc_EWs = L_two_Lya_Ws + L_three_plus_Lya_Ws
        
        
        print
        print 'len(all_assoc_azimuths): ',all_assoc_azimuths
        print
        print 'len(all_not_assoc_azimuths): ',all_not_assoc_azimuths
        print
        
        
        # mean L_isolated_associated
        
#         plot1 = ax.scatter(all_assoc_azimuths,
#                             all_assoc_EWs,
#                             c='black',
#                             s=50,
#                             label=r'$\rm All~Assoc.$',
#                             marker=symbol_diamond,
#                             alpha=alpha_isolated)
        
        bin_means,edges,binNumber = stats.binned_statistic(array(all_assoc_azimuths), 
                                                            array(all_assoc_EWs),
                                                            statistic='mean',
                                                            bins=bins)
                                                            
        left,right = edges[:-1],edges[1:]
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='solid',color='black',lw=2.0,alpha=alpha+0.1,label=r'$\rm All~Assoc.$')
        
        
        # mean L_associated
#         plot1 = ax.scatter(all_not_assoc_azimuths,
#                             all_not_assoc_EWs,
#                             c=color_orange,
#                             s=50,
#                             label=r'$\rm Two+$',
#                             marker=symbol_circle,
#                             alpha=alpha_assoc)
        
        bin_means,edges,binNumber = stats.binned_statistic(array(all_not_assoc_azimuths), 
                                                            array(all_not_assoc_EWs),
                                                            statistic='mean',
                                                            bins=bins)
                                                            
        left,right = edges[:-1],edges[1:]
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        plt.plot(X,Y, ls='dashed',color=color_orange,lw=2.0,alpha=alpha+0.1,label=r'$\rm Two+$')
        
        
        
#         bin_means,edges,binNumber = stats.binned_statistic(array(all_assoc_azimuths), 
#                                                             array(all_assoc_EWs),
#                                                             statistic='mean',
#                                                             bins=bins)
#                                                             
#         left,right = edges[:-1],edges[1:]
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plt.plot(X,Y, ls='solid',color='black',lw=2.0,alpha=alpha+0.1,label=r'$\rm All ~Assoc.$')
        
        
        


        # mean L_two
#         bin_means,edges,binNumber = stats.binned_statistic(array(all_not_assoc_azimuths), 
#                                                             array(all_not_assoc_EWs),
#                                                             statistic='mean',
#                                                             bins=bins)
#                                                             
#         left,right = edges[:-1],edges[1:]
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plt.plot(X,Y, ls='solid',color=color_red,lw=2.0,alpha=alpha+0.1,label=r'$\rm Two+$')
        
        
        
        # 90% percentile
#         bin_means,edges,binNumber = stats.binned_statistic(array(allAz), array(allW), \
#         statistic=lambda y: perc90(y), bins=bins)
#         left,right = edges[:-1],edges[1:]        
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plt.plot(X,Y, ls='dashed',color='dimgrey',lw=2.0,alpha=alpha+0.1,label=r'$\rm 90th\% ~EW$')
        
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
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
        
#         title('W(azimuth) for red vs blue shifted absorption')
        xlabel(r'$\rm Azimuth ~[deg]$')
        ylabel(r'$\rm Equivalent ~Width ~[m\AA]$')
        
#         legend(scatterpoints=1,fancybox=True,prop={'size':15},loc=2)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0,90)
        tight_layout()
    
    
        print
        print 'len(azimuths): ',len(azimuths)
        print
        print 'max(Lya_Ws): ',max(Lya_Ws)

        if plot_EW_az_assoc_vs_not_save:
            savefig('{0}/W(azimuth)_assoc_vs_not.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()
            


###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    