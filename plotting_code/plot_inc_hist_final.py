#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_inc_hist_final.py, v7 06/06/18



Taken from:
$Id:  plot_inc_hist_all2.py, v6 01/03/18

Plot histograms of the inclinations of the associated galaxies (and all galaxies too)

    - this now combines plotCosIncDifHist_full2.py, plotFancyIncHist_full2.py, 
    plotFancyCosIncDifHist_full2.py, plotCosIncHist_full2.py, plotIncHist_full2.py,
    and plotFancyCosIncHist_full2.py all into this one program (split into sections below)

v5: updated to work with the new, automatically updated LG_correlation_combined5.csv
    (12/04/15) - original updates to the individual files
    
v5.1: updated for LG_correlation_combined5_8_edit2.csv
    (2/17/2016)
    
v5.2: remake plots with v_hel instead of vcorr (4/22/16)

v5.3: remake plots with new large galaxy sample (7/13/16) -> /plots4/

v5.4: include ability to limit results based on environment number (7/14/16)
    also likelihood limits

v6: Updated for AAS_2017 results (01/03/18)

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

    
def main():
    
    
    # plotting options:
    # This one plots two histograms. The top one has red and blue shifted absorbers
    # overlaid, with a combined black outline histogram. The bottom plot is the full 
    # galaxy table
    plot_adjustedInc_redblue_overlaid_all = True
    plot_adjustedInc_redblue_overlaid_all_save = True
    
    # This one plots a CDF for red, blue and all galaxy inclinations all in the same frame
    plot_CDF_red_blue_all = False
    plot_CDF_red_blue_all_save = False
    
    # this plots histograms of inclination for the L_associated_isolated, L_associated,
    # and L_group all at once
    plot_adjustedInc_all_associated_isolated_group = True
    plot_adjustedInc_all_associated_isolated_group_save = True
    
    # 
    plot_CDF_adjustedInc_all_associated_isolated_group = False
    plot_CDF_adjustedInc_all_associated_isolated_group_save = False
    
    
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
    

    
    Lya_vs = L_associated_isolated['Lya_vs']
    e_Lya_vs = L_associated_isolated['e_Lya_vs']
    Lya_Ws = L_associated_isolated['Lya_Ws']
    e_Lya_Ws = L_associated_isolated['e_Lya_Ws']
    Nas = L_associated_isolated['Nas']
    e_Nas = L_associated_isolated['e_Nas']
    bs = L_associated_isolated['bs']
    e_bs = L_associated_isolated['e_bs']
    Ws = L_associated_isolated['Ws']
    e_Ws = L_associated_isolated['e_Ws']
    targets = L_associated_isolated['targets']
    z_targets = L_associated_isolated['z_targets']
    RA_targets = L_associated_isolated['RA_targets']
    Dec_targets = L_associated_isolated['Dec_targets']
    Names = L_associated_isolated['Names']
    RA_galaxies = L_associated_isolated['RA_galaxies']
    Dec_galaxies = L_associated_isolated['Dec_galaxies']
    impacts = L_associated_isolated['impacts']
    azimuths = L_associated_isolated['azimuths']
    PAs = L_associated_isolated['PAs']
    incs = L_associated_isolated['incs']
    adjustedIncs = L_associated_isolated['adjustedIncs']
    ls = L_associated_isolated['ls']
    l_cuss = L_associated_isolated['l_cuss']
    R_virs = L_associated_isolated['R_virs']
    cuss = L_associated_isolated['cuss']
    MajDiams = L_associated_isolated['MajDiams']
    MTypes = L_associated_isolated['MTypes']
    Vhels = L_associated_isolated['Vhels']
    vcorrs = L_associated_isolated['vcorrs']
    bestDists = L_associated_isolated['bestDists']
    e_bestDists = L_associated_isolated['e_bestDists']
    group_nums = L_associated_isolated['group_nums']
    group_mems = L_associated_isolated['group_mems']
    group_dists = L_associated_isolated['group_dists']
    Lstar_meds = L_associated_isolated['Lstar_meds']
    e_Lstar_meds = L_associated_isolated['e_Lstar_meds']
    Bmags = L_associated_isolated['Bmags']


        
    # lists for the full galaxy dataset    
#     allPA = gtDict['allPA']
#     allInclinations = gtDict['allInclinations']
#     allCosInclinations = gtDict['allCosInclinations']
#     allFancyInclinations = gtDict['allFancyInclinations']
#     allCosFancyInclinations = gtDict['allCosFancyInclinations']
#     allDiameter = gtDict['allDiameters']


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
    

            
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
    # fancyInclination histograms for redshifted vs blueshifted distributions of absorbers
    # plotted together and separately all in one, then the whole table as a subplot
    #
    
    if plot_adjustedInc_redblue_overlaid_all:
    
        fig = figure(figsize=(7.7, 5.7))
        subplots_adjust(hspace=0.200)

#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        bins = arange(0,100,10)
        blue = []
        red = []
        blue_assoc = []
        red_assoc = []
        
        alpha_red = 0.6
        alpha_blue = 0.7
        alpha_both = 0.7
        alpha_all = 0.7
                
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        adjustedIncs_associated = L_associated['adjustedIncs']
        Lya_vs_associated = L_associated['Lya_vs']
        Vhels_associated = L_associated['Vhels']

                
        for Lya_v, e_Lya_v, i, w, Vhel in zip(Lya_vs, e_Lya_vs, adjustedIncs, Lya_Ws, Vhels):
            vel_dif = Lya_v - Vhel
            
            if vel_dif >=0:
                # blue shifted galaxy, but absorber is REDSHIFTED
#                 print 'd: ',d
                red.append(i)
                redLya.append(w)
                redLyaErr.append(e_Lya_v)
            else:
                # red shifted absorber, but galaxy is BLUESHIFTED
                blue.append(i)
                blueLya.append(w)
                blueLyaErr.append(e_Lya_v)
                
                
        for Lya_v, i, Vhel in zip(Lya_vs_associated, adjustedIncs_associated, Vhels_associated):
            vel_dif = Lya_v - Vhel
            
            if vel_dif >=0:
                # blue shifted galaxy, but absorber is REDSHIFTED
#                 print 'd: ',d
                red_assoc.append(i)
            else:
                # red shifted absorber, but galaxy is BLUESHIFTED
                blue_assoc.append(i)                
                
                
        # just associated
        ax = fig.add_subplot(211)
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(5)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(2.5)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

#         hist(red,bins=bins,histtype='bar',color='red',hatch='\\',lw=1.5,alpha = alpha,label=r'$\rm Redshifted$')        
#         hist(blue,bins=bins,histtype='bar',color='Blue',lw=1.5,alpha = alpha,label=r'$\rm Blueshifted$')

        # blue
        hist(blue+blue_assoc,
        bins=bins,
        histtype='bar',
        color=color_blue,
        lw=1.5,
        alpha=alpha_blue,
        edgecolor='black',
        label=r'$\rm v_{{Ly\alpha}} < v_{{sys}}$')

        # red
        hist(red+red_assoc,
        bins=bins,
        histtype='bar',
        color=color_red,
        hatch='/',
        lw=1.5,
        alpha=alpha_red,
        edgecolor='black',
        label=r'$\rm v_{{Ly\alpha}} \ge v_{{sys}}$')

        # together
        hist(adjustedIncs+adjustedIncs_associated,
        bins=bins,
        histtype='step',
        color='black',
        lw=1.5,
        alpha=alpha_both,
        edgecolor='black',
        label=r'$\rm All ~ Associated$')
        
#         hist(blue, bins=bins,histtype='bar',color=color_blue, hatch='\\',lw=1.7,alpha = alpha+0.25,label=r'$\rm v_{{Ly\alpha}} < v_{{sys}}$')
#         hist(red, bins=bins,histtype='bar',color=color_red, lw=1.7,alpha = alpha,label=r'$\rm v_{{Ly\alpha}} \ge v_{{sys}}$')        
#         hist(adjustedIncs, bins=bins,histtype='step',color='Black',lw=2.1,alpha = 0.9,label=r'$\rm All ~ Associated$')
        
        legend(scatterpoints=1,prop={'size':15},loc=2,fancybox=True)
        ylabel(r'$\rm Number$')
        

        # full table
        ax = fig.add_subplot(212)
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(5000)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(2500)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        
        # full galaxy table
        hist(allAdjustedIncs,
        bins=bins,
        histtype='bar',
        color='green',
        lw=1.7,
        alpha=alpha_all,
        edgecolor='black',
        label=r'$\rm All ~ Galaxies$')
        
#         hist(allAdjustedIncs, bins=bins, histtype='bar', lw=1.5,color = 'green',alpha=0.9,label=r'$\rm All$')

        legend(scatterpoints=1,prop={'size':12},loc=2,fancybox=True)
        xlabel(r'$\rm Galaxy ~ Inclination ~ [deg]$')
        ylabel(r'$\rm Number$')
#         tight_layout()

        if plot_adjustedInc_redblue_overlaid_all_save:
            savefig('{0}/hist(adjustedInc)_redblue_overlaid_allassoc.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()
            

            
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
    # fancyInclination CDF for redshifted vs blueshifted distributions of absorbers
    # plotted together and separately all in one, then the whole table as a subplot
    #
    # dv = v_absorber - v_galaxy -> positive is redshifted absorber
    
    if plot_CDF_red_blue_all:
    
        fig = figure(figsize=(7.7, 7.7))
        subplots_adjust(hspace=0.200)
        
        alpha_red = 0.65
        alpha_blue = 0.75
        alpha_both = 0.6
        alpha_all = 0.8

        bins = arange(0,90,0.5)
        blue = []
        red = []
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        for Lya_v, e_Lya_v, i, w, Vhel in zip(Lya_vs, e_Lya_vs, adjustedIncs, Lya_Ws, Vhels):
            vel_dif = Lya_v - Vhel
            
            if vel_dif >=0:
                # blue shifted galaxy, but absorber is REDSHIFTED
#                 print 'd: ',d
                red.append(i)
                redLya.append(w)
                redLyaErr.append(e_Lya_v)
            else:
                # red shifted absorber, but galaxy is BLUESHIFTED
                blue.append(i)
                blueLya.append(w)
                blueLyaErr.append(e_Lya_v)
                
                
        # just associated
        ax = fig.add_subplot(111)


        # blue
        hist(blue,
        bins=bins,
        histtype='step',
        cumulative=True,
        normed=1,
        color=color_blue,
        lw=1.0,
        alpha=alpha_blue,
        label=r'$\rm Blueshifted ~CDF$')
        
#         n_blues, bins_blues, patches_blues = hist(blue, bins, normed=1, histtype="step",\
#         cumulative=True,lw=1,label=r'$\rm Blueshifted ~CDF$',color='blue')

        # red
        hist(red,
        bins=bins,
        histtype='step',
        cumulative=True,
        normed=1,
        color=color_red,
        lw=1.0,
        alpha=alpha_red,
        label=r'$\rm Redshifted ~CDF$')
        
#         n_reds, bins_reds, patches_reds = hist(red, bins, normed=1, histtype="step",\
#         cumulative=True, color='red',lw=1,label=r'$\rm Redshifted ~CDF$')

        # together
        hist(adjustedIncs,
        bins=bins,
        histtype='step',
        cumulative=True,
        normed=1,
        color='black',
        lw=1.0,
        alpha=alpha_both,
        label=r'$\rm All ~ Associated$')

        
        # full galaxy table
        hist(allAdjustedIncs,
        bins=bins,
        histtype='step',
        cumulative=True,
        normed=1,
        color='grey',
        lw=1.0,
        alpha=alpha_both,
        label=r'$\rm All ~ Galaxies$')
        
#         n_all, bins_all, patches_blues = hist(allAdjustedIncs, bins=bins, normed=1, \
#         histtype='step',cumulative=True, color='Black',lw=1.5,alpha = 0.9,label=r'$\rm All$')
    
    
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(0.1)
        majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
        minorLocator   = MultipleLocator(0.05)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        xlim(0,90)
        ylim(0,1)
    
        legend(scatterpoints=1,prop={'size':12},loc=2,fancybox=True)
        xlabel(r'$\rm Galaxy ~ Inclination ~ [deg]$')
        ylabel(r'$\rm CDF $')
        tight_layout()

        if plot_CDF_red_blue_all_save:
            savefig('{0}/CDF(fancy_inclination)_red_blue_full_all_{0}.pdf'.format(saveDirectory, ),format='pdf',bbox_inches='tight')
        else:
            show()

            
            
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################

    # fancyInclination histograms for redshifted vs blueshifted distributions of absorbers
    # plotted together and separately all in one, then the whole table as a subplot
    #
    
    if plot_adjustedInc_all_associated_isolated_group:
    
        fig = figure(figsize=(7.7, 5.7))
        subplots_adjust(hspace=0.200)

#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        bins = arange(0,100,10)
        blue = []
        red = []
        
        alpha_L_associated_isolated = 0.7
        alpha_L_associated = 0.7
        alpha_L_two = 0.7
        alpha_L_two_plus = 0.7
        alpha_L_group = 0.7
        alpha_all = 0.8
                
        color_L_associated_isolated = '#e41a1c'
        color_L_associated = '#377eb8'
        color_L_two = '#4daf4a'
        color_L_two_plus = '#984ea3'
        color_L_group = '#ff7f00'
        
        color_L_associated_isolated = color_blue
        color_L_associated = color_red
        color_L_group = 'black'
                
        adjustedIncs_associated_isolated = L_associated_isolated['adjustedIncs']
        adjustedIncs_associated = L_associated['adjustedIncs']
        adjustedIncs_group = L_group['adjustedIncs']
        group_mems = L_group['group_mems']
        
        group_incs = []
        for inc, group_mem in zip(adjustedIncs_group, group_mems):
            if float(group_mem) >=2:
                group_incs.append(inc)

#         L_nonassociated
#         L_two
#         L_two_plus
        
        # just associated
        ax = fig.add_subplot(211)
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(5)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(2.5)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

#         hist(red,bins=bins,histtype='bar',color='red',hatch='\\',lw=1.5,alpha = alpha,label=r'$\rm Redshifted$')        
#         hist(blue,bins=bins,histtype='bar',color='Blue',lw=1.5,alpha = alpha,label=r'$\rm Blueshifted$')

        # adjustedIncs_associated_isolated
        hist(adjustedIncs_associated_isolated,
        bins=bins,
        histtype='step',
        color=color_L_associated_isolated,
        lw=1.5,
        alpha=alpha_L_associated_isolated,
        label=r'$\rm Associated-isolated$')

        # adjustedIncs_associated
        hist(adjustedIncs_associated,
        bins=bins,
        histtype='step',
        color=color_L_associated,
        lw=1.5,
        alpha=alpha_L_associated,
        label=r'$\rm Associated$')

        # together
        hist(group_incs,
        bins=bins,
        histtype='step',
        color=color_L_group,
        lw=1.5,
        alpha=alpha_L_group,
        label=r'$\rm Group$')
        
#         hist(blue, bins=bins,histtype='bar',color=color_blue, hatch='\\',lw=1.7,alpha = alpha+0.25,label=r'$\rm v_{{Ly\alpha}} < v_{{sys}}$')
#         hist(red, bins=bins,histtype='bar',color=color_red, lw=1.7,alpha = alpha,label=r'$\rm v_{{Ly\alpha}} \ge v_{{sys}}$')        
#         hist(adjustedIncs, bins=bins,histtype='step',color='Black',lw=2.1,alpha = 0.9,label=r'$\rm All ~ Associated$')
        
        legend(scatterpoints=1,prop={'size':12},loc=2,fancybox=True)
        ylabel(r'$\rm Number$')
        

        # full table
        ax = fig.add_subplot(212)
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(5000)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(2500)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        
        # full galaxy table
        hist(allAdjustedIncs,
        bins=bins,
        histtype='bar',
        color='green',
        lw=1.7,
        alpha=alpha_all,
        edgecolor='black',
        label=r'$\rm All ~ Galaxies$')
        
#         hist(allAdjustedIncs, bins=bins, histtype='bar', lw=1.5,color = 'green',alpha=0.9,label=r'$\rm All$')

        legend(scatterpoints=1,prop={'size':12},loc=2,fancybox=True)
        xlabel(r'$\rm Galaxy ~ Inclination ~ [deg]$')
        ylabel(r'$\rm Number$')
#         tight_layout()

        if plot_adjustedInc_all_associated_isolated_group_save:
            savefig('{0}/hist(adjustedInc)_associtated_isolated_group_all_overlaid.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()
            
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################

    # fancyInclination histograms for redshifted vs blueshifted distributions of absorbers
    # plotted together and separately all in one, then the whole table as a subplot
    #
    
    if plot_CDF_adjustedInc_all_associated_isolated_group:
    
        fig = figure(figsize=(7.7, 5.7))
        subplots_adjust(hspace=0.200)

#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
#         bins = arange(0,100,10)
        bins = arange(0,100,5)

        blue = []
        red = []
        
        alpha_L_associated_isolated = 0.7
        alpha_L_associated = 0.7
        alpha_L_two = 0.7
        alpha_L_two_plus = 0.7
        alpha_L_group = 0.7
        alpha_all = 0.8
                
        color_L_associated_isolated = '#e41a1c'
        color_L_associated = '#377eb8'
        color_L_two = '#4daf4a'
        color_L_two_plus = '#984ea3'
        color_L_group = '#ff7f00'
        
        color_L_associated_isolated = color_blue
        color_L_associated = color_red
        color_L_group = 'black'
                
        adjustedIncs_associated_isolated = L_associated_isolated['adjustedIncs']
        adjustedIncs_associated = L_associated['adjustedIncs']
        adjustedIncs_group = L_group['adjustedIncs']

#         L_nonassociated
#         L_two
#         L_two_plus
        
        # just associated
        ax = fig.add_subplot(111)
        
#         hist(red,bins=bins,histtype='bar',color='red',hatch='\\',lw=1.5,alpha = alpha,label=r'$\rm Redshifted$')        
#         hist(blue,bins=bins,histtype='bar',color='Blue',lw=1.5,alpha = alpha,label=r'$\rm Blueshifted$')

        # adjustedIncs_associated_isolated
        hist(adjustedIncs_associated_isolated,
        bins=bins,
        histtype='step',
        cumulative=True,
        normed=1,
        color=color_L_associated_isolated,
        lw=1.5,
        alpha=alpha_L_associated_isolated,
        label=r'$\rm Associated-isolated$')

        # adjustedIncs_associated
        hist(adjustedIncs_associated,
        bins=bins,
        histtype='step',
        cumulative=True,
        normed=1,
        color=color_L_associated,
        lw=1.5,
        alpha=alpha_L_associated,
        label=r'$\rm Associated$')

        # together
        hist(adjustedIncs_group,
        bins=bins,
        histtype='step',
        cumulative=True,
        normed=1,
        color=color_L_group,
        lw=1.5,
        alpha=alpha_L_group,
        label=r'$\rm Group$')
        
#         hist(blue, bins=bins,histtype='bar',color=color_blue, hatch='\\',lw=1.7,alpha = alpha+0.25,label=r'$\rm v_{{Ly\alpha}} < v_{{sys}}$')
#         hist(red, bins=bins,histtype='bar',color=color_red, lw=1.7,alpha = alpha,label=r'$\rm v_{{Ly\alpha}} \ge v_{{sys}}$')        
#         hist(adjustedIncs, bins=bins,histtype='step',color='Black',lw=2.1,alpha = 0.9,label=r'$\rm All ~ Associated$')
        
        # full galaxy table
        hist(allAdjustedIncs,
        bins=bins,
        histtype='step',
        cumulative=True,
        normed=1,
        color='green',
        lw=1.7,
        alpha=alpha_all,
        label=r'$\rm All ~ Galaxies$')
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(0.1)
        majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
        minorLocator   = MultipleLocator(0.05)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        xlim(0,90)
        ylim(0,1)
        
        
#         hist(allAdjustedIncs, bins=bins, histtype='bar', lw=1.5,color = 'green',alpha=0.9,label=r'$\rm All$')

        legend(scatterpoints=1,prop={'size':12},loc=2,fancybox=True)
        xlabel(r'$\rm Galaxy ~ Inclination ~ [deg]$')
        ylabel(r'$\rm CDF $')
#         tight_layout()

        if plot_CDF_adjustedInc_all_associated_isolated_group_save:
            savefig('{0}/CDF(adjustedInc)_associated_isolated_group_all_overlaid.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


if __name__=="__main__":
    # do the work
    main()
    