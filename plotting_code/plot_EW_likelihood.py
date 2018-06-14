#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_EW_likelihood.py, v 6.0 06/13/18

Plot EW as a function of likelihood and likelihood_cus


Based on:
$Id:  plotW_Diameter2_final.py, v 5.8 01/05/18

Plot equivalent width, NaV and Doppler parameter as a function of galaxy diameter and R_vir


This is the plotW_Diameter bit from histograms3.py. Now is separated, and loads in a pickle
file of the relevant data, as created by "buildDataLists.py"


Previous (from histograms3.py):
    Plot some stuff for the 100largest initial results

    Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

    Updated for the pilot paper (05/06/15)

v5: updated for pilot paper.
    (1/5/16)
    
v5.1: updated for LG_correlation_combined5_8_edit2.csv with l_min = 0.001 (02/24/2016)

v5.2: remake plots with v_hel instead of vcorr (4/22/16)

v5.3: remake plots with new large galaxy sample (7/13/16) -> /plots4/

v5.4: add the ability to limit results based on 'environment' number (7/14/16)

v5.5: minor formatting updates for pilot paper (8/08/16)

v5.6: update for LG_correlation_combined5_11_25cut_edit4.csv -> /plots5/ (9/26/16)

v5.7: update for first paper revision (12/16/2016)
    - changed median bin calculation to stats.binned_statistic, increased marker size,
        and changed points with imp/R_vir <=1 to open symbols
    
v5.8: update for AAS_2017 (01/05/18)
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


def errors(a):
    # return the standard error in the mean for the input array
    return stats.sem(a)

    

def main():
    # plot EW as a function of likelihood, include median
    # EW histograms
    plot_EW_likelihood_median = False
    plot_EW_likelihood_median_save = False

    # plot EW as a function of likelihood_cus, include median
    # EW histograms
    plot_EW_likelihood_cus_median = False
    plot_EW_likelihood_cus_median_save = False
    
    # plot EW as a function of the sum of all likelihoods around that line, include
    # median histograms
    plot_EW_summed_median = True
    plot_EW_summed_median_save = True
    
    

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
        
        isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/isolated3.p'
        L_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_isolated3.p'
        L_associated_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated3.p'
        L_associated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated3.p'
        L_nonassociated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated3.p'
        L_two_filename = '/Users/frenchd/Research/inclination/git_inclination/L_two3.p'
        L_two_plus_filename = '/Users/frenchd/Research/inclination/git_inclination/L_two_plus3.p'
        L_group_filename = '/Users/frenchd/Research/inclination/git_inclination/L_group3.p'
        L_summed_filename = '/Users/frenchd/Research/inclination/git_inclination/L_summed3.p'


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
    L_summed_file = open(L_summed_filename,'r')


    # unload the data from them
    isolated = pickle.load(isolated_file)
    L_isolated = pickle.load(L_isolated_file)
    L_associated_isolated = pickle.load(L_associated_isolated_file)
    L_associated = pickle.load(L_associated_file)
    L_nonassociated = pickle.load(L_nonassociated_file)
    L_two = pickle.load(L_two_file)
    L_two_plus = pickle.load(L_two_plus_file)
    L_group = pickle.load(L_group_file)
    L_summed = pickle.load(L_summed_file)

    
    # close the files
    isolated_file.close()
    L_isolated_file.close()
    L_associated_isolated_file.close()
    L_associated_file.close()
    L_nonassociated_file.close()
    L_two_file.close()
    L_two_plus_file.close()
    L_group_file.close()
    L_summed_file.close()

    
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
    
##########################################################################################
##########################################################################################
    
    if plot_EW_likelihood_median:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        
        second = 'Associated'
        
        color_green = '#1b9e77'
        color_purple = '#7570b3'
        color_orange = '#d95f02'
        color_purple2 = '#984ea3'

        alpha_isolated = 0.45
        alpha_associated = 0.35
        alpha_two = 0.35
        alpha_three = 0.35
        alpha_group = 0.35
        alpha_second = 0.35
        alpha_bins = 0.94
        markerSize = 60
        
        binSize = 0.2
        bins = arange(0, 1.2, binSize)
        
        label_isolated = r'$\rm Isolated$'
        label_associated = r'$\rm Associated$'
        label_two = r'$\rm Two$'
        label_three = r'$\rm Three+$'
        label_group = r'$\rm Group$'
        label_second = r'$\rm {0}$'.format(second)

        symbol_isolated = 'D'
        symbol_associated = 'o'
        symbol_two = 'o'
        symbol_three = 'o'
        symbol_group = 'o'
        symbol_second = 'o'


        color_isolated = 'black'
        color_associated = color_red
        color_two = color_red
        color_three = color_purple
        color_group = color_orange
        color_second = color_orange
        
        maxW = 1500.

        # define the x and y data for the isolated set
        isolated_xs = np.array(ls)
        isolated_ys = Lya_Ws

        # grab the associated data and define the x and y data
        associated_Lya_Ws = L_associated['Lya_Ws']
        associated_R_virs = L_associated['R_virs']
        associated_impacts = L_associated['impacts']
        associated_MajDiams = L_associated['MajDiams']
        associated_ls = L_associated['ls']
        associated_xs = np.array(associated_ls)
        associated_ys = associated_Lya_Ws
        
        
        # grab the two data and define the x and y data
        two_Lya_Ws = L_two['Lya_Ws']
        two_R_virs = L_two['R_virs']
        two_impacts = L_two['impacts']
        two_MajDiams = L_two['MajDiams']
        two_ls = L_two['ls']
        two_xs = np.array(two_ls)
        two_ys = np.array(two_Lya_Ws)
        
        
        # grab the two_plus data and define the x and y data
        three_Lya_Ws = L_two_plus['Lya_Ws']
        three_R_virs = L_two_plus['R_virs']
        three_impacts = L_two_plus['impacts']
        three_MajDiams = L_two_plus['MajDiams']
        three_ls = L_two_plus['ls']
        three_xs = np.array(three_ls)
        three_ys = np.array(three_Lya_Ws)
        
        
        # grab the group data and define the x and y data
        group_Lya_Ws = L_group['Lya_Ws']
        group_R_virs = L_group['R_virs']
        group_impacts = L_group['impacts']
        group_MajDiams = L_group['MajDiams']
        group_ls = L_group['ls']
        group_xs = np.array(group_ls)
        group_ys = np.array(group_Lya_Ws)
        
        if second == 'Group':
            second_xs = group_xs
            second_ys = group_ys

        if second == 'Associated':
            second_xs = associated_xs
            second_ys = associated_ys
            
        
        # isolated
        plot1 = scatter(isolated_xs,
                        isolated_ys,
                        marker=symbol_isolated,
                        c=color_isolated,
                        s=markerSize,
                        edgecolor='black',
                        alpha=alpha_isolated,
                        label=label_isolated)
        
        if second:                
            # second data set
            plot1 = scatter(second_xs,
                            second_ys,
                            marker=symbol_second,
                            c=color_second,
                            s=markerSize,
                            edgecolor='black',
                            alpha=alpha_second,
                            label=label_second)

    
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
        
        
        if second:
            # histogram two
            bin_means, edges, binNumber = stats.binned_statistic(second_xs,
                                                                second_ys,
                                                                statistic='median',
                                                                bins=bins)
            left,right = edges[:-1],edges[1:]        
            X = array([left,right]).T.flatten()
            Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
            plot(X,
                Y,
                ls='dashed',
                color=color_second,
                lw=2.0,
                alpha=alpha_bins,
                label=r'$\rm {0} ~Median ~EW$'.format(second))
        
    
        
        # x-axis
        majorLocator   = MultipleLocator(0.25)
        majorFormatter = FormatStrFormatter(r'$\rm %.2f$')
        minorLocator   = MultipleLocator(0.05)
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
        
        
        xlabel(r'$\rm \mathcal{L} $')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        leg = ax.legend(scatterpoints=1,prop={'size':14},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0, 1)

        if plot_EW_likelihood_median_save:
            savefig('{0}/W(likelihood)_median_plus_{1}_binSize{2}.pdf'.format(saveDirectory, second, binSize),format='pdf',bbox_inches='tight')
        else:
            show()


##########################################################################################
##########################################################################################
    
##########################################################################################
##########################################################################################
    
    if plot_EW_likelihood_cus_median:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        
        second = 'Associated'
        
        color_green = '#1b9e77'
        color_purple = '#7570b3'
        color_orange = '#d95f02'
        color_purple2 = '#984ea3'

        alpha_isolated = 0.45
        alpha_associated = 0.35
        alpha_two = 0.35
        alpha_three = 0.35
        alpha_group = 0.35
        alpha_second = 0.35
        alpha_bins = 0.94
        markerSize = 60
        
        binSize = 0.2
        bins = arange(0, 1.2, binSize)
        
        label_isolated = r'$\rm Isolated$'
        label_associated = r'$\rm Associated$'
        label_two = r'$\rm Two$'
        label_three = r'$\rm Three+$'
        label_group = r'$\rm Group$'
        label_second = r'$\rm {0}$'.format(second)

        symbol_isolated = 'D'
        symbol_associated = 'o'
        symbol_two = 'o'
        symbol_three = 'o'
        symbol_group = 'o'
        symbol_second = 'o'


        color_isolated = 'black'
        color_associated = color_red
        color_two = color_red
        color_three = color_purple
        color_group = color_orange
        color_second = color_orange
        
        maxW = 1500.

        # define the x and y data for the isolated set
        isolated_xs = np.array(l_cuss)
        isolated_ys = Lya_Ws

        # grab the associated data and define the x and y data
        associated_Lya_Ws = L_associated['Lya_Ws']
        associated_R_virs = L_associated['R_virs']
        associated_impacts = L_associated['impacts']
        associated_MajDiams = L_associated['MajDiams']
        associated_ls = L_associated['ls']
        associated_l_cuss = L_associated['l_cuss']
        associated_xs = np.array(associated_l_cuss)
        associated_ys = associated_Lya_Ws
        
        
        # grab the two data and define the x and y data
        two_Lya_Ws = L_two['Lya_Ws']
        two_R_virs = L_two['R_virs']
        two_impacts = L_two['impacts']
        two_MajDiams = L_two['MajDiams']
        two_ls = L_two['ls']
        two_l_cuss = L_two['l_cuss']
        two_xs = np.array(two_l_cuss)
        two_ys = np.array(two_Lya_Ws)
        
        
        # grab the two_plus data and define the x and y data
        three_Lya_Ws = L_two_plus['Lya_Ws']
        three_R_virs = L_two_plus['R_virs']
        three_impacts = L_two_plus['impacts']
        three_MajDiams = L_two_plus['MajDiams']
        three_ls = L_two_plus['ls']
        three_l_cuss = L_two_plus['l_cuss']
        three_xs = np.array(three_l_cuss)
        three_ys = np.array(three_Lya_Ws)
        
        
        # grab the group data and define the x and y data
        group_Lya_Ws = L_group['Lya_Ws']
        group_R_virs = L_group['R_virs']
        group_impacts = L_group['impacts']
        group_MajDiams = L_group['MajDiams']
        group_ls = L_group['ls']
        group_l_cuss = L_group['l_cuss']
        group_xs = np.array(group_l_cuss)
        group_ys = np.array(group_Lya_Ws)
        
        if second == 'Group':
            second_xs = group_xs
            second_ys = group_ys

        if second == 'Associated':
            second_xs = associated_xs
            second_ys = associated_ys
            
        
        # isolated
        plot1 = scatter(isolated_xs,
                        isolated_ys,
                        marker=symbol_isolated,
                        c=color_isolated,
                        s=markerSize,
                        edgecolor='black',
                        alpha=alpha_isolated,
                        label=label_isolated)
        
        if second:                
            # second data set
            plot1 = scatter(second_xs,
                            second_ys,
                            marker=symbol_second,
                            c=color_second,
                            s=markerSize,
                            edgecolor='black',
                            alpha=alpha_second,
                            label=label_second)

    
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
        
        
        if second:
            # histogram two
            bin_means, edges, binNumber = stats.binned_statistic(second_xs,
                                                                second_ys,
                                                                statistic='median',
                                                                bins=bins)
            left,right = edges[:-1],edges[1:]        
            X = array([left,right]).T.flatten()
            Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
            plot(X,
                Y,
                ls='dashed',
                color=color_second,
                lw=2.0,
                alpha=alpha_bins,
                label=r'$\rm {0} ~Median ~EW$'.format(second))
        
    
        
        # x-axis
        majorLocator   = MultipleLocator(0.25)
        majorFormatter = FormatStrFormatter(r'$\rm %.2f$')
        minorLocator   = MultipleLocator(0.05)
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
        
        
        xlabel(r'$\rm \mathcal{L}_{{1.5}} $')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        leg = ax.legend(scatterpoints=1,prop={'size':14},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0, 1)

        if plot_EW_likelihood_cus_median_save:
            savefig('{0}/W(likelihood_cus)_median_plus_{1}_binSize{2}.pdf'.format(saveDirectory, second, binSize),format='pdf',bbox_inches='tight')
        else:
            show()


  
##########################################################################################
##########################################################################################
    
    if plot_EW_summed_median:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        
        second = False
        
        color_green = '#1b9e77'
        color_purple = '#7570b3'
        color_orange = '#d95f02'
        color_purple2 = '#984ea3'

        alpha_isolated = 0.45
        alpha_associated = 0.35
        alpha_two = 0.35
        alpha_three = 0.35
        alpha_group = 0.35
        alpha_second = 0.35
        alpha_bins = 0.94
        markerSize = 60
        
        binSize = 0.2
        bins = arange(0, 1.2, binSize)
        
        label_isolated = r'$\rm Isolated$'
        label_associated = r'$\rm Associated$'
        label_two = r'$\rm Two$'
        label_three = r'$\rm Three+$'
        label_group = r'$\rm Group$'
        label_second = r'$\rm {0}$'.format(second)

        symbol_isolated = 'D'
        symbol_associated = 'o'
        symbol_two = 'o'
        symbol_three = 'o'
        symbol_group = 'o'
        symbol_second = 'o'


        color_isolated = 'black'
        color_associated = color_red
        color_two = color_red
        color_three = color_purple
        color_group = color_orange
        color_second = color_orange
        
        maxW = 1500.

        # define the x and y data for the summed set
        summed_ls = L_summed['summed_ls']
        summed_Lya_Ws = L_summed['Lya_Ws']
        
        main_xs = np.sum(np.array(summed_ls))
        main_ys = Lya_Ws

        # grab the associated data and define the x and y data
        associated_Lya_Ws = L_associated['Lya_Ws']
        associated_R_virs = L_associated['R_virs']
        associated_impacts = L_associated['impacts']
        associated_MajDiams = L_associated['MajDiams']
        associated_ls = L_associated['ls']
        associated_l_cuss = L_associated['l_cuss']
        associated_xs = np.array(associated_l_cuss)
        associated_ys = associated_Lya_Ws
        
        
        # grab the two data and define the x and y data
        two_Lya_Ws = L_two['Lya_Ws']
        two_R_virs = L_two['R_virs']
        two_impacts = L_two['impacts']
        two_MajDiams = L_two['MajDiams']
        two_ls = L_two['ls']
        two_l_cuss = L_two['l_cuss']
        two_xs = np.array(two_l_cuss)
        two_ys = np.array(two_Lya_Ws)
        
        
        # grab the two_plus data and define the x and y data
        three_Lya_Ws = L_two_plus['Lya_Ws']
        three_R_virs = L_two_plus['R_virs']
        three_impacts = L_two_plus['impacts']
        three_MajDiams = L_two_plus['MajDiams']
        three_ls = L_two_plus['ls']
        three_l_cuss = L_two_plus['l_cuss']
        three_xs = np.array(three_l_cuss)
        three_ys = np.array(three_Lya_Ws)
        
        
        # grab the group data and define the x and y data
        group_Lya_Ws = L_group['Lya_Ws']
        group_R_virs = L_group['R_virs']
        group_impacts = L_group['impacts']
        group_MajDiams = L_group['MajDiams']
        group_ls = L_group['ls']
        group_l_cuss = L_group['l_cuss']
        group_xs = np.array(group_l_cuss)
        group_ys = np.array(group_Lya_Ws)
        
        
        
        
        if second == 'Group':
            second_xs = group_xs
            second_ys = group_ys

        if second == 'Associated':
            second_xs = associated_xs
            second_ys = associated_ys
            
        
            
        
        # isolated
        plot1 = scatter(main_xs,
                        main_ys,
                        marker=symbol_isolated,
                        c=color_isolated,
                        s=markerSize,
                        edgecolor='black',
                        alpha=alpha_isolated,
                        label=label_isolated)
        
        if second:                
            # second data set
            plot1 = scatter(second_xs,
                            second_ys,
                            marker=symbol_second,
                            c=color_second,
                            s=markerSize,
                            edgecolor='black',
                            alpha=alpha_second,
                            label=label_second)

    
        # histogram isolated
        bin_means, edges, binNumber = stats.binned_statistic(main_xs,
                                                            main_ys,
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
        
        
        if second:
            # histogram two
            bin_means, edges, binNumber = stats.binned_statistic(second_xs,
                                                                second_ys,
                                                                statistic='median',
                                                                bins=bins)
            left,right = edges[:-1],edges[1:]        
            X = array([left,right]).T.flatten()
            Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
            plot(X,
                Y,
                ls='dashed',
                color=color_second,
                lw=2.0,
                alpha=alpha_bins,
                label=r'$\rm {0} ~Median ~EW$'.format(second))
        
    
        
        # x-axis
        majorLocator   = MultipleLocator(0.25)
        majorFormatter = FormatStrFormatter(r'$\rm %.2f$')
        minorLocator   = MultipleLocator(0.05)
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
        
        
        xlabel(r'$\rm \Sigma \mathcal{L} $')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        leg = ax.legend(scatterpoints=1,prop={'size':14},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0, 1)

        if plot_EW_summed_median_save:
            savefig('{0}/W(summed_ls)_median_plus_{1}_binSize{2}.pdf'.format(saveDirectory, second, binSize),format='pdf',bbox_inches='tight')
        else:
            show()

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

if __name__=="__main__":
    # do the work
    main()
    