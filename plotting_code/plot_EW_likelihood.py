#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_EW_likelihood.py, v 6.0 06/13/18

Plot EW as a function of likelihood and likelihood_cus

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
    plot_EW_likelihood_mean_plus = True
    plot_EW_likelihood_mean_plus_save = True

    # plot EW as a function of likelihood_cus, include median
    # EW histograms
    plot_EW_summed_likelihood_mean_plus = False
    plot_EW_summed_likelihood_mean_plus_save = False

    # plot_number = 1 for just the isolated sample, =2 adds the associated, =3 adds two+
    # =4 adds groups with 2 or more members
    plot_number = 2
    
    # which data set to use? Options are '', '_min001', '_cus','_min001_cus',
    # '_min001_double', '_min005_v150', '_min005_v250'
    data_set = '_min001'
    
    min_EW = 0
    max_EW = 1500
    

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
    L_summed_file = open(L_summed_filename,'r')


    # unload the data from them
    isolated = pickle.load(isolated_file)
    L_isolated = pickle.load(L_isolated_file)
    L_associated_isolated = pickle.load(L_associated_isolated_file)
    L_associated = pickle.load(L_associated_file)
    L_nonassociated = pickle.load(L_nonassociated_file)
    L_two = pickle.load(L_two_file)
    L_three_plus = pickle.load(L_three_plus_file)
    L_group = pickle.load(L_group_file)
    L_summed = pickle.load(L_summed_file)

    
    # close the files
    isolated_file.close()
    L_isolated_file.close()
    L_associated_isolated_file.close()
    L_associated_file.close()
    L_nonassociated_file.close()
    L_two_file.close()
    L_three_plus_file.close()
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
    
    if plot_EW_likelihood_mean_plus:
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
        

        alpha_isolated = 0.6
        alpha_assoc = 0.6
        alpha_two = 0.6
        alpha_three = 0.6
        alpha_group = 0.6
        alpha_bins = 0.99
        markerSize = 20
        
#         binSize = 50
#         bins = arange(0, 550, binSize)
        binSize = 0.25
        bins = arange(0, 1.25, binSize)
        
        bins = np.logspace(-3, 0, num=6, endpoint=True)

        
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
#         color_two = color_purple2
        color_two = color_red
        color_group = color_blue
        

        # define the x and y data for the isolated set
        
        Lya_Ws2 = []
        R_virs2 = []
        impacts2 = []
        bs2 = []
        ls2 = []
        for w, r, i, b, l in zip(Lya_Ws, R_virs, impacts, bs, ls):
            if w <= max_EW and w >= min_EW:
                Lya_Ws2.append(w)
                R_virs2.append(r)
                impacts2.append(i)
                bs2.append(b)
                ls2.append(l)
        
        isolated_xs = np.array(ls2)
        isolated_ys = np.array(Lya_Ws2)
        
        

        # grab the associated data and define the x and y data
        associated_Lya_Ws = L_associated['Lya_Ws']
        associated_R_virs = L_associated['R_virs']
        associated_impacts = L_associated['impacts']
        associated_bs = L_associated['bs']
        associated_ls = L_associated['ls']
        
        associated_Lya_Ws2 = []
        associated_R_virs2 = []
        associated_impacts2 = []
        associated_bs2 = []
        associated_ls2 = []
        for w, r, i, b, l in zip(associated_Lya_Ws, associated_R_virs, associated_impacts, associated_bs, associated_ls):
            if w <= max_EW and w >= min_EW:
                associated_Lya_Ws2.append(w)
                associated_R_virs2.append(r)
                associated_impacts2.append(i)
                associated_bs2.append(b)
                associated_ls2.append(l)
        
        associated_xs = np.array(associated_ls2)
        associated_ys = np.array(associated_Lya_Ws2)
        
        
        # grab the two data and define the x and y data
        two_Lya_Ws = L_two['Lya_Ws']
        two_R_virs = L_two['R_virs']
        two_impacts = L_two['impacts']
        two_bs = L_two['bs']
        two_ls = L_two['ls']


        two_Lya_Ws2 = []
        two_R_virs2 = []
        two_impacts2 = []
        two_bs2 = []
        two_ls2 = []
        for w, r, i, b, l in zip(two_Lya_Ws, two_R_virs, two_impacts, two_bs, two_ls):
            if w <= max_EW and w >= min_EW:
                two_Lya_Ws2.append(w)
                two_R_virs2.append(r)
                two_impacts2.append(i)
                two_bs2.append(b)
                two_ls2.append(l)

        two_xs = np.array(two_ls2)
        two_ys = np.array(two_Lya_Ws2)
        
        
        # grab the three_plus data and define the x and y data
        three_Lya_Ws = L_three_plus['Lya_Ws']
        three_R_virs = L_three_plus['R_virs']
        three_impacts = L_three_plus['impacts']
        three_bs = L_three_plus['bs']
        three_ls = L_three_plus['ls']

        three_Lya_Ws2 = []
        three_R_virs2 = []
        three_impacts2 = []
        three_bs2 = []
        three_ls2 = []
        for w, r, i, b, l in zip(three_Lya_Ws, three_R_virs, three_impacts, three_bs, three_ls):
            if w <= max_EW and w >= min_EW:
                three_Lya_Ws2.append(w)
                three_R_virs2.append(r)
                three_impacts2.append(i)
                three_bs2.append(b)
                three_ls2.append(l)
        
        three_xs = np.array(three_ls2)
        three_ys = np.array(three_Lya_Ws2)
        
        
        
        # grab the group data and define the x and y data
        group_Lya_Ws = L_group['Lya_Ws']
        group_R_virs = L_group['R_virs']
        group_impacts = L_group['impacts']
        group_mems = L_group['group_mems']
        group_bs = L_group['bs']
        group_ls = L_group['ls']

        group_Lya_Ws2 = []
        group_R_virs2 = []
        group_impacts2 = []
        group_bs2 = []
        group_ls2 = []
        for w, r, i, group, b, l in zip(group_Lya_Ws, group_R_virs, group_impacts, group_mems, group_bs, group_ls):
            if float(group) >= 2 and w <= max_EW and w >= min_EW:
                group_Lya_Ws2.append(w)
                group_R_virs2.append(r)
                group_impacts2.append(i)
                group_bs2.append(b)
                group_ls2.append(l)

        group_xs = np.array(group_ls2)
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
            label=r'$\rm Isolated~ Mean ~EW$')
            
        ax.set_xscale("log")

            
        # associated
        if plot_number >=2:
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
                label=r'$\rm Assoc. ~Mean ~EW$')

           
        # two
        if plot_number >=3:
            two_plus_xs = np.array(list(two_xs) + list(three_xs))
            two_plus_ys = np.array(list(two_ys) + list(three_ys))
            
            plot1 = scatter(two_plus_xs,
                            two_plus_ys,
                            marker=symbol_two,
                            c=color_two,
                            s=markerSize,
                            edgecolor='black',
                            alpha=alpha_two,
                            label=label_two)
        
            # histogram two
            bin_means, edges, binNumber = stats.binned_statistic(two_plus_xs,
                                                                two_plus_ys,
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
                label=r'$\rm Two+ ~Mean ~EW$')
           
                        
#         # group
        if plot_number ==4:
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
#         majorLocator   = MultipleLocator(0.2)
#         majorFormatter = FormatStrFormatter(r'$\rm %.2f$')
#         minorLocator   = MultipleLocator(0.1)
#         ax.xaxis.set_major_locator(majorLocator)
#         ax.xaxis.set_major_formatter(majorFormatter)
#         ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        
        xlabel(r'$\rm \mathcal{L}$')
        ylabel(r'$\rm EW~ [m \AA]$')
        leg = ax.legend(scatterpoints=1,prop={'size':12},loc='upper left',fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        # grid
        ax.grid(b=None,which='major',axis='both')
        ax.set_axisbelow(True)
        ax.yaxis.grid(color='gray', linestyle='dashed')
        
        # axis limits
        ylim(0,1500)
#         xlim(0, 2.5)

        if plot_EW_likelihood_mean_plus_save:
            savefig('{0}/W(likelihood)_binSize{1}_cut{2}-{3}_plus{4}_dataset{5}.pdf'.format(saveDirectory, binSize, min_EW, max_EW, plot_number, data_set),format='pdf',bbox_inches='tight')
        else:
            show()

##########################################################################################
##########################################################################################


##########################################################################################
##########################################################################################
    
    if plot_EW_summed_likelihood_mean_plus:
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
        bins = arange(0, 3.5, binSize)

        
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
        bs2 = []
        ls2 = []
        for w, r, i, b, l in zip(Lya_Ws, R_virs, impacts, bs, ls):
            if float(w) <= max_EW and float(w) >= min_EW:
                Lya_Ws2.append(w)
                R_virs2.append(r)
                impacts2.append(i)
                bs2.append(b)
                ls2.append(l)
        
        isolated_xs = np.array(ls2)
        isolated_ys = np.array(Lya_Ws2)
        
        

        # grab the associated data and define the x and y data
        associated_Lya_Ws = L_associated['Lya_Ws']
        associated_R_virs = L_associated['R_virs']
        associated_impacts = L_associated['impacts']
        associated_bs = L_associated['bs']
        associated_ls = L_associated['ls']
        
        associated_Lya_Ws2 = []
        associated_R_virs2 = []
        associated_impacts2 = []
        associated_bs2 = []
        associated_ls2 = []
        for w, r, i, b, l in zip(associated_Lya_Ws, associated_R_virs, associated_impacts, associated_bs, associated_ls):
            if float(w) <= max_EW and float(w) >= min_EW:
                associated_Lya_Ws2.append(w)
                associated_R_virs2.append(r)
                associated_impacts2.append(i)
                associated_bs2.append(b)
                associated_ls2.append(l)
        
        associated_xs = np.array(associated_ls2)
        associated_ys = np.array(associated_Lya_Ws2)
        
        
        # grab the two data and define the x and y data
        two_Lya_Ws = L_two['Lya_Ws']
        two_R_virs = L_two['R_virs']
        two_impacts = L_two['impacts']
        two_bs = L_two['bs']
        two_ls = L_two['ls']


        two_Lya_Ws2 = []
        two_R_virs2 = []
        two_impacts2 = []
        two_bs2 = []
        two_ls2 = []
        for w, r, i, b, l in zip(two_Lya_Ws, two_R_virs, two_impacts, two_bs, two_ls):
            if float(w) <= max_EW and float(w) >= min_EW:
                two_Lya_Ws2.append(w)
                two_R_virs2.append(r)
                two_impacts2.append(i)
                two_bs2.append(b)
                two_ls2.append(l)

        two_xs = np.array(two_ls2)
        two_ys = np.array(two_Lya_Ws2)
        
        
        # grab the three_plus data and define the x and y data
        three_Lya_Ws = L_three_plus['Lya_Ws']
        three_R_virs = L_three_plus['R_virs']
        three_impacts = L_three_plus['impacts']
        three_bs = L_three_plus['bs']
        three_ls = L_three_plus['ls']

        three_Lya_Ws2 = []
        three_R_virs2 = []
        three_impacts2 = []
        three_bs2 = []
        three_ls2 = []
        for w, r, i, b, l in zip(three_Lya_Ws, three_R_virs, three_impacts, three_bs, three_ls):
            if float(w) <= max_EW and float(w) >= min_EW:
                three_Lya_Ws2.append(w)
                three_R_virs2.append(r)
                three_impacts2.append(i)
                three_bs2.append(b)
                three_ls2.append(l)
        
        three_xs = np.array(three_ls2)
        three_ys = np.array(three_Lya_Ws2)
        
        
        
        # grab the group data and define the x and y data
        group_Lya_Ws = L_group['Lya_Ws']
        group_R_virs = L_group['R_virs']
        group_impacts = L_group['impacts']
        group_mems = L_group['group_mems']
        group_bs = L_group['bs']
        group_ls = L_group['ls']

        group_Lya_Ws2 = []
        group_R_virs2 = []
        group_impacts2 = []
        group_bs2 = []
        group_ls2 = []
        for w, r, i, group, b, l in zip(group_Lya_Ws, group_R_virs, group_impacts, group_mems, group_bs, group_ls):
            if float(group) >= 2 and float(w) <= max_EW and float(w) >= min_EW:
                group_Lya_Ws2.append(w)
                group_R_virs2.append(r)
                group_impacts2.append(i)
                group_bs2.append(b)
                group_bs2.append(l)

        group_xs = np.array(group_ls2)
        group_ys = np.array(group_Lya_Ws2)

                
        
        # define the x and y data for the summed set
        summed_ls = L_summed['summed_ls']
        summed_bs = L_summed['bs']
        summed_Lya_Ws = L_summed['Lya_Ws']
        summed_l_cuss = L_summed['summed_l_cuss']
        
        summed_ls2 = []
        summed_bs2 = []
        summed_Lya_Ws2 = []
        summed_l_cuss2 = []
        for b, l, w , l_cuss in zip(summed_bs, summed_ls, summed_Lya_Ws, summed_l_cuss):
            if float(w) <= max_EW and float(w) >= min_EW:
                summed_ls2.append(np.sum(np.array(l)))
                summed_bs2.append(b)
                summed_Lya_Ws2.append(w)
                summed_l_cuss2.append(np.sum(np.array(l_cuss)))

        main_xs = np.array(summed_l_cuss2)
        main_ys = np.array(summed_Lya_Ws2)


##########################################################################################
        # do the plotting 
        
        # summed
        plot1 = scatter(main_xs,
                        main_ys,
                        marker=symbol_isolated,
                        c=color_isolated,
                        s=markerSize,
                        edgecolor='black',
                        alpha=alpha_isolated)
                        
                        
        # histogram summed
        bin_means, edges, binNumber = stats.binned_statistic(main_xs,
                                                            main_ys,
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
            label=r'$\rm Mean ~EW$')
            
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
        
        
        xlabel(r'$\rm \Sigma \small \mathcal{L_{{cus}}}$')
        ylabel(r'$\rm EW~ [m \AA]$')
#         leg = ax.legend(scatterpoints=1,prop={'size':12},loc=1,fancybox=True)

        # grid lines
        ax.grid(b=None,which='major',axis='both')
        ax.set_axisbelow(True)
        ax.yaxis.grid(color='gray', linestyle='dashed')

        # axis limits
        ylim(0,1000)
#         xlim(0, 1.5)

        if plot_EW_summed_likelihood_mean_plus_save:
            savefig('{0}/W(summed_likelihood_cus)_binSize{1}_cut{2}-{3}_dataset{4}.pdf'.format(saveDirectory, binSize, min_EW, max_EW, data_set),format='pdf',bbox_inches='tight')
        else:
            show()

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

if __name__=="__main__":
    # do the work
    main()
    