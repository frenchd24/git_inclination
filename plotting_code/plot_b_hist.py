#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_b_hist.py, v 1.0 06/20/18

Plot b histograms and CDFs


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


def errors(a):
    # return the standard error in the mean for the input array
    return stats.sem(a)

    

def main():
    # plot CDF of bs
    plot_b_cdf = True
    plot_b_cdf_save = True

    # plot CDF of bs - old way of doing this
    plot_b_cdf_old = False
    plot_b_cdf_old_save = False

    # plot histogram of b
    plot_b_hist = False
    plot_b_hist_save = False
    
    # plot histogram of b in 3 bins of azimuth
    plot_b_hist_az = False
    plot_b_hist_az_save = False
    
    # plot b histograms as a function of number of group members
    plot_b_hist_group = False
    plot_b_hist_group_save = False
    
    # plot b histograms as a function of MType
    plot_b_hist_MType = False
    plot_b_hist_MType_save = False

    # which data set to use?
    data_set = '_double'
    
    plot_errors = True

    min_EW = 50
    max_EW = 10000
    
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
    all_file = open(all_filename,'r')


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
    all = pickle.load(all_file)

    
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
    all_file.close()

    
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
    
    if plot_b_cdf:
        
        countb = 0
        countr = 0
        count = -1
        
        second = 'Associated'
        
        color_green = '#1b9e77'
        color_purple = '#7570b3'
        color_orange = '#d95f02'
        color_purple2 = '#984ea3'
        
        
        color_green = '#1b9e77'
        color_orange = '#d95f02'
        color_purple3 = '#7570b3'
        color_pink = '#e7298a'
        color_lime = '#66a61e'
        color_yellow = '#e6ab02'
        color_brown = '#a6761d'
        color_coal = '#666666'
        
        color_blue = '#377eb8'

        alpha_green = 0.99
        alpha_orange = 0.99
        alpha_purple = 0.99
        alpha_pink = 0.99
        alpha_lime = 0.99
        alpha_yellow = 0.99
        alpha_brown = 0.99
        alpha_coal = 0.99
        alpha_blue = 0.99
        alpha_red = 0.99
        alpha_black = 0.99
        
        alpha_error = 0.4
        
        binSize = 1
        bins = arange(0, 120+binSize, binSize)
#         bins = 100
        
        
        label_isolated = r'$Isolated$'
        label_associated = r'$\rm Associated$'
        label_two = r'$\rm Two$'
        label_three = r'$\rm Three+$'
        label_group = r'$\rm Group$'
        label_second = r'$\rm {0}$'.format(second)


        color_isolated = 'black'
        color_associated = color_red
        color_two = color_red
        color_three = color_purple
        color_group = color_orange
        color_second = color_orange
        
        
        # grab the full set
        all_bs = all['bs']
        all_e_bs = all['e_bs']
        all_Lya_Ws = all['Lya_Ws']
        all_e_Lya_Ws = all['e_Lya_Ws']
        
        all_bs_cut = []
        all_e_bs_cut = []
        for w, e, b, e_b in zip(all_Lya_Ws, all_e_Lya_Ws, all_bs, all_e_bs):
            if w <= max_EW and w >= min_EW:
                all_bs_cut.append(b)
                all_e_bs_cut.append(e_b)

                
        # isolated
        isolated_bs = isolated['bs']
        isolated_e_bs = isolated['e_bs']
        isolated_Lya_Ws = isolated['Lya_Ws']
        isolated_e_Lya_Ws = isolated['e_Lya_Ws']
        
        isolated_bs_cut = []
        isolated_e_bs_cut = []
        for w, e, b, e_b in zip(isolated_Lya_Ws, isolated_e_Lya_Ws, isolated_bs, isolated_e_bs):
            if w <= max_EW and w >= min_EW:
                isolated_bs_cut.append(b)
                isolated_e_bs_cut.append(e_b)


        # L_isolated
        L_isolated_bs = L_isolated['bs']
        L_isolated_e_bs = L_isolated['e_bs']
        L_isolated_Lya_Ws = L_isolated['Lya_Ws']
        L_isolated_e_Lya_Ws = L_isolated['e_Lya_Ws']
        
        L_isolated_bs_cut = []
        L_isolated_e_bs_cut = []
        for w, e, b, e_b in zip(L_isolated_Lya_Ws, L_isolated_e_Lya_Ws, L_isolated_bs, L_isolated_e_bs):
            if w <= max_EW and w >= min_EW:
                L_isolated_bs_cut.append(b)
                L_isolated_e_bs_cut.append(e_b)


        # L_nonassociated
        L_nonassociated_bs = L_nonassociated['bs']
        L_nonassociated_e_bs = L_nonassociated['e_bs']
        L_nonassociated_Lya_Ws = L_nonassociated['Lya_Ws']
        L_nonassociated_e_Lya_Ws = L_nonassociated['e_Lya_Ws']
        
        
        L_nonassociated_bs_cut = []
        L_nonassociated_e_bs_cut = []
        for w, e, b, e_b in zip(L_nonassociated_Lya_Ws, L_nonassociated_e_Lya_Ws, L_nonassociated_bs, L_nonassociated_e_bs):
            if w <= max_EW and w >= min_EW:
                L_nonassociated_bs_cut.append(b)
                L_nonassociated_e_bs_cut.append(e_b)


        # L_isolate_associated
        L_isolated_associated_bs_cut = []
        L_isolated_associated_e_bs_cut = []
        for w, e, b, e_b in zip(Lya_Ws, e_Lya_Ws, bs, e_bs):
            if w <= max_EW and w >= min_EW:
                L_isolated_associated_bs_cut.append(b)
                L_isolated_associated_e_bs_cut.append(e_b)


        # grab the associated data 
        L_associated_bs = L_associated['bs']
        L_associated_e_bs = L_associated['e_bs']
        L_associated_Lya_Ws = L_associated['Lya_Ws']
        L_associated_e_Lya_Ws = L_associated['e_Lya_Ws']
        L_associated_R_virs = L_associated['R_virs']
        L_associated_impacts = L_associated['impacts']
        L_associated_MajDiams = L_associated['MajDiams']
        L_associated_ls = L_associated['ls']
        
        L_associated_bs_cut = []
        L_associated_e_bs_cut = []
        for w, e, b, e_b in zip(L_associated_Lya_Ws, L_associated_e_Lya_Ws, L_associated_bs, L_associated_e_bs):
            if w <= max_EW and w >= min_EW:
                L_associated_bs_cut.append(b)
                L_associated_e_bs_cut.append(e_b)


        # grab the two data 
        L_two_bs = L_two['bs']
        L_two_e_bs = L_two['e_bs']
        L_two_Lya_Ws = L_two['Lya_Ws']
        L_two_e_Lya_Ws = L_two['e_Lya_Ws']
        L_two_R_virs = L_two['R_virs']
        L_two_impacts = L_two['impacts']
        L_two_MajDiams = L_two['MajDiams']
        L_two_ls = L_two['ls']
        
        L_two_bs_cut = []
        L_two_e_bs_cut = []
        for w, e, b, e_b in zip(L_two_Lya_Ws, L_two_e_Lya_Ws, L_two_bs, L_two_e_bs):
            if w <= max_EW and w >= min_EW:
                L_two_bs_cut.append(b)
                L_two_e_bs_cut.append(e_b)

        
        # grab the three_plus data 
        L_three_bs = L_three_plus['bs']
        L_three_e_bs = L_three_plus['e_bs']
        L_three_Lya_Ws = L_three_plus['Lya_Ws']
        L_three_e_Lya_Ws = L_three_plus['e_Lya_Ws']
        L_three_R_virs = L_three_plus['R_virs']
        L_three_impacts = L_three_plus['impacts']
        L_three_MajDiams = L_three_plus['MajDiams']
        L_three_ls = L_three_plus['ls']
        
        L_three_bs_cut = []
        L_three_e_bs_cut = []
        for w, e, b, e_b in zip(L_three_Lya_Ws, L_three_e_Lya_Ws, L_three_bs, L_three_e_bs):
            if w <= max_EW and w >= min_EW:
                L_three_bs_cut.append(b)
                L_three_e_bs_cut.append(e_b)

        
        # grab the group data and define the x and y data
        L_group_bs = L_group['bs']
        L_group_e_bs = L_group['e_bs']
        L_group_Lya_Ws = L_group['Lya_Ws']
        L_group_e_Lya_Ws = L_group['e_Lya_Ws']
        L_group_R_virs = L_group['R_virs']
        L_group_impacts = L_group['impacts']
        L_group_MajDiams = L_group['MajDiams']
        L_group_mems = L_group['group_mems']
        L_group_ls = L_group['ls']
        
        L_group_bs_cut = []
        L_group_e_bs_cut = []
        for w, m, e, b, e_b in zip(L_group_Lya_Ws, L_group_e_Lya_Ws, L_group_bs, L_group_mems, L_group_e_bs):
            if w <= max_EW and w >= min_EW and int(m) >= 2:
                L_group_bs_cut.append(b)
                L_group_e_bs_cut.append(e_b)

##########################################################################################
        bins_right = bins[1:]                           

        count_list = [1,2,3,4,5,6]
        
        for count in count_list:
            print 'Plotting number {0}'.format(count)
            
            fig = figure(figsize=(7.7,5.7))
            ax = fig.add_subplot(111)
        
        
            # all EWs
            hist(all_bs_cut,
            bins=bins,
            histtype='step',
            cumulative=True,
            normed=1,
            color='black',
            lw=1.5,
            alpha=alpha_black,
            label=r'$ All$')
         
            if plot_errors:
                min_n, min_bins, min_p = hist(np.array(all_bs_cut) - np.array(all_e_bs_cut),
                                                bins=bins,
                                                histtype='step',
                                                cumulative=True,
                                                normed=1,
                                                color='black',
                                                lw=1.5,
                                                alpha=0)
                                        
                max_n, max_bins, max_p = hist(np.array(all_bs_cut) + np.array(all_e_bs_cut),
                                                bins=bins,
                                                histtype='step',
                                                cumulative=True,
                                                normed=1,
                                                color='black',
                                                lw=1.5,
                                                alpha=0)           

                fill_between(bins_right, min_n, max_n, facecolor='black', alpha=alpha_error)
        
        
        
            # isolated
            if count >=2:
                hist(isolated_bs_cut,
                bins=bins,
                histtype='step',
                cumulative=True,
                normed=1,
                color='grey',
                lw=1.5,
                alpha=alpha_black,
                label=r'$ Isolated$')

                if plot_errors:
                    min_n, min_bins, min_p = hist(np.array(isolated_bs_cut) - np.array(isolated_e_bs_cut),
                                                    bins=bins,
                                                    histtype='step',
                                                    cumulative=True,
                                                    normed=1,
                                                    color='grey',
                                                    lw=1.5,
                                                    alpha=0)
                                        
                    max_n, max_bins, max_p = hist(np.array(isolated_bs_cut) + np.array(isolated_e_bs_cut),
                                                    bins=bins,
                                                    histtype='step',
                                                    cumulative=True,
                                                    normed=1,
                                                    color='grey',
                                                    lw=1.5,
                                                    alpha=0)           

                    fill_between(bins_right, min_n, max_n, facecolor='grey', alpha=alpha_error)
        
        
        
            # L_isolated
            if count >=3:
                hist(L_isolated_bs_cut,
                bins=bins,
                histtype='step',
                cumulative=True,
                normed=1,
                color=color_brown,
                lw=1.5,
                alpha=alpha_brown,
                label=r'$ \mathcal{L}-isolated$')
        
        
                if plot_errors:
                    min_n, min_bins, min_p = hist(np.array(L_isolated_bs_cut) - np.array(L_isolated_e_bs_cut),
                                                    bins=bins,
                                                    histtype='step',
                                                    cumulative=True,
                                                    normed=1,
                                                    color=color_brown,
                                                    lw=1.5,
                                                    alpha=0)
                                        
                    max_n, max_bins, max_p = hist(np.array(L_isolated_bs_cut) + np.array(L_isolated_e_bs_cut),
                                                    bins=bins,
                                                    histtype='step',
                                                    cumulative=True,
                                                    normed=1,
                                                    color=color_brown,
                                                    lw=1.5,
                                                    alpha=0)           

                    fill_between(bins_right, min_n, max_n, facecolor=color_brown, alpha=alpha_error)
        
        
        

            # isolated associated
            if count >=4:
                hist(L_isolated_associated_bs_cut,
                bins=bins,
                histtype='step',
                cumulative=True,
                normed=1,
                color=color_green,
                lw=1.5,
                alpha=alpha_green,
                label=r'$ \mathcal{L}-associated-isolated$')
        
                if plot_errors:
                    min_n, min_bins, min_p = hist(np.array(L_isolated_associated_bs_cut) - np.array(L_isolated_associated_e_bs_cut),
                                                    bins=bins,
                                                    histtype='step',
                                                    cumulative=True,
                                                    normed=1,
                                                    color=color_green,
                                                    lw=1.5,
                                                    alpha=0)
                                        
                    max_n, max_bins, max_p = hist(np.array(L_isolated_associated_bs_cut) + np.array(L_isolated_associated_e_bs_cut),
                                                    bins=bins,
                                                    histtype='step',
                                                    cumulative=True,
                                                    normed=1,
                                                    color=color_green,
                                                    lw=1.5,
                                                    alpha=0)           

                    fill_between(bins_right, min_n, max_n, facecolor=color_green, alpha=alpha_error)
        
        
        
        
            # associated
            if count >=5:
                hist(L_associated_bs_cut,
                bins=bins,
                histtype='step',
                cumulative=True,
                normed=1,
                color=color_orange,
                lw=1.5,
                alpha=alpha_orange,
                label=r'$ \mathcal{L}-associated$')
                
                if plot_errors:
                    min_n, min_bins, min_p = hist(np.array(L_associated_bs_cut) - np.array(L_associated_e_bs_cut),
                                                    bins=bins,
                                                    histtype='step',
                                                    cumulative=True,
                                                    normed=1,
                                                    color=color_orange,
                                                    lw=1.5,
                                                    alpha=0)
                                        
                    max_n, max_bins, max_p = hist(np.array(L_associated_bs_cut) + np.array(L_associated_e_bs_cut),
                                                    bins=bins,
                                                    histtype='step',
                                                    cumulative=True,
                                                    normed=1,
                                                    color=color_orange,
                                                    lw=1.5,
                                                    alpha=0)           

                    fill_between(bins_right, min_n, max_n, facecolor=color_orange, alpha=alpha_error)
        


        # L_two+
#         L_twoplus_bs_cut = L_two_bs_cut + L_three_bs
#         L_twoplus_e_bs_cut = L_two_e_bs_cut + L_three_e_bs
# 
#         hist(L_twoplus_bs_cut,
#         bins=bins,
#         histtype='step',
#         cumulative=True,
#         normed=1,
#         color=color_purple3,
#         lw=1.5,
#         alpha=alpha_purple,
#         label=r'$\rm L\_two+$')
#         
#         min_n, min_bins, min_p = hist(np.array(L_twoplus_bs_cut) - np.array(L_twoplus_e_bs_cut),
#                                         bins=bins,
#                                         histtype='step',
#                                         cumulative=True,
#                                         normed=1,
#                                         color=color_purple3,
#                                         lw=1.5,
#                                         alpha=0)
#                                         
#         max_n, max_bins, max_p = hist(np.array(L_twoplus_bs_cut) + np.array(L_twoplus_e_bs_cut),
#                                         bins=bins,
#                                         histtype='step',
#                                         cumulative=True,
#                                         normed=1,
#                                         color=color_purple3,
#                                         lw=1.5,
#                                         alpha=0)           
# 
#         fill_between(bins_right, min_n, max_n, facecolor=color_purple3, alpha=0.5)
        
        
            # two
            if count >=6:
                two_plus_bs_cut = np.array(list(L_two_bs_cut) + list(L_three_bs_cut))
                two_plus_e_bs_cut = np.array(list(L_two_e_bs_cut) + list(L_three_e_bs_cut))

                hist(two_plus_bs_cut,
                bins=bins,
                histtype='step',
                cumulative=True,
                normed=1,
                color=color_purple3,
                lw=1.5,
                alpha=alpha_purple,
                label=r'$ \mathcal{L}-two+$')
        
                if plot_errors:
                    min_n, min_bins, min_p = hist(np.array(two_plus_bs_cut) - np.array(two_plus_e_bs_cut),
                                                    bins=bins,
                                                    histtype='step',
                                                    cumulative=True,
                                                    normed=1,
                                                    color=color_purple3,
                                                    lw=1.5,
                                                    alpha=0)
                                        
                    max_n, max_bins, max_p = hist(np.array(two_plus_bs_cut) + np.array(two_plus_e_bs_cut),
                                                    bins=bins,
                                                    histtype='step',
                                                    cumulative=True,
                                                    normed=1,
                                                    color=color_purple3,
                                                    lw=1.5,
                                                    alpha=0)

                    fill_between(bins_right, min_n, max_n, facecolor=color_purple3, alpha=alpha_error)
        
        
            # L_group
            if count >=7:
                hist(L_group_bs_cut,
                bins=bins,
                histtype='step',
                cumulative=True,
                normed=1,
                color=color_blue,
                lw=1.5,
                alpha=alpha_blue,
                label=r'$ \mathcal{L}-group$')
        
                if plot_errors:
                    min_n, min_bins, min_p = hist(np.array(L_group_bs_cut) - np.array(L_group_e_bs_cut),
                                                    bins=bins,
                                                    histtype='step',
                                                    cumulative=True,
                                                    normed=1,
                                                    color=color_blue,
                                                    lw=1.5,
                                                    alpha=0)
                                        
                    max_n, max_bins, max_p = hist(np.array(L_group_bs_cut) + np.array(L_group_e_bs_cut),
                                                    bins=bins,
                                                    histtype='step',
                                                    cumulative=True,
                                                    normed=1,
                                                    color=color_blue,
                                                    lw=1.5,
                                                    alpha=0)           

                    fill_between(bins_right, min_n, max_n, facecolor=color_blue, alpha=alpha_error)
#         

        
            # L_nonassociated
    #         hist(L_nonassociated_bs_cut,
    #         bins=bins,
    #         histtype='step',
    #         cumulative=True,
    #         normed=1,
    #         color=color_pink,
    #         lw=1.5,
    #         alpha=alpha_pink,
    #         label=r'$\rm L\_nonassoc.$')
        
        
            # x-axis
            majorLocator   = MultipleLocator(20)
            majorFormatter = FormatStrFormatter(r'$\rm %d$')
            minorLocator   = MultipleLocator(10)
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
        
            ylim(0,1)

            legend(scatterpoints=1, prop={'size':14}, loc='lower right', fancybox=True)
            xlabel(r'$\rm b ~[km s^{{-1}}]$')
            ylabel(r'$\rm CDF$')

            # grid
            ax.grid(b=None,which='major',axis='both')
            ax.set_axisbelow(True)
            ax.yaxis.grid(color='gray', linestyle='solid',alpha=0.5)
            ax.xaxis.grid(color='gray', linestyle='solid',alpha=0.5)


            xlim(0, 120)

            if plot_b_cdf_save:
                savefig('{0}/hist(b)_all{1}_bins{2}_EWcut{3}-{4}_err{5}_dataset{6}.pdf'.format(saveDirectory, count, binSize, min_EW, max_EW, plot_errors, data_set),format='pdf',bbox_inches='tight')
            else:
                show()
    
    
##########################################################################################
##########################################################################################
    
    if plot_b_cdf_old:
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
        
        color_green = '#1b9e77'
        color_orange = '#d95f02'
        color_purple3 = '#7570b3'
        color_pink = '#e7298a'
        color_lime = '#66a61e'
        color_yellow = '#e6ab02'
        color_brown = '#a6761d'
        color_coal = '#666666'


        alpha_green = 0.99
        alpha_orange = 0.99
        alpha_purple = 0.99
        alpha_pink = 0.99
        alpha_lime = 0.99
        alpha_yellow = 0.99
        alpha_brown = 0.99
        alpha_coal = 0.99
        alpha_blue = 0.99
        alpha_red = 0.99
        alpha_black = 0.99
        
        binSize = 1
        bins = arange(0, 121, binSize)
#         bins = 100
        
        
        label_isolated = r'$\rm Isolated$'
        label_associated = r'$\rm Associated$'
        label_two = r'$\rm Two$'
        label_three = r'$\rm Three+$'
        label_group = r'$\rm Group$'
        label_second = r'$\rm {0}$'.format(second)


        color_isolated = 'black'
        color_associated = color_red
        color_two = color_red
        color_three = color_purple
        color_group = color_orange
        color_second = color_orange
        
        
        # grab the full set
        all_Lya_Ws = all['Lya_Ws']
        all_bs = all['bs']

        all_Lya_Ws_cut = []
        all_bs_cut = []
        for w, b in zip(all_Lya_Ws, all_bs):
            if w <= max_EW and w >= min_EW:
                all_Lya_Ws_cut.append(w)
                all_bs_cut.append(b)
                
        # isolated
        isolated_Lya_Ws = isolated['Lya_Ws']
        isolated_bs = isolated['bs']

        isolated_Lya_Ws_cut = []
        isolated_bs_cut = []
        for w, b in zip(isolated_Lya_Ws, isolated_bs):
            if w <= max_EW and w >= min_EW:
                isolated_Lya_Ws_cut.append(w)
                isolated_bs_cut.append(b)

                
        # L_isolated
        L_isolated_Lya_Ws = L_isolated['Lya_Ws']
        L_isolated_bs = L_isolated['bs']

        L_isolated_Lya_Ws_cut = []
        L_isolated_bs_cut = []
        for w, b in zip(L_isolated_Lya_Ws, L_isolated_bs):
            if w <= max_EW and w >= min_EW:
                L_isolated_Lya_Ws_cut.append(w)
                L_isolated_bs_cut.append(b)


        # L_nonassociated
        L_nonassociated_Lya_Ws = L_nonassociated['Lya_Ws']
        L_nonassociated_bs = L_nonassociated['bs']

        L_nonassociated_Lya_Ws_cut = []
        L_nonassociated_bs_cut = []
        for w, b in zip(L_nonassociated_Lya_Ws, L_nonassociated_bs):
            if w <= max_EW and w >= min_EW:
                L_nonassociated_Lya_Ws_cut.append(w)
                L_nonassociated_bs_cut.append(b)


        # L_isolate_associated
        L_isolated_associated_Lya_Ws_cut = []
        L_isolated_associated_bs_cut = []
        for w, b in zip(Lya_Ws, bs):
            if w <= max_EW and w >= min_EW:
                L_isolated_associated_Lya_Ws_cut.append(w)
                L_isolated_associated_bs_cut.append(b)

        # grab the associated data 
        L_associated_Lya_Ws = L_associated['Lya_Ws']
        L_associated_R_virs = L_associated['R_virs']
        L_associated_impacts = L_associated['impacts']
        L_associated_MajDiams = L_associated['MajDiams']
        L_associated_ls = L_associated['ls']
        L_associated_bs = L_associated['bs']
        
        L_associated_Lya_Ws_cut = []
        L_associated_bs_cut = []
        for w, b in zip(L_associated_Lya_Ws, L_associated_bs):
            if w <= max_EW and w >= min_EW:
                L_associated_Lya_Ws_cut.append(w)
                L_associated_bs_cut.append(b)
        
        
        # grab the two data 
        L_two_Lya_Ws = L_two['Lya_Ws']
        L_two_R_virs = L_two['R_virs']
        L_two_impacts = L_two['impacts']
        L_two_MajDiams = L_two['MajDiams']
        L_two_ls = L_two['ls']
        L_two_bs = L_two['bs']
        
        L_two_Lya_Ws_cut = []
        L_two_bs_cut = []
        for w, b in zip(L_two_Lya_Ws, L_two_bs):
            if w <= max_EW and w >= min_EW:
                L_two_Lya_Ws_cut.append(w)
                L_two_bs_cut.append(b)

        
        # grab the three_plus data 
        L_three_Lya_Ws = L_three_plus['Lya_Ws']
        L_three_R_virs = L_three_plus['R_virs']
        L_three_impacts = L_three_plus['impacts']
        L_three_MajDiams = L_three_plus['MajDiams']
        L_three_ls = L_three_plus['ls']
        L_three_bs = L_three_plus['bs']
        
        L_three_Lya_Ws_cut = []
        L_three_bs_cut = []
        for w, b in zip(L_three_Lya_Ws, L_three_bs):
            if w <= max_EW and w >= min_EW:
                L_three_Lya_Ws_cut.append(w)
                L_three_bs_cut.append(b)


        # grab the group data and define the x and y data
        L_group_Lya_Ws = L_group['Lya_Ws']
        L_group_R_virs = L_group['R_virs']
        L_group_impacts = L_group['impacts']
        L_group_MajDiams = L_group['MajDiams']
        L_group_mems = L_group['group_mems']
        L_group_ls = L_group['ls']
        L_group_bs = L_group['bs']
        
        L_group_Lya_Ws_cut = []
        L_group_bs_cut = []
        for w, m, b in zip(L_group_Lya_Ws, L_group_mems, L_group_bs):
            if w <= max_EW and w >= min_EW and int(m) >= 2:
                L_group_Lya_Ws_cut.append(w)
                L_group_bs_cut.append(b)

                
        # all EWs
        hist(all_bs_cut,
        bins=bins,
        histtype='step',
        cumulative=True,
        normed=1,
        color='black',
        lw=1.5,
        alpha=alpha_black,
        label=r'$\rm All$')
        
        
        # isolated 
#         hist(isolated_bs_cut,
#         bins=bins,
#         histtype='step',
#         cumulative=True,
#         normed=1,
#         color='grey',
#         lw=1.5,
#         alpha=alpha_black,
#         label=r'$\rm Isolated$')
        
        
        # L_isolated 
#         hist(L_isolated_bs_cut,
#         bins=bins,
#         histtype='step',
#         cumulative=True,
#         normed=1,
#         color=color_brown,
#         lw=1.5,
#         alpha=alpha_brown,
#         label=r'$\rm L\_isolated$')
        

        # isolated associated
#         hist(L_isolated_associated_bs_cut,
#         bins=bins,
#         histtype='step',
#         cumulative=True,
#         normed=1,
#         color=color_green,
#         lw=1.5,
#         alpha=alpha_green,
#         label=r'$\rm L\_isolated\_assoc.$')
        
        
        # associated
#         hist(L_associated_bs_cut,
#         bins=bins,
#         histtype='step',
#         cumulative=True,
#         normed=1,
#         color=color_orange,
#         lw=1.5,
#         alpha=alpha_orange,
#         label=r'$\rm L\_assoc.$')


        # L_two
#         hist(L_two_bs_cut,
#         bins=bins,
#         histtype='step',
#         cumulative=True,
#         normed=1,
#         color=color_purple3,
#         lw=1.5,
#         alpha=alpha_purple,
#         label=r'$\rm L\_two$')
        
        
        
        # L_nonassociated
#         hist(L_nonassociated_Lya_Ws_cut,
#         bins=bins,
#         histtype='step',
#         cumulative=True,
#         normed=1,
#         color=color_pink,
#         lw=1.5,
#         alpha=alpha_pink,
#         label=r'$\rm L\_nonassoc.$')
        
        
        # L_group
#         hist(L_group_Lya_Ws_cut,
#         bins=bins,
#         histtype='step',
#         cumulative=True,
#         normed=1,
#         color=color_blue,
#         lw=1.5,
#         alpha=alpha_yellow,
#         label=r'$\rm L\_group$')
        
        
#         color_green = '#1b9e77'
#         color_orange = '#d95f02'
#         color_purple3 = '#7570b3'
#         color_pink = '#e7298a'
#         color_lime = '#66a61e'
#         color_yellow = '#e6ab02'
#         color_brown = '#a6761d'
#         color_coal = '#666666'
        
        # x-axis
        majorLocator   = MultipleLocator(20)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(10)
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
        
        ylim(0,1)
        

        legend(scatterpoints=1, prop={'size':14}, loc='lower right', fancybox=True)
        xlabel(r'$\rm b ~[km s^{{-1}}]$')
        ylabel(r'$\rm CDF$')

        ax.grid(b=None,which='major',axis='both')
#         ylim(0,1300)
        xlim(20, 120)

        if plot_b_cdf_old_save:
            savefig('{0}/hist(b)_all6_bins{1}_minmaxEW_{2}-{3}_dataset{4}.pdf'.format(saveDirectory, binSize, min_EW, max_EW, data_set),format='pdf',bbox_inches='tight')
        else:
            show()


##########################################################################################
##########################################################################################

    
##########################################################################################
##########################################################################################
    
    if plot_b_hist:
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
        
        
        color_green = '#1b9e77'
        color_orange = '#d95f02'
        color_purple3 = '#7570b3'
        color_pink = '#e7298a'
        color_lime = '#66a61e'
        color_yellow = '#e6ab02'
        color_brown = '#a6761d'
        color_coal = '#666666'


        alpha_green = 0.7
        alpha_orange = 0.7
        alpha_purple = 0.7
        alpha_pink = 0.7
        alpha_lime = 0.7
        alpha_yellow = 0.7
        alpha_brown = 0.7
        alpha_coal = 0.7
        alpha_blue = 0.7
        alpha_red = 0.7
        alpha_black = 0.7
        
        binSize = 10
        bins = arange(0, 200+binSize, binSize)
        
        label_isolated = r'$\rm Isolated$'
        label_associated = r'$\rm Associated$'
        label_two = r'$\rm Two$'
        label_three = r'$\rm Three+$'
        label_group = r'$\rm Group$'
        label_second = r'$\rm {0}$'.format(second)


        color_isolated = 'black'
        color_associated = color_red
        color_two = color_red
        color_three = color_purple
        color_group = color_orange
        color_second = color_orange
        
        # grab the full set
        all_Lya_Ws = all['Lya_Ws']
        all_bs = all['bs']

        all_Lya_Ws_cut = []
        all_bs_cut = []
        for w, b in zip(all_Lya_Ws, all_bs):
            if w <= max_EW and w >= min_EW:
                all_Lya_Ws_cut.append(w)
                all_bs_cut.append(b)
                
        # isolated
        isolated_Lya_Ws = isolated['Lya_Ws']
        isolated_bs = isolated['bs']

        isolated_Lya_Ws_cut = []
        isolated_bs_cut = []
        for w, b in zip(isolated_Lya_Ws, isolated_bs):
            if w <= max_EW and w >= min_EW:
                isolated_Lya_Ws_cut.append(w)
                isolated_bs_cut.append(b)

                
        # L_isolated
        L_isolated_Lya_Ws = L_isolated['Lya_Ws']
        L_isolated_bs = L_isolated['bs']

        L_isolated_Lya_Ws_cut = []
        L_isolated_bs_cut = []
        for w, b in zip(L_isolated_Lya_Ws, L_isolated_bs):
            if w <= max_EW and w >= min_EW:
                L_isolated_Lya_Ws_cut.append(w)
                L_isolated_bs_cut.append(b)


        # L_nonassociated
        L_nonassociated_Lya_Ws = L_nonassociated['Lya_Ws']
        L_nonassociated_bs = L_nonassociated['bs']

        L_nonassociated_Lya_Ws_cut = []
        L_nonassociated_bs_cut = []
        for w, b in zip(L_nonassociated_Lya_Ws, L_nonassociated_bs):
            if w <= max_EW and w >= min_EW:
                L_nonassociated_Lya_Ws_cut.append(w)
                L_nonassociated_bs_cut.append(b)


        # L_isolate_associated
        L_isolated_associated_Lya_Ws_cut = []
        L_isolated_associated_bs_cut = []
        for w, b in zip(Lya_Ws, bs):
            if w <= max_EW and w >= min_EW:
                L_isolated_associated_Lya_Ws_cut.append(w)
                L_isolated_associated_bs_cut.append(b)

        # grab the associated data 
        L_associated_Lya_Ws = L_associated['Lya_Ws']
        L_associated_R_virs = L_associated['R_virs']
        L_associated_impacts = L_associated['impacts']
        L_associated_MajDiams = L_associated['MajDiams']
        L_associated_ls = L_associated['ls']
        L_associated_bs = L_associated['bs']
        
        L_associated_Lya_Ws_cut = []
        L_associated_bs_cut = []
        for w, b in zip(L_associated_Lya_Ws, L_associated_bs):
            if w <= max_EW and w >= min_EW:
                L_associated_Lya_Ws_cut.append(w)
                L_associated_bs_cut.append(b)
        
        
        # grab the two data 
        L_two_Lya_Ws = L_two['Lya_Ws']
        L_two_R_virs = L_two['R_virs']
        L_two_impacts = L_two['impacts']
        L_two_MajDiams = L_two['MajDiams']
        L_two_ls = L_two['ls']
        L_two_bs = L_two['bs']
        
        L_two_Lya_Ws_cut = []
        L_two_bs_cut = []
        for w, b in zip(L_two_Lya_Ws, L_two_bs):
            if w <= max_EW and w >= min_EW:
                L_two_Lya_Ws_cut.append(w)
                L_two_bs_cut.append(b)

        
        # grab the three_plus data 
        L_three_Lya_Ws = L_three_plus['Lya_Ws']
        L_three_R_virs = L_three_plus['R_virs']
        L_three_impacts = L_three_plus['impacts']
        L_three_MajDiams = L_three_plus['MajDiams']
        L_three_ls = L_three_plus['ls']
        L_three_bs = L_three_plus['bs']
        
        L_three_Lya_Ws_cut = []
        L_three_bs_cut = []
        for w, b in zip(L_three_Lya_Ws, L_three_bs):
            if w <= max_EW and w >= min_EW:
                L_three_Lya_Ws_cut.append(w)
                L_three_bs_cut.append(b)


        # grab the group data and define the x and y data
        L_group_Lya_Ws = L_group['Lya_Ws']
        L_group_R_virs = L_group['R_virs']
        L_group_impacts = L_group['impacts']
        L_group_MajDiams = L_group['MajDiams']
        L_group_mems = L_group['group_mems']
        L_group_ls = L_group['ls']
        L_group_bs = L_group['bs']
        
        L_group_Lya_Ws_cut = []
        L_group_bs_cut = []
        for w, m, b in zip(L_group_Lya_Ws, L_group_mems, L_group_bs):
            if w <= max_EW and w >= min_EW and int(m) >= 3:
                L_group_Lya_Ws_cut.append(w)
                L_group_bs_cut.append(b)
                
                
        # all EWs
#         hist(all_Lya_Ws_cut,
#         bins=bins,
#         histtype='step',
#         cumulative=True,
#         normed=1,
#         color='black',
#         lw=1.5,
#         alpha=alpha_black,
#         label=r'$\rm All$')
        
        
        # isolated 
#         hist(isolated_Lya_Ws_cut,
#         bins=bins,
#         histtype='bar',
#         color='black',
#         lw=1.5,
#         alpha=alpha_black,
#         label=r'$\rm Isolated$')
        
        
        # L_isolated 
#         hist(L_isolated_Lya_Ws_cut,
#         bins=bins,
#         histtype='bar',
#         color=color_brown,
#         lw=1.5,
#         alpha=alpha_brown,
#         label=r'$\rm L\_isolated$')
        

        # isolated associated
        hist(L_isolated_associated_bs_cut,
        bins=bins,
        histtype='bar',
        color=color_green,
        lw=1.5,
        alpha=alpha_green,
        label=r'$\rm L\_isolated\_assoc.$')
        
        
        # associated
        hist(L_associated_bs_cut,
        bins=bins,
        histtype='step',
        color=color_orange,
        lw=2.5,
        alpha=alpha_orange+0.1,
        label=r'$\rm L\_assoc.$')


        # L_two
        L_two_plus = np.array(list(L_two_bs_cut) + list(L_three_bs_cut))
        hist(L_two_bs_cut,
        bins=bins,
        histtype='bar',
        color=color_purple3,
        lw=1.5,
        alpha=alpha_purple,
        label=r'$\rm L\_two+$')

        # three
#         hist(L_three_Lya_Ws_cut,
#         bins=bins,
#         histtype='bar',
#         color=color_red,
#         lw=1.5,
#         alpha=alpha_red,
#         label=r'$\rm L\_three$')
        
        
        # L_nonassociated
#         hist(L_nonassociated_Lya_Ws_cut,
#         bins=bins,
#         histtype='bar',
#         color=color_pink,
#         lw=1.5,
#         alpha=alpha_pink,
#         label=r'$\rm L\_nonassoc.$')
        
        
        # L_group
#         hist(L_group_bs_cut,
#         bins=bins,
#         histtype='bar',
#         color=color_blue,
#         lw=1.5,
#         alpha=alpha_yellow,
#         label=r'$\rm L\_group$')
        
        
#         color_green = '#1b9e77'
#         color_orange = '#d95f02'
#         color_purple3 = '#7570b3'
#         color_pink = '#e7298a'
#         color_lime = '#66a61e'
#         color_yellow = '#e6ab02'
#         color_brown = '#a6761d'
#         color_coal = '#666666'

        
        # x-axis
        majorLocator   = MultipleLocator(20)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(10)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)


        legend(scatterpoints=1, prop={'size':14}, loc='upper right', fancybox=True)
        xlabel(r'$\rm b ~[km s^{{-1}}]$')
        ylabel(r'$\rm Number$')

        ax.grid(b=None,which='major',axis='both')
#         ylim(0,1300)
        xlim(0, 200)

        if plot_b_hist_save:
            savefig('{0}/hist(b)_bins{1}_EWcut{2}-{3}_bothassoc_twoplus_dataset{4}.pdf'.format(saveDirectory, binSize, min_EW, max_EW, data_set),format='pdf',bbox_inches='tight')
        else:
            show()


##########################################################################################
##########################################################################################



##########################################################################################
##########################################################################################
    
    if plot_b_hist_az:
        fig = figure(figsize=(7.7,6.7))
        ax1 = fig.add_subplot(311)
        
        
        countb = 0
        countr = 0
        count = -1
        
        # azimuth bin right edges
        az1 = 30
        az2 = 60
        az3 = 90
        
        second = 'Associated'
        
        color_green = '#1b9e77'
        color_purple = '#7570b3'
        color_orange = '#d95f02'
        color_purple2 = '#984ea3'
        
        
        color_green = '#1b9e77'
        color_orange = '#d95f02'
        color_purple3 = '#7570b3'
        color_pink = '#e7298a'
        color_lime = '#66a61e'
        color_yellow = '#e6ab02'
        color_brown = '#a6761d'
        color_coal = '#666666'


        alpha_green = 0.7
        alpha_orange = 0.7
        alpha_purple = 0.7
        alpha_pink = 0.7
        alpha_lime = 0.7
        alpha_yellow = 0.7
        alpha_brown = 0.7
        alpha_coal = 0.7
        alpha_blue = 0.7
        alpha_red = 0.7
        alpha_black = 0.7
        
        binSize = 5
        bins = arange(0, 130+binSize, binSize)
        
        label_isolated = r'$\rm Isolated$'
        label_associated = r'$\rm Associated$'
        label_two = r'$\rm Two$'
        label_three = r'$\rm Three+$'
        label_group = r'$\rm Group$'
        label_second = r'$\rm {0}$'.format(second)


        color_isolated = 'black'
        color_associated = color_red
        color_two = color_red
        color_three = color_purple
        color_group = color_orange
        color_second = color_orange
        
        # grab the full set
        all_Lya_Ws = all['Lya_Ws']
        all_bs = all['bs']

        all_Lya_Ws_cut = []
        all_bs_cut = []
        for w, b in zip(all_Lya_Ws, all_bs):
            if w <= max_EW and w >= min_EW:
                all_Lya_Ws_cut.append(w)
                all_bs_cut.append(b)
                
        # isolated
        isolated_Lya_Ws = isolated['Lya_Ws']
        isolated_bs = isolated['bs']

        isolated_Lya_Ws_cut = []
        isolated_bs_cut = []
        for w, b in zip(isolated_Lya_Ws, isolated_bs):
            if w <= max_EW and w >= min_EW:
                isolated_Lya_Ws_cut.append(w)
                isolated_bs_cut.append(b)

                
        # L_isolated
        L_isolated_Lya_Ws = L_isolated['Lya_Ws']
        L_isolated_bs = L_isolated['bs']

        L_isolated_Lya_Ws_cut = []
        L_isolated_bs_cut = []
        for w, b in zip(L_isolated_Lya_Ws, L_isolated_bs):
            if w <= max_EW and w >= min_EW:
                L_isolated_Lya_Ws_cut.append(w)
                L_isolated_bs_cut.append(b)


        # L_nonassociated
        L_nonassociated_Lya_Ws = L_nonassociated['Lya_Ws']
        L_nonassociated_bs = L_nonassociated['bs']

        L_nonassociated_Lya_Ws_cut = []
        L_nonassociated_bs_cut = []
        for w, b in zip(L_nonassociated_Lya_Ws, L_nonassociated_bs):
            if w <= max_EW and w >= min_EW:
                L_nonassociated_Lya_Ws_cut.append(w)
                L_nonassociated_bs_cut.append(b)


        # L_isolated_associated
        L_isolated_associated_Lya_Ws_cut = []
        L_isolated_associated_bs_cut_az1 = []
        L_isolated_associated_bs_cut_az2 = []
        L_isolated_associated_bs_cut_az3 = []
        
        for w, b, az in zip(Lya_Ws, bs, azimuths):
            if w <= max_EW and w >= min_EW:
                L_isolated_associated_Lya_Ws_cut.append(w)
                if az <= az1:
                    L_isolated_associated_bs_cut_az1.append(b)
                elif az > az1 and az<= az2:
                    L_isolated_associated_bs_cut_az2.append(b)
                else:
                    L_isolated_associated_bs_cut_az3.append(b)
                    

        # grab the associated data 
        L_associated_Lya_Ws = L_associated['Lya_Ws']
        L_associated_R_virs = L_associated['R_virs']
        L_associated_impacts = L_associated['impacts']
        L_associated_MajDiams = L_associated['MajDiams']
        L_associated_ls = L_associated['ls']
        L_associated_bs = L_associated['bs']
        L_associated_azimuths = L_associated['azimuths']

        
        L_associated_Lya_Ws_cut = []
        L_associated_bs_cut_az1 = []
        L_associated_bs_cut_az2 = []
        L_associated_bs_cut_az3 = []
        for w, b, az in zip(L_associated_Lya_Ws, L_associated_bs, L_associated_azimuths):
            if w <= max_EW and w >= min_EW:
                L_associated_Lya_Ws_cut.append(w)
                if az <= az1:
                    L_associated_bs_cut_az1.append(b)
                elif az > az1 and az<= az2:
                    L_associated_bs_cut_az2.append(b)
                else:
                    L_associated_bs_cut_az3.append(b)
        
        # grab the two data 
        L_two_Lya_Ws = L_two['Lya_Ws']
        L_two_R_virs = L_two['R_virs']
        L_two_impacts = L_two['impacts']
        L_two_MajDiams = L_two['MajDiams']
        L_two_ls = L_two['ls']
        L_two_bs = L_two['bs']
        L_two_azimuths = L_two['azimuths']

        L_two_Lya_Ws_cut = []
        L_two_bs_cut_az1 = []
        L_two_bs_cut_az2 = []
        L_two_bs_cut_az3 = []
        for w, b, az in zip(L_two_Lya_Ws, L_two_bs, L_two_azimuths):
            if w <= max_EW and w >= min_EW:
                L_two_Lya_Ws_cut.append(w)
                if az <= az1:
                    L_two_bs_cut_az1.append(b)
                elif az > az1 and az<= az2:
                    L_two_bs_cut_az2.append(b)
                else:
                    L_two_bs_cut_az3.append(b)
        
        # grab the three_plus data 
        L_three_Lya_Ws = L_three_plus['Lya_Ws']
        L_three_R_virs = L_three_plus['R_virs']
        L_three_impacts = L_three_plus['impacts']
        L_three_MajDiams = L_three_plus['MajDiams']
        L_three_ls = L_three_plus['ls']
        L_three_bs = L_three_plus['bs']
        L_three_azimuths = L_three_plus['azimuths']

        L_three_Lya_Ws_cut = []
        L_three_bs_cut_az1 = []
        L_three_bs_cut_az2 = []
        L_three_bs_cut_az3 = []
        for w, b, az in zip(L_three_Lya_Ws, L_three_bs, L_three_azimuths):
            if w <= max_EW and w >= min_EW:
                L_three_Lya_Ws_cut.append(w)
                if az <= az1:
                    L_three_bs_cut_az1.append(b)
                elif az > az1 and az<= az2:
                    L_three_bs_cut_az2.append(b)
                else:
                    L_three_bs_cut_az3.append(b)

        # grab the group data and define the x and y data
        L_group_Lya_Ws = L_group['Lya_Ws']
        L_group_R_virs = L_group['R_virs']
        L_group_impacts = L_group['impacts']
        L_group_MajDiams = L_group['MajDiams']
        L_group_mems = L_group['group_mems']
        L_group_ls = L_group['ls']
        L_group_bs = L_group['bs']
        
        L_group_Lya_Ws_cut = []
        L_group_bs_cut = []
        for w, m, b in zip(L_group_Lya_Ws, L_group_mems, L_group_bs):
            if w <= max_EW and w >= min_EW and int(m) >= 3:
                L_group_Lya_Ws_cut.append(w)
                L_group_bs_cut.append(b)
                
                
        # all EWs
#         hist(all_Lya_Ws_cut,
#         bins=bins,
#         histtype='step',
#         cumulative=True,
#         normed=1,
#         color='black',
#         lw=1.5,
#         alpha=alpha_black,
#         label=r'$\rm All$')
        
        
        # isolated 
#         hist(isolated_Lya_Ws_cut,
#         bins=bins,
#         histtype='bar',
#         color='black',
#         lw=1.5,
#         alpha=alpha_black,
#         label=r'$\rm Isolated$')
        
        
        # L_isolated 
#         hist(L_isolated_Lya_Ws_cut,
#         bins=bins,
#         histtype='bar',
#         color=color_brown,
#         lw=1.5,
#         alpha=alpha_brown,
#         label=r'$\rm L\_isolated$')
        

        all_assoc_az1 = np.array(list(L_isolated_associated_bs_cut_az1) + list(L_associated_bs_cut_az1))

        # isolated associated
#         hist(L_isolated_associated_bs_cut_az1,
#         bins=bins,
#         histtype='bar',
#         color=color_green,
#         lw=1.5,
#         alpha=alpha_green,
#         label=r'$\rm L\_isolated\_assoc.$')
#         
#         # associated
#         hist(L_associated_bs_cut_az1,
#         bins=bins,
#         histtype='step',
#         color=color_orange,
#         lw=2.5,
#         alpha=alpha_orange+0.1,
#         label=r'$\rm L\_assoc.$')
        
        hist(all_assoc_az1,
        bins=bins,
        histtype='step',
        color=color_orange,
        lw=2.5,
        alpha=alpha_orange+0.1,
        label=r'$\rm All~Assoc.$')


        # L_two
        L_two_plus_az1 = np.array(list(L_two_bs_cut_az1) + list(L_three_bs_cut_az1))

        hist(L_two_plus_az1,
        bins=bins,
        histtype='bar',
        color=color_purple3,
        lw=1.5,
        alpha=alpha_purple,
        label=r'$\rm L\_two+$')

        # x-axis
        majorLocator   = MultipleLocator(20)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(10)
        ax1.xaxis.set_major_locator(majorLocator)
        ax1.xaxis.set_major_formatter(majorFormatter)
        ax1.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax1.yaxis.set_major_locator(majorLocator)
        ax1.yaxis.set_major_formatter(majorFormatter)
        ax1.yaxis.set_minor_locator(minorLocator)

        ax1.set_title(r'$\rm Az \leq 30^{{\circ}}$')

        ############################################################################
        ax2 = fig.add_subplot(312)
        
        # isolated associated
        
        all_assoc_az2 = np.array(list(L_isolated_associated_bs_cut_az2) + list(L_associated_bs_cut_az2))
        
#         hist(L_isolated_associated_bs_cut_az2,
#         bins=bins,
#         histtype='bar',
#         color=color_green,
#         lw=1.5,
#         alpha=alpha_green,
#         label=r'$\rm L\_isolated\_assoc.$')
#         
#         # associated
#         hist(L_associated_bs_cut_az2,
#         bins=bins,
#         histtype='step',
#         color=color_orange,
#         lw=2.5,
#         alpha=alpha_orange+0.1,
#         label=r'$\rm L\_assoc.$')

        hist(all_assoc_az2,
        bins=bins,
        histtype='step',
        color=color_orange,
        lw=2.5,
        alpha=alpha_orange+0.1,
        label=r'$\rm All~Assoc.$')


        # L_two
        L_two_plus_az2 = np.array(list(L_two_bs_cut_az2) + list(L_three_bs_cut_az2))

        hist(L_two_plus_az2,
        bins=bins,
        histtype='bar',
        color=color_purple3,
        lw=1.5,
        alpha=alpha_purple,
        label=r'$\rm L\_two+$')
        
        # x-axis
        majorLocator   = MultipleLocator(20)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(10)
        ax2.xaxis.set_major_locator(majorLocator)
        ax2.xaxis.set_major_formatter(majorFormatter)
        ax2.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax2.yaxis.set_major_locator(majorLocator)
        ax2.yaxis.set_major_formatter(majorFormatter)
        ax2.yaxis.set_minor_locator(minorLocator)
        
        ax2.set_title(r'$\rm 30 < Az \leq 60^{{\circ}}$')

        ############################################################################
        ax3 = fig.add_subplot(313)
        
        all_assoc_az3 = np.array(list(L_isolated_associated_bs_cut_az3) + list(L_associated_bs_cut_az3))

        
        # isolated associated
#         hist(L_isolated_associated_bs_cut_az3,
#         bins=bins,
#         histtype='bar',
#         color=color_green,
#         lw=1.5,
#         alpha=alpha_green,
#         label=r'$\rm L\_isolated\_assoc.$')
#         
#         # associated
#         hist(L_associated_bs_cut_az3,
#         bins=bins,
#         histtype='step',
#         color=color_orange,
#         lw=2.5,
#         alpha=alpha_orange+0.1,
#         label=r'$\rm L\_assoc.$')


        hist(all_assoc_az3,
        bins=bins,
        histtype='step',
        color=color_orange,
        lw=2.5,
        alpha=alpha_orange+0.1,
        label=r'$\rm All~Assoc.$')


        # L_two
        L_two_plus_az3 = np.array(list(L_two_bs_cut_az3) + list(L_three_bs_cut_az3))

        hist(L_two_plus_az3,
        bins=bins,
        histtype='bar',
        color=color_purple3,
        lw=1.5,
        alpha=alpha_purple,
        label=r'$\rm L\_two+$')

        ax3.set_title(r'$\rm 60 < Az \leq 90^{{\circ}}$')

        # x-axis
        majorLocator   = MultipleLocator(20)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(10)
        ax3.xaxis.set_major_locator(majorLocator)
        ax3.xaxis.set_major_formatter(majorFormatter)
        ax3.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(4)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(2)
        ax3.yaxis.set_major_locator(majorLocator)
        ax3.yaxis.set_major_formatter(majorFormatter)
        ax3.yaxis.set_minor_locator(minorLocator)


        legend(scatterpoints=1, prop={'size':12}, loc='upper right', fancybox=True)
        xlabel(r'$\rm b ~[km s^{{-1}}]$')
        ylabel(r'$\rm Number$')

        ax1.grid(b=None,which='major',axis='both')
        ax2.grid(b=None,which='major',axis='both')
        ax3.grid(b=None,which='major',axis='both')
        
        tight_layout()

        xlim(0, 140)

        if plot_b_hist_az_save:
            savefig('{0}/hist(b-az)_bins{1}_EWcut{2}-{3}_allAssoc_twoplus_dataset{4}.pdf'.format(saveDirectory, binSize, min_EW, max_EW, data_set),format='pdf',bbox_inches='tight')
        else:
            show()


##########################################################################################
##########################################################################################



##########################################################################################
##########################################################################################
    
    if plot_b_hist_group:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        ax.grid(b=None,which='major',axis='both')
        
        countb = 0
        countr = 0
        count = -1
        
        second = 'Associated'
        
        color_green = '#1b9e77'
        color_purple = '#7570b3'
        color_orange = '#d95f02'
        color_purple2 = '#984ea3'
        
        
        color_green = '#1b9e77'
        color_orange = '#d95f02'
        color_purple3 = '#7570b3'
        color_pink = '#e7298a'
        color_lime = '#66a61e'
        color_yellow = '#e6ab02'
        color_brown = '#a6761d'
        color_coal = '#666666'


        alpha_green = 0.8
        alpha_orange = 0.8
        alpha_purple = 0.8
        alpha_pink = 0.8
        alpha_lime = 0.8
        alpha_yellow = 0.99
        alpha_brown = 0.8
        alpha_coal = 0.8
        alpha_blue = 0.8
        alpha_red = 0.8
        alpha_black = 0.8
        
        binSize = 10
        bins = arange(0, 200, binSize)
        
        label_isolated = r'$\rm Isolated$'
        label_associated = r'$\rm Associated$'
        label_two = r'$\rm Two$'
        label_three = r'$\rm Three+$'
        label_group = r'$\rm Group$'
        label_second = r'$\rm {0}$'.format(second)


        color_isolated = 'black'
        color_associated = color_red
        color_two = color_red
        color_three = color_purple
        color_group = color_orange
        color_second = color_orange
        
        
        # grab the full set
        all_Lya_Ws = all['Lya_Ws']
        
        all_Lya_Ws_cut = []
        for w in all_Lya_Ws:
            if w <= max_EW:
                all_Lya_Ws_cut.append(w)
                
                
        # isolated
        isolated_Lya_Ws = isolated['Lya_Ws']
        isolated_Lya_Ws_cut = []
        for w in isolated_Lya_Ws:
            if w <= max_EW:
                isolated_Lya_Ws_cut.append(w)
                
        # L_isolated
        L_isolated_Lya_Ws = L_isolated['Lya_Ws']
        L_isolated_Lya_Ws_cut = []
        for w in L_isolated_Lya_Ws:
            if w <= max_EW:
                L_isolated_Lya_Ws_cut.append(w)
                

        # L_nonassociated
        L_nonassociated_Lya_Ws = L_nonassociated['Lya_Ws']
        L_nonassociated_Lya_Ws_cut = []
        for w in L_nonassociated_Lya_Ws:
            if w <= max_EW:
                L_nonassociated_Lya_Ws_cut.append(w)


        # L_isolate_associated
        L_isolated_associated_Lya_Ws_cut = []
        for w in Lya_Ws:
            if w <= max_EW:
                L_isolated_associated_Lya_Ws_cut.append(w)
        
        # grab the associated data 
        L_associated_Lya_Ws = L_associated['Lya_Ws']
        L_associated_R_virs = L_associated['R_virs']
        L_associated_impacts = L_associated['impacts']
        L_associated_MajDiams = L_associated['MajDiams']
        L_associated_ls = L_associated['ls']
        
        L_associated_Lya_Ws_cut = []
        for w in L_associated_Lya_Ws:
            if w <= max_EW:
                L_associated_Lya_Ws_cut.append(w)
        
        # grab the two data 
        L_two_Lya_Ws = L_two['Lya_Ws']
        L_two_R_virs = L_two['R_virs']
        L_two_impacts = L_two['impacts']
        L_two_MajDiams = L_two['MajDiams']
        L_two_ls = L_two['ls']
        
        L_two_Lya_Ws_cut = []
        for w in L_two_Lya_Ws:
            if w <= max_EW:
                L_two_Lya_Ws_cut.append(w)
        
        
        # grab the three_plus data 
        L_three_Lya_Ws = L_three_plus['Lya_Ws']
        L_three_R_virs = L_three_plus['R_virs']
        L_three_impacts = L_three_plus['impacts']
        L_three_MajDiams = L_three_plus['MajDiams']
        L_three_ls = L_three_plus['ls']
        
        L_three_Lya_Ws_cut = []
        for w in L_three_Lya_Ws:
            if w <= max_EW:
                L_three_Lya_Ws_cut.append(w)
        
        # grab the group data and define the x and y data
        L_group_Lya_Ws = L_group['Lya_Ws']
        L_group_R_virs = L_group['R_virs']
        L_group_impacts = L_group['impacts']
        L_group_MajDiams = L_group['MajDiams']
        L_group_mems = L_group['group_mems']
        L_group_ls = L_group['ls']
        L_group_bs = L_group['bs']

        L_group_Lya_Ws_cut24 = []
        L_group_Lya_Ws_cut57 = []
        L_group_Lya_Ws_cut810 = []
        L_group_Lya_Ws_cut11plus = []
        
        L_group_bs_cut24 = []
        L_group_bs_cut57 = []
        L_group_bs_cut810 = []
        L_group_bs_cut11plus = []
        for w, m, b in zip(L_group_Lya_Ws, L_group_mems, L_group_bs):
            if w <= max_EW:
                if int(m) >=2 and int(m) <=4:
                    L_group_Lya_Ws_cut24.append(w)
                    L_group_bs_cut24.append(b)
                    
                if int(m) >4 and int(m) <=7:
                    L_group_Lya_Ws_cut57.append(w)
                    L_group_bs_cut57.append(b)
                    
                if int(m) >7 and int(m) <=10:
                    L_group_Lya_Ws_cut810.append(w)
                    L_group_bs_cut810.append(b)
                    
                if int(m) >10:
                    L_group_Lya_Ws_cut11plus.append(w)
                    L_group_bs_cut11plus.append(b)

        # L_group
        hist(L_group_bs_cut24,
        bins=bins,
        histtype='step',
        color=color_blue,
        lw=2.5,
        alpha=alpha_yellow,
        label=r'$\rm L\_group~2-4$')
        
        # L_group
        hist(L_group_bs_cut57,
        bins=bins,
        histtype='step',
        color=color_green,
        lw=2.5,
        alpha=alpha_yellow,
        label=r'$\rm L\_group~5-7$')
        
        # L_group
        hist(L_group_bs_cut810,
        bins=bins,
        histtype='step',
        color=color_red,
        lw=2.5,
        alpha=alpha_yellow,
        label=r'$\rm L\_group~8-10$')
        
        # L_group
        hist(L_group_bs_cut11plus,
        bins=bins,
        histtype='step',
        color='black',
        lw=2.5,
        alpha=alpha_yellow,
        label=r'$\rm L\_group~11+$')

        
#         color_green = '#1b9e77'
#         color_orange = '#d95f02'
#         color_purple3 = '#7570b3'
#         color_pink = '#e7298a'
#         color_lime = '#66a61e'
#         color_yellow = '#e6ab02'
#         color_brown = '#a6761d'
#         color_coal = '#666666'

        
        # x-axis
        majorLocator   = MultipleLocator(20)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(10)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)


        legend(scatterpoints=1, prop={'size':14}, loc='upper right', fancybox=True)
        xlabel(r'$\rm b ~[km s^{{-1}}]$')
        ylabel(r'$\rm Number$')

#         ylim(0,1300)
        xlim(0, 200)

        if plot_b_hist_group_save:
            savefig('{0}/hist(b_group)_bins{1}.pdf'.format(saveDirectory, binSize),format='pdf',bbox_inches='tight')
        else:
            show()



##########################################################################################
##########################################################################################
    
    if plot_b_hist_MType:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(211)
        
        countb = 0
        countr = 0
        count = -1
        
        second = 'Associated'
        
        color_green = '#1b9e77'
        color_purple = '#7570b3'
        color_orange = '#d95f02'
        color_purple2 = '#984ea3'
        
        
        color_green = '#1b9e77'
        color_orange = '#d95f02'
        color_purple3 = '#7570b3'
        color_pink = '#e7298a'
        color_lime = '#66a61e'
        color_yellow = '#e6ab02'
        color_brown = '#a6761d'
        color_coal = '#666666'


        alpha_green = 0.8
        alpha_orange = 0.8
        alpha_purple = 0.8
        alpha_pink = 0.8
        alpha_lime = 0.8
        alpha_yellow = 0.99
        alpha_brown = 0.8
        alpha_coal = 0.8
        alpha_blue = 0.8
        alpha_red = 0.8
        alpha_black = 0.8
        
        binSize = 10
        bins = arange(0, 200, binSize)
        
        label_isolated = r'$\rm Isolated$'
        label_associated = r'$\rm Associated$'
        label_two = r'$\rm Two$'
        label_three = r'$\rm Three+$'
        label_group = r'$\rm Group$'
        label_second = r'$\rm {0}$'.format(second)


        color_isolated = 'black'
        color_associated = color_red
        color_two = color_red
        color_three = color_purple
        color_group = color_orange
        color_second = color_orange


        # L_nonassociated
        L_nonassociated_Lya_Ws = L_nonassociated['Lya_Ws']
        L_nonassociated_MTypes = L_nonassociated['MTypes']
        L_nonassociated_bs = L_nonassociated['bs']
        
        L_nonassociated_Lya_Ws_cut = []
        L_nonassociated_MTypes_cut = []
        L_nonassociated_bs_cut = []
        for w, type, b in zip(L_nonassociated_Lya_Ws, L_nonassociated_MTypes, L_nonassociated_bs):
            if w <= max_EW:
                L_nonassociated_Lya_Ws_cut.append(w)
                L_nonassociated_MTypes_cut.append(type)
                L_nonassociated_bs_cut.append(b)


        # L_isolate_associated
        L_isolated_associated_MTypes = MTypes
        L_isolated_associated_bs = bs
        
        L_isolated_associated_Lya_Ws_cut = []
        L_isolated_associated_MTypes_cut = []
        L_isolated_associated_bs_cut = []
        for w, type, b in zip(Lya_Ws, L_isolated_associated_MTypes, L_isolated_associated_bs):
            if w <= max_EW:
                L_isolated_associated_Lya_Ws_cut.append(w)
                L_isolated_associated_MTypes_cut.append(type)
                L_isolated_associated_bs_cut.append(b)


        # grab the associated data 
        L_associated_Lya_Ws = L_associated['Lya_Ws']
        L_associated_R_virs = L_associated['R_virs']
        L_associated_impacts = L_associated['impacts']
        L_associated_MajDiams = L_associated['MajDiams']
        L_associated_ls = L_associated['ls']
        L_associated_MTypes = L_associated['MTypes']
        L_associated_bs = L_associated['bs']

        L_associated_Lya_Ws_cut_e = []
        L_associated_Lya_Ws_cut_s = []
        
        L_associated_bs_cut_e = []
        L_associated_bs_cut_s = []


        L_associated_MTypes_cut_e = []
        L_associated_MTypes_cut_s = []

        for w, type, b in zip(L_associated_Lya_Ws, L_associated_MTypes, L_associated_bs):
            if w <= max_EW:                
                if type[0] == 'E':
                    L_associated_MTypes_cut_e.append(type)
                    L_associated_Lya_Ws_cut_e.append(w)
                    L_associated_bs_cut_e.append(b)

                else:
                    L_associated_MTypes_cut_s.append(type)
                    L_associated_Lya_Ws_cut_s.append(w)
                    L_associated_bs_cut_s.append(b)

        
        # grab the two data 
        L_two_Lya_Ws = L_two['Lya_Ws']
        L_two_R_virs = L_two['R_virs']
        L_two_impacts = L_two['impacts']
        L_two_MajDiams = L_two['MajDiams']
        L_two_ls = L_two['ls']
        L_two_MTypes = L_two['MTypes']
        
        L_two_Lya_Ws_cut = []
        L_two_MTypes_cut = []
        for w, type in zip(L_two_Lya_Ws, L_two_MTypes):
            if w <= max_EW:
                L_two_Lya_Ws_cut.append(w)
                L_two_MTypes_cut.append(type)

        
        # grab the three_plus data 
        L_three_Lya_Ws = L_three_plus['Lya_Ws']
        L_three_R_virs = L_three_plus['R_virs']
        L_three_impacts = L_three_plus['impacts']
        L_three_MajDiams = L_three_plus['MajDiams']
        L_three_ls = L_three_plus['ls']
        
        L_three_Lya_Ws_cut = []
        for w in L_three_Lya_Ws:
            if w <= max_EW:
                L_three_Lya_Ws_cut.append(w)
        
        # grab the group data and define the x and y data
        L_group_Lya_Ws = L_group['Lya_Ws']
        L_group_R_virs = L_group['R_virs']
        L_group_impacts = L_group['impacts']
        L_group_MajDiams = L_group['MajDiams']
        L_group_mems = L_group['group_mems']
        L_group_ls = L_group['ls']
        
        L_group_Lya_Ws_cut24 = []
        L_group_Lya_Ws_cut57 = []
        L_group_Lya_Ws_cut810 = []
        L_group_Lya_Ws_cut11plus = []
        for w, m in zip(L_group_Lya_Ws, L_group_mems):
            if w <= max_EW:
                if int(m) >=2 and int(m) <=4:
                    L_group_Lya_Ws_cut24.append(w)
                    
        # L_group
        hist(L_associated_bs_cut_e,
        bins=bins,
        histtype='bar',
        color=color_red,
        lw=2.5,
        alpha=alpha_red,
        edgecolor='black',
        label=r'$\rm E-Type$')
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        legend(scatterpoints=1, prop={'size':14}, loc='upper right', fancybox=True)
        xlabel(r'$\rm b ~[km s^{{-1}}]$')
        ylabel(r'$\rm Number$')
        
#         ylim(0,1300)
#         xlim(0, 1000)
#         ax.grid(b=None,which='major',axis='both')
#         ax.set_axisbelow(True)
#         ax.yaxis.grid(color='gray', linestyle='dashed')

##########################################################################################
        ax = fig.add_subplot(212)

        # L_group
        hist(L_associated_bs_cut_s,
        bins=bins,
        histtype='bar',
        color=color_blue,
        lw=2.5,
        alpha=alpha_blue,
        edgecolor='black',
        label=r'$\rm S-type$')
        
        
        # x-axis
        majorLocator   = MultipleLocator(20)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(10)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(4)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(2)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        legend(scatterpoints=1, prop={'size':14}, loc='upper right', fancybox=True)
        xlabel(r'$\rm b ~[km s^{{-1}}]$')
        ylabel(r'$\rm Number$')
        
#         ylim(0,1300)
#         xlim(0, 1000)
        
#         ax.grid(b=None,which='major',axis='both')
#         ax.set_axisbelow(True)
#         ax.yaxis.grid(color='gray', linestyle='dashed')

        if plot_b_hist_MType_save:
            savefig('{0}/hist(b)_MType_bins{1}.pdf'.format(saveDirectory, binSize),format='pdf',bbox_inches='tight')
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
    