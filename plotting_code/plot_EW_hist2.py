#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_EW_hist2.py, v 2.0 06/20/18

Plot EW histograms and CDFs

v2: include errors
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
    # plot CDF of EWs
    plot_EW_cdf = True
    plot_EW_cdf_save = True

    # plot histogram of EWs
    plot_EW_hist = False
    plot_EW_hist_save = False
    
    # plot EW histograms as a function of number of group members
    plot_EW_hist_group = False
    plot_EW_hist_group_save = False
    
    # plot EW histograms as a function of number of group members
    plot_EW_hist_MType = False
    plot_EW_hist_MType_save = False



    # plot_number = 1 for just the isolated sample, =2 adds the associated, =3 adds two+
    # =4 adds groups with 2 or more members
    plot_number = 4
    
    # which data set to use? Options are '', '_min001', '_cus','_min001_cus',
    # '_min001_double', '_min005_v150', '_min005_v250'
    data_set = ''
    
    min_EW = 50
    max_EW = 15000

    
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
    
    if plot_EW_cdf:
        
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
        
        alpha_error = 0.5
        
        binSize = 10
        bins = arange(0, 1000+binSize, binSize)
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
        all_e_Lya_Ws = all['e_Lya_Ws']
        
        all_Lya_Ws_cut = []
        all_e_Lya_Ws_cut = []
        for w, e in zip(all_Lya_Ws, all_e_Lya_Ws):
            if w <= max_EW and w >= min_EW:
                all_Lya_Ws_cut.append(w)
                all_e_Lya_Ws_cut.append(e)

                
        # isolated
        isolated_Lya_Ws = isolated['Lya_Ws']
        isolated_e_Lya_Ws = isolated['e_Lya_Ws']
        
        isolated_Lya_Ws_cut = []
        isolated_e_Lya_Ws_cut = []
        for w, e in zip(isolated_Lya_Ws, isolated_e_Lya_Ws):
            if w <= max_EW and w >= min_EW:
                isolated_Lya_Ws_cut.append(w)
                isolated_e_Lya_Ws_cut.append(e)


        # L_isolated
        L_isolated_Lya_Ws = L_isolated['Lya_Ws']
        L_isolated_e_Lya_Ws = L_isolated['e_Lya_Ws']
        
        L_isolated_Lya_Ws_cut = []
        L_isolated_e_Lya_Ws_cut = []
        for w, e in zip(L_isolated_Lya_Ws, L_isolated_e_Lya_Ws):
            if w <= max_EW and w >= min_EW:
                L_isolated_Lya_Ws_cut.append(w)
                L_isolated_e_Lya_Ws_cut.append(e)


        # L_nonassociated
        L_nonassociated_Lya_Ws = L_nonassociated['Lya_Ws']
        L_nonassociated_e_Lya_Ws = L_nonassociated['e_Lya_Ws']
        
        L_nonassociated_Lya_Ws_cut = []
        L_nonassociated_e_Lya_Ws_cut = []
        for w, e in zip(L_nonassociated_Lya_Ws, L_nonassociated_e_Lya_Ws):
            if w <= max_EW and w >= min_EW:
                L_nonassociated_Lya_Ws_cut.append(w)
                L_nonassociated_e_Lya_Ws_cut.append(e)


        # L_isolate_associated
        L_isolated_associated_Lya_Ws_cut = []
        L_isolated_associated_e_Lya_Ws_cut = []
        for w, e in zip(Lya_Ws, e_Lya_Ws):
            if w <= max_EW and w >= min_EW:
                L_isolated_associated_Lya_Ws_cut.append(w)
                L_isolated_associated_e_Lya_Ws_cut.append(e)


        # grab the associated data 
        L_associated_Lya_Ws = L_associated['Lya_Ws']
        L_associated_e_Lya_Ws = L_associated['e_Lya_Ws']
        L_associated_R_virs = L_associated['R_virs']
        L_associated_impacts = L_associated['impacts']
        L_associated_MajDiams = L_associated['MajDiams']
        L_associated_ls = L_associated['ls']
        
        L_associated_Lya_Ws_cut = []
        L_associated_e_Lya_Ws_cut = []
        for w, e in zip(L_associated_Lya_Ws, L_associated_e_Lya_Ws):
            if w <= max_EW and w >= min_EW:
                L_associated_Lya_Ws_cut.append(w)
                L_associated_e_Lya_Ws_cut.append(e)


        # grab the two data 
        L_two_Lya_Ws = L_two['Lya_Ws']
        L_two_e_Lya_Ws = L_two['e_Lya_Ws']
        L_two_R_virs = L_two['R_virs']
        L_two_impacts = L_two['impacts']
        L_two_MajDiams = L_two['MajDiams']
        L_two_ls = L_two['ls']
        
        L_two_Lya_Ws_cut = []
        L_two_e_Lya_Ws_cut = []
        for w, e in zip(L_two_Lya_Ws, L_two_e_Lya_Ws):
            if w <= max_EW and w >= min_EW:
                L_two_Lya_Ws_cut.append(w)
                L_two_e_Lya_Ws_cut.append(e)

        
        # grab the three_plus data 
        L_three_Lya_Ws = L_three_plus['Lya_Ws']
        L_three_e_Lya_Ws = L_three_plus['e_Lya_Ws']
        L_three_R_virs = L_three_plus['R_virs']
        L_three_impacts = L_three_plus['impacts']
        L_three_MajDiams = L_three_plus['MajDiams']
        L_three_ls = L_three_plus['ls']
        
        L_three_Lya_Ws_cut = []
        L_three_e_Lya_Ws_cut = []
        for w, e in zip(L_three_Lya_Ws, L_three_e_Lya_Ws):
            if w <= max_EW and w >= min_EW:
                L_three_Lya_Ws_cut.append(w)
                L_three_e_Lya_Ws_cut.append(e)

        
        # grab the group data and define the x and y data
        L_group_Lya_Ws = L_group['Lya_Ws']
        L_group_e_Lya_Ws = L_group['e_Lya_Ws']
        L_group_R_virs = L_group['R_virs']
        L_group_impacts = L_group['impacts']
        L_group_MajDiams = L_group['MajDiams']
        L_group_mems = L_group['group_mems']
        L_group_ls = L_group['ls']
        
        L_group_Lya_Ws_cut = []
        L_group_e_Lya_Ws_cut = []
        for w, m, e in zip(L_group_Lya_Ws, L_group_mems, L_group_e_Lya_Ws):
            if w <= max_EW and w >= min_EW and int(m) >= 2:
                L_group_Lya_Ws_cut.append(w)
                L_group_e_Lya_Ws_cut.append(e)

##########################################################################################
        bins_right = bins[1:]                           

        count_list = [1,2,3,4,5,6]
        
        for count in count_list:
            print 'Plotting number {0}'.format(count)
            
            fig = figure(figsize=(7.7,5.7))
            ax = fig.add_subplot(111)
        
        
            # all EWs
            hist(all_Lya_Ws_cut,
            bins=bins,
            histtype='step',
            cumulative=True,
            normed=1,
            color='black',
            lw=1.5,
            alpha=alpha_black,
            label=r'$\rm All$')
         
            min_n, min_bins, min_p = hist(np.array(all_Lya_Ws_cut) - np.array(all_e_Lya_Ws_cut),
                                            bins=bins,
                                            histtype='step',
                                            cumulative=True,
                                            normed=1,
                                            color='black',
                                            lw=1.5,
                                            alpha=0)
                                        
            max_n, max_bins, max_p = hist(np.array(all_Lya_Ws_cut) + np.array(all_e_Lya_Ws_cut),
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
                hist(isolated_Lya_Ws_cut,
                bins=bins,
                histtype='step',
                cumulative=True,
                normed=1,
                color='grey',
                lw=1.5,
                alpha=alpha_black,
                label=r'$\rm Isolated$')

                min_n, min_bins, min_p = hist(np.array(isolated_Lya_Ws_cut) - np.array(isolated_e_Lya_Ws_cut),
                                                bins=bins,
                                                histtype='step',
                                                cumulative=True,
                                                normed=1,
                                                color='grey',
                                                lw=1.5,
                                                alpha=0)
                                        
                max_n, max_bins, max_p = hist(np.array(isolated_Lya_Ws_cut) + np.array(isolated_e_Lya_Ws_cut),
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
                hist(L_isolated_Lya_Ws_cut,
                bins=bins,
                histtype='step',
                cumulative=True,
                normed=1,
                color=color_brown,
                lw=1.5,
                alpha=alpha_brown,
                label=r'$\rm \mathcal{L}-Isolated$')
        
                min_n, min_bins, min_p = hist(np.array(L_isolated_Lya_Ws_cut) - np.array(L_isolated_e_Lya_Ws_cut),
                                                bins=bins,
                                                histtype='step',
                                                cumulative=True,
                                                normed=1,
                                                color=color_brown,
                                                lw=1.5,
                                                alpha=0)
                                        
                max_n, max_bins, max_p = hist(np.array(L_isolated_Lya_Ws_cut) + np.array(L_isolated_e_Lya_Ws_cut),
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
                hist(L_isolated_associated_Lya_Ws_cut,
                bins=bins,
                histtype='step',
                cumulative=True,
                normed=1,
                color=color_green,
                lw=1.5,
                alpha=alpha_green,
                label=r'$\rm \mathcal{L}-Isolated-Assoc.$')
        
                min_n, min_bins, min_p = hist(np.array(L_isolated_associated_Lya_Ws_cut) - np.array(L_isolated_associated_e_Lya_Ws_cut),
                                                bins=bins,
                                                histtype='step',
                                                cumulative=True,
                                                normed=1,
                                                color=color_green,
                                                lw=1.5,
                                                alpha=0)
                                        
                max_n, max_bins, max_p = hist(np.array(L_isolated_associated_Lya_Ws_cut) + np.array(L_isolated_associated_e_Lya_Ws_cut),
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
                hist(L_associated_Lya_Ws_cut,
                bins=bins,
                histtype='step',
                cumulative=True,
                normed=1,
                color=color_orange,
                lw=1.5,
                alpha=alpha_orange,
                label=r'$\rm \mathcal{L}-Assoc.$')
        
                min_n, min_bins, min_p = hist(np.array(L_associated_Lya_Ws_cut) - np.array(L_associated_e_Lya_Ws_cut),
                                                bins=bins,
                                                histtype='step',
                                                cumulative=True,
                                                normed=1,
                                                color=color_orange,
                                                lw=1.5,
                                                alpha=0)
                                        
                max_n, max_bins, max_p = hist(np.array(L_associated_Lya_Ws_cut) + np.array(L_associated_e_Lya_Ws_cut),
                                                bins=bins,
                                                histtype='step',
                                                cumulative=True,
                                                normed=1,
                                                color=color_orange,
                                                lw=1.5,
                                                alpha=0)           

                fill_between(bins_right, min_n, max_n, facecolor=color_orange, alpha=alpha_error)
        


        # L_two+
#         L_twoplus_Lya_Ws_cut = L_two_Lya_Ws_cut + L_three_Lya_Ws
#         L_twoplus_e_Lya_Ws_cut = L_two_e_Lya_Ws_cut + L_three_e_Lya_Ws
# 
#         hist(L_twoplus_Lya_Ws_cut,
#         bins=bins,
#         histtype='step',
#         cumulative=True,
#         normed=1,
#         color=color_purple3,
#         lw=1.5,
#         alpha=alpha_purple,
#         label=r'$\rm L\_two+$')
#         
#         min_n, min_bins, min_p = hist(np.array(L_twoplus_Lya_Ws_cut) - np.array(L_twoplus_e_Lya_Ws_cut),
#                                         bins=bins,
#                                         histtype='step',
#                                         cumulative=True,
#                                         normed=1,
#                                         color=color_purple3,
#                                         lw=1.5,
#                                         alpha=0)
#                                         
#         max_n, max_bins, max_p = hist(np.array(L_twoplus_Lya_Ws_cut) + np.array(L_twoplus_e_Lya_Ws_cut),
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
                two_plus_Lya_Ws_cut = np.array(list(L_two_Lya_Ws_cut) + list(L_three_Lya_Ws_cut))
                two_plus_e_Lya_Ws_cut = np.array(list(L_two_e_Lya_Ws_cut) + list(L_three_e_Lya_Ws_cut))

                hist(two_plus_Lya_Ws_cut,
                bins=bins,
                histtype='step',
                cumulative=True,
                normed=1,
                color=color_purple3,
                lw=1.5,
                alpha=alpha_purple,
                label=r'$\rm \mathcal{L}-Two+$')
        
                min_n, min_bins, min_p = hist(np.array(two_plus_Lya_Ws_cut) - np.array(two_plus_e_Lya_Ws_cut),
                                                bins=bins,
                                                histtype='step',
                                                cumulative=True,
                                                normed=1,
                                                color=color_purple3,
                                                lw=1.5,
                                                alpha=0)
                                        
                max_n, max_bins, max_p = hist(np.array(two_plus_Lya_Ws_cut) + np.array(two_plus_e_Lya_Ws_cut),
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
                hist(L_group_Lya_Ws_cut,
                bins=bins,
                histtype='step',
                cumulative=True,
                normed=1,
                color=color_blue,
                lw=1.5,
                alpha=alpha_blue,
                label=r'$\rm L\_group$')
        
       
                min_n, min_bins, min_p = hist(np.array(L_group_Lya_Ws_cut) - np.array(L_group_e_Lya_Ws_cut),
                                                bins=bins,
                                                histtype='step',
                                                cumulative=True,
                                                normed=1,
                                                color=color_blue,
                                                lw=1.5,
                                                alpha=0)
                                        
                max_n, max_bins, max_p = hist(np.array(L_group_Lya_Ws_cut) + np.array(L_group_e_Lya_Ws_cut),
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
    #         hist(L_nonassociated_Lya_Ws_cut,
    #         bins=bins,
    #         histtype='step',
    #         cumulative=True,
    #         normed=1,
    #         color=color_pink,
    #         lw=1.5,
    #         alpha=alpha_pink,
    #         label=r'$\rm L\_nonassoc.$')
        
        
            # x-axis
            majorLocator   = MultipleLocator(200)
            majorFormatter = FormatStrFormatter(r'$\rm %d$')
            minorLocator   = MultipleLocator(100)
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
            xlabel(r'$\rm Equivalent~Width ~[m\AA]$')
            ylabel(r'$\rm CDF$')


            # grid
            ax.grid(b=None,which='major',axis='both')
            ax.set_axisbelow(True)
            ax.yaxis.grid(color='gray', linestyle='solid',alpha=0.5)
            ax.xaxis.grid(color='gray', linestyle='solid',alpha=0.5)


    #         ylim(0,1300)
            xlim(0, 1000)

            if plot_EW_cdf_save:
                savefig('{0}/hist(EW)_bins{1}_{2}_EWcut{3}-{4}_err_dataset{5}.pdf'.format(saveDirectory, binSize, count, min_EW, max_EW, data_set),format='pdf',bbox_inches='tight')
            else:
                show()


##########################################################################################
##########################################################################################

    
##########################################################################################
##########################################################################################
    
    if plot_EW_hist:
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


        alpha_green = 0.8
        alpha_orange = 0.8
        alpha_purple = 0.8
        alpha_pink = 0.8
        alpha_lime = 0.8
        alpha_yellow = 0.8
        alpha_brown = 0.8
        alpha_coal = 0.8
        alpha_blue = 0.8
        alpha_red = 0.8
        alpha_black = 0.8
        
        binSize = 50
        bins = arange(0, 1050, binSize)
        
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
            if w <= max_EW and w >= min_EW:
                all_Lya_Ws_cut.append(w)
                
                
        # isolated
        isolated_Lya_Ws = isolated['Lya_Ws']
        isolated_Lya_Ws_cut = []
        for w in isolated_Lya_Ws:
            if w <= max_EW and w >= min_EW:
                isolated_Lya_Ws_cut.append(w)
                
        # L_isolated
        L_isolated_Lya_Ws = L_isolated['Lya_Ws']
        L_isolated_Lya_Ws_cut = []
        for w in L_isolated_Lya_Ws:
            if w <= max_EW and w >= min_EW:
                L_isolated_Lya_Ws_cut.append(w)
                

        # L_nonassociated
        L_nonassociated_Lya_Ws = L_nonassociated['Lya_Ws']
        L_nonassociated_Lya_Ws_cut = []
        for w in L_nonassociated_Lya_Ws:
            if w <= max_EW and w >= min_EW:
                L_nonassociated_Lya_Ws_cut.append(w)


        # L_isolate_associated
        L_isolated_associated_Lya_Ws_cut = []
        for w in Lya_Ws:
            if w <= max_EW and w >= min_EW:
                L_isolated_associated_Lya_Ws_cut.append(w)
        
        # grab the associated data 
        L_associated_Lya_Ws = L_associated['Lya_Ws']
        L_associated_R_virs = L_associated['R_virs']
        L_associated_impacts = L_associated['impacts']
        L_associated_MajDiams = L_associated['MajDiams']
        L_associated_ls = L_associated['ls']
        
        L_associated_Lya_Ws_cut = []
        for w in L_associated_Lya_Ws:
            if w <= max_EW and w >= min_EW:
                L_associated_Lya_Ws_cut.append(w)
        
        # grab the two data 
        L_two_Lya_Ws = L_two['Lya_Ws']
        L_two_R_virs = L_two['R_virs']
        L_two_impacts = L_two['impacts']
        L_two_MajDiams = L_two['MajDiams']
        L_two_ls = L_two['ls']
        
        L_two_Lya_Ws_cut = []
        for w in L_two_Lya_Ws:
            if w <= max_EW and w >= min_EW:
                L_two_Lya_Ws_cut.append(w)
        
        
        # grab the three_plus data 
        L_three_Lya_Ws = L_three_plus['Lya_Ws']
        L_three_R_virs = L_three_plus['R_virs']
        L_three_impacts = L_three_plus['impacts']
        L_three_MajDiams = L_three_plus['MajDiams']
        L_three_ls = L_three_plus['ls']
        
        L_three_Lya_Ws_cut = []
        for w in L_three_Lya_Ws:
            if w <= max_EW and w >= min_EW:
                L_three_Lya_Ws_cut.append(w)
        
        # grab the group data and define the x and y data
        L_group_Lya_Ws = L_group['Lya_Ws']
        L_group_R_virs = L_group['R_virs']
        L_group_impacts = L_group['impacts']
        L_group_MajDiams = L_group['MajDiams']
        L_group_mems = L_group['group_mems']
        L_group_ls = L_group['ls']
        
        L_group_Lya_Ws_cut = []
        for w, m in zip(L_group_Lya_Ws, L_group_mems):
            if w <= max_EW and w >= min_EW and int(m) >= 9:
                print 'm: ',m
                L_group_Lya_Ws_cut.append(w)
                
                
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
#         hist(L_isolated_associated_Lya_Ws_cut,
#         bins=bins,
#         histtype='bar',
#         color=color_green,
#         lw=1.5,
#         alpha=alpha_green,
#         label=r'$\rm L\_isolated\_assoc.$')
        
        
        # associated
#         hist(L_associated_Lya_Ws_cut,
#         bins=bins,
#         histtype='bar',
#         color=color_orange,
#         lw=1.5,
#         alpha=alpha_orange,
#         label=r'$\rm L\_assoc.$')


        # L_two
#         hist(L_two_Lya_Ws_cut,
#         bins=bins,
#         histtype='bar',
#         color=color_purple3,
#         lw=1.5,
#         alpha=alpha_purple,
#         label=r'$\rm L\_two$')

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
        hist(L_group_Lya_Ws_cut,
        bins=bins,
        histtype='bar',
        color=color_blue,
        lw=1.5,
        alpha=alpha_yellow,
        label=r'$\rm L\_group$')
        
        
#         color_green = '#1b9e77'
#         color_orange = '#d95f02'
#         color_purple3 = '#7570b3'
#         color_pink = '#e7298a'
#         color_lime = '#66a61e'
#         color_yellow = '#e6ab02'
#         color_brown = '#a6761d'
#         color_coal = '#666666'

        
        # x-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
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
        xlabel(r'$\rm Equivalent~Width ~[m\AA]$')
        ylabel(r'$\rm Number$')

        ax.grid(b=None,which='major',axis='both')
#         ylim(0,1300)
        xlim(0, 1000)

        if plot_EW_hist_save:
            savefig('{0}/hist(EW)_bins{1}_L_group9.pdf'.format(saveDirectory, binSize),format='pdf',bbox_inches='tight')
        else:
            show()


##########################################################################################
##########################################################################################


##########################################################################################
##########################################################################################
    
    if plot_EW_hist_group:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        ax.grid(b=None,which='major',axis='both')

        EW_cut = 1000
        
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
        
        binSize = 100
        bins = arange(0, 1100, binSize)
        
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
            if w <= EW_cut:
                all_Lya_Ws_cut.append(w)
                
                
        # isolated
        isolated_Lya_Ws = isolated['Lya_Ws']
        isolated_Lya_Ws_cut = []
        for w in isolated_Lya_Ws:
            if w <= EW_cut:
                isolated_Lya_Ws_cut.append(w)
                
        # L_isolated
        L_isolated_Lya_Ws = L_isolated['Lya_Ws']
        L_isolated_Lya_Ws_cut = []
        for w in L_isolated_Lya_Ws:
            if w <= EW_cut:
                L_isolated_Lya_Ws_cut.append(w)
                

        # L_nonassociated
        L_nonassociated_Lya_Ws = L_nonassociated['Lya_Ws']
        L_nonassociated_Lya_Ws_cut = []
        for w in L_nonassociated_Lya_Ws:
            if w <= EW_cut:
                L_nonassociated_Lya_Ws_cut.append(w)


        # L_isolate_associated
        L_isolated_associated_Lya_Ws_cut = []
        for w in Lya_Ws:
            if w <= EW_cut:
                L_isolated_associated_Lya_Ws_cut.append(w)
        
        # grab the associated data 
        L_associated_Lya_Ws = L_associated['Lya_Ws']
        L_associated_R_virs = L_associated['R_virs']
        L_associated_impacts = L_associated['impacts']
        L_associated_MajDiams = L_associated['MajDiams']
        L_associated_ls = L_associated['ls']
        
        L_associated_Lya_Ws_cut = []
        for w in L_associated_Lya_Ws:
            if w <= EW_cut:
                L_associated_Lya_Ws_cut.append(w)
        
        # grab the two data 
        L_two_Lya_Ws = L_two['Lya_Ws']
        L_two_R_virs = L_two['R_virs']
        L_two_impacts = L_two['impacts']
        L_two_MajDiams = L_two['MajDiams']
        L_two_ls = L_two['ls']
        
        L_two_Lya_Ws_cut = []
        for w in L_two_Lya_Ws:
            if w <= 1200:
                L_two_Lya_Ws_cut.append(w)
        
        
        # grab the three_plus data 
        L_three_Lya_Ws = L_three_plus['Lya_Ws']
        L_three_R_virs = L_three_plus['R_virs']
        L_three_impacts = L_three_plus['impacts']
        L_three_MajDiams = L_three_plus['MajDiams']
        L_three_ls = L_three_plus['ls']
        
        L_three_Lya_Ws_cut = []
        for w in L_three_Lya_Ws:
            if w <= EW_cut:
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
            if w <= EW_cut:
                if int(m) >=2 and int(m) <=4:
                    L_group_Lya_Ws_cut24.append(w)
                    
                if int(m) >4 and int(m) <=7:
                    L_group_Lya_Ws_cut57.append(w)
                    
                if int(m) >7 and int(m) <=10:
                    L_group_Lya_Ws_cut810.append(w)
                    
                if int(m) >10:
                    L_group_Lya_Ws_cut11plus.append(w)
                    
        # L_group
        hist(L_group_Lya_Ws_cut24,
        bins=bins,
        histtype='step',
        color=color_blue,
        lw=2.5,
        alpha=alpha_yellow,
        label=r'$\rm L\_group~2-4$')
        
        # L_group
        hist(L_group_Lya_Ws_cut57,
        bins=bins,
        histtype='step',
        color=color_green,
        lw=2.5,
        alpha=alpha_yellow,
        label=r'$\rm L\_group~5-7$')
        
        # L_group
        hist(L_group_Lya_Ws_cut810,
        bins=bins,
        histtype='step',
        color=color_red,
        lw=2.5,
        alpha=alpha_yellow,
        label=r'$\rm L\_group~8-10$')
        
        # L_group
        hist(L_group_Lya_Ws_cut11plus,
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
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
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
        xlabel(r'$\rm Equivalent~Width ~[m\AA]$')
        ylabel(r'$\rm Number$')

#         ylim(0,1300)
        xlim(0, 1000)

        if plot_EW_hist_group_save:
            savefig('{0}/hist(EW_group)_bins{1}.pdf'.format(saveDirectory, binSize),format='pdf',bbox_inches='tight')
        else:
            show()



##########################################################################################
##########################################################################################
    
    if plot_EW_hist_MType:
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
        
        binSize = 100
        bins = arange(0, 1100, binSize)
        
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
        
        L_nonassociated_Lya_Ws_cut = []
        L_nonassociated_MTypes_cut = []
        for w, type in zip(L_nonassociated_Lya_Ws, L_nonassociated_MTypes):
            if w <= max_EW:
                L_nonassociated_Lya_Ws_cut.append(w)
                L_nonassociated_MTypes_cut.append(type)


        # L_isolate_associated
        L_isolated_associated_MTypes = MTypes
        
        L_isolated_associated_Lya_Ws_cut = []
        L_isolated_associated_MTypes_cut = []
        for w, type in zip(Lya_Ws, L_isolated_associated_MTypes):
            if w <= max_EW:
                L_isolated_associated_Lya_Ws_cut.append(w)
                L_isolated_associated_MTypes_cut.append(type)


        # grab the associated data 
        L_associated_Lya_Ws = L_associated['Lya_Ws']
        L_associated_R_virs = L_associated['R_virs']
        L_associated_impacts = L_associated['impacts']
        L_associated_MajDiams = L_associated['MajDiams']
        L_associated_ls = L_associated['ls']
        L_associated_MTypes = L_associated['MTypes']
        
        L_associated_Lya_Ws_cut_e = []
        L_associated_Lya_Ws_cut_s = []

        L_associated_MTypes_cut_e = []
        L_associated_MTypes_cut_s = []

        for w, type in zip(L_associated_Lya_Ws, L_associated_MTypes):
            if w <= max_EW:                
                if type[0] == 'E':
                    L_associated_MTypes_cut_e.append(type)
                    L_associated_Lya_Ws_cut_e.append(w)

                else:
                    L_associated_MTypes_cut_s.append(type)
                    L_associated_Lya_Ws_cut_s.append(w)


        
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
        hist(L_associated_Lya_Ws_cut_e,
        bins=bins,
        histtype='bar',
        color=color_red,
        lw=2.5,
        alpha=alpha_red,
        edgecolor='black',
        label=r'$\rm E-Type$')
        
        # x-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
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
        xlabel(r'$\rm Equivalent~Width ~[m\AA]$')
        ylabel(r'$\rm Number$')
        
#         ylim(0,1300)
#         xlim(0, 1000)
#         ax.grid(b=None,which='major',axis='both')
#         ax.set_axisbelow(True)
#         ax.yaxis.grid(color='gray', linestyle='dashed')

##########################################################################################
        ax = fig.add_subplot(212)

        # L_group
        hist(L_associated_Lya_Ws_cut_s,
        bins=bins,
        histtype='bar',
        color=color_blue,
        lw=2.5,
        alpha=alpha_blue,
        edgecolor='black',
        label=r'$\rm S-type$')
        
        
        # x-axis
        majorLocator   = MultipleLocator(200)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(100)
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
        xlabel(r'$\rm Equivalent~Width ~[m\AA]$')
        ylabel(r'$\rm Number$')
        
#         ylim(0,1300)
#         xlim(0, 1000)
        
#         ax.grid(b=None,which='major',axis='both')
#         ax.set_axisbelow(True)
#         ax.yaxis.grid(color='gray', linestyle='dashed')

        if plot_EW_hist_MType_save:
            savefig('{0}/hist(EW)_MType_bins{1}.pdf'.format(saveDirectory, binSize),format='pdf',bbox_inches='tight')
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
    