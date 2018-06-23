#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_EW_group_mems.py, v 1.0 06/19/18

Plot EW as a function of group membership of nearby galaxy


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
    # plot EW histograms as a function of number of group members
    plot_EW_vs_group_mems = True
    plot_EW_vs_group_mems_save = True

    
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
        L_summed_filename = '/Users/frenchd/Research/inclination/git_inclination/L_summed6.p'
        all_filename = '/Users/frenchd/Research/inclination/git_inclination/all6.p'


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
    all_file = open(all_filename,'r')


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
    all = pickle.load(all_file)

    
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
    
    if plot_EW_vs_group_mems:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)

        EW_cut = 1000
        markerSize = 20
        
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
        
        binSize = 1
        bins = arange(0, 30, binSize)
        
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
        
        
        # grab the two_plus data 
        L_three_Lya_Ws = L_two_plus['Lya_Ws']
        L_three_R_virs = L_two_plus['R_virs']
        L_three_impacts = L_two_plus['impacts']
        L_three_MajDiams = L_two_plus['MajDiams']
        L_three_ls = L_two_plus['ls']
        
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
        
        L_group_Lya_Ws_cut = []
        L_group_group_mems_cut = []

        for w, m in zip(L_group_Lya_Ws, L_group_mems):
            if w <= EW_cut:
                L_group_Lya_Ws_cut.append(w)
                L_group_group_mems_cut.append(m)
                    
        # L_group
        # associated
        plot1 = scatter(L_group_group_mems_cut,
                        L_group_Lya_Ws_cut,
                        marker='D',
                        c=color_blue,
                        s=markerSize,
                        edgecolor='black',
                        alpha=0.8)
        
        # histogram associated
#         bin_means, edges, binNumber = stats.binned_statistic(L_group_group_mems_cut,
#                                                             L_group_Lya_Ws_cut,
#                                                             statistic='mean',
#                                                             bins=bins)
#         left,right = edges[:-1],edges[1:]
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plot(X,
#             Y,
#             ls='solid',
#             color='black',
#             lw=2.0,
#             alpha=0.99,
#             label=r'$\rm Mean ~EW$')
        
        
#         color_green = '#1b9e77'
#         color_orange = '#d95f02'
#         color_purple3 = '#7570b3'
#         color_pink = '#e7298a'
#         color_lime = '#66a61e'
#         color_yellow = '#e6ab02'
#         color_brown = '#a6761d'
#         color_coal = '#666666'

        
        # x-axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
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


        legend(scatterpoints=1, prop={'size':14}, loc='upper right', fancybox=True)
        ylabel(r'$\rm Equivalent~Width ~[m\AA]$')
        xlabel(r'$\rm Number of Group Members$')

#         ylim(0,1300)
#         xlim(0, 1000)
        ax.grid(b=None,which='major',axis='both')

        ax.set_axisbelow(True)
        ax.yaxis.grid(color='gray', linestyle='dashed')

        if plot_EW_vs_group_mems_save:
            savefig('{0}/W(group_mems)_bins{1}.pdf'.format(saveDirectory, binSize),format='pdf',bbox_inches='tight')
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
    