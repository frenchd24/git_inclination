#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_b_inc.py, v 1.0 06/18/18

Plot dopplar b-parameter as a function of inclination

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
    # plot dopplar b parameter as a function of diameter w/o splitting, add other sets if you want
    plot_b_inc_mean_plus = True
    plot_b_inc_mean_plus_save = True

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



#########################################################################################
#########################################################################################


##########################################################################################
##########################################################################################
    
    if plot_b_inc_mean_plus:
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
        binSize = 15
        bins = arange(0, 100, binSize)

        
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
        
        maxEW = 15000.

        # define the x and y data for the isolated set
        
        Lya_Ws2 = []
        R_virs2 = []
        impacts2 = []
        bs2 = []
        ls2 = []
        MajDiams2 = []
        azimuths2 = []
        for w, r, i, b, l, maj, az in zip(Lya_Ws, R_virs, impacts, bs, ls, MajDiams, azimuths):
            if float(w) <= maxEW:
                Lya_Ws2.append(w)
                R_virs2.append(r)
                impacts2.append(i)
                bs2.append(b)
                ls2.append(l)
                MajDiams2.append(maj)
                azimuths2.append(az)

        isolated_xs = np.array(adjustedIncs)
        isolated_ys = np.array(bs2)
        

        # grab the associated data and define the x and y data
        associated_Lya_Ws = L_associated['Lya_Ws']
        associated_R_virs = L_associated['R_virs']
        associated_impacts = L_associated['impacts']
        associated_bs = L_associated['bs']
        associated_ls = L_associated['ls']
        associated_MajDiams = L_associated['MajDiams']
        associated_azimuths = L_associated['azimuths']
        associated_adjustedIncs = L_associated['adjustedIncs']

        
        associated_Lya_Ws2 = []
        associated_R_virs2 = []
        associated_impacts2 = []
        associated_bs2 = []
        associated_ls2 = []
        associated_MajDiams2 = []
        associated_azimuths2 = []
        for w, r, i, b, l, maj, az in zip(associated_Lya_Ws, associated_R_virs, associated_impacts, associated_bs, associated_ls, associated_MajDiams, associated_azimuths):
            if float(w) <= maxEW:
                associated_Lya_Ws2.append(w)
                associated_R_virs2.append(r)
                associated_impacts2.append(i)
                associated_bs2.append(b)
                associated_ls2.append(l)
                associated_MajDiams2.append(maj)
                associated_azimuths2.append(az)
        
        associated_xs = np.array(associated_adjustedIncs)
        associated_ys = np.array(associated_bs2)
        
        
        # grab the two data and define the x and y data
        two_Lya_Ws = L_two['Lya_Ws']
        two_R_virs = L_two['R_virs']
        two_impacts = L_two['impacts']
        two_bs = L_two['bs']
        two_ls = L_two['ls']
        two_MajDiams = L_two['MajDiams']
        two_azimuths = L_two['azimuths']
        two_adjustedIncs = L_two['adjustedIncs']

        two_Lya_Ws2 = []
        two_R_virs2 = []
        two_impacts2 = []
        two_bs2 = []
        two_ls2 = []
        two_MajDiams2 = []
        two_azimuths2 = []
        for w, r, i, b, l, maj, az in zip(two_Lya_Ws, two_R_virs, two_impacts, two_bs, two_ls, two_MajDiams, two_azimuths):
            if float(w) <= maxEW:
                two_Lya_Ws2.append(w)
                two_R_virs2.append(r)
                two_impacts2.append(i)
                two_bs2.append(b)
                two_ls2.append(l)
                two_MajDiams2.append(maj)
                two_azimuths2.append(az)

        two_xs = np.array(two_adjustedIncs)
        two_ys = np.array(two_bs2)
        
        
        # grab the two_plus data and define the x and y data
        three_Lya_Ws = L_two_plus['Lya_Ws']
        three_R_virs = L_two_plus['R_virs']
        three_impacts = L_two_plus['impacts']
        three_bs = L_two_plus['bs']
        three_ls = L_two_plus['ls']
        three_MajDiams = L_two_plus['MajDiams']
        three_azimuths = L_two_plus['azimuths']
        three_adjustedIncs = L_two_plus['adjustedIncs']

        three_Lya_Ws2 = []
        three_R_virs2 = []
        three_impacts2 = []
        three_bs2 = []
        three_ls2 = []
        three_MajDiams2 = []
        three_azimuths2 = []
        for w, r, i, b, l, maj, az in zip(three_Lya_Ws, three_R_virs, three_impacts, three_bs, three_ls, three_MajDiams, three_azimuths):
            if float(w) <= maxEW:
                three_Lya_Ws2.append(w)
                three_R_virs2.append(r)
                three_impacts2.append(i)
                three_bs2.append(b)
                three_ls2.append(l)
                three_MajDiams2.append(maj)
                three_azimuths2.append(az)

        three_xs = np.array(three_adjustedIncs)
        three_ys = np.array(three_bs2)
        
        
        
        # grab the group data and define the x and y data
        group_Lya_Ws = L_group['Lya_Ws']
        group_R_virs = L_group['R_virs']
        group_impacts = L_group['impacts']
        group_mems = L_group['group_mems']
        group_bs = L_group['bs']
        group_ls = L_group['ls']
        group_MajDiams = L_group['MajDiams']
        group_azimuths = L_group['azimuths']
        group_adjustedIncs = L_group['adjustedIncs']

        group_Lya_Ws2 = []
        group_R_virs2 = []
        group_impacts2 = []
        group_bs2 = []
        group_ls2 = []
        group_MajDiams2 = []
        group_azimuths2 = []
        group_adjustedIncs2 = []
        for w, r, i, group, b, l, maj, az, inc in zip(group_Lya_Ws, group_R_virs, group_impacts, group_mems, group_bs, group_ls, group_MajDiams, group_azimuths, group_adjustedIncs):
            if float(group) >= 2 and float(w) <= maxEW:
                group_Lya_Ws2.append(w)
                group_R_virs2.append(r)
                group_impacts2.append(i)
                group_bs2.append(b)
                group_ls2.append(l)
                group_MajDiams2.append(maj)
                group_azimuths2.append(az)
                group_adjustedIncs2.append(inc)

        group_xs = np.array(group_adjustedIncs2)
        group_ys = np.array(group_bs2)
        
        
        
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
            label=r'$\rm Isolated~ Mean ~b$')
            
            
        # associated
#         plot1 = scatter(associated_xs,
#                         associated_ys,
#                         marker=symbol_assoc,
#                         c=color_assoc,
#                         s=markerSize,
#                         edgecolor='black',
#                         alpha=alpha_assoc,
#                         label=label_assoc)
#         
#         # histogram associated
#         bin_means, edges, binNumber = stats.binned_statistic(associated_xs,
#                                                             associated_ys,
#                                                             statistic='mean',
#                                                             bins=bins)
#         left,right = edges[:-1],edges[1:]        
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plot(X,
#             Y,
#             ls='solid',
#             color=color_assoc,
#             lw=2.0,
#             alpha=alpha_bins,
#             label=r'$\rm Assoc. ~Mean ~b$')

           
        # two
#         plot1 = scatter(two_xs,
#                         two_ys,
#                         marker=symbol_two,
#                         c=color_two,
#                         s=markerSize,
#                         edgecolor='black',
#                         alpha=alpha_two,
#                         label=label_two)
#         
#         # histogram group
#         bin_means, edges, binNumber = stats.binned_statistic(two_xs,
#                                                             two_ys,
#                                                             statistic='mean',
#                                                             bins=bins)
#         left,right = edges[:-1],edges[1:]
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plot(X,
#             Y,
#             ls='solid',
#             color=color_two,
#             lw=2.0,
#             alpha=alpha_bins,
#             label=r'$\rm Two ~Mean ~b$')
           
                        
#         # group
#         plot1 = scatter(group_xs,
#                         group_ys,
#                         marker=symbol_group,
#                         c=color_group,
#                         s=markerSize,
#                         edgecolor='black',
#                         alpha=alpha_group,
#                         label=label_group)
#         
#         # histogram group
#         bin_means, edges, binNumber = stats.binned_statistic(group_xs,
#                                                             group_ys,
#                                                             statistic='mean',
#                                                             bins=bins)
#         left,right = edges[:-1],edges[1:]
#         X = array([left,right]).T.flatten()
#         Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
#         plot(X,
#             Y,
#             ls='solid',
#             color=color_group,
#             lw=2.0,
#             alpha=alpha_bins,
#             label=r'$\rm Group ~Mean ~b$')
        
    
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(20)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(10)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        
        xlabel(r'$\rm Inc.~[Deg]$')
        ylabel(r'$\rm b~ [km s^{{-1}}]$')
        leg = ax.legend(scatterpoints=1,prop={'size':12},loc='upper left',fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,200)
#         xlim(0, 2.5)

        if plot_b_inc_mean_plus_save:
            savefig('{0}/b(inc)_mean_binSize{1}_plus1.pdf'.format(saveDirectory, binSize),format='pdf',bbox_inches='tight')
        else:
            show()

##########################################################################################
##########################################################################################


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


if __name__=="__main__":
    # do the work
    main()
    