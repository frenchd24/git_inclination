#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_b_impact2.py, v 6.0 06/11/18



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
    
    # plot b as a function of impact/r_vir
    plot_b_impact_vir_MType = True
    plot_b_impact_vir_MType_save = True
    
    # plot b as a function of impact 
    plot_b_impact_MType = True
    plot_b_impact_MType_save = True
    
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


#########################################################################################
#########################################################################################
##########################################################################################
##########################################################################################


#########################################################################################
#########################################################################################
##########################################################################################
##########################################################################################
    
    if plot_b_impact_vir_MType:
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
        bs2_S = []

        Lya_Ws2_S0 = []
        R_virs2_S0 = []
        impacts2_S0 = []
        bs2_S0 = []

        Lya_Ws2_E = []
        R_virs2_E = []
        impacts2_E = []
        bs2_E = []

        Lya_Ws2_I = []
        R_virs2_I = []
        impacts2_I = []
        bs2_I = []

        Lya_Ws2_other = []
        R_virs2_other = []
        impacts2_other = []
        bs2_other = []
        for w, r, i, name, b in zip(Lya_Ws, R_virs, impacts, Names, bs):
            mtype = morph_dict[name]
            
            if float(w) <= max_EW and float(w) >= min_EW:
                if mtype == 'e':
                    Lya_Ws2_E.append(w)
                    R_virs2_E.append(r)
                    impacts2_E.append(i)
                    bs2_E.append(b)
                    print 'e: ',name
                
                elif mtype == 's0':
                    Lya_Ws2_S0.append(w)
                    R_virs2_S0.append(r)
                    impacts2_S0.append(i)
                    bs2_S0.append(b)
                    print 's0: ',name
                    
                elif mtype == 'i':
                    Lya_Ws2_I.append(w)
                    R_virs2_I.append(r)
                    impacts2_I.append(i)
                    bs2_I.append(b)

                elif mtype == 'sa' or mtype == 'sb':
                    # combine sa, sb here
                    Lya_Ws2_S.append(w)
                    R_virs2_S.append(r)
                    impacts2_S.append(i)
                    bs2_S.append(b)

                else:
                    # unknown type
                    Lya_Ws2_other.append(w)
                    R_virs2_other.append(r)
                    impacts2_other.append(i)
                    bs2_other.append(b)

        
        isolated_xs_E = np.array(impacts2_E)/np.array(R_virs2_E)
        isolated_ys_E = np.array(bs2_E)
        
        isolated_xs_S = np.array(impacts2_S)/np.array(R_virs2_S)
        isolated_ys_S = np.array(bs2_S)
        
        isolated_xs_S0 = np.array(impacts2_S0)/np.array(R_virs2_S0)
        isolated_ys_S0 = np.array(bs2_S0)
        
        isolated_xs_I = np.array(impacts2_I)/np.array(R_virs2_I)
        isolated_ys_I = np.array(bs2_I)

        isolated_xs_other = np.array(impacts2_other)/np.array(R_virs2_other)
        isolated_ys_other = np.array(bs2_other)

        # grab the associated data and define the x and y data
        associated_Lya_Ws = L_associated['Lya_Ws']
        associated_R_virs = L_associated['R_virs']
        associated_impacts = L_associated['impacts']
        associated_MTypes = L_associated['MTypes']
        associated_Names = L_associated['Names']
        associated_bs = L_associated['bs']

        associated_Lya_Ws2_S = []
        associated_R_virs2_S = []
        associated_impacts2_S = []
        associated_bs2_S = []

        associated_Lya_Ws2_S0 = []
        associated_R_virs2_S0 = []
        associated_impacts2_S0 = []
        associated_bs2_S0 = []

        associated_Lya_Ws2_E = []
        associated_R_virs2_E = []
        associated_impacts2_E = []
        associated_bs2_E = []

        associated_Lya_Ws2_I = []
        associated_R_virs2_I = []
        associated_impacts2_I = []
        associated_bs2_I = []

        associated_Lya_Ws2_other = []
        associated_R_virs2_other = []
        associated_impacts2_other = []
        associated_bs2_other = []
        for w, r, i, name, b in zip(associated_Lya_Ws, associated_R_virs, associated_impacts, associated_Names, associated_bs):
            mtype = morph_dict[name]
            
            if float(w) <= max_EW and float(w) >= min_EW:            
                if mtype == 'e':
                    associated_Lya_Ws2_E.append(w)
                    associated_R_virs2_E.append(r)
                    associated_impacts2_E.append(i)
                    associated_bs2_E.append(b)

                    print 'e: ', name
        
                elif mtype == 's0':
                    associated_Lya_Ws2_S0.append(w)
                    associated_R_virs2_S0.append(r)
                    associated_impacts2_S0.append(i)
                    associated_bs2_S0.append(b)

                    print 's0: ', name
        
                elif mtype == 'i':
                    associated_Lya_Ws2_I.append(w)
                    associated_R_virs2_I.append(r)
                    associated_impacts2_I.append(i)
                    associated_bs2_I.append(b)

                elif mtype == 'sa' or mtype == 'sb':
                    # combine sa, sb here
                    associated_Lya_Ws2_S.append(w)
                    associated_R_virs2_S.append(r)
                    associated_impacts2_S.append(i)
                    associated_bs2_S.append(b)

                else:
                    # unknown type
                    associated_Lya_Ws2_other.append(w)
                    associated_R_virs2_other.append(r)
                    associated_impacts2_other.append(i)
                    associated_bs2_other.append(b)

        associated_xs_E = np.array(associated_impacts2_E)/np.array(associated_R_virs2_E)
        associated_ys_E = np.array(associated_bs2_E)
        
        associated_xs_S = np.array(associated_impacts2_S)/np.array(associated_R_virs2_S)
        associated_ys_S = np.array(associated_bs2_S)
        
        associated_xs_S0 = np.array(associated_impacts2_S0)/np.array(associated_R_virs2_S0)
        associated_ys_S0 = np.array(associated_bs2_S0)
        
        associated_xs_I = np.array(associated_impacts2_I)/np.array(associated_R_virs2_I)
        associated_ys_I = np.array(associated_bs2_I)
        
        associated_xs_other = np.array(associated_impacts2_other)/np.array(associated_R_virs2_other)
        associated_ys_other = np.array(associated_bs2_other)
        
        
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
                label=r'$\rm S-type ~\emph{b}$')
                
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
                label=r'$\rm E-type ~\emph{b}$')
            
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
                label=r'$\rm S0-type ~\emph{b}$')
                
                
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
                label=r'$\rm I-type ~\emph{b}$')
                
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
                label=r'$\rm ?-type ~\emph{b}$')
        
        
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
            print "rho/rvir linregress(a, b) = ",linregress(x_all,y_all)
            print
    
        
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %.2f$')
        minorLocator   = MultipleLocator(0.25)
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
        
        
        xlabel(r'$\rm \rho / R_{{vir}}$')
        ylabel(r'$\rm \emph{b} ~ [km s^{{-1}}]$')
        leg = ax.legend(scatterpoints=1,prop={'size':12},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,140)
        xlim(0, 2.5)

        if plot_b_impact_vir_MType_save:
            savefig('{0}/b(impact_vir)_MType_binSize{1}_EWcut{2}-{3}_dataset{4}.pdf'.format(saveDirectory, binSize, min_EW, max_EW, data_set),format='pdf',bbox_inches='tight')
        else:
            show()
            
            
##########################################################################################
##########################################################################################
    
    if plot_b_impact_MType:
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
        bs2_S = []

        Lya_Ws2_S0 = []
        R_virs2_S0 = []
        impacts2_S0 = []
        bs2_S0 = []

        Lya_Ws2_E = []
        R_virs2_E = []
        impacts2_E = []
        bs2_E = []


        Lya_Ws2_I = []
        R_virs2_I = []
        impacts2_I = []
        bs2_I = []

        Lya_Ws2_other = []
        R_virs2_other = []
        impacts2_other = []
        bs2_other = []
        for w, r, i, name, b in zip(Lya_Ws, R_virs, impacts, Names, bs):
            mtype = morph_dict[name]
            
            if float(w) <= max_EW and float(w) >= min_EW:
                
                if mtype == 'e':
                    Lya_Ws2_E.append(w)
                    R_virs2_E.append(r)
                    impacts2_E.append(i)
                    bs2_E.append(b)

                elif mtype == 's0':
                    Lya_Ws2_S0.append(w)
                    R_virs2_S0.append(r)
                    impacts2_S0.append(i)
                    bs2_S0.append(b)

                elif mtype == 'i':
                    Lya_Ws2_I.append(w)
                    R_virs2_I.append(r)
                    impacts2_I.append(i)
                    bs2_I.append(b)

                elif mtype == 'sa' or mtype == 'sb':
                    # combine sa, sb here
                    Lya_Ws2_S.append(w)
                    R_virs2_S.append(r)
                    impacts2_S.append(i)
                    bs2_S.append(b)

                else:
                    # unknown type
                    Lya_Ws2_other.append(w)
                    R_virs2_other.append(r)
                    impacts2_other.append(i)
                    bs2_other.append(b)

        
        isolated_xs_E = np.array(impacts2_E)
        isolated_ys_E = np.array(bs2_E)
        
        isolated_xs_S = np.array(impacts2_S)
        isolated_ys_S = np.array(bs2_S)
        
        isolated_xs_S0 = np.array(impacts2_S0)
        isolated_ys_S0 = np.array(bs2_S0)
        
        isolated_xs_I = np.array(impacts2_I)
        isolated_ys_I = np.array(bs2_I)
        
        isolated_xs_other = np.array(impacts2_other)
        isolated_ys_other = np.array(bs2_other)

        # grab the associated data and define the x and y data
        associated_Lya_Ws = L_associated['Lya_Ws']
        associated_R_virs = L_associated['R_virs']
        associated_impacts = L_associated['impacts']
        associated_MTypes = L_associated['MTypes']
        associated_Names = L_associated['Names']
        associated_bs = L_associated['bs']

        associated_Lya_Ws2_S = []
        associated_R_virs2_S = []
        associated_impacts2_S = []
        associated_bs2_S = []

        associated_Lya_Ws2_S0 = []
        associated_R_virs2_S0 = []
        associated_impacts2_S0 = []
        associated_bs2_S0 = []

        associated_Lya_Ws2_E = []
        associated_R_virs2_E = []
        associated_impacts2_E = []
        associated_bs2_E = []

        associated_Lya_Ws2_I = []
        associated_R_virs2_I = []
        associated_impacts2_I = []
        associated_bs2_I = []

        associated_Lya_Ws2_other = []
        associated_R_virs2_other = []
        associated_impacts2_other = []
        associated_bs2_other = []

        for w, r, i, name, b in zip(associated_Lya_Ws, associated_R_virs, associated_impacts, associated_Names, associated_bs):
            mtype = morph_dict[name]
            
            if float(w) <= max_EW and float(w) >= min_EW:            
                if mtype == 'e':
                    associated_Lya_Ws2_E.append(w)
                    associated_R_virs2_E.append(r)
                    associated_impacts2_E.append(i)
                    associated_bs2_E.append(b)

                elif mtype == 's0':
                    associated_Lya_Ws2_S0.append(w)
                    associated_R_virs2_S0.append(r)
                    associated_impacts2_S0.append(i)
                    associated_bs2_S0.append(b)

                elif mtype == 'i':
                    associated_Lya_Ws2_I.append(w)
                    associated_R_virs2_I.append(r)
                    associated_impacts2_I.append(i)
                    associated_bs2_I.append(b)

                elif mtype == 'sa' or mtype == 'sb':
                    # combine sa, sb here
                    associated_Lya_Ws2_S.append(w)
                    associated_R_virs2_S.append(r)
                    associated_impacts2_S.append(i)
                    associated_bs2_S.append(b)

                else:
                    # unknown types
                    associated_Lya_Ws2_other.append(w)
                    associated_R_virs2_other.append(r)
                    associated_impacts2_other.append(i)
                    associated_bs2_other.append(b)

        
        associated_xs_E = np.array(associated_impacts2_E)
        associated_ys_E = np.array(associated_bs2_E)
        
        associated_xs_S = np.array(associated_impacts2_S)
        associated_ys_S = np.array(associated_bs2_S)
        
        associated_xs_S0 = np.array(associated_impacts2_S0)
        associated_ys_S0 = np.array(associated_bs2_S0)
        
        associated_xs_I = np.array(associated_impacts2_I)
        associated_ys_I = np.array(associated_bs2_I)
        
        associated_xs_other = np.array(associated_impacts2_other)
        associated_ys_other = np.array(associated_bs2_other)
        
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
                label=r'$\rm S-type ~\emph{b}$')
            
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
                label=r'$\rm E-type ~\emph{b}$')
            
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
                label=r'$\rm S0-type ~\emph{b}$')


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
                label=r'$\rm I-type ~\emph{b}$')
        
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
                label=r'$\rm ?-type ~\emph{b}$')



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

        # x-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)
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
        
        
        xlabel(r'$\rm \rho ~[kpc]$')
        ylabel(r'$\rm \emph{b}~[km s^{{-1}}]$')
        leg = ax.legend(scatterpoints=1,prop={'size':12},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,140)
        xlim(0, 500)

        if plot_b_impact_MType_save:
            savefig('{0}/b(impact)_MType_binSize{1}_EWcut{2}-{3}_dataset{4}.pdf'.format(saveDirectory, binSize, min_EW, max_EW, data_set),format='pdf',bbox_inches='tight')
        else:
            show()
            

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


if __name__=="__main__":
    # do the work
    main()
    