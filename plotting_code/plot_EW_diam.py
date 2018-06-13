#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

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
    # plot EW as a function of diameter, include median
    # EW histograms
    plot_EW_diam_median = False
    plot_EW_diam_median_save = False

    # plot EW as a function of R_vir, include median
    # EW histograms
    plot_EW_vir_median = True
    plot_EW_vir_median_save = True

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
    # plot equivalent width as a function of virial radius for red vs blue
    # shifted absorption, include average histograms
    #
    
    plotW_vir_avg = False
    save = False
    
    if plotW_vir_avg:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        alpha = 0.7
        binSize = 50
        markerSize = 60

        
        bSymbol = 'D'
        rSymbol = 'o'
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        
        placeArrayr = zeros(7)
        placeCountr = zeros(7)
        placeArrayb = zeros(7)
        placeCountb = zeros(7)
        redW = []
        blueW = []
        redVir = []
        blueVir = []
        
        virList = [int(round(float(v),0)) for v in virList]
        
        
        for d,v,w,m in zip(difList,virList,lyaWList,majList):
            # check if all the values are good
            if isNumber(d) and isNumber(v) and isNumber(w) and isNumber(m):
                if d!=-99 and v!=-99 and w!=-99 and m!=-99:
                    w = float(w)
                    m = float(m)
                    d = float(d)
                    if d>0:
                        # galaxy is behind absorber, so gas is blue shifted
                        color = 'Blue'
                        symbol = bSymbol
                        
                        blueW.append(w)
                        blueVir.append(v)
                        
                        # which bin does it belong too?
                        place = v/binSize
                        print 'place: ',place
                        placeArrayb[place] += float(w)
                        print 'placeArrayb: ',placeArrayb
                        placeCountb[place] +=1.
                        print 'placecountb: ',placeCountb
                        
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(v,w,marker=symbol,alpha=alpha,c='Blue',\
                            s=markerSize)

                    if d<0:
                        # gas is red shifted compared to galaxy
                        color = 'Red'
                        symbol = rSymbol
                        
                        redW.append(w)
                        redVir.append(v)
                        
                        # which bin does it belong too?
                        place = v/binSize
                        placeArrayr[place] += float(w)
                        placeCountr[place] +=1.
                        
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(v,w,marker=symbol,alpha=alpha,c='Red',\
                            s=markerSize)
            
                    plot1 = scatter(v,w,marker=symbol,alpha=alpha,c=color,s=markerSize)
                    
        rHist = placeArrayr/placeCountr
        print 'rHist: ',rHist
        bHist = placeArrayb/placeCountb
        print 'bHist: ',bHist
        
        totalrHist = []
        totalrVir = []
        totalbHist = []
        totalbVir = []
        
        for r,v in zip(rHist,arange(0,max(virList),binSize)):
            if not isNumber(r):
                r = 0
            
            totalrHist.append(r)
            totalrHist.append(r)

            totalrVir.append(v)
            totalrVir.append(v+binSize)
            
        for b,v in zip(bHist,arange(0,max(virList),binSize)):
            if not isNumber(b):
                b = 0
            totalbHist.append(b)
            totalbHist.append(b)

            totalbVir.append(v)
            totalbVir.append(v+binSize)
        
        print 'totalrVir: ',totalrVir
        print 'totalrHist: ',totalrHist
        print
        print 'totalbVir: ',totalbVir
        print 'totalbHist: ',totalbHist
        print
        
        
        # x-axis
        majorLocator   = MultipleLocator(50)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(25)
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
        
        # bins are in Rvir
        bins = arange(150,400,binSize)
           
#         bin_means,edges,binNumber = stats.binned_statistic(array(redVir), array(redW), statistic='mean', bins=bins)
#         left,right = edges[:-1],edges[1:]
#         X = array([left,right]).T.flatten()
#         Y = array([bin_means,bin_means]).T.flatten()
#         plot(X,Y, c='red',ls='dotted',lw=2.5,alpha=alpha,label=r'\rm $Mean ~Redshifted ~EW$')
#             
#         bin_means,edges,binNumber = stats.binned_statistic(array(blueVir), array(blueW), statistic='mean', bins=bins)
#         left,right = edges[:-1],edges[1:]
#         X = array([left,right]).T.flatten()
#         Y = array([bin_means,bin_means]).T.flatten()
#         plot(X,Y, c='blue',ls='dashed',lw=1.5,alpha=alpha,label=r'\rm $Mean ~Blueshifted ~EW$')
        
        plot2 = ax.plot(totalrVir,totalrHist,c='Red',lw=2.5,ls='dotted',label=r'$\rm Mean ~ Redshifted ~ EW$')
        plot3 = ax.plot(totalbVir,totalbHist,c='Blue',lw=1.5,ls='dashed',label=r'$\rm Mean ~ Blueshifted ~ EW$')
        
        xlabel(r'$\rm R_{vir} ~ [kpc]$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        ax.legend(scatterpoints=1,prop={'size':15},loc=2,fancybox=True)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1200)
        xlim(150,350)

        if save:
            savefig('{0}/W(vir)_avgHistograms2_maxEnv{1}.pdf'.format(saveDirectory,maxEnv),format='pdf',bbox_inches='tight')
        else:
            show()



##########################################################################################
##########################################################################################
    
##########################################################################################
##########################################################################################
    
    if plot_EW_diam_median:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        
        second = 'Group'
        
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
        
        binSize = 20
        bins = arange(0, 100, binSize)
        
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
        isolated_xs = np.array(MajDiams)
        isolated_ys = Lya_Ws

        # grab the associated data and define the x and y data
        associated_Lya_Ws = L_associated['Lya_Ws']
        associated_R_virs = L_associated['R_virs']
        associated_impacts = L_associated['impacts']
        associated_MajDiams = L_associated['MajDiams']
        associated_xs = np.array(associated_MajDiams)
        associated_ys = associated_Lya_Ws
        
        
        # grab the two data and define the x and y data
        two_Lya_Ws = L_two['Lya_Ws']
        two_R_virs = L_two['R_virs']
        two_impacts = L_two['impacts']
        two_MajDiams = L_two['MajDiams']
        two_xs = np.array(two_MajDiams)
        two_ys = np.array(two_Lya_Ws)
        
        
        # grab the two_plus data and define the x and y data
        three_Lya_Ws = L_two_plus['Lya_Ws']
        three_R_virs = L_two_plus['R_virs']
        three_impacts = L_two_plus['impacts']
        three_MajDiams = L_two_plus['MajDiams']
        three_xs = np.array(three_MajDiams)
        three_ys = np.array(three_Lya_Ws)
        
        
        # grab the group data and define the x and y data
        group_Lya_Ws = L_group['Lya_Ws']
        group_R_virs = L_group['R_virs']
        group_impacts = L_group['impacts']
        group_MajDiams = L_group['MajDiams']
        group_xs = np.array(group_MajDiams)
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
        
        
        xlabel(r'$\rm D~[kpc]$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        leg = ax.legend(scatterpoints=1,prop={'size':14},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0, 100)

        if plot_EW_diam_median_save:
            savefig('{0}/W(diam)_median_plus_{1}_binSize{2}.pdf'.format(saveDirectory, second, binSize),format='pdf',bbox_inches='tight')
        else:
            show()


##########################################################################################
##########################################################################################
    
    if plot_EW_vir_median:
        fig = figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(111)
        
        countb = 0
        countr = 0
        count = -1
        
        second = 'Group'
        
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
        
        binSize = 20
        bins = arange(0, 100, binSize)
        
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
        isolated_xs = np.array(MajDiams)
        isolated_ys = Lya_Ws

        # grab the associated data and define the x and y data
        associated_Lya_Ws = L_associated['Lya_Ws']
        associated_R_virs = L_associated['R_virs']
        associated_impacts = L_associated['impacts']
        associated_MajDiams = L_associated['MajDiams']
        associated_xs = np.array(associated_R_virs)
        associated_ys = associated_Lya_Ws
        
        
        # grab the two data and define the x and y data
        two_Lya_Ws = L_two['Lya_Ws']
        two_R_virs = L_two['R_virs']
        two_impacts = L_two['impacts']
        two_MajDiams = L_two['MajDiams']
        two_xs = np.array(two_R_virs)
        two_ys = np.array(two_Lya_Ws)
        
        
        # grab the two_plus data and define the x and y data
        three_Lya_Ws = L_two_plus['Lya_Ws']
        three_R_virs = L_two_plus['R_virs']
        three_impacts = L_two_plus['impacts']
        three_MajDiams = L_two_plus['MajDiams']
        three_xs = np.array(three_R_virs)
        three_ys = np.array(three_Lya_Ws)
        
        
        # grab the group data and define the x and y data
        group_Lya_Ws = L_group['Lya_Ws']
        group_R_virs = L_group['R_virs']
        group_impacts = L_group['impacts']
        group_MajDiams = L_group['MajDiams']
        group_xs = np.array(group_R_virs)
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
        
        
        xlabel(r'$\rm R_{{vir}}~[kpc]$')
        ylabel(r'$\rm Equivalent ~ Width ~ [m\AA]$')
        leg = ax.legend(scatterpoints=1,prop={'size':14},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(0,1300)
        xlim(0, 100)

        if plot_EW_vir_median_save:
            savefig('{0}/W(vir)_median_plus_{1}_binSize{2}.pdf'.format(saveDirectory, second, binSize),format='pdf',bbox_inches='tight')
        else:
            show()



##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

if __name__=="__main__":
    # do the work
    main()
    