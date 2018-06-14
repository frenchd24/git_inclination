#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_az_hist2.py, v 5.7 06/10/18



Comes from:
$Id:  plotAzHist2.py, v 5.6 01/05/18

This is the plotAzHist bit from histograms3.py. Now is separated, and loads in a pickle
file of the relevant data, as created by "buildDataLists.py"

Plots a histogram of azimuth angles for red and blue shifted absorbers

Previous (from histograms3.py):
    Plot some stuff for the 100largest initial results

    Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

    Updated for the pilot paper (05/06/15)

v5: updated to work with the new, automatically updated LG_correlation_combined5.csv
    (12/04/15) more on (12/28/2015)
    
v5.1: updated to work with LG_correlation_combined5_8_edit2.csv and l_min = 0.001
    (02/22/2016)

v5.2: remake plots with v_hel instead of vcorr (4/21/16)

v5.3: remake plots with newest large galaxy sample (07/13/16) -> /plots4/

v5.4: include ability to limit results based on environment number (7/14/16)

v5.5: update for LG_correlation_combined5_11_25cut_edit4.csv -> /plots5/ (9/26/16)

v5.6: updated for AAS_2017 (01/05/18)
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

# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# ## for Palatino and other serif fonts use:
# #rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)

from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})

# plt.clf()
# plt.rcdefaults()
# fontScale = 18
# params = {'axes.labelsize': fontScale,
# 'axes.titlesize': fontScale,
# 'font.size': fontScale,
# 'xtick.labelsize': fontScale,
# 'ytick.labelsize': fontScale,
# 'font.weight': 450,
# 'axes.labelweight': 450,
# 'xtick.major.size': 4,
# 'xtick.major.width': 1.2,
# 'xtick.minor.size': 3,
# 'xtick.minor.width': 1.2,
# 'ytick.major.size': 4,
# 'ytick.major.width': 1.2,
# 'ytick.minor.size': 3,
# 'ytick.minor.width': 1.2
# }


fontScale = 24
rc('text', usetex=True)
rc('font', size=24, family='serif', weight='normal')
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
    # Plot the azimuth distribution with red and blue shifted absorbers separated but
    # overlaid, and the combined sample plotted over that as a step hist
    plot_az_hist_all_overlaid = False
    plot_az_hist_all_overlaid_save = False
    
    # Plot azimuth distribution with red and blue shifted absorbers separated but overlaid
    plot_az_hist_dif = False
    plot_az_hist_dif_save = False
    
    # plot all azimuth distributions for each subset (e.g., L_associated, L_two, etc)
    plot_az_all_subsets = True
    plot_az_all_subsets_save = True
    

    color_blue = '#436bad'      # french blue
    color_red = '#ec2d01'     # tomato red
    

    if getpass.getuser() == 'frenchd':

#         pickleFilename = '/Users/frenchd/Research/inclination/git_inclination/picklePilot_plusSALT_14.p'
#         gtPickleFilename = '/Users/frenchd/Research/inclination/git_inclination/pickleGT.p'
#         saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/plotting_code/figs'

#         gtPickleFilename = '/Users/frenchd/Research/inclination/git_inclination/pickleGT_filteredAll.p'
        gtPickleFilename = '/Users/frenchd/Research/GT_update2/pickleGT_filteredAll.p'

        saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/plotting_code/figs/'
        
        isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/isolated4.p'
        L_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_isolated4.p'
        L_associated_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated4.p'
        L_associated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated4.p'
        L_nonassociated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated4.p'
        L_two_filename = '/Users/frenchd/Research/inclination/git_inclination/L_two4.p'
        L_two_plus_filename = '/Users/frenchd/Research/inclination/git_inclination/L_two_plus4.p'
        L_group_filename = '/Users/frenchd/Research/inclination/git_inclination/L_group4.p'


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
    
    
    # save each plot?
    save = False
    
#     results = open(resultsFilename,'rU')
#     reader = csv.DictReader(results)
    
#     WS = open(WS09data,'rU')
#     WSreader = csv.DictReader(WS,delimiter=';')
    
    virInclude = False
    cusInclude = False
    finalInclude = 1
    
    maxEnv = 3000
    minL = 0.001
    
    # if match, then the includes in the file have to MATCH the includes above. e.g., if 
    # virInclude = False, cusInclude = True, finalInclude = False, then only systems
    # matching those three would be included. Otherwise, all cusInclude = True would be included
    # regardless of the others
    match = False
    
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
    if plot_az_hist_dif:
    
        fig = figure(figsize=(10,6))
        subplots_adjust(hspace=0.200)

#         bins = [0,15,30,45,60,75,90]
#         bins = arange(0,90,10)
        binsize = 10
        bins = arange(0,100,binsize)
        blue = []
        red = []
        alpha_red = 0.70
        alpha_blue = 0.75
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        for Lya_v, Vhel, a, w, e in zip(Lya_vs, Vhels, azimuths, Lya_Ws, e_Lya_Ws):
            vel_dif = Lya_v - Vhel

            if vel_dif >=0:
                # blue shifted galaxy, but absorber is REDSHIFTED
#                 print 'd: ',d
                red.append(a)
                redLya.append(w)
                redLyaErr.append(e)
            else:
                # red shifted absorber, but galaxy is BLUESHIFTED
                blue.append(a)
                blueLya.append(w)
                blueLyaErr.append(e)
        
        print 'stats-----'
        print 'max red: ',max(red)
        print 'min red: ',min(red)
        print 'max blue:' ,max(blue)
        print 'min blue: ',min(blue)
        
        ax = fig.add_subplot(211)        
        hist(red,
        bins=bins,
        histtype='bar',
        edgecolor='black',
        color=color_red,
        alpha=alpha_red,
        label='Redshifted absorbers')
        
        ylabel("Number")
        ylim(0,7)
        xlim(0,90)
        legend(scatterpoints=1,prop={'size':14},loc=2,fancybox=True)

        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        

        ax = fig.add_subplot(212)
        hist(blue,
        bins=bins,
        histtype='bar',
        edgecolor='black',
        color=color_blue,
        alpha=alpha_blue,
        label='Blueshifted absorbers')
        
        ylabel('Number')
        xlabel("Azimuth (deg)")
        xlim(0,90)
        ylim(0,7)
        
        legend(scatterpoints=1,prop={'size':14},loc=2,fancybox=True)
        
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        

        tight_layout()

        if plot_az_hist_dif_save:
            savefig('{0}/hist(azimuth)_dif_{1}.pdf'.format(saveDirectory, binsize),format='pdf')
        else:
            show()

            
#########################################################################################
#########################################################################################
    
    if plot_az_hist_all_overlaid:
    
        fig = figure(figsize=(10,5))
#         subplots_adjust(hspace=0.200)
        
        alpha_red = 0.65
        alpha_blue = 0.7
        alpha_both = 0.8

        bins = arange(0,100,10)
        blue = []
        red = []
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []
        
        for Lya_v, Vhel, a, w, e in zip(Lya_vs, Vhels, azimuths, Lya_Ws, e_Lya_Ws):
            vel_dif = Lya_v - Vhel
            
            print 'vel_dif: ',vel_dif

            if vel_dif >=0:
                # blue shifted galaxy, but absorber is REDSHIFTED
#                 print 'd: ',d
                red.append(a)
                redLya.append(w)
                redLyaErr.append(e)
            else:
                # red shifted absorber, but galaxy is BLUESHIFTED
                blue.append(a)
                blueLya.append(w)
                blueLyaErr.append(e)
        
        
        ax = fig.add_subplot(111)


        # blue
        hist(blue,
        bins=bins,
        histtype='bar',
        color=color_blue,
        hatch='/',
        lw=1.0,
        alpha=alpha_blue,
        edgecolor='black',
        label=r'$\rm Blueshifted$')

        # red
        hist(red,
        bins=bins,
        histtype='bar',
        color=color_red,
        lw=1.0,
        alpha=alpha_red,
        edgecolor='black',
        label=r'$\rm Redshifted$')

        # together
        hist(azimuths,
        bins=bins,
        histtype='step',
        color='black',
        lw=1.0,
        alpha=alpha_both,
        edgecolor='black',
        label=r'$\rm All ~ Associated$')


        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y axis
        majorLocator   = MultipleLocator(2)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        ylabel(r'$\rm Number$')
        xlabel(r'$\rm Azimuth ~ [deg]$')
        xlim(0,90)
#         ylim(0,20)
        legend(scatterpoints=1,prop={'size':14},loc=2,fancybox=True)


#         tight_layout()

        if plot_az_hist_all_overlaid_save:
            savefig('{0}/hist(azimuth)_overlaid_all.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()

#########################################################################################
#########################################################################################    
    if plot_az_all_subsets:
    
        fig = figure(figsize=(10,6))
        subplots_adjust(hspace=0.200)
        bins = arange(0,100,10)
        
        L_associated_azimuths = L_associated['azimuths']
        L_nonassociated_azimuths = L_nonassociated['azimuths']
        L_two_azimuths = L_two['azimuths']
        L_two_plus_azimuths = L_two_plus['azimuths']
        L_group_azimuths = L_group['azimuths']
        
        alpha_L_associated_isolated = 1.0
        alpha_L_associated = 0.8
        alpha_L_two = 1.0
        alpha_L_two_plus = 1.0
        alpha_L_group = 1.0
        
        lw = 2.0
        
        color_blue = '#377eb8' # blue
        color_black = 'black'
        color_orange = '#ff7f00' # orange
#         color_L_two_plus = '#984ea3'
        color_grey = 'grey'

#         color_L_associated_isolated = '#e41a1c' # red
#         color_L_associated = '#377eb8' # blue
#         color_L_two = '#4daf4a' # green
#         color_L_two_plus = '#984ea3' # purple
#         color_L_group = '#ff7f00' # orange
        

        ax = fig.add_subplot(111)
        
        all_associated_azimuths = azimuths + L_associated_azimuths

        # L_associated_isolated
        hist(all_associated_azimuths,
        bins=bins,
        histtype='step',
        color=color_black,
        lw=lw,
        alpha=alpha_L_associated_isolated,
        label=r'$\rm All~ Assoc.$')


        # L_associated
#         hist(L_associated_azimuths,
#         bins=bins,
#         histtype='step',
#         color=color_L_associated,
#         lw=lw,
#         alpha=alpha_L_associated,
#         label=r'$\rm L\_associated$')


        # L_two
        hist(L_two_azimuths + L_two_plus_azimuths,
        bins=bins,
        histtype='step',
        color=color_orange,
        lw=lw,
        alpha=alpha_L_two,
        label=r'$\rm L\_two+$')
        
        # L_two_plus_azimuths
#         hist(L_two_plus_azimuths,
#         bins=bins,
#         histtype='step',
#         color=color_L_two_plus,
#         lw=lw,
#         alpha=alpha_L_two_plus,
#         label=r'$\rm L\_two\_plus$')
        
        # L_group_azimuths
#         hist(L_group_azimuths,
#         bins=bins,
#         histtype='step',
#         color=color_L_group,
#         lw=lw,
#         alpha=alpha_L_group,
#         label=r'$\rm L\_group$')


        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y axis
        majorLocator   = MultipleLocator(4)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(2)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        ylabel(r'$\rm Number$')
        xlabel(r'$\rm Azimuth ~ [deg]$')
        xlim(0,90)
        ylim(0,40)
        legend(scatterpoints=1,prop={'size':14},loc=2,fancybox=True)



        if plot_az_all_subsets_save:
            savefig('{0}/hist(azimuth)_all_assoc_subsets.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


if __name__=="__main__":
    # do the work
    main()
    