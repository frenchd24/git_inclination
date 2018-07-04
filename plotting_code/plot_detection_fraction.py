#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_detection_fraction.py, v 1.0 07/02/18

Plot detection fraction as a function of impact parameter and some other shit

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
    # plot dopplar b parameter as a function of impact w/o splitting, add other sets if you want
    plot_detection_fraction_both = False
    plot_detection_fraction_both_save = False
    
    # plot dopplar b parameter as a function of impact w/o splitting, add other sets if you want
    plot_detection_fraction_impact = True
    plot_detection_fraction_impact_save = True
    
    # plot dopplar b parameter as a function of impact w/o splitting, add other sets if you want
    plot_detection_fraction_likelihood = True
    plot_detection_fraction_likelihood_save = True
    
    # plot_number = 1 for just the isolated sample, =2 adds the associated, =3 adds two+
    # =4 adds groups with 2 or more members
    Lstar_min = 0.5

    # some colors
    color_blue = '#436bad'      # french blue
    color_red = '#ec2d01'     # tomato red

    if getpass.getuser() == 'frenchd':

#         gtPickleFilename = '/Users/frenchd/Research/GT_update2/pickleGT_filteredAll.p'

        saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/plotting_code/figs/'
        
#         detection_fraction_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction3.p'
        detection_fraction_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_lstarmin_05.p'

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # pickle file for the whole galaxy table:
#     gtPickleFile = open(gtPickleFilename,'rU')
#     gtDict = pickle.load(gtPickleFile)
#     gtPickleFile.close()
    
    
    # open all the pickle files
    detection_fraction_file = open(detection_fraction_filename,'r')

    # unload the data from them
    full_dict = pickle.load(detection_fraction_file)

    # close the files
    detection_fraction_file.close()
    
    
    # unload all the stuff
 
    # \Delta v lists for impact parameter
    
    dv_1000 = full_dict['dv_1000']
    dv_750 = full_dict['dv_750']
    dv_500 = full_dict['dv_500']
    dv_400 = full_dict['dv_400']
    dv_300 = full_dict['dv_300']
    dv_200 = full_dict['dv_200']
    dv_100 = full_dict['dv_100']
    dv_50 = full_dict['dv_50']

    # now impact parameter detection counts
    dv400_imp1000_det = full_dict['dv400_imp1000_det']
    dv400_imp750_det = full_dict['dv400_imp750_det']
    dv400_imp500_det = full_dict['dv400_imp500_det']
    dv400_imp400_det = full_dict['dv400_imp400_det']
    dv400_imp300_det = full_dict['dv400_imp300_det']
    dv400_imp200_det = full_dict['dv400_imp200_det']
    dv400_imp100_det = full_dict['dv400_imp100_det']
    dv400_imp50_det = full_dict['dv400_imp50_det']
    
    
    # now impact parameter non-detection counts
    dv400_imp1000_non = full_dict['dv400_imp1000_non']
    dv400_imp750_non = full_dict['dv400_imp750_non']
    dv400_imp500_non = full_dict['dv400_imp500_non']
    dv400_imp400_non = full_dict['dv400_imp400_non']
    dv400_imp300_non = full_dict['dv400_imp300_non']
    dv400_imp200_non = full_dict['dv400_imp200_non']
    dv400_imp100_non = full_dict['dv400_imp100_non']
    dv400_imp50_non = full_dict['dv400_imp50_non']


    # now for likelihood thresholds
    dv_l001 = full_dict['dv_l001']
    dv_l005 = full_dict['dv_l005']
    dv_l01 = full_dict['dv_l01']
    dv_l05 = full_dict['dv_l05']
    dv_l1 = full_dict['dv_l1']
    dv_l5 = full_dict['dv_l5']
    dv_l75 = full_dict['dv_l75']

    
    # now for likelihood detections
    dv400_l001_det = full_dict['dv400_l001_det']
    dv400_l005_det = full_dict['dv400_l005_det']
    dv400_l01_det = full_dict['dv400_l01_det']
    dv400_l05_det = full_dict['dv400_l05_det']
    dv400_l1_det = full_dict['dv400_l1_det']
    dv400_l5_det = full_dict['dv400_l5_det']
    dv400_l75_det = full_dict['dv400_l75_det']

    
    # now for likelihood non-detections
    dv400_l001_non = full_dict['dv400_l001_non']
    dv400_l005_non = full_dict['dv400_l005_non']
    dv400_l01_non = full_dict['dv400_l01_non']
    dv400_l05_non = full_dict['dv400_l05_non']
    dv400_l1_non = full_dict['dv400_l1_non']
    dv400_l5_non = full_dict['dv400_l5_non']
    dv400_l75_non = full_dict['dv400_l75_non']

    print 'Detection fraction for 1000 kpc: ', float(dv400_imp1000_det) / (dv400_imp1000_det + dv400_imp1000_non)
    print 'Detection fraction for 750 kpc: ', float(dv400_imp750_det) / (dv400_imp750_det + dv400_imp750_non)
    print 'Detection fraction for 500 kpc: ', float(dv400_imp500_det) / (dv400_imp500_det + dv400_imp500_non)
    print 'Detection fraction for 400 kpc: ', float(dv400_imp400_det) / (dv400_imp400_det + dv400_imp400_non)
    print 'Detection fraction for 300 kpc: ', float(dv400_imp300_det) / (dv400_imp300_det + dv400_imp300_non)
    print 'Detection fraction for 200 kpc: ', float(dv400_imp200_det) / (dv400_imp200_det + dv400_imp200_non)
    print 'Detection fraction for 100 kpc: ', float(dv400_imp100_det) / (dv400_imp100_det + dv400_imp100_non)


#########################################################################################
#########################################################################################

##########################################################################################
##########################################################################################
    
    if plot_detection_fraction_impact:
        fig = figure(figsize=(7.7,5.7))
        ax1 = fig.add_subplot(111)
        
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
        

        alpha_likelihood = 0.8
        alpha_impact = 0.8
        markerSize = 15
        lw = 2.
        
        label_likelihood = r'$\rm \mathcal{L}-Detection~Fraction$'
        label_impact = r'$\rm \rho - Detection~Fraction$'

        symbol_likelihood = 'D'
        symbol_impact = 'o'
        
        color_likelihood = color_blue
        color_impact = color_coal

        maxEW = 15000.

##########################################################################################
        # do the plotting 
        frac_imp1000 = float(dv400_imp1000_det) / float(dv400_imp1000_det + dv400_imp1000_non)
        frac_imp750 = float(dv400_imp750_det) / float(dv400_imp750_det + dv400_imp750_non)
        frac_imp500 = float(dv400_imp500_det) / float(dv400_imp500_det + dv400_imp500_non)
        frac_imp400 = float(dv400_imp400_det) / float(dv400_imp400_det + dv400_imp400_non)
        frac_imp300 = float(dv400_imp300_det) / float(dv400_imp300_det + dv400_imp300_non)
        frac_imp200 = float(dv400_imp200_det) / float(dv400_imp200_det + dv400_imp200_non)
        
        try:
            frac_imp100 = float(dv400_imp100_det) / float(dv400_imp100_det + dv400_imp100_non)
        except Exception, e:
            print 'error: ',e
            frac_imp100 = 0.
            
        try:
            frac_imp50  = float(dv400_imp50_det)  / float(dv400_imp50_det  + dv400_imp50_non)
        except Exception, e:
            print 'error: ',e
            frac_imp50 = 0.

        impact_x = [50, 100, 200, 300, 400, 500, 750, 1000]
        impact_y = [frac_imp50, frac_imp100, frac_imp200, frac_imp300, frac_imp400, frac_imp500, frac_imp750, frac_imp1000]


        # impact detection fraction
        ax1.plot(impact_x,
                impact_y,
                marker=symbol_impact,
                c=color_impact,
                ms=markerSize,
                markeredgecolor='black',
                lw = lw,
                alpha=alpha_impact,
                label=label_impact)


        # x-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)
        ax1.xaxis.set_major_locator(majorLocator)
        ax1.xaxis.set_major_formatter(majorFormatter)
        ax1.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
        majorLocator   = MultipleLocator(0.2)
        majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
        minorLocator   = MultipleLocator(0.1)
        ax1.yaxis.set_major_locator(majorLocator)
        ax1.yaxis.set_major_formatter(majorFormatter)
        ax1.yaxis.set_minor_locator(minorLocator)

        ax1.set_xlabel(r'$\rm \rho ~[kpc]$')
        
        ax1.set_ylabel(r'$\rm Detection~Fraction$')

        leg = ax1.legend(scatterpoints=1,prop={'size':12},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax1.grid(b=None,which='major',axis='both')
        ylim(0, 1.)
        xlim(0, 1000.)

        if plot_detection_fraction_impact_save:
            savefig('{0}/detection_fraction_impact_Lstarmin{1}.pdf'.format(saveDirectory, Lstar_min),format='pdf',bbox_inches='tight')
        else:
            show()

##########################################################################################
##########################################################################################
    
    if plot_detection_fraction_likelihood:
        fig = figure(figsize=(7.7,5.7))
        ax1 = fig.add_subplot(111)
        
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
        

        alpha_likelihood = 0.8
        alpha_impact = 0.8
        markerSize = 15
        lw = 2.
        
        binSize = 100
        bins = arange(0, 600, binSize)

        
        label_likelihood = r'$\rm \mathcal{L}-Detection~Fraction$'
        label_impact = r'$\rm \rho - Detection~Fraction$'

        symbol_likelihood = 'D'
        symbol_impact = 'o'
        
        color_likelihood = color_blue
        color_impact = color_coal

        maxEW = 15000.

##########################################################################################
        # do the plotting
        
        frac_l001 = float(dv400_l001_det) / float(dv400_l001_det + dv400_l001_non)
        frac_l005 = float(dv400_l005_det) / float(dv400_l005_det + dv400_l005_non)
        frac_l01  = float(dv400_l01_det)  / float(dv400_l01_det  + dv400_l01_non)
        frac_l05  = float(dv400_l05_det)  / float(dv400_l05_det  + dv400_l05_non)
        frac_l1   = float(dv400_l1_det)   / float(dv400_l1_det   + dv400_l1_non)
        
        try:
            frac_l5   = float(dv400_l5_det)   / float(dv400_l5_det   + dv400_l5_non)
        except Exception,e:
            print 'error: ',e
            frac_l5 = 0.
            
        try:
            frac_l75   = float(dv400_l75_det)   / float(dv400_l75_det   + dv400_l75_non)
        except Exception,e:
            print 'error: ',e
            frac_l75 = 0.

        likelihood_x = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 0.75]
        likelihood_y = [frac_l001, frac_l005, frac_l01, frac_l05, frac_l1, frac_l5, frac_l75]


        # likelihood detection fraction
        ax1.plot(likelihood_x,
                likelihood_y,
                marker=symbol_likelihood,
                c=color_likelihood,
                ms=markerSize,
                markeredgecolor='black',
                lw = lw,
                alpha=alpha_likelihood,
                label=label_likelihood)

        ax1.set_xlabel(r'$\rm \mathcal{L}$')
        ax1.set_xscale("log")

        
        # x-axis
#         majorLocator   = MultipleLocator(0.01)
#         majorFormatter = FormatStrFormatter(r'$\rm %d$')
#         minorLocator   = MultipleLocator()
#         ax1.xaxis.set_major_locator(majorLocator)
#         ax1.xaxis.set_major_formatter(majorFormatter)
#         ax1.xaxis.set_minor_locator(minorLocator)

        # y-axis
        majorLocator   = MultipleLocator(0.2)
        majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
        minorLocator   = MultipleLocator(0.1)
        ax1.yaxis.set_major_locator(majorLocator)
        ax1.yaxis.set_major_formatter(majorFormatter)
        ax1.yaxis.set_minor_locator(minorLocator)
        
        ax1.set_ylabel(r'$\rm Detection~Fraction$')

        leg = ax1.legend(scatterpoints=1,prop={'size':12},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax1.grid(b=None,which='major',axis='both')
        ylim(0., 1.)
        xlim(0.0001, 1.)

        if plot_detection_fraction_likelihood_save:
            savefig('{0}/detection_fraction_likelihood_Lstarmin{1}.pdf'.format(saveDirectory, Lstar_min),format='pdf',bbox_inches='tight')
        else:
            show()

##########################################################################################
##########################################################################################



##########################################################################################
##########################################################################################
    
    if plot_detection_fraction_both:
        fig = figure(figsize=(7.7,5.7))
        ax1 = fig.add_subplot(111)
        
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
        

        alpha_likelihood = 0.8
        alpha_impact = 0.8
        markerSize = 15
        lw = 2.

        
        label_likelihood = r'$\rm \mathcal{L}-Detection~Fraction$'
        label_impact = r'$\rm \rho - Detection~Fraction$'

        symbol_likelihood = 'D'
        symbol_impact = 'o'
        
        color_likelihood = color_blue
        color_impact = color_coal

        maxEW = 15000.

##########################################################################################
        # do the plotting 
        
        frac_imp500 = float(dv400_imp500_det) / float(dv400_imp500_det + dv400_imp500_non)
        frac_imp400 = float(dv400_imp400_det) / float(dv400_imp400_det + dv400_imp400_non)
        frac_imp300 = float(dv400_imp300_det) / float(dv400_imp300_det + dv400_imp300_non)
        frac_imp200 = float(dv400_imp200_det) / float(dv400_imp200_det + dv400_imp200_non)
        frac_imp100 = float(dv400_imp100_det) / float(dv400_imp100_det + dv400_imp100_non)
        frac_imp50  = float(dv400_imp50_det)  / float(dv400_imp50_det  + dv400_imp50_non)

        impact_x = [50, 100, 200, 300, 400, 500]
        impact_y = [frac_imp50, frac_imp100, frac_imp200, frac_imp300, frac_imp400, frac_imp500]
        
        
        frac_l001 = float(dv400_l001_det) / float(dv400_l001_det + dv400_l001_non)
        frac_l005 = float(dv400_l005_det) / float(dv400_l005_det + dv400_l005_non)
        frac_l01  = float(dv400_l01_det)  / float(dv400_l01_det  + dv400_l01_non)
        frac_l05  = float(dv400_l05_det)  / float(dv400_l05_det  + dv400_l05_non)
        frac_l1   = float(dv400_l1_det)   / float(dv400_l1_det   + dv400_l1_non)
        frac_l5   = float(dv400_l5_det)   / float(dv400_l5_det   + dv400_l5_non)

        likelihood_x = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5]
        likelihood_y = [frac_l001, frac_l005, frac_l01, frac_l05, frac_l1, frac_l5]


        # impact detection fraction
        ax1.plot(impact_x,
                impact_y,
                marker=symbol_impact,
                c=color_impact,
                ms=markerSize,
                markeredgecolor='black',
                lw = lw,
                alpha=alpha_impact,
                label=label_impact)


        # x-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)
        ax1.xaxis.set_major_locator(majorLocator)
        ax1.xaxis.set_major_formatter(majorFormatter)
        ax1.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
#         majorLocator   = MultipleLocator(0.2)
#         majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
#         minorLocator   = MultipleLocator(0.1)
#         ax1.yaxis.set_major_locator(majorLocator)
#         ax1.yaxis.set_major_formatter(majorFormatter)
#         ax1.yaxis.set_minor_locator(minorLocator)

        ax1.set_xlabel(r'$\rm \rho ~[kpc]$')

        # share a y-axis, have different top and bottom x-scales
        ax2 = ax1.twiny()

        # likelihood detection fraction
        ax2.plot(likelihood_x,
                likelihood_y,
                marker=symbol_likelihood,
                c=color_likelihood,
                ms=markerSize,
                markeredgecolor='black',
                lw = lw,
                alpha=alpha_likelihood,
                label=label_likelihood)

        ax2.set_xlabel(r'$\rm \mathcal{L}$')
    
        
        # x-axis
#         majorLocator   = MultipleLocator(0.01)
#         majorFormatter = FormatStrFormatter(r'$\rm %d$')
#         minorLocator   = MultipleLocator()
#         ax.xaxis.set_major_locator(majorLocator)
#         ax.xaxis.set_major_formatter(majorFormatter)
#         ax.xaxis.set_minor_locator(minorLocator)
        
        # y-axis
#         majorLocator   = MultipleLocator(10)
#         majorFormatter = FormatStrFormatter(r'$\rm %d$')
#         minorLocator   = MultipleLocator(5)
#         ax.yaxis.set_major_locator(majorLocator)
#         ax.yaxis.set_major_formatter(majorFormatter)
#         ax.yaxis.set_minor_locator(minorLocator)
        
        
        majorLocator   = MultipleLocator(0.2)
        majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
        minorLocator   = MultipleLocator(0.1)
        ax1.yaxis.set_major_locator(majorLocator)
        ax1.yaxis.set_major_formatter(majorFormatter)
        ax1.yaxis.set_minor_locator(minorLocator)
        
        ax1.set_ylabel(r'$\rm Detection~Fraction$')

        leg = ax1.legend(scatterpoints=1,prop={'size':12},loc=1,fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax1.grid(b=None,which='major',axis='both')
        ylim(0.0001, 1.)
#         xlim(0, 2.5)

        if plot_detection_fraction_both_save:
            savefig('{0}/detection_fraction.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()



#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


if __name__=="__main__":
    # do the work
    main()
    