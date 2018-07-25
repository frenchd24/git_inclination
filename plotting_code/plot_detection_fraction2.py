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

fontScale = 14
rc('text', usetex=True)
rc('font', size=14, family='serif', weight='normal')
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
    
    
def bmean(a):
    # compute the mean of 'a' but just return -99.99 if the array is empty

    try:
        return np.mean(a)
    except Exception, e:
        print 'Error: ',e
        return -99.99
        
        
def bmedian(a):
    # compute the median of 'a' but just return -99.99 if the array is empty

    try:
        return np.median(a)
    except Exception, e:
        print 'Error: ',e
        return -99.99
        
    
    
def main():
    # plot detection fraction as a function of both impact parameter and likelihood - 
    # does not work right now
    plot_detection_fraction_both = False
    plot_detection_fraction_both_save = False
    
    # plot detection fraction as a function of impact parameter
    plot_detection_fraction_impact = False
    plot_detection_fraction_impact_save = False
    
    # plot detection fraction as a function of likelihood
    plot_detection_fraction_likelihood = False
    plot_detection_fraction_likelihood_save = False
    
    # plot detection fraction as a function of likelihood
    plot_detection_fraction_likelihood_inc = False
    plot_detection_fraction_likelihood_inc_save = False

    # plot detection fraction as a function of likelihood
    plot_detection_fraction_impact_inc = False
    plot_detection_fraction_impact_inc_save = False

    
    # plot_number = 1 for just the isolated sample, =2 adds the associated, =3 adds two+
    # =4 adds groups with 2 or more members
    Lstar_min = 0.5

    # some colors
    color_blue = '#436bad'   # french blue
    color_red = '#ec2d01'    # tomato red

    if getpass.getuser() == 'frenchd':

#         gtPickleFilename = '/Users/frenchd/Research/GT_update2/pickleGT_filteredAll.p'

        saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/plotting_code/figs/'
        
#         detection_fraction_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction3.p'
        detection_fraction_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_lstarcut_all.p'

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
    dv400_imp25_det = full_dict['dv400_imp25_det']

    
    # now impact parameter non-detection counts
    dv400_imp1000_non = full_dict['dv400_imp1000_non']
    dv400_imp750_non = full_dict['dv400_imp750_non']
    dv400_imp500_non = full_dict['dv400_imp500_non']
    dv400_imp400_non = full_dict['dv400_imp400_non']
    dv400_imp300_non = full_dict['dv400_imp300_non']
    dv400_imp200_non = full_dict['dv400_imp200_non']
    dv400_imp100_non = full_dict['dv400_imp100_non']
    dv400_imp50_non = full_dict['dv400_imp50_non']
    dv400_imp25_non = full_dict['dv400_imp25_non']
    
    
    # impact parameter detection inclinations
    dv400_imp1000_det_inc = np.array(full_dict['dv400_imp1000_det_inc'])
    dv400_imp750_det_inc = np.array(full_dict['dv400_imp750_det_inc'])
    dv400_imp500_det_inc = np.array(full_dict['dv400_imp500_det_inc'])
    dv400_imp400_det_inc = np.array(full_dict['dv400_imp400_det_inc'])
    dv400_imp300_det_inc = np.array(full_dict['dv400_imp300_det_inc'])
    dv400_imp200_det_inc = np.array(full_dict['dv400_imp200_det_inc'])
    dv400_imp100_det_inc = np.array(full_dict['dv400_imp100_det_inc'])
    dv400_imp50_det_inc = np.array(full_dict['dv400_imp50_det_inc'])
    dv400_imp25_det_inc = np.array(full_dict['dv400_imp25_det_inc'])

    
    # now impact parameter non-detection counts
    dv400_imp1000_non_inc = np.array(full_dict['dv400_imp1000_non_inc'])
    dv400_imp750_non_inc = np.array(full_dict['dv400_imp750_non_inc'])
    dv400_imp500_non_inc = np.array(full_dict['dv400_imp500_non_inc'])
    dv400_imp400_non_inc = np.array(full_dict['dv400_imp400_non_inc'])
    dv400_imp300_non_inc = np.array(full_dict['dv400_imp300_non_inc'])
    dv400_imp200_non_inc = np.array(full_dict['dv400_imp200_non_inc'])
    dv400_imp100_non_inc = np.array(full_dict['dv400_imp100_non_inc'])
    dv400_imp50_non_inc = np.array(full_dict['dv400_imp50_non_inc'])
    dv400_imp25_non_inc = np.array(full_dict['dv400_imp25_non_inc'])




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

    
    # now for likelihood detections
    dv400_l001_det_inc = np.array(full_dict['dv400_l001_det_inc'])
    dv400_l005_det_inc = np.array(full_dict['dv400_l005_det_inc'])
    dv400_l01_det_inc = np.array(full_dict['dv400_l01_det_inc'])
    dv400_l05_det_inc = np.array(full_dict['dv400_l05_det_inc'])
    dv400_l1_det_inc = np.array(full_dict['dv400_l1_det_inc'])
    dv400_l5_det_inc = np.array(full_dict['dv400_l5_det_inc'])
    dv400_l75_det_inc = np.array(full_dict['dv400_l75_det_inc'])

    
    # now for likelihood non-detections
    dv400_l001_non_inc = np.array(full_dict['dv400_l001_non_inc'])
    dv400_l005_non_inc = np.array(full_dict['dv400_l005_non_inc'])
    dv400_l01_non_inc = np.array(full_dict['dv400_l01_non_inc'])
    dv400_l05_non_inc = np.array(full_dict['dv400_l05_non_inc'])
    dv400_l1_non_inc = np.array(full_dict['dv400_l1_non_inc'])
    dv400_l5_non_inc = np.array(full_dict['dv400_l5_non_inc'])
    dv400_l75_non_inc = np.array(full_dict['dv400_l75_non_inc'])
    
    
    
    
    print 'dv400_l001_non_inc: ',dv400_l001_non_inc
    print 'dv400_l001_det_inc: ',dv400_l001_det_inc
    
    dv400_imp1000_det_inc = dv400_imp1000_det_inc.astype(np.float)
    dv400_imp750_det_inc = dv400_imp750_det_inc.astype(np.float)
    dv400_imp500_det_inc = dv400_imp500_det_inc.astype(np.float)
    dv400_imp400_det_inc = dv400_imp400_det_inc.astype(np.float)
    dv400_imp300_det_inc = dv400_imp300_det_inc.astype(np.float)
    dv400_imp200_det_inc = dv400_imp200_det_inc.astype(np.float)
    dv400_imp100_det_inc = dv400_imp100_det_inc.astype(np.float)
    dv400_imp50_det_inc = dv400_imp50_det_inc.astype(np.float)
    dv400_imp25_det_inc = dv400_imp25_det_inc.astype(np.float)

    dv400_imp1000_non_inc = dv400_imp1000_non_inc.astype(np.float)
    dv400_imp750_non_inc = dv400_imp750_non_inc.astype(np.float)
    dv400_imp500_non_inc = dv400_imp500_non_inc.astype(np.float)
    dv400_imp400_non_inc = dv400_imp400_non_inc.astype(np.float)
    dv400_imp300_non_inc = dv400_imp300_non_inc.astype(np.float)
    dv400_imp200_non_inc = dv400_imp200_non_inc.astype(np.float)
    dv400_imp100_non_inc = dv400_imp100_non_inc.astype(np.float)
    dv400_imp50_non_inc = dv400_imp50_non_inc.astype(np.float)
    dv400_imp25_non_inc = dv400_imp25_non_inc.astype(np.float)
    
    dv400_l001_det_inc = dv400_l001_det_inc.astype(np.float)
    dv400_l005_det_inc = dv400_l005_det_inc.astype(np.float)
    dv400_l01_det_inc = dv400_l01_det_inc.astype(np.float)
    dv400_l05_det_inc = dv400_l05_det_inc.astype(np.float)
    dv400_l1_det_inc = dv400_l1_det_inc.astype(np.float)
    dv400_l5_det_inc = dv400_l5_det_inc.astype(np.float)
    dv400_l75_det_inc = dv400_l75_det_inc.astype(np.float)
    
    dv400_l001_non_inc = dv400_l001_non_inc.astype(np.float)
    dv400_l005_non_inc = dv400_l005_non_inc.astype(np.float)
    dv400_l01_non_inc = dv400_l01_non_inc.astype(np.float)
    dv400_l05_non_inc = dv400_l05_non_inc.astype(np.float)
    dv400_l1_non_inc = dv400_l1_non_inc.astype(np.float)
    dv400_l5_non_inc = dv400_l5_non_inc.astype(np.float)
    dv400_l75_non_inc = dv400_l75_non_inc.astype(np.float)
    
    # remove '-99' values from the likelihood inclination lists
    dv400_l0001_det_inc = list(filter(lambda x: x!= -99, dv400_l0001_det_inc))
    dv400_l0005_det_inc = list(filter(lambda x: x!= -99, dv400_l0005_det_inc))
    dv400_l001_det_inc = list(filter(lambda x: x!= -99, dv400_l001_det_inc))
    dv400_l005_det_inc = list(filter(lambda x: x!= -99, dv400_l005_det_inc))
    dv400_l01_det_inc = list(filter(lambda x: x!= -99, dv400_l01_det_inc))
    dv400_l05_det_inc = list(filter(lambda x: x!= -99, dv400_l05_det_inc))
    dv400_l1_det_inc = list(filter(lambda x: x!= -99, dv400_l1_det_inc))
    dv400_l5_det_inc = list(filter(lambda x: x!= -99, dv400_l5_det_inc))
    dv400_l75_det_inc = list(filter(lambda x: x!= -99, dv400_l75_det_inc))

    dv400_l0001_non_inc = list(filter(lambda x: x!= -99, dv400_l0001_non_inc))
    dv400_l0005_non_inc = list(filter(lambda x: x!= -99, dv400_l0005_non_inc))
    dv400_l001_non_inc = list(filter(lambda x: x!= -99, dv400_l001_non_inc))
    dv400_l005_non_inc = list(filter(lambda x: x!= -99, dv400_l005_non_inc))
    dv400_l01_non_inc = list(filter(lambda x: x!= -99, dv400_l01_non_inc))
    dv400_l05_non_inc = list(filter(lambda x: x!= -99, dv400_l05_non_inc))
    dv400_l1_non_inc = list(filter(lambda x: x!= -99, dv400_l1_non_inc))
    dv400_l5_non_inc = list(filter(lambda x: x!= -99, dv400_l5_non_inc))
    dv400_l75_non_inc = list(filter(lambda x: x!= -99, dv400_l75_non_inc))
    
    # remove '-99' values from the impact parameter inclination lists
    dv400_imp1000_det_inc = list(filter(lambda x: x!= -99, dv400_imp1000_det_inc))
    dv400_imp750_det_inc = list(filter(lambda x: x!= -99, dv400_imp750_det_inc))
    dv400_imp500_det_inc = list(filter(lambda x: x!= -99, dv400_imp500_det_inc))
    dv400_imp400_det_inc = list(filter(lambda x: x!= -99, dv400_imp400_det_inc))
    dv400_imp300_det_inc = list(filter(lambda x: x!= -99, dv400_imp300_det_inc))
    dv400_imp200_det_inc = list(filter(lambda x: x!= -99, dv400_imp200_det_inc))
    dv400_imp100_det_inc = list(filter(lambda x: x!= -99, dv400_imp100_det_inc))
    dv400_imp50_det_inc = list(filter(lambda x: x!= -99, dv400_imp50_det_inc))
    dv400_imp25_det_inc = list(filter(lambda x: x!= -99, dv400_imp25_det_inc))

    dv400_imp1000_non_inc = list(filter(lambda x: x!= -99, dv400_imp1000_non_inc))
    dv400_imp750_non_inc = list(filter(lambda x: x!= -99, dv400_imp750_non_inc))
    dv400_imp500_non_inc = list(filter(lambda x: x!= -99, dv400_imp500_non_inc))
    dv400_imp400_non_inc = list(filter(lambda x: x!= -99, dv400_imp400_non_inc))
    dv400_imp300_non_inc = list(filter(lambda x: x!= -99, dv400_imp300_non_inc))
    dv400_imp200_non_inc = list(filter(lambda x: x!= -99, dv400_imp200_non_inc))
    dv400_imp100_non_inc = list(filter(lambda x: x!= -99, dv400_imp100_non_inc))
    dv400_imp50_non_inc = list(filter(lambda x: x!= -99, dv400_imp50_non_inc))
    dv400_imp25_non_inc = list(filter(lambda x: x!= -99, dv400_imp25_non_inc))
    
    
    
    
    print 'dv400_imp1000_det_inc: ',dv400_imp1000_det_inc

    print
    print 'Detection fraction for 1000 kpc: ', float(dv400_imp1000_det) / (dv400_imp1000_det + dv400_imp1000_non)
    print 'Detection fraction for 750 kpc: ', float(dv400_imp750_det) / (dv400_imp750_det + dv400_imp750_non)
    print 'Detection fraction for 500 kpc: ', float(dv400_imp500_det) / (dv400_imp500_det + dv400_imp500_non)
    print 'Detection fraction for 400 kpc: ', float(dv400_imp400_det) / (dv400_imp400_det + dv400_imp400_non)
    print 'Detection fraction for 300 kpc: ', float(dv400_imp300_det) / (dv400_imp300_det + dv400_imp300_non)
    print 'Detection fraction for 200 kpc: ', float(dv400_imp200_det) / (dv400_imp200_det + dv400_imp200_non)
    print 'Detection fraction for 100 kpc: ', float(dv400_imp100_det) / (dv400_imp100_det + dv400_imp100_non)
    print 'Detection fraction for 50 kpc: ', float(dv400_imp50_det) / (dv400_imp50_det + dv400_imp50_non)
    print 'Detection fraction for 25 kpc: ', float(dv400_imp25_det) / (dv400_imp25_det + dv400_imp25_non)
    print
    print 'dv400_imp50_det: ',dv400_imp50_det
    print
    print 'dv400_imp50_non: ',dv400_imp50_non
    print
    print 'dv400_imp25_det: ',dv400_imp25_det
    print
    print 'dv400_imp25_non: ',dv400_imp25_non
    print
    print
    print
    
    try:
        print 'Detection fraction for 25 kpc: ', float(dv400_imp25_det) / (dv400_imp25_det + dv400_imp25_non)
    except Exception, e:
        print 'Error: ',e
    print
    print
    print 'Detection fraction for 1000 kpc mean inc: ', bmean(np.array(dv400_imp1000_det_inc)), ', ',bmean(np.array(dv400_imp1000_non_inc))
    print 'Detection fraction for 750 kpc mean inc: ', bmean(np.array(dv400_imp750_det_inc)), ', ',bmean(np.array(dv400_imp750_non_inc))
    print 'Detection fraction for 500 kpc mean inc: ', bmean(np.array(dv400_imp500_det_inc)), ', ',bmean(np.array(dv400_imp500_non_inc))
    print 'Detection fraction for 400 kpc mean inc: ', bmean(np.array(dv400_imp400_det_inc)), ', ',bmean(np.array(dv400_imp400_non_inc))
    print 'Detection fraction for 300 kpc mean inc: ', bmean(np.array(dv400_imp300_det_inc)), ', ',bmean(np.array(dv400_imp300_non_inc))
    print 'Detection fraction for 200 kpc mean inc: ', bmean(np.array(dv400_imp200_det_inc)), ', ',bmean(np.array(dv400_imp200_non_inc))
    print 'Detection fraction for 100 kpc mean inc: ', bmean(np.array(dv400_imp100_det_inc)), ', ',bmean(np.array(dv400_imp100_non_inc))
    print 'Detection fraction for 50 kpc mean inc: ', bmean(np.array(dv400_imp50_det_inc)), ', ',bmean(np.array(dv400_imp50_non_inc))
    print 'Detection fraction for 25 kpc mean inc: ', bmean(np.array(dv400_imp25_det_inc)), ', ',bmean(np.array(dv400_imp25_non_inc))
    print
    print
    print
    print
    print 'Detection fraction for 0.001 L mean inc: ', bmean(np.array(dv400_l001_det_inc)), ', ',bmean(np.array(dv400_l001_non_inc))
    print 'Detection fraction for 0.00 L mean inc: ', bmean(np.array(dv400_l005_det_inc)), ', ',bmean(np.array(dv400_l005_non_inc))
    print 'Detection fraction for 0.01 L mean inc: ', bmean(np.array(dv400_l01_det_inc)), ', ',bmean(np.array(dv400_l01_non_inc))
    print 'Detection fraction for 0.05 L mean inc: ', bmean(np.array(dv400_l05_det_inc)), ', ',bmean(np.array(dv400_l05_non_inc))
    print 'Detection fraction for 0.1 L mean inc: ', bmean(np.array(dv400_l1_det_inc)), ', ',bmean(np.array(dv400_l1_non_inc))
    print 'Detection fraction for 0.5 L mean inc: ', bmean(np.array(dv400_l5_det_inc)), ', ',bmean(np.array(dv400_l5_non_inc))
    print 'Detection fraction for 0.75 L mean inc: ', bmean(np.array(dv400_l75_det_inc)), ', ',bmean(np.array(dv400_l75_non_inc))



#########################################################################################
#########################################################################################
    # bootstrap this shit
    
#     reps = 1000
#     xb = np.random.choice(x, (n, reps), replace=True)
#     yb = 1/np.arange(1, n+1)[:, None] * np.cumsum(xb, axis=0)
#     upper, lower = np.percentile(yb, [2.5, 97.5], axis=1)


    print 'dv400_l001_det_inc: ',dv400_l001_det_inc
    print
    print 'dv400_l005_det_inc: ',dv400_l005_det_inc
    print
    print 'dv400_l01_det_inc: ',dv400_l01_det_inc
    print
    print 'dv400_l05_det_inc: ',dv400_l05_det_inc
    print
    print 'dv400_l1_det_inc: ',dv400_l1_det_inc
    print
    print 'dv400_l5_det_inc: ',dv400_l5_det_inc
    print
    print 'dv400_l75_det_inc: ',dv400_l75_det_inc
    print
    print 'dv400_l75_non_inc: ',dv400_l75_non_inc
    print



    x001 = np.array(dv400_l001_det_inc)
    x005 = np.array(dv400_l005_det_inc)
    x01 = np.array(dv400_l01_det_inc)
    x05 = np.array(dv400_l05_det_inc)
    x1 = np.array(dv400_l1_det_inc)
    x5 = np.array(dv400_l5_det_inc)
    x75 = np.array(dv400_l75_det_inc)

    n001 = len(x001)
    n005 = len(x005)
    n01 = len(x01)
    n05 = len(x05)
    n1 = len(x1)
    n5 = len(x5)
    n75 = len(x75)

    reps = 1000
    xb001 = np.random.choice(x001, (reps, n001), replace=True)
    xb005 = np.random.choice(x005, (reps, n005), replace=True)
    xb01 = np.random.choice(x01, (reps, n01), replace=True)
    xb05 = np.random.choice(x05, (reps, n05), replace=True)
    xb1 = np.random.choice(x1, (reps, n1), replace=True)
    xb5 = np.random.choice(x5, (reps, n5), replace=True)
    xb75 = np.random.choice(x75, (reps, n75), replace=True)

    meds001 = []
    for i in xb001:
        med = np.mean(i)
        meds001.append(med)

    meds005 = []
    for i in xb005:
        med = np.mean(i)
        meds005.append(med)
        
    meds01 = []
    for i in xb01:
        med = np.mean(i)
        meds01.append(med)
        
    meds05 = []
    for i in xb05:
        med = np.mean(i)
        meds05.append(med)
        
    meds1 = []
    for i in xb1:
        med = np.mean(i)
        meds1.append(med)
        
    meds5 = []
    for i in xb001:
        med = np.mean(i)
        meds5.append(med)
        
    meds75 = []
    for i in xb75:
        med = np.mean(i)
        meds75.append(med)
    
    print
    print 'dv400_l001_det_inc - median(meds, std(meds): ',mean(meds001), ', ',std(meds001)
    print 'dv400_l005_det_inc - median(meds, std(meds): ',mean(meds005), ', ',std(meds005)
    print 'dv400_l01_det_inc - median(meds, std(meds): ',mean(meds01), ', ',std(meds01)
    print 'dv400_l05_det_inc - median(meds, std(meds): ',mean(meds05), ', ',std(meds05)
    print 'dv400_l1_det_inc - median(meds, std(meds): ',mean(meds1), ', ',std(meds1)
    print 'dv400_l5_det_inc - median(meds, std(meds): ',mean(meds5), ', ',std(meds5)
    print 'dv400_l75_det_inc - median(meds, std(meds): ',mean(meds75), ', ',std(meds75)
    print
    print

    x001 = np.array(dv400_l001_non_inc)
    x005 = np.array(dv400_l005_non_inc)
    x01 = np.array(dv400_l01_non_inc)
    x05 = np.array(dv400_l05_non_inc)
    x1 = np.array(dv400_l1_non_inc)
    x5 = np.array(dv400_l5_non_inc)
    x75 = np.array(dv400_l75_non_inc)

    n001 = len(x001)
    n005 = len(x005)
    n01 = len(x01)
    n05 = len(x05)
    n1 = len(x1)
    n5 = len(x5)
    n75 = len(x75)

    reps = 1000
    xb001 = np.random.choice(x001, (reps, n001), replace=True)
    xb005 = np.random.choice(x005, (reps, n005), replace=True)
    xb01 = np.random.choice(x01, (reps, n01), replace=True)
    xb05 = np.random.choice(x05, (reps, n05), replace=True)
    xb1 = np.random.choice(x1, (reps, n1), replace=True)
    xb5 = np.random.choice(x5, (reps, n5), replace=True)
    xb75 = np.random.choice(x75, (reps, n75), replace=True)


    meds001 = []
    for i in xb001:
        med = np.mean(i)
        meds001.append(med)

    meds005 = []
    for i in xb005:
        med = np.mean(i)
        meds005.append(med)
        
    meds01 = []
    for i in xb01:
        med = np.mean(i)
        meds01.append(med)
        
    meds05 = []
    for i in xb05:
        med = np.mean(i)
        meds05.append(med)
        
    meds1 = []
    for i in xb1:
        med = np.mean(i)
        meds1.append(med)
        
    meds5 = []
    for i in xb001:
        med = np.mean(i)
        meds5.append(med)
        
    meds75 = []
    for i in xb75:
        med = np.mean(i)
        meds75.append(med)
    
    print
    print 'dv400_l001_non_inc - median(meds, std(meds): ',mean(meds001), ', ',std(meds001)
    print 'dv400_l005_non_inc - median(meds, std(meds): ',mean(meds005), ', ',std(meds005)
    print 'dv400_l01_non_inc - median(meds, std(meds): ',mean(meds01), ', ',std(meds01)
    print 'dv400_l05_non_inc - median(meds, std(meds): ',mean(meds05), ', ',std(meds05)
    print 'dv400_l1_non_inc - median(meds, std(meds): ',mean(meds1), ', ',std(meds1)
    print 'dv400_l5_non_inc - median(meds, std(meds): ',mean(meds5), ', ',std(meds5)
    print 'dv400_l75_non_inc - median(meds, std(meds): ',mean(meds75), ', ',std(meds75)
    print
    print






##########################################################################################
##########################################################################################
    
    if plot_detection_fraction_likelihood_inc:
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
        

        alpha_det = 0.8
        alpha_non = 0.8
        markerSize = 10
        lw = 2.
        
        binSize = 100
        bins = arange(0, 600, binSize)

        
        label_det = r'$\rm Detections$'
        label_non = r'$\rm Non-Detections$'

        symbol_det = 'D'
        symbol_non = 'X'
        
        color_det = color_blue
        color_non = color_red

        maxEW = 15000.
        use_mean = False

##########################################################################################
        # do the plotting
        
        print 'len(dv400_l5_det_inc) : ',len(dv400_l1_det_inc)
        print 'len(dv400_l5_det_inc) : ',len(dv400_l5_det_inc)
        print 'len(dv400_l75_det_inc) : ',len(dv400_l75_det_inc)
        print
        print
        print 'len(dv400_l5_non_inc) : ',len(dv400_l1_non_inc)
        print 'len(dv400_l5_non_inc) : ',len(dv400_l5_non_inc)
        print 'len(dv400_l75_non_inc) : ',len(dv400_l75_non_inc)
        print
            
            

        x = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 0.75]
        if use_mean:
            y_det = [bmean(dv400_l001_det_inc),
                    bmean(dv400_l005_det_inc), 
                    bmean(dv400_l01_det_inc),
                    bmean(dv400_l05_det_inc),
                    bmean(dv400_l1_det_inc),
                    bmean(dv400_l5_det_inc),
                    bmean(dv400_l75_det_inc)]
                
                
            y_non = [bmean(dv400_l001_non_inc),
                    bmean(dv400_l005_non_inc), 
                    bmean(dv400_l01_non_inc),
                    bmean(dv400_l05_non_inc),
                    bmean(dv400_l1_non_inc),
                    bmean(dv400_l5_non_inc),
                    bmean(dv400_l75_non_inc)]
                    
        else:
            y_det = [bmedian(dv400_l001_det_inc),
                    bmedian(dv400_l005_det_inc), 
                    bmedian(dv400_l01_det_inc),
                    bmedian(dv400_l05_det_inc),
                    bmedian(dv400_l1_det_inc),
                    bmedian(dv400_l5_det_inc),
                    bmedian(dv400_l75_det_inc)]
                
                
            y_non = [bmedian(dv400_l001_non_inc),
                    bmedian(dv400_l005_non_inc), 
                    bmedian(dv400_l01_non_inc),
                    bmedian(dv400_l05_non_inc),
                    bmedian(dv400_l1_non_inc),
                    bmedian(dv400_l5_non_inc),
                    bmedian(dv400_l75_non_inc)]



        # likelihood detection fraction mean inc
        ax1.plot(x,
                y_det,
                marker=symbol_det,
                c=color_det,
                ms=markerSize,
                markeredgecolor='black',
                lw = lw,
                alpha=alpha_det,
                label=label_det)

        # likelihood non-detection fraction mean inc
        ax1.plot(x,
                y_non,
                marker=symbol_non,
                c=color_non,
                ms=markerSize,
                markeredgecolor='black',
                lw = lw,
                alpha=alpha_non,
                label=label_non)

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
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax1.yaxis.set_major_locator(majorLocator)
        ax1.yaxis.set_major_formatter(majorFormatter)
        ax1.yaxis.set_minor_locator(minorLocator)
        
        ax1.set_ylabel(r'$\rm Inclination~[deg]$')
        
        leg = ax1.legend(scatterpoints=1,prop={'size':12},loc='lower left',fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax1.grid(b=None,which='major',axis='both')
        ylim(0., 90.)
        xlim(0.001, 1.)
        ax1.invert_xaxis()


        if plot_detection_fraction_likelihood_inc_save:
            savefig('{0}/detection_fraction_likelihood_inc_median_Lstarcut01-2.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()

##########################################################################################
##########################################################################################


##########################################################################################
##########################################################################################
    
    if plot_detection_fraction_impact_inc:
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
        

        alpha_det = 0.8
        alpha_non = 0.8
        markerSize = 10
        lw = 2.
        
        binSize = 100
        bins = arange(0, 600, binSize)

        
        label_det = r'$\rm Detections$'
        label_non = r'$\rm Non-Detections$'

        symbol_det = 'D'
        symbol_non = 'X'
        
        color_det = color_blue
        color_non = color_red

        maxEW = 15000.
        use_mean = False

##########################################################################################
        # do the plotting
        
        print 'len(dv400_imp100_det_inc) : ',len(dv400_imp100_det_inc)
        print 'len(dv400_imp50_det_inc) : ',len(dv400_imp50_det_inc)
        print 'len(dv400_imp25_det_inc) : ',len(dv400_imp25_det_inc)
        print
        print
        print 'len(dv400_imp100_non_inc) : ',len(dv400_imp100_non_inc)
        print 'len(dv400_imp50_non_inc) : ',len(dv400_imp50_non_inc)
        print 'len(dv400_imp25_non_inc) : ',len(dv400_imp25_non_inc)
        print

        dv400_imp1000_det_inc = list(filter(lambda x: x!= -99, dv400_imp1000_det_inc))
        dv400_imp750_det_inc = list(filter(lambda x: x!= -99, dv400_imp750_det_inc))
        dv400_imp500_det_inc = list(filter(lambda x: x!= -99, dv400_imp500_det_inc))
        dv400_imp400_det_inc = list(filter(lambda x: x!= -99, dv400_imp400_det_inc))
        dv400_imp300_det_inc = list(filter(lambda x: x!= -99, dv400_imp300_det_inc))
        dv400_imp200_det_inc = list(filter(lambda x: x!= -99, dv400_imp200_det_inc))
        dv400_imp100_det_inc = list(filter(lambda x: x!= -99, dv400_imp100_det_inc))
        dv400_imp50_det_inc = list(filter(lambda x: x!= -99, dv400_imp50_det_inc))
        dv400_imp25_det_inc = list(filter(lambda x: x!= -99, dv400_imp25_det_inc))


        dv400_imp1000_non_inc = list(filter(lambda x: x!= -99, dv400_imp1000_non_inc))
        dv400_imp750_non_inc = list(filter(lambda x: x!= -99, dv400_imp750_non_inc))
        dv400_imp500_non_inc = list(filter(lambda x: x!= -99, dv400_imp500_non_inc))
        dv400_imp400_non_inc = list(filter(lambda x: x!= -99, dv400_imp400_non_inc))
        dv400_imp300_non_inc = list(filter(lambda x: x!= -99, dv400_imp300_non_inc))
        dv400_imp200_non_inc = list(filter(lambda x: x!= -99, dv400_imp200_non_inc))
        dv400_imp100_non_inc = list(filter(lambda x: x!= -99, dv400_imp100_non_inc))
        dv400_imp50_non_inc = list(filter(lambda x: x!= -99, dv400_imp50_non_inc))
        dv400_imp25_non_inc = list(filter(lambda x: x!= -99, dv400_imp25_non_inc))


        x = [1000, 750, 500, 400, 300, 200, 100, 50, 25]
        if use_mean:
            y_det = [bmean(dv400_imp1000_det_inc),
                    bmean(dv400_imp750_det_inc),
                    bmean(dv400_imp500_det_inc),
                    bmean(dv400_imp400_det_inc), 
                    bmean(dv400_imp300_det_inc),
                    bmean(dv400_imp200_det_inc),
                    bmean(dv400_imp100_det_inc),
                    bmean(dv400_imp50_det_inc),
                    bmean(dv400_imp25_det_inc)]
                
                
            y_non = [bmean(dv400_imp1000_non_inc),
                    bmean(dv400_imp750_non_inc),
                    bmean(dv400_imp500_non_inc),
                    bmean(dv400_imp400_non_inc), 
                    bmean(dv400_imp300_non_inc),
                    bmean(dv400_imp200_non_inc),
                    bmean(dv400_imp100_non_inc),
                    bmean(dv400_imp50_non_inc),
                    bmean(dv400_imp25_non_inc)]
                    
        else:
            y_det = [bmedian(dv400_imp1000_det_inc),
                    bmedian(dv400_imp750_det_inc),
                    bmedian(dv400_imp500_det_inc),
                    bmedian(dv400_imp400_det_inc), 
                    bmedian(dv400_imp300_det_inc),
                    bmedian(dv400_imp200_det_inc),
                    bmedian(dv400_imp100_det_inc),
                    bmedian(dv400_imp50_det_inc),
                    bmedian(dv400_imp25_det_inc)]
                
                
            y_non = [bmedian(dv400_imp1000_non_inc),
                    bmedian(dv400_imp750_non_inc),
                    bmedian(dv400_imp500_non_inc),
                    bmedian(dv400_imp400_non_inc), 
                    bmedian(dv400_imp300_non_inc),
                    bmedian(dv400_imp200_non_inc),
                    bmedian(dv400_imp100_non_inc),
                    bmedian(dv400_imp50_non_inc),
                    bmedian(dv400_imp25_non_inc)]



        # likelihood detection fraction mean inc
        ax1.plot(x,
                y_det,
                marker=symbol_det,
                c=color_det,
                ms=markerSize,
                markeredgecolor='black',
                lw = lw,
                alpha=alpha_det,
                label=label_det)

        # likelihood non-detection fraction mean inc
        ax1.plot(x,
                y_non,
                marker=symbol_non,
                c=color_non,
                ms=markerSize,
                markeredgecolor='black',
                lw = lw,
                alpha=alpha_non,
                label=label_non)
                

        ax1.set_xlabel(r'$\rm \rho~[kpc]$')
        
        # x-axis
#         majorLocator   = MultipleLocator(0.01)
#         majorFormatter = FormatStrFormatter(r'$\rm %d$')
#         minorLocator   = MultipleLocator()
#         ax1.xaxis.set_major_locator(majorLocator)
#         ax1.xaxis.set_major_formatter(majorFormatter)
#         ax1.xaxis.set_minor_locator(minorLocator)

        # y-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax1.yaxis.set_major_locator(majorLocator)
        ax1.yaxis.set_major_formatter(majorFormatter)
        ax1.yaxis.set_minor_locator(minorLocator)
        
        ax1.set_ylabel(r'$\rm Inclination~[deg]$')

        leg = ax1.legend(scatterpoints=1,prop={'size':12},loc='lower left',fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax1.grid(b=None,which='major',axis='both')
        ylim(0., 90.)
        xlim(0, 1000)

        if plot_detection_fraction_impact_inc_save:
            savefig('{0}/detection_fraction_impact_inc_median_Lstarcut01-2.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
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
    