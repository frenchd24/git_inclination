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
        
        
        
def return_bootstrap_errors(x, reps):
    # takes in x as an array of values, calculates bootstrap mean and median errors based 
    # on 'reps' number of iterations
    
    x = np.array(x)
    lenx = len(x)

    x_resample = np.random.choice(x, (reps, lenx), replace=True)

    resampled_means = []
    for i in x_resample:
        m = np.mean(i)
        resampled_means.append(m)
    
    mean_err = std(resampled_means)


    resampled_medians = []
    for i in x_resample:
        med = np.median(i)
        resampled_medians.append(med)
    
    median_err = std(resampled_medians)
        

    return mean_err, median_err
    


def main():
    # plot detection fraction as a function of both impact parameter and likelihood - 
    # does not work right now
    plot_detection_fraction_both = False
    plot_detection_fraction_both_save = False
    
    # plot detection fraction as a function of both impact parameter and likelihood - 
    # does not work right now
    plot_detection_fraction_velcut_both = True
    plot_detection_fraction_velcut_both_save = True
    
    
    # plot_number = 1 for just the isolated sample, =2 adds the associated, =3 adds two+
    # =4 adds groups with 2 or more members
    Lstar_min = 0.5
    
    # which lstar cut subset to use?
    lstar_cut  = 'include3'

    # some colors
    color_blue = '#436bad'  # french blue
    color_red = '#ec2d01'   # tomato red

    if getpass.getuser() == 'frenchd':

#         gtPickleFilename = '/Users/frenchd/Research/GT_update2/pickleGT_filteredAll.p'

        saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/plotting_code/detection_fraction_figs/'
        
        detection_fraction_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_lstarcut_{0}.p'.format(lstar_cut)
#         detection_fraction_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_minEW0_closestonly.p'
#         detection_fraction_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_strictbins.p'

#         detection_fraction_vcut_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_vcut2500.p'
#         detection_fraction_vcut_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_lstarcut01_vcut2500.p'
#         detection_fraction_vcut_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_dmin75_vcut2500_minEW50.p'
        detection_fraction_min50_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_minEW50_2.p'
        detection_fraction_min100_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_minEW100_closestonly.p'
        detection_fraction_min200_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_minEW200_closestonly.p'
        detection_fraction_min300_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_minEW300_closestonly.p'

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # pickle file for the whole galaxy table:
#     gtPickleFile = open(gtPickleFilename,'rU')
#     gtDict = pickle.load(gtPickleFile)
#     gtPickleFile.close()
    
    
    # open all the pickle files
    detection_fraction_file = open(detection_fraction_filename,'r')
    detection_fraction_min50_file = open(detection_fraction_min50_filename,'r')
    detection_fraction_min100_file = open(detection_fraction_min100_filename,'r')
    detection_fraction_min200_file = open(detection_fraction_min200_filename,'r')
    detection_fraction_min300_file = open(detection_fraction_min300_filename,'r')

    # unload the data from them
    full_dict = pickle.load(detection_fraction_file)
    full_dict_min50 = pickle.load(detection_fraction_min50_file)
    full_dict_min100 = pickle.load(detection_fraction_min100_file)
    full_dict_min200 = pickle.load(detection_fraction_min200_file)
    full_dict_min300 = pickle.load(detection_fraction_min300_file)

    # close the files
    detection_fraction_file.close()
    detection_fraction_min50_file.close()
    detection_fraction_min100_file.close()
    detection_fraction_min200_file.close()
    detection_fraction_min300_file.close()

    
    # unload all the stuff
     
    ##############
    # minEW 50 results
    # now impact parameter detection counts
    dv400_imp1000_min50_det = full_dict_min50['dv400_imp1000_det']
    dv400_imp750_min50_det = full_dict_min50['dv400_imp750_det']
    dv400_imp500_min50_det = full_dict_min50['dv400_imp500_det']
    dv400_imp400_min50_det = full_dict_min50['dv400_imp400_det']
    dv400_imp300_min50_det = full_dict_min50['dv400_imp300_det']
    dv400_imp200_min50_det = full_dict_min50['dv400_imp200_det']
    dv400_imp100_min50_det = full_dict_min50['dv400_imp100_det']
    dv400_imp50_min50_det = full_dict_min50['dv400_imp50_det']
    dv400_imp25_min50_det = full_dict_min50['dv400_imp25_det']

    # now impact parameter non-detection counts
    dv400_imp1000_min50_non = full_dict_min50['dv400_imp1000_non']
    dv400_imp750_min50_non = full_dict_min50['dv400_imp750_non']
    dv400_imp500_min50_non = full_dict_min50['dv400_imp500_non']
    dv400_imp400_min50_non = full_dict_min50['dv400_imp400_non']
    dv400_imp300_min50_non = full_dict_min50['dv400_imp300_non']
    dv400_imp200_min50_non = full_dict_min50['dv400_imp200_non']
    dv400_imp100_min50_non = full_dict_min50['dv400_imp100_non']
    dv400_imp50_min50_non = full_dict_min50['dv400_imp50_non']
    dv400_imp25_min50_non = full_dict_min50['dv400_imp25_non']
    
    
    # now for likelihood detections
    dv400_l0001_min50_det = full_dict_min50['dv400_l0001_det']
    dv400_l0005_min50_det = full_dict_min50['dv400_l0005_det']
    dv400_l001_min50_det = full_dict_min50['dv400_l001_det']
    dv400_l005_min50_det = full_dict_min50['dv400_l005_det']
    dv400_l01_min50_det = full_dict_min50['dv400_l01_det']
    dv400_l05_min50_det = full_dict_min50['dv400_l05_det']
    dv400_l1_min50_det = full_dict_min50['dv400_l1_det']
    dv400_l5_min50_det = full_dict_min50['dv400_l5_det']
    dv400_l75_min50_det = full_dict_min50['dv400_l75_det']
    
    # now for likelihood non-detections
    dv400_l0001_min50_non = full_dict_min50['dv400_l0001_non']
    dv400_l0005_min50_non = full_dict_min50['dv400_l0005_non']
    dv400_l001_min50_non = full_dict_min50['dv400_l001_non']
    dv400_l005_min50_non = full_dict_min50['dv400_l005_non']
    dv400_l01_min50_non = full_dict_min50['dv400_l01_non']
    dv400_l05_min50_non = full_dict_min50['dv400_l05_non']
    dv400_l1_min50_non = full_dict_min50['dv400_l1_non']
    dv400_l5_min50_non = full_dict_min50['dv400_l5_non']
    dv400_l75_min50_non = full_dict_min50['dv400_l75_non']
    
    
    ##############
    # minEW 100 results
    # now impact parameter detection counts
    dv400_imp1000_min100_det = full_dict_min100['dv400_imp1000_det']
    dv400_imp750_min100_det = full_dict_min100['dv400_imp750_det']
    dv400_imp500_min100_det = full_dict_min100['dv400_imp500_det']
    dv400_imp400_min100_det = full_dict_min100['dv400_imp400_det']
    dv400_imp300_min100_det = full_dict_min100['dv400_imp300_det']
    dv400_imp200_min100_det = full_dict_min100['dv400_imp200_det']
    dv400_imp100_min100_det = full_dict_min100['dv400_imp100_det']
    dv400_imp50_min100_det = full_dict_min100['dv400_imp50_det']
    dv400_imp25_min100_det = full_dict_min100['dv400_imp25_det']

    # now impact parameter non-detection counts
    dv400_imp1000_min100_non = full_dict_min100['dv400_imp1000_non']
    dv400_imp750_min100_non = full_dict_min100['dv400_imp750_non']
    dv400_imp500_min100_non = full_dict_min100['dv400_imp500_non']
    dv400_imp400_min100_non = full_dict_min100['dv400_imp400_non']
    dv400_imp300_min100_non = full_dict_min100['dv400_imp300_non']
    dv400_imp200_min100_non = full_dict_min100['dv400_imp200_non']
    dv400_imp100_min100_non = full_dict_min100['dv400_imp100_non']
    dv400_imp50_min100_non = full_dict_min100['dv400_imp50_non']
    dv400_imp25_min100_non = full_dict_min100['dv400_imp25_non']
    
    
    # now for likelihood detections
    dv400_l0001_min100_det = full_dict_min100['dv400_l0001_det']
    dv400_l0005_min100_det = full_dict_min100['dv400_l0005_det']
    dv400_l001_min100_det = full_dict_min100['dv400_l001_det']
    dv400_l005_min100_det = full_dict_min100['dv400_l005_det']
    dv400_l01_min100_det = full_dict_min100['dv400_l01_det']
    dv400_l05_min100_det = full_dict_min100['dv400_l05_det']
    dv400_l1_min100_det = full_dict_min100['dv400_l1_det']
    dv400_l5_min100_det = full_dict_min100['dv400_l5_det']
    dv400_l75_min100_det = full_dict_min100['dv400_l75_det']
    
    # now for likelihood non-detections
    dv400_l0001_min100_non = full_dict_min100['dv400_l0001_non']
    dv400_l0005_min100_non = full_dict_min100['dv400_l0005_non']
    dv400_l001_min100_non = full_dict_min100['dv400_l001_non']
    dv400_l005_min100_non = full_dict_min100['dv400_l005_non']
    dv400_l01_min100_non = full_dict_min100['dv400_l01_non']
    dv400_l05_min100_non = full_dict_min100['dv400_l05_non']
    dv400_l1_min100_non = full_dict_min100['dv400_l1_non']
    dv400_l5_min100_non = full_dict_min100['dv400_l5_non']
    dv400_l75_min100_non = full_dict_min100['dv400_l75_non']
    
    
    
    ##############
    # minEW 200 results
    # now impact parameter detection counts
    dv400_imp1000_min200_det = full_dict_min200['dv400_imp1000_det']
    dv400_imp750_min200_det = full_dict_min200['dv400_imp750_det']
    dv400_imp500_min200_det = full_dict_min200['dv400_imp500_det']
    dv400_imp400_min200_det = full_dict_min200['dv400_imp400_det']
    dv400_imp300_min200_det = full_dict_min200['dv400_imp300_det']
    dv400_imp200_min200_det = full_dict_min200['dv400_imp200_det']
    dv400_imp100_min200_det = full_dict_min200['dv400_imp100_det']
    dv400_imp50_min200_det = full_dict_min200['dv400_imp50_det']
    dv400_imp25_min200_det = full_dict_min200['dv400_imp25_det']

    # now impact parameter non-detection counts
    dv400_imp1000_min200_non = full_dict_min200['dv400_imp1000_non']
    dv400_imp750_min200_non = full_dict_min200['dv400_imp750_non']
    dv400_imp500_min200_non = full_dict_min200['dv400_imp500_non']
    dv400_imp400_min200_non = full_dict_min200['dv400_imp400_non']
    dv400_imp300_min200_non = full_dict_min200['dv400_imp300_non']
    dv400_imp200_min200_non = full_dict_min200['dv400_imp200_non']
    dv400_imp100_min200_non = full_dict_min200['dv400_imp100_non']
    dv400_imp50_min200_non = full_dict_min200['dv400_imp50_non']
    dv400_imp25_min200_non = full_dict_min200['dv400_imp25_non']
    
    
    # now for likelihood detections
    dv400_l0001_min200_det = full_dict_min200['dv400_l0001_det']
    dv400_l0005_min200_det = full_dict_min200['dv400_l0005_det']
    dv400_l001_min200_det = full_dict_min200['dv400_l001_det']
    dv400_l005_min200_det = full_dict_min200['dv400_l005_det']
    dv400_l01_min200_det = full_dict_min200['dv400_l01_det']
    dv400_l05_min200_det = full_dict_min200['dv400_l05_det']
    dv400_l1_min200_det = full_dict_min200['dv400_l1_det']
    dv400_l5_min200_det = full_dict_min200['dv400_l5_det']
    dv400_l75_min200_det = full_dict_min200['dv400_l75_det']
    
    # now for likelihood non-detections
    dv400_l0001_min200_non = full_dict_min200['dv400_l0001_non']
    dv400_l0005_min200_non = full_dict_min200['dv400_l0005_non']
    dv400_l001_min200_non = full_dict_min200['dv400_l001_non']
    dv400_l005_min200_non = full_dict_min200['dv400_l005_non']
    dv400_l01_min200_non = full_dict_min200['dv400_l01_non']
    dv400_l05_min200_non = full_dict_min200['dv400_l05_non']
    dv400_l1_min200_non = full_dict_min200['dv400_l1_non']
    dv400_l5_min200_non = full_dict_min200['dv400_l5_non']
    dv400_l75_min200_non = full_dict_min200['dv400_l75_non']
    
    
    
    ##############
    # minEW 300 results
    # now impact parameter detection counts
    dv400_imp1000_min300_det = full_dict_min300['dv400_imp1000_det']
    dv400_imp750_min300_det = full_dict_min300['dv400_imp750_det']
    dv400_imp500_min300_det = full_dict_min300['dv400_imp500_det']
    dv400_imp400_min300_det = full_dict_min300['dv400_imp400_det']
    dv400_imp300_min300_det = full_dict_min300['dv400_imp300_det']
    dv400_imp200_min300_det = full_dict_min300['dv400_imp200_det']
    dv400_imp100_min300_det = full_dict_min300['dv400_imp100_det']
    dv400_imp50_min300_det = full_dict_min300['dv400_imp50_det']
    dv400_imp25_min300_det = full_dict_min300['dv400_imp25_det']

    # now impact parameter non-detection counts
    dv400_imp1000_min300_non = full_dict_min300['dv400_imp1000_non']
    dv400_imp750_min300_non = full_dict_min300['dv400_imp750_non']
    dv400_imp500_min300_non = full_dict_min300['dv400_imp500_non']
    dv400_imp400_min300_non = full_dict_min300['dv400_imp400_non']
    dv400_imp300_min300_non = full_dict_min300['dv400_imp300_non']
    dv400_imp200_min300_non = full_dict_min300['dv400_imp200_non']
    dv400_imp100_min300_non = full_dict_min300['dv400_imp100_non']
    dv400_imp50_min300_non = full_dict_min300['dv400_imp50_non']
    dv400_imp25_min300_non = full_dict_min300['dv400_imp25_non']
    
    
    # now for likelihood detections
    dv400_l0001_min300_det = full_dict_min300['dv400_l0001_det']
    dv400_l0005_min300_det = full_dict_min300['dv400_l0005_det']
    dv400_l001_min300_det = full_dict_min300['dv400_l001_det']
    dv400_l005_min300_det = full_dict_min300['dv400_l005_det']
    dv400_l01_min300_det = full_dict_min300['dv400_l01_det']
    dv400_l05_min300_det = full_dict_min300['dv400_l05_det']
    dv400_l1_min300_det = full_dict_min300['dv400_l1_det']
    dv400_l5_min300_det = full_dict_min300['dv400_l5_det']
    dv400_l75_min300_det = full_dict_min300['dv400_l75_det']
    
    # now for likelihood non-detections
    dv400_l0001_min300_non = full_dict_min300['dv400_l0001_non']
    dv400_l0005_min300_non = full_dict_min300['dv400_l0005_non']
    dv400_l001_min300_non = full_dict_min300['dv400_l001_non']
    dv400_l005_min300_non = full_dict_min300['dv400_l005_non']
    dv400_l01_min300_non = full_dict_min300['dv400_l01_non']
    dv400_l05_min300_non = full_dict_min300['dv400_l05_non']
    dv400_l1_min300_non = full_dict_min300['dv400_l1_non']
    dv400_l5_min300_non = full_dict_min300['dv400_l5_non']
    dv400_l75_min300_non = full_dict_min300['dv400_l75_non']
    
##########################################################################################
    
    
##########################################################################################
    
##########################################################################################
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


##########################################################################################
    # now for likelihood thresholds
    dv_l0001 = full_dict['dv_l0001']
    dv_l0005 = full_dict['dv_l0005']
    dv_l001 = full_dict['dv_l001']
    dv_l005 = full_dict['dv_l005']
    dv_l01 = full_dict['dv_l01']
    dv_l05 = full_dict['dv_l05']
    dv_l1 = full_dict['dv_l1']
    dv_l5 = full_dict['dv_l5']
    dv_l75 = full_dict['dv_l75']

    # now for likelihood detections
    dv400_l0001_det = full_dict['dv400_l0001_det']
    dv400_l0005_det = full_dict['dv400_l0005_det']
    dv400_l001_det = full_dict['dv400_l001_det']
    dv400_l005_det = full_dict['dv400_l005_det']
    dv400_l01_det = full_dict['dv400_l01_det']
    dv400_l05_det = full_dict['dv400_l05_det']
    dv400_l1_det = full_dict['dv400_l1_det']
    dv400_l5_det = full_dict['dv400_l5_det']
    dv400_l75_det = full_dict['dv400_l75_det']
    
    # now for likelihood non-detections
    dv400_l0001_non = full_dict['dv400_l0001_non']
    dv400_l0005_non = full_dict['dv400_l0005_non']
    dv400_l001_non = full_dict['dv400_l001_non']
    dv400_l005_non = full_dict['dv400_l005_non']
    dv400_l01_non = full_dict['dv400_l01_non']
    dv400_l05_non = full_dict['dv400_l05_non']
    dv400_l1_non = full_dict['dv400_l1_non']
    dv400_l5_non = full_dict['dv400_l5_non']
    dv400_l75_non = full_dict['dv400_l75_non']

    print 'Detection for 1000 kpc: ', float(dv400_imp1000_det)
    print 'Detection for 750 kpc: ', float(dv400_imp750_det)
    print 'Detection for 500 kpc: ', float(dv400_imp500_det)
    print 'Detection for 400 kpc: ', float(dv400_imp400_det)
    print 'Detection for 300 kpc: ', float(dv400_imp300_det)
    print 'Detection for 200 kpc: ', float(dv400_imp200_det)
    print 'Detection for 100 kpc: ', float(dv400_imp100_det)
    print 'Detection for 50 kpc: ', float(dv400_imp50_det)
    print 'Detection for 25 kpc: ', float(dv400_imp25_det)


    print 'Detection fraction for 1000 kpc: ', float(dv400_imp1000_det) / (dv400_imp1000_det + dv400_imp1000_non)
    print 'Detection fraction for 750 kpc: ', float(dv400_imp750_det) / (dv400_imp750_det + dv400_imp750_non)
    print 'Detection fraction for 500 kpc: ', float(dv400_imp500_det) / (dv400_imp500_det + dv400_imp500_non)
    print 'Detection fraction for 400 kpc: ', float(dv400_imp400_det) / (dv400_imp400_det + dv400_imp400_non)
#     print 'Detection fraction for 300 kpc: ', float(dv400_imp300_det) / (dv400_imp300_det + dv400_imp300_non)
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



#########################################################################################
#########################################################################################

#########################################################################################
#########################################################################################
    # prepare the data
    frac_imp1000 = float(dv400_imp1000_det) / float(dv400_imp1000_det + dv400_imp1000_non)
    frac_imp750 = float(dv400_imp750_det) / float(dv400_imp750_det + dv400_imp750_non)
    frac_imp500 = float(dv400_imp500_det) / float(dv400_imp500_det + dv400_imp500_non)
    frac_imp400 = float(dv400_imp400_det) / float(dv400_imp400_det + dv400_imp400_non)
    frac_imp300 = float(dv400_imp300_det) / float(dv400_imp300_det + dv400_imp300_non)
    frac_imp200 = float(dv400_imp200_det) / float(dv400_imp200_det + dv400_imp200_non)
    frac_imp100 = float(dv400_imp100_det) / float(dv400_imp100_det + dv400_imp100_non)
    frac_imp50  = float(dv400_imp50_det)  / float(dv400_imp50_det  + dv400_imp50_non)
    
    frac_imp1000_err = np.sqrt(dv400_imp1000_det) / float(dv400_imp1000_det + dv400_imp1000_non)
    frac_imp750_err = np.sqrt(dv400_imp750_det) / float(dv400_imp750_det + dv400_imp750_non)
    frac_imp500_err = np.sqrt(dv400_imp500_det) / float(dv400_imp500_det + dv400_imp500_non)
    frac_imp400_err = np.sqrt(dv400_imp400_det) / float(dv400_imp400_det + dv400_imp400_non)
    frac_imp300_err = np.sqrt(dv400_imp300_det) / float(dv400_imp300_det + dv400_imp300_non)
    frac_imp200_err = np.sqrt(dv400_imp200_det) / float(dv400_imp200_det + dv400_imp200_non)
    frac_imp100_err = np.sqrt(dv400_imp100_det) / float(dv400_imp100_det + dv400_imp100_non)
    frac_imp50_err  = np.sqrt(dv400_imp50_det)  / float(dv400_imp50_det  + dv400_imp50_non)
    
    
    try:
        frac_imp25  = float(dv400_imp25_det)  / float(dv400_imp25_det  + dv400_imp25_non)
        frac_imp25_err  = np.sqrt(dv400_imp25_det)  / float(dv400_imp25_det  + dv400_imp25_non)

    except Exception,e:
        print 'Error: ',e
        frac_imp25 = 0
        frac_imp25_err = 0

    impact_x = [25, 50, 100, 200, 300, 400, 500, 750, 1000]
    impact_y = [frac_imp25, frac_imp50, frac_imp100, frac_imp200, frac_imp300, frac_imp400, frac_imp500, frac_imp750, frac_imp1000]
    impact_y_err = [frac_imp25_err,
                    frac_imp50_err,
                    frac_imp100_err,
                    frac_imp200_err,
                    frac_imp300_err,
                    frac_imp400_err,
                    frac_imp500_err,
                    frac_imp750_err,
                    frac_imp1000_err]

    frac_l0001 = float(dv400_l0001_det) / float(dv400_l0001_det + dv400_l0001_non)
    frac_l0005 = float(dv400_l0005_det) / float(dv400_l0005_det + dv400_l0005_non)
    frac_l001 = float(dv400_l001_det) / float(dv400_l001_det + dv400_l001_non)
    frac_l005 = float(dv400_l005_det) / float(dv400_l005_det + dv400_l005_non)
    frac_l01  = float(dv400_l01_det)  / float(dv400_l01_det  + dv400_l01_non)
    frac_l05  = float(dv400_l05_det)  / float(dv400_l05_det  + dv400_l05_non)
    frac_l1   = float(dv400_l1_det)   / float(dv400_l1_det   + dv400_l1_non)
    frac_l5   = float(dv400_l5_det)   / float(dv400_l5_det   + dv400_l5_non)
    frac_l75   = float(dv400_l75_det)   / float(dv400_l75_det   + dv400_l75_non)
    
    frac_l0001_err = np.sqrt(dv400_l0001_det) / float(dv400_l0001_det + dv400_l0001_non)
    frac_l0005_err = np.sqrt(dv400_l0005_det) / float(dv400_l0005_det + dv400_l0005_non)
    frac_l001_err = np.sqrt(dv400_l001_det) / float(dv400_l001_det + dv400_l001_non)
    frac_l005_err = np.sqrt(dv400_l005_det) / float(dv400_l005_det + dv400_l005_non)
    frac_l01_err  = np.sqrt(dv400_l01_det)  / float(dv400_l01_det  + dv400_l01_non)
    frac_l05_err  = np.sqrt(dv400_l05_det)  / float(dv400_l05_det  + dv400_l05_non)
    frac_l1_err   = np.sqrt(dv400_l1_det)   / float(dv400_l1_det   + dv400_l1_non)
    frac_l5_err   = np.sqrt(dv400_l5_det)   / float(dv400_l5_det   + dv400_l5_non)
    frac_l75_err   = np.sqrt(dv400_l75_det)   / float(dv400_l75_det   + dv400_l75_non)
    

    likelihood_x = [0.75, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001]
    likelihood_y = [frac_l75, frac_l5, frac_l1, frac_l05, frac_l01, frac_l005, frac_l001, frac_l0005, frac_l0001]
    likelihood_y_err = [frac_l75_err,
                        frac_l5_err,
                        frac_l1_err,
                        frac_l05_err,
                        frac_l01_err,
                        frac_l005_err,
                        frac_l001_err,
                        frac_l0005_err,
                        frac_l0001_err]
    


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
        

        alpha_likelihood = 0.9
        alpha_impact = 0.9
        markerSize = 12
        lw = 2.5
        legend_size = 12
        legend_font = 12

        
        label_likelihood = r'$\rm \mathcal{L}-Detection~Fraction$'
        label_impact = r'$\rm \rho - Detection~Fraction$'

        symbol_likelihood = 'D'
        symbol_impact = 'o'
        
        color_likelihood = color_blue
        color_impact = color_coal
        
        ls_likelihood = 'dashed'
        ls_impact = 'solid'
        

        maxEW = 15000.

##########################################################################################
        print
        print
        print '---------------'
        print
        print 'dv400_imp1000_det: ', dv400_imp1000_det, dv400_imp1000_non
        print 'dv400_imp750_det: ', dv400_imp750_det, dv400_imp750_non
        print 'dv400_imp500_det: ', dv400_imp500_det, dv400_imp500_non
        print 'dv400_imp400_det: ', dv400_imp400_det, dv400_imp400_non
        print 'dv400_imp300_det: ', dv400_imp300_det, dv400_imp300_non
        print 'dv400_imp200_det: ', dv400_imp200_det, dv400_imp200_non
        print 'dv400_imp100_det: ', dv400_imp100_det, dv400_imp100_non
        print 'dv400_imp50_det: ', dv400_imp50_det, dv400_imp50_non
        print 'dv400_imp25_det: ', dv400_imp25_det, dv400_imp25_non
        print
        print
        print
        print '---------------'
        print
        print 'dv400_l0001_det: ', dv400_l0001_det, dv400_l0001_non
        print 'dv400_l0005_det: ', dv400_l0005_det, dv400_l0005_non
        print 'dv400_l001_det: ', dv400_l001_det, dv400_l001_non
        print 'dv400_l005_det: ', dv400_l005_det, dv400_l005_non
        print 'dv400_l01_det: ', dv400_l01_det, dv400_l01_non
        print 'dv400_l05_det: ', dv400_l05_det, dv400_l05_non
        print 'dv400_l1_det: ', dv400_l1_det, dv400_l1_non
        print 'dv400_l5_det: ', dv400_l5_det, dv400_l5_non
        print 'dv400_l75_det: ', dv400_l75_det, dv400_l75_non
        print
        print
        print

        # do the plotting 


        # impact detection fraction
#         plot1 = ax1.plot(impact_x,
#                 impact_y,
#                 marker=symbol_impact,
#                 c=color_impact,
#                 ms=markerSize,
#                 markeredgecolor='black',
#                 lw = lw,
#                 alpha=alpha_impact,
#                 label=label_impact)
                
                
        plot1 = ax1.errorbar(impact_x,
                    impact_y,
                    yerr=impact_y_err,
                    marker=symbol_impact,
                    c=color_impact,
                    ms=markerSize,
                    markeredgecolor='black',
                    ls = ls_impact,
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
#         plot2 = ax2.plot(likelihood_x,
#                 likelihood_y,
#                 marker=symbol_likelihood,
#                 c=color_likelihood,
#                 ms=markerSize,
#                 markeredgecolor='black',
#                 lw = lw,
#                 alpha=alpha_likelihood,
#                 label=label_likelihood)
                
        plot2 = ax2.errorbar(likelihood_x,
                    likelihood_y,
                    yerr=likelihood_y_err,
                    marker=symbol_likelihood,
                    c=color_likelihood,
                    ms=markerSize,
                    markeredgecolor='black',
                    ls = ls_likelihood,
                    lw = lw,
                    alpha=alpha_likelihood,
                    label=label_likelihood)
                
                

        ax2.set_xlabel(r'$\rm \mathcal{L}$')
        ax2.set_xscale("log")
        ax2.invert_xaxis()
        
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
        
        
        majorLocator   = MultipleLocator(0.1)
        majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
        minorLocator   = MultipleLocator(0.05)
        ax1.yaxis.set_major_locator(majorLocator)
        ax1.yaxis.set_major_formatter(majorFormatter)
        ax1.yaxis.set_minor_locator(minorLocator)
        
        ax1.set_ylabel(r'$\rm Detection~Fraction$')

        # deal with legends for multiple axes
#         plots = plot1 + plot2
#         labs = [p.get_label() for p in plots]
#         ax1.legend(plots, labs, prop={'size':12},loc='lower left',fancybox=True)


        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines
                              
        impact = mlines.Line2D([], [], color=color_impact, marker=symbol_impact,lw=2.5,
                                  markersize=legend_size, markeredgecolor='black', label=label_impact)
                              
        likelihood = mlines.Line2D([], [], color=color_likelihood, marker=symbol_likelihood,lw=2.5,
                                markeredgecolor='black', markersize=legend_size, label=label_likelihood)
                              
        plt.legend(handles=[impact, likelihood],loc='upper right', 
                                borderpad=0.8, fontsize=legend_font, fancybox=True)


#         ax1.grid(b=None,which='major',axis='both')
        ylim(0.0, 1.)
#         xlim(0, 2.5)

        if plot_detection_fraction_both_save:
            savefig('{0}/detection_fraction_strictbins.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()


##########################################################################################
##########################################################################################
##########################################################################################





##########################################################################################
##########################################################################################
    
    if plot_detection_fraction_velcut_both:
        fig = figure(figsize=(7.7,6.7))
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
        
        
        alpha_likelihood = 0.99
        alpha_likelihood_min50 = 0.5
        alpha_likelihood_min100 = 0.5
        alpha_likelihood_min200 = 0.5
        alpha_likelihood_min300 = 0.5

        alpha_impact = 0.99
        alpha_impact_min50 = 0.5
        alpha_impact_min100 = 0.5
        alpha_impact_min200 = 0.5
        alpha_impact_min300 = 0.5

        
        legend_size = 10
        legend_font = 10

        markerSize = 12
        markerSize_min50 = 6
        markerSize_min100 = 6
        markerSize_min200 = 6
        markerSize_min300 = 6

        lw = 2.5
        lw_min50 = 0.8
        lw_min100 = 0.8
        lw_min200 = 0.8
        lw_min300 = 0.8

        
        label_likelihood = r'$\rm \mathcal{L}-Detection~Fraction$'
        label_likelihood_min50 = r'$\rm \mathcal{L}-Detection~Fraction~(EW \ge 50~m\AA)$'
        label_likelihood_min100 = r'$\rm \mathcal{L}-Detection~Fraction~(EW \ge= 100~m\AA)$'
        label_likelihood_min200 = r'$\rm \mathcal{L}-Detection~Fraction~(EW \ge= 200~m\AA)$'
        label_likelihood_min300 = r'$\rm \mathcal{L}-Detection~Fraction~(EW \ge= 300~m\AA)$'

        label_impact = r'$\rm \rho - Detection~Fraction$'
        label_impact_min50 = r'$\rm \rho - Detection~Fraction~(EW \ge 50~m\AA)$'
        label_impact_min100 = r'$\rm \rho - Detection~Fraction~(EW \ge 100~m\AA)$'
        label_impact_min200 = r'$\rm \rho - Detection~Fraction~(EW \ge 200~m\AA)$'
        label_impact_min300 = r'$\rm \rho - Detection~Fraction~(EW \ge 300~m\AA)$'

        label_min50 = r'$\rm EW \ge 50~m\AA$'
        label_min100 = r'$\rm EW \ge 100~m\AA$'
        label_min200 = r'$\rm EW \ge 200~m\AA$'
        label_min300 = r'$\rm EW \ge 300~m\AA$'


        symbol_likelihood = 'D'
        symbol_likelihood_min50 = 'D'
        symbol_likelihood_min100 = 'D'
        symbol_likelihood_min200 = 'D'
        symbol_likelihood_min300 = 'D'

        symbol_impact = 'o'
        symbol_impact_min50 = 'o'
        symbol_impact_min100 = 'o'
        symbol_impact_min200 = 'o'
        symbol_impact_min300 = 'o'

        color_likelihood = color_blue
        color_likelihood_min50 = color_purple2
        color_likelihood_min100 = color_green
        color_likelihood_min200 = color_yellow
        color_likelihood_min300 = color_red

        color_impact = color_coal
        color_impact_min50 = color_purple2
        color_impact_min100 = color_green
        color_impact_min200 = color_yellow
        color_impact_min300 = color_red


        ls_likelihood = 'dashed'
        ls_likelihood_min50 = 'dashed'
        ls_likelihood_min100 = 'dashed'
        ls_likelihood_min200 = 'dashed'
        ls_likelihood_min300 = 'dashed'

        ls_impact = 'solid'
        ls_impact_min50 = 'solid'
        ls_impact_min100 = 'solid'
        ls_impact_min200 = 'solid'
        ls_impact_min300 = 'solid'

        
        maxEW = 15000.


##########################################################################################
        # prepare the data 
        frac_imp1000_min50 = float(dv400_imp1000_min50_det) / float(dv400_imp1000_min50_det + dv400_imp1000_min50_non)
        frac_imp750_min50 = float(dv400_imp750_min50_det) / float(dv400_imp750_min50_det + dv400_imp750_min50_non)
        frac_imp500_min50 = float(dv400_imp500_min50_det) / float(dv400_imp500_min50_det + dv400_imp500_min50_non)
        frac_imp400_min50 = float(dv400_imp400_min50_det) / float(dv400_imp400_min50_det + dv400_imp400_min50_non)
        frac_imp300_min50 = float(dv400_imp300_min50_det) / float(dv400_imp300_min50_det + dv400_imp300_min50_non)
        frac_imp200_min50 = float(dv400_imp200_min50_det) / float(dv400_imp200_min50_det + dv400_imp200_min50_non)
        frac_imp100_min50 = float(dv400_imp100_min50_det) / float(dv400_imp100_min50_det + dv400_imp100_min50_non)
        frac_imp50_min50  = float(dv400_imp50_min50_det)  / float(dv400_imp50_min50_det  + dv400_imp50_min50_non)
        
        frac_imp1000_min50_err = np.sqrt(dv400_imp1000_min50_det) / float(dv400_imp1000_min50_det + dv400_imp1000_min50_non)
        frac_imp750_min50_err = np.sqrt(dv400_imp750_min50_det) / float(dv400_imp750_min50_det + dv400_imp750_min50_non)
        frac_imp500_min50_err = np.sqrt(dv400_imp500_min50_det) / float(dv400_imp500_min50_det + dv400_imp500_min50_non)
        frac_imp400_min50_err = np.sqrt(dv400_imp400_min50_det) / float(dv400_imp400_min50_det + dv400_imp400_min50_non)
        frac_imp300_min50_err = np.sqrt(dv400_imp300_min50_det) / float(dv400_imp300_min50_det + dv400_imp300_min50_non)
        frac_imp200_min50_err = np.sqrt(dv400_imp200_min50_det) / float(dv400_imp200_min50_det + dv400_imp200_min50_non)
        frac_imp100_min50_err = np.sqrt(dv400_imp100_min50_det) / float(dv400_imp100_min50_det + dv400_imp100_min50_non)
        frac_imp50_min50_err  = np.sqrt(dv400_imp50_min50_det)  / float(dv400_imp50_min50_det  + dv400_imp50_min50_non)
      
        try:
            frac_imp25_min50  = float(dv400_imp25_min50_det)  / float(dv400_imp25_min50_det  + dv400_imp25_min50_non)
            frac_imp25_min50_err  = np.sqrt(dv400_imp25_min50_det)  / float(dv400_imp25_min50_det  + dv400_imp25_min50_non)

        except Exception,e:
            print 'Error: ',e
            frac_imp25_min50 = 0
            frac_imp25_min50_err = 0

        # minEW = 100
        frac_imp1000_min100 = float(dv400_imp1000_min100_det) / float(dv400_imp1000_min100_det + dv400_imp1000_min100_non)
        frac_imp750_min100 = float(dv400_imp750_min100_det) / float(dv400_imp750_min100_det + dv400_imp750_min100_non)
        frac_imp500_min100 = float(dv400_imp500_min100_det) / float(dv400_imp500_min100_det + dv400_imp500_min100_non)
        frac_imp400_min100 = float(dv400_imp400_min100_det) / float(dv400_imp400_min100_det + dv400_imp400_min100_non)
        frac_imp300_min100 = float(dv400_imp300_min100_det) / float(dv400_imp300_min100_det + dv400_imp300_min100_non)
        frac_imp200_min100 = float(dv400_imp200_min100_det) / float(dv400_imp200_min100_det + dv400_imp200_min100_non)
        frac_imp100_min100 = float(dv400_imp100_min100_det) / float(dv400_imp100_min100_det + dv400_imp100_min100_non)
        frac_imp50_min100  = float(dv400_imp50_min100_det)  / float(dv400_imp50_min100_det  + dv400_imp50_min100_non)
        
        frac_imp1000_min100_err = np.sqrt(dv400_imp1000_min100_det) / float(dv400_imp1000_min100_det + dv400_imp1000_min100_non)
        frac_imp750_min100_err = np.sqrt(dv400_imp750_min100_det) / float(dv400_imp750_min100_det + dv400_imp750_min100_non)
        frac_imp500_min100_err = np.sqrt(dv400_imp500_min100_det) / float(dv400_imp500_min100_det + dv400_imp500_min100_non)
        frac_imp400_min100_err = np.sqrt(dv400_imp400_min100_det) / float(dv400_imp400_min100_det + dv400_imp400_min100_non)
        frac_imp300_min100_err = np.sqrt(dv400_imp300_min100_det) / float(dv400_imp300_min100_det + dv400_imp300_min100_non)
        frac_imp200_min100_err = np.sqrt(dv400_imp200_min100_det) / float(dv400_imp200_min100_det + dv400_imp200_min100_non)
        frac_imp100_min100_err = np.sqrt(dv400_imp100_min100_det) / float(dv400_imp100_min100_det + dv400_imp100_min100_non)
        frac_imp50_min100_err  = np.sqrt(dv400_imp50_min100_det)  / float(dv400_imp50_min100_det  + dv400_imp50_min100_non)
      
        try:
            frac_imp25_min100  = float(dv400_imp25_min100_det)  / float(dv400_imp25_min100_det  + dv400_imp25_min100_non)
            frac_imp25_min100_err  = np.sqrt(dv400_imp25_min100_det)  / float(dv400_imp25_min100_det  + dv400_imp25_min100_non)

        except Exception,e:
            print 'Error: ',e
            frac_imp25_min100 = 0
            frac_imp25_min100_err = 0
      
      
        # minEW = 200
        frac_imp1000_min200 = float(dv400_imp1000_min200_det) / float(dv400_imp1000_min200_det + dv400_imp1000_min200_non)
        frac_imp750_min200 = float(dv400_imp750_min200_det) / float(dv400_imp750_min200_det + dv400_imp750_min200_non)
        frac_imp500_min200 = float(dv400_imp500_min200_det) / float(dv400_imp500_min200_det + dv400_imp500_min200_non)
        frac_imp400_min200 = float(dv400_imp400_min200_det) / float(dv400_imp400_min200_det + dv400_imp400_min200_non)
        frac_imp300_min200 = float(dv400_imp300_min200_det) / float(dv400_imp300_min200_det + dv400_imp300_min200_non)
        frac_imp200_min200 = float(dv400_imp200_min200_det) / float(dv400_imp200_min200_det + dv400_imp200_min200_non)
        frac_imp100_min200 = float(dv400_imp100_min200_det) / float(dv400_imp100_min200_det + dv400_imp100_min200_non)
        frac_imp50_min200  = float(dv400_imp50_min200_det)  / float(dv400_imp50_min200_det  + dv400_imp50_min200_non)
        
        frac_imp1000_min200_err = np.sqrt(dv400_imp1000_min200_det) / float(dv400_imp1000_min200_det + dv400_imp1000_min200_non)
        frac_imp750_min200_err = np.sqrt(dv400_imp750_min200_det) / float(dv400_imp750_min200_det + dv400_imp750_min200_non)
        frac_imp500_min200_err = np.sqrt(dv400_imp500_min200_det) / float(dv400_imp500_min200_det + dv400_imp500_min200_non)
        frac_imp400_min200_err = np.sqrt(dv400_imp400_min200_det) / float(dv400_imp400_min200_det + dv400_imp400_min200_non)
        frac_imp300_min200_err = np.sqrt(dv400_imp300_min200_det) / float(dv400_imp300_min200_det + dv400_imp300_min200_non)
        frac_imp200_min200_err = np.sqrt(dv400_imp200_min200_det) / float(dv400_imp200_min200_det + dv400_imp200_min200_non)
        frac_imp100_min200_err = np.sqrt(dv400_imp100_min200_det) / float(dv400_imp100_min200_det + dv400_imp100_min200_non)
        frac_imp50_min200_err  = np.sqrt(dv400_imp50_min200_det)  / float(dv400_imp50_min200_det  + dv400_imp50_min200_non)
            
        try:
            frac_imp25_min200  = float(dv400_imp25_min200_det)  / float(dv400_imp25_min200_det  + dv400_imp25_min200_non)
            frac_imp25_min200_err  = np.sqrt(dv400_imp25_min200_det)  / float(dv400_imp25_min200_det  + dv400_imp25_min200_non)

        except Exception,e:
            print 'Error: ',e
            frac_imp25_min200 = 0
            frac_imp25_min200_err = 0
            
        # minEW = 300
        frac_imp1000_min300 = float(dv400_imp1000_min300_det) / float(dv400_imp1000_min300_det + dv400_imp1000_min300_non)
        frac_imp750_min300 = float(dv400_imp750_min300_det) / float(dv400_imp750_min300_det + dv400_imp750_min300_non)
        frac_imp500_min300 = float(dv400_imp500_min300_det) / float(dv400_imp500_min300_det + dv400_imp500_min300_non)
        frac_imp400_min300 = float(dv400_imp400_min300_det) / float(dv400_imp400_min300_det + dv400_imp400_min300_non)
        frac_imp300_min300 = float(dv400_imp300_min300_det) / float(dv400_imp300_min300_det + dv400_imp300_min300_non)
        frac_imp200_min300 = float(dv400_imp200_min300_det) / float(dv400_imp200_min300_det + dv400_imp200_min300_non)
        frac_imp100_min300 = float(dv400_imp100_min300_det) / float(dv400_imp100_min300_det + dv400_imp100_min300_non)
        frac_imp50_min300  = float(dv400_imp50_min300_det)  / float(dv400_imp50_min300_det  + dv400_imp50_min300_non)
        
        frac_imp1000_min300_err = np.sqrt(dv400_imp1000_min300_det) / float(dv400_imp1000_min300_det + dv400_imp1000_min300_non)
        frac_imp750_min300_err = np.sqrt(dv400_imp750_min300_det) / float(dv400_imp750_min300_det + dv400_imp750_min300_non)
        frac_imp500_min300_err = np.sqrt(dv400_imp500_min300_det) / float(dv400_imp500_min300_det + dv400_imp500_min300_non)
        frac_imp400_min300_err = np.sqrt(dv400_imp400_min300_det) / float(dv400_imp400_min300_det + dv400_imp400_min300_non)
        frac_imp300_min300_err = np.sqrt(dv400_imp300_min300_det) / float(dv400_imp300_min300_det + dv400_imp300_min300_non)
        frac_imp200_min300_err = np.sqrt(dv400_imp200_min300_det) / float(dv400_imp200_min300_det + dv400_imp200_min300_non)
        frac_imp100_min300_err = np.sqrt(dv400_imp100_min300_det) / float(dv400_imp100_min300_det + dv400_imp100_min300_non)
        frac_imp50_min300_err  = np.sqrt(dv400_imp50_min300_det)  / float(dv400_imp50_min300_det  + dv400_imp50_min300_non)
            
        try:
            frac_imp25_min300  = float(dv400_imp25_min300_det)  / float(dv400_imp25_min300_det  + dv400_imp25_min300_non)
            frac_imp25_min300_err  = np.sqrt(dv400_imp25_min300_det)  / float(dv400_imp25_min300_det  + dv400_imp25_min300_non)

        except Exception,e:
            print 'Error: ',e
            frac_imp25_min300 = 0
            frac_imp25_min300_err = 0

        impact_x = [25, 50, 100, 200, 300, 400, 500, 750, 1000]
        impact_y_min50 = [frac_imp25_min50,
                    frac_imp50_min50,
                    frac_imp100_min50,
                    frac_imp200_min50,
                    frac_imp300_min50,
                    frac_imp400_min50,
                    frac_imp500_min50,
                    frac_imp750_min50,
                    frac_imp1000_min50]
                    
        impact_y_min50_err = [frac_imp25_min50_err,
                        frac_imp50_min50_err,
                        frac_imp100_min50_err,
                        frac_imp200_min50_err,
                        frac_imp300_min50_err,
                        frac_imp400_min50_err,
                        frac_imp500_min50_err,
                        frac_imp750_min50_err,
                        frac_imp1000_min50_err]

        impact_y_min100 = [frac_imp25_min100,
                    frac_imp50_min100,
                    frac_imp100_min100,
                    frac_imp200_min100,
                    frac_imp300_min100,
                    frac_imp400_min100,
                    frac_imp500_min100,
                    frac_imp750_min100,
                    frac_imp1000_min100]
                    
        impact_y_min100_err = [frac_imp25_min100_err,
                        frac_imp50_min100_err,
                        frac_imp100_min100_err,
                        frac_imp200_min100_err,
                        frac_imp300_min100_err,
                        frac_imp400_min100_err,
                        frac_imp500_min100_err,
                        frac_imp750_min100_err,
                        frac_imp1000_min100_err]

        impact_y_min200 = [frac_imp25_min200,
                    frac_imp50_min200,
                    frac_imp100_min200,
                    frac_imp200_min200,
                    frac_imp300_min200,
                    frac_imp400_min200,
                    frac_imp500_min200,
                    frac_imp750_min200,
                    frac_imp1000_min200]
                    
        impact_y_min200_err = [frac_imp25_min200_err,
                        frac_imp50_min200_err,
                        frac_imp100_min200_err,
                        frac_imp200_min200_err,
                        frac_imp300_min200_err,
                        frac_imp400_min200_err,
                        frac_imp500_min200_err,
                        frac_imp750_min200_err,
                        frac_imp1000_min200_err]

        impact_y_min300 = [frac_imp25_min300,
                    frac_imp50_min300,
                    frac_imp100_min300,
                    frac_imp200_min300,
                    frac_imp300_min300,
                    frac_imp400_min300,
                    frac_imp500_min300,
                    frac_imp750_min300,
                    frac_imp1000_min300]
                    
        impact_y_min300_err = [frac_imp25_min300_err,
                        frac_imp50_min300_err,
                        frac_imp100_min300_err,
                        frac_imp200_min300_err,
                        frac_imp300_min300_err,
                        frac_imp400_min300_err,
                        frac_imp500_min300_err,
                        frac_imp750_min300_err,
                        frac_imp1000_min300_err]


###################
        # now likelihoods; minEW=50 first
        frac_l0001_min50 = float(dv400_l0001_min50_det) / float(dv400_l0001_min50_det + dv400_l0001_min50_non)
        frac_l0005_min50 = float(dv400_l0005_min50_det) / float(dv400_l0005_min50_det + dv400_l0005_min50_non)
        frac_l001_min50 = float(dv400_l001_min50_det) / float(dv400_l001_min50_det + dv400_l001_min50_non)
        frac_l005_min50 = float(dv400_l005_min50_det) / float(dv400_l005_min50_det + dv400_l005_min50_non)
        frac_l01_min50 = float(dv400_l01_min50_det)  / float(dv400_l01_min50_det  + dv400_l01_min50_non)
        frac_l05_min50 = float(dv400_l05_min50_det)  / float(dv400_l05_min50_det  + dv400_l05_min50_non)
        frac_l1_min50 = float(dv400_l1_min50_det)   / float(dv400_l1_min50_det   + dv400_l1_min50_non)
        frac_l5_min50 = float(dv400_l5_min50_det)   / float(dv400_l5_min50_det   + dv400_l5_min50_non)
        frac_l75_min50 = float(dv400_l75_min50_det)   / float(dv400_l75_min50_det   + dv400_l75_min50_non)

        frac_l0001_min50_err = np.sqrt(dv400_l0001_min50_det) / float(dv400_l0001_min50_det + dv400_l0001_min50_non)
        frac_l0005_min50_err = np.sqrt(dv400_l0005_min50_det) / float(dv400_l0005_min50_det + dv400_l0005_min50_non)
        frac_l001_min50_err = np.sqrt(dv400_l001_min50_det) / float(dv400_l001_min50_det + dv400_l001_min50_non)
        frac_l005_min50_err = np.sqrt(dv400_l005_min50_det) / float(dv400_l005_min50_det + dv400_l005_min50_non)
        frac_l01_min50_err = np.sqrt(dv400_l01_min50_det)  / float(dv400_l01_min50_det  + dv400_l01_min50_non)
        frac_l05_min50_err = np.sqrt(dv400_l05_min50_det)  / float(dv400_l05_min50_det  + dv400_l05_min50_non)
        frac_l1_min50_err = np.sqrt(dv400_l1_min50_det)   / float(dv400_l1_min50_det   + dv400_l1_min50_non)
        frac_l5_min50_err = np.sqrt(dv400_l5_min50_det)   / float(dv400_l5_min50_det   + dv400_l5_min50_non)
        frac_l75_min50_err = np.sqrt(dv400_l75_min50_det)   / float(dv400_l75_min50_det   + dv400_l75_min50_non)

        # minEW = 100
        frac_l0001_min100 = float(dv400_l0001_min100_det) / float(dv400_l0001_min100_det + dv400_l0001_min100_non)
        frac_l0005_min100 = float(dv400_l0005_min100_det) / float(dv400_l0005_min100_det + dv400_l0005_min100_non)
        frac_l001_min100 = float(dv400_l001_min100_det) / float(dv400_l001_min100_det + dv400_l001_min100_non)
        frac_l005_min100 = float(dv400_l005_min100_det) / float(dv400_l005_min100_det + dv400_l005_min100_non)
        frac_l01_min100 = float(dv400_l01_min100_det)  / float(dv400_l01_min100_det  + dv400_l01_min100_non)
        frac_l05_min100 = float(dv400_l05_min100_det)  / float(dv400_l05_min100_det  + dv400_l05_min100_non)
        frac_l1_min100 = float(dv400_l1_min100_det)   / float(dv400_l1_min100_det   + dv400_l1_min100_non)
        frac_l5_min100 = float(dv400_l5_min100_det)   / float(dv400_l5_min100_det   + dv400_l5_min100_non)
        frac_l75_min100 = float(dv400_l75_min100_det)   / float(dv400_l75_min100_det   + dv400_l75_min100_non)

        frac_l0001_min100_err = np.sqrt(dv400_l0001_min100_det) / float(dv400_l0001_min100_det + dv400_l0001_min100_non)
        frac_l0005_min100_err = np.sqrt(dv400_l0005_min100_det) / float(dv400_l0005_min100_det + dv400_l0005_min100_non)
        frac_l001_min100_err = np.sqrt(dv400_l001_min100_det) / float(dv400_l001_min100_det + dv400_l001_min100_non)
        frac_l005_min100_err = np.sqrt(dv400_l005_min100_det) / float(dv400_l005_min100_det + dv400_l005_min100_non)
        frac_l01_min100_err = np.sqrt(dv400_l01_min100_det)  / float(dv400_l01_min100_det  + dv400_l01_min100_non)
        frac_l05_min100_err = np.sqrt(dv400_l05_min100_det)  / float(dv400_l05_min100_det  + dv400_l05_min100_non)
        frac_l1_min100_err = np.sqrt(dv400_l1_min100_det)   / float(dv400_l1_min100_det   + dv400_l1_min100_non)
        frac_l5_min100_err = np.sqrt(dv400_l5_min100_det)   / float(dv400_l5_min100_det   + dv400_l5_min100_non)
        frac_l75_min100_err = np.sqrt(dv400_l75_min100_det)   / float(dv400_l75_min100_det   + dv400_l75_min100_non)


        # minEW = 200
        frac_l0001_min200 = float(dv400_l0001_min200_det) / float(dv400_l0001_min200_det + dv400_l0001_min200_non)
        frac_l0005_min200 = float(dv400_l0005_min200_det) / float(dv400_l0005_min200_det + dv400_l0005_min200_non)
        frac_l001_min200 = float(dv400_l001_min200_det) / float(dv400_l001_min200_det + dv400_l001_min200_non)
        frac_l005_min200 = float(dv400_l005_min200_det) / float(dv400_l005_min200_det + dv400_l005_min200_non)
        frac_l01_min200 = float(dv400_l01_min200_det)  / float(dv400_l01_min200_det  + dv400_l01_min200_non)
        frac_l05_min200 = float(dv400_l05_min200_det)  / float(dv400_l05_min200_det  + dv400_l05_min200_non)
        frac_l1_min200 = float(dv400_l1_min200_det)   / float(dv400_l1_min200_det   + dv400_l1_min200_non)
        frac_l5_min200 = float(dv400_l5_min200_det)   / float(dv400_l5_min200_det   + dv400_l5_min200_non)
        frac_l75_min200 = float(dv400_l75_min200_det)   / float(dv400_l75_min200_det   + dv400_l75_min200_non)

        frac_l0001_min200_err = np.sqrt(dv400_l0001_min200_det) / float(dv400_l0001_min200_det + dv400_l0001_min200_non)
        frac_l0005_min200_err = np.sqrt(dv400_l0005_min200_det) / float(dv400_l0005_min200_det + dv400_l0005_min200_non)
        frac_l001_min200_err = np.sqrt(dv400_l001_min200_det) / float(dv400_l001_min200_det + dv400_l001_min200_non)
        frac_l005_min200_err = np.sqrt(dv400_l005_min200_det) / float(dv400_l005_min200_det + dv400_l005_min200_non)
        frac_l01_min200_err = np.sqrt(dv400_l01_min200_det)  / float(dv400_l01_min200_det  + dv400_l01_min200_non)
        frac_l05_min200_err = np.sqrt(dv400_l05_min200_det)  / float(dv400_l05_min200_det  + dv400_l05_min200_non)
        frac_l1_min200_err = np.sqrt(dv400_l1_min200_det)   / float(dv400_l1_min200_det   + dv400_l1_min200_non)
        frac_l5_min200_err = np.sqrt(dv400_l5_min200_det)   / float(dv400_l5_min200_det   + dv400_l5_min200_non)
        frac_l75_min200_err = np.sqrt(dv400_l75_min200_det)   / float(dv400_l75_min200_det   + dv400_l75_min200_non)


       # minEW = 300
        frac_l0001_min300 = float(dv400_l0001_min300_det) / float(dv400_l0001_min300_det + dv400_l0001_min300_non)
        frac_l0005_min300 = float(dv400_l0005_min300_det) / float(dv400_l0005_min300_det + dv400_l0005_min300_non)
        frac_l001_min300 = float(dv400_l001_min300_det) / float(dv400_l001_min300_det + dv400_l001_min300_non)
        frac_l005_min300 = float(dv400_l005_min300_det) / float(dv400_l005_min300_det + dv400_l005_min300_non)
        frac_l01_min300 = float(dv400_l01_min300_det)  / float(dv400_l01_min300_det  + dv400_l01_min300_non)
        frac_l05_min300 = float(dv400_l05_min300_det)  / float(dv400_l05_min300_det  + dv400_l05_min300_non)
        frac_l1_min300 = float(dv400_l1_min300_det)   / float(dv400_l1_min300_det   + dv400_l1_min300_non)
        frac_l5_min300 = float(dv400_l5_min300_det)   / float(dv400_l5_min300_det   + dv400_l5_min300_non)
        frac_l75_min300 = float(dv400_l75_min300_det)   / float(dv400_l75_min300_det   + dv400_l75_min300_non)

        frac_l0001_min300_err = np.sqrt(dv400_l0001_min300_det) / float(dv400_l0001_min300_det + dv400_l0001_min300_non)
        frac_l0005_min300_err = np.sqrt(dv400_l0005_min300_det) / float(dv400_l0005_min300_det + dv400_l0005_min300_non)
        frac_l001_min300_err = np.sqrt(dv400_l001_min300_det) / float(dv400_l001_min300_det + dv400_l001_min300_non)
        frac_l005_min300_err = np.sqrt(dv400_l005_min300_det) / float(dv400_l005_min300_det + dv400_l005_min300_non)
        frac_l01_min300_err = np.sqrt(dv400_l01_min300_det)  / float(dv400_l01_min300_det  + dv400_l01_min300_non)
        frac_l05_min300_err = np.sqrt(dv400_l05_min300_det)  / float(dv400_l05_min300_det  + dv400_l05_min300_non)
        frac_l1_min300_err = np.sqrt(dv400_l1_min300_det)   / float(dv400_l1_min300_det   + dv400_l1_min300_non)
        frac_l5_min300_err = np.sqrt(dv400_l5_min300_det)   / float(dv400_l5_min300_det   + dv400_l5_min300_non)
        frac_l75_min300_err = np.sqrt(dv400_l75_min300_det)   / float(dv400_l75_min300_det   + dv400_l75_min300_non)


        likelihood_x = [0.75, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001]
        likelihood_y_min50 = [frac_l75_min50,
                        frac_l5_min50,
                        frac_l1_min50,
                        frac_l05_min50,
                        frac_l01_min50,
                        frac_l005_min50,
                        frac_l001_min50,
                        frac_l0005_min50,
                        frac_l0001_min50]
                        
        likelihood_y_min50_err = [frac_l75_min50_err,
                            frac_l5_min50_err,
                            frac_l1_min50_err,
                            frac_l05_min50_err,
                            frac_l01_min50_err,
                            frac_l005_min50_err,
                            frac_l001_min50_err,
                            frac_l0005_min50_err,
                            frac_l0001_min50_err]

        likelihood_y_min100 = [frac_l75_min100,
                        frac_l5_min100,
                        frac_l1_min100,
                        frac_l05_min100,
                        frac_l01_min100,
                        frac_l005_min100,
                        frac_l001_min100,
                        frac_l0005_min100,
                        frac_l0001_min100]
                        
        likelihood_y_min100_err = [frac_l75_min100_err,
                            frac_l5_min100_err,
                            frac_l1_min100_err,
                            frac_l05_min100_err,
                            frac_l01_min100_err,
                            frac_l005_min100_err,
                            frac_l001_min100_err,
                            frac_l0005_min100_err,
                            frac_l0001_min100_err]

        likelihood_y_min200 = [frac_l75_min200,
                        frac_l5_min200,
                        frac_l1_min200,
                        frac_l05_min200,
                        frac_l01_min200,
                        frac_l005_min200,
                        frac_l001_min200,
                        frac_l0005_min200,
                        frac_l0001_min200]
                        
        likelihood_y_min200_err = [frac_l75_min200_err,
                            frac_l5_min200_err,
                            frac_l1_min200_err,
                            frac_l05_min200_err,
                            frac_l01_min200_err,
                            frac_l005_min200_err,
                            frac_l001_min200_err,
                            frac_l0005_min200_err,
                            frac_l0001_min200_err]

        likelihood_y_min300 = [frac_l75_min300,
                        frac_l5_min300,
                        frac_l1_min300,
                        frac_l05_min300,
                        frac_l01_min300,
                        frac_l005_min300,
                        frac_l001_min300,
                        frac_l0005_min300,
                        frac_l0001_min300]
                        
        likelihood_y_min300_err = [frac_l75_min300_err,
                            frac_l5_min300_err,
                            frac_l1_min300_err,
                            frac_l05_min300_err,
                            frac_l01_min300_err,
                            frac_l005_min300_err,
                            frac_l001_min300_err,
                            frac_l0005_min300_err,
                            frac_l0001_min300_err]

        # impact detection fraction
#         plot1 = ax1.plot(impact_x,
#                 impact_y,
#                 marker=symbol_impact,
#                 c=color_impact,
#                 ms=markerSize,
#                 markeredgecolor='black',
#                 lw = lw,
#                 alpha=alpha_impact,
#                 label=label_impact)
                
        # first all the stuff
        plot1 = ax1.errorbar(impact_x,
                    impact_y,
                    yerr=impact_y_err,
                    marker=symbol_impact,
                    c=color_impact,
                    ms=markerSize,
                    markeredgecolor='black',
                    ls = ls_impact,
                    lw = lw,
                    alpha=alpha_impact,
                    label=label_impact)
                
        plot1_min50 = ax1.errorbar(impact_x,
                    impact_y_min50,
                    yerr=impact_y_min50_err,
                    marker=symbol_impact_min50,
                    c=color_impact_min50,
                    ms=markerSize_min50,
                    markeredgecolor='black',
                    ls=ls_impact_min50,
                    lw=lw_min50,
                    alpha=alpha_impact_min50,
                    label=label_impact_min50)

        plot1_min100 = ax1.errorbar(impact_x,
                    impact_y_min100,
                    yerr=impact_y_min100_err,
                    marker=symbol_impact_min100,
                    c=color_impact_min100,
                    ms=markerSize_min100,
                    markeredgecolor='black',
                    ls=ls_impact_min100,
                    lw=lw_min100,
                    alpha=alpha_impact_min100,
                    label=label_impact_min100)

        plot1_min200 = ax1.errorbar(impact_x,
                    impact_y_min200,
                    yerr=impact_y_min200_err,
                    marker=symbol_impact_min200,
                    c=color_impact_min200,
                    ms=markerSize_min200,
                    markeredgecolor='black',
                    ls=ls_impact_min200,
                    lw=lw_min200,
                    alpha=alpha_impact_min200,
                    label=label_impact_min200)

        plot1_min300 = ax1.errorbar(impact_x,
                    impact_y_min300,
                    yerr=impact_y_min300_err,
                    marker=symbol_impact_min300,
                    c=color_impact_min300,
                    ms=markerSize_min300,
                    markeredgecolor='black',
                    ls=ls_impact_min300,
                    lw=lw_min300,
                    alpha=alpha_impact_min300,
                    label=label_impact_min300)

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
#         plot2 = ax2.plot(likelihood_x,
#                 likelihood_y,
#                 marker=symbol_likelihood,
#                 c=color_likelihood,
#                 ms=markerSize,
#                 markeredgecolor='black',
#                 lw = lw,
#                 alpha=alpha_likelihood,
#                 label=label_likelihood)
                
        # first all of it
        plot2 = ax2.errorbar(likelihood_x,
                    likelihood_y,
                    yerr=likelihood_y_err,
                    marker=symbol_likelihood,
                    c=color_likelihood,
                    ms=markerSize,
                    markeredgecolor='black',
                    ls = ls_likelihood,
                    lw = lw,
                    alpha=alpha_likelihood,
                    label=label_likelihood)
                
        plot2_min50 = ax2.errorbar(likelihood_x,
                    likelihood_y_min50,
                    yerr=likelihood_y_min50_err,
                    marker=symbol_likelihood_min50,
                    c=color_likelihood_min50,
                    ms=markerSize_min50,
                    markeredgecolor='black',
                    ls = ls_likelihood_min50,
                    lw = lw_min50,
                    alpha=alpha_likelihood_min50,
                    label=label_likelihood_min50)

        plot2_min100 = ax2.errorbar(likelihood_x,
                    likelihood_y_min100,
                    yerr=likelihood_y_min100_err,
                    marker=symbol_likelihood_min100,
                    c=color_likelihood_min100,
                    ms=markerSize_min100,
                    markeredgecolor='black',
                    ls = ls_likelihood_min100,
                    lw = lw_min100,
                    alpha=alpha_likelihood_min100,
                    label=label_likelihood_min100)

        plot2_min200 = ax2.errorbar(likelihood_x,
                    likelihood_y_min200,
                    yerr=likelihood_y_min200_err,
                    marker=symbol_likelihood_min200,
                    c=color_likelihood_min200,
                    ms=markerSize_min200,
                    markeredgecolor='black',
                    ls = ls_likelihood_min200,
                    lw = lw_min200,
                    alpha=alpha_likelihood_min200,
                    label=label_likelihood_min200)

        plot2_min300 = ax2.errorbar(likelihood_x,
                    likelihood_y_min300,
                    yerr=likelihood_y_min300_err,
                    marker=symbol_likelihood_min300,
                    c=color_likelihood_min300,
                    ms=markerSize_min300,
                    markeredgecolor='black',
                    ls = ls_likelihood_min300,
                    lw = lw_min300,
                    alpha=alpha_likelihood_min300,
                    label=label_likelihood_min300)

        ax2.set_xlabel(r'$\rm \mathcal{L}$')
        ax2.set_xscale("log")
        ax2.invert_xaxis()
        
        # y-axis 
        majorLocator   = MultipleLocator(0.1)
        majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
        minorLocator   = MultipleLocator(0.05)
        ax1.yaxis.set_major_locator(majorLocator)
        ax1.yaxis.set_major_formatter(majorFormatter)
        ax1.yaxis.set_minor_locator(minorLocator)
        
        ax1.set_ylabel(r'$\rm Detection~Fraction$')

        # deal with legends for multiple axes
#         plots = plot1 + plot2
#         labs = [p.get_label() for p in plots]
#         ax1.legend(plots, labs, prop={'size':12},loc='lower left',fancybox=True)
        
        
        import matplotlib.patches as mpatches
        import matplotlib.lines as mlines
                              
        impact = mlines.Line2D([], [], color=color_impact, marker=symbol_impact,lw=2.5,
                                  markersize=legend_size, markeredgecolor='black', label=label_impact)
                              
        likelihood = mlines.Line2D([], [], color=color_likelihood, marker=symbol_likelihood,lw=2.5,
                                markeredgecolor='black', markersize=legend_size, label=label_likelihood)
                                
        min50 = mlines.Line2D([], [], color=color_impact_min50, marker=None,lw=lw_min50,
                                  markersize=legend_size, markeredgecolor='black', label=label_min50)

        min100 = mlines.Line2D([], [], color=color_impact_min100, marker=None,lw=lw_min100,
                                  markersize=legend_size, markeredgecolor='black', label=label_min100)

        min200 = mlines.Line2D([], [], color=color_impact_min200, marker=None,lw=lw_min200,
                                  markersize=legend_size, markeredgecolor='black', label=label_min200)

        min300 = mlines.Line2D([], [], color=color_impact_min300, marker=None,lw=lw_min300,
                                  markersize=legend_size, markeredgecolor='black', label=label_min300)
                                  
        plt.legend(handles=[impact, likelihood, min50, min100, min200, min300],loc='lower left',
                                borderpad=0.8, fontsize=legend_font, fancybox=True)
        
        
        
        
#         import matplotlib.patches as mpatches
#         import matplotlib.lines as mlines
#                               
#         impact = mlines.Line2D([], [], color=color_impact, marker=symbol_impact,lw=2.5,
#                                   markersize=legend_size, markeredgecolor='black', label=label_impact)
#                               
#         likelihood = mlines.Line2D([], [], color=color_likelihood, marker=symbol_likelihood,lw=2.5,
#                                 markeredgecolor='black', markersize=legend_size, label=label_likelihood)
#                                 
#         impact_min50 = mlines.Line2D([], [], color=color_impact_min50, marker=symbol_impact_min50,lw=lw_min50,
#                                   markersize=legend_size, markeredgecolor='black', label=label_impact_min50)
#                               
#         likelihood_min50 = mlines.Line2D([], [], color=color_likelihood_min50, marker=symbol_likelihood_min50,lw=lw_min50,
#                                 markeredgecolor='black', markersize=legend_size, label=label_likelihood_min50)
# 
#         impact_min100 = mlines.Line2D([], [], color=color_impact_min100, marker=symbol_impact_min100,lw=lw_min100,
#                                   markersize=legend_size, markeredgecolor='black', label=label_impact_min100)
#                               
#         likelihood_min100 = mlines.Line2D([], [], color=color_likelihood_min100, marker=symbol_likelihood_min100,lw=lw_min100,
#                                 markeredgecolor='black', markersize=legend_size, label=label_likelihood_min100)
# 
#         impact_min200 = mlines.Line2D([], [], color=color_impact_min200, marker=symbol_impact_min200,lw=lw_min200,
#                                   markersize=legend_size, markeredgecolor='black', label=label_impact_min200)
#                               
#         likelihood_min200 = mlines.Line2D([], [], color=color_likelihood_min200, marker=symbol_likelihood_min200,lw=lw_min200,
#                                 markeredgecolor='black', markersize=legend_size, label=label_likelihood_min200)
#                                 
#         plt.legend(handles=[impact, likelihood, impact_min50, likelihood_min50, impact_min100, likelihood_min100, impact_min200, likelihood_min200],loc='lower left', 
#                                 borderpad=0.8, fontsize=legend_font, fancybox=True)
        
        
        
        
#         from matplotlib.legend_handler import HandlerBase
# 
#         class AnyObjectHandler(HandlerBase):
#             def create_artists(self, legend, orig_handle,
#                                x0, y0, width, height, fontsize, trans):
#                 l1 = plt.Line2D([x0,y0+width], [0.7*height,0.7*height],
#                                    linestyle=orig_handle[1], color='k')
#                 l2 = plt.Line2D([x0,y0+width], [0.3*height,0.3*height], 
#                                    color=orig_handle[0])
#                 return [l1, l2]


#         x = np.linspace(0, 3)
#         fig, axL = plt.subplots(figsize=(4,3))
#         axR = axL.twinx()
# 
#         axL.plot(x, np.sin(x), color='k', linestyle='--')
#         axR.plot(x, 100*np.cos(x), color='r')
# 
#         axL.plot(x, .3*np.sin(x), color='k', linestyle=':')
#         axR.plot(x, 20*np.cos(x), color='limegreen')
# 
#         axL.set_ylabel('sin(x)', color='k')
#         axR.set_ylabel('100 cos(x)', color='r')
#         axR.tick_params('y', colors='r')
# 
#         plt.legend([("r","--"), ("limegreen",":")], ['label', "label2"],
#                    handler_map={tuple: AnyObjectHandler()})
        
        
#         plt.legend([(color_likelihood_min50, ls_likelihood_min50),
#                     (color_likelihood_min100, ls_likelihood_min100),
#                     (color_likelihood_min200, ls_likelihood_min200)], 
#                     [label_likelihood_min50, label_likelihood_min100, label_likelihood_min200],
#                    handler_map={tuple: AnyObjectHandler()},
#                    loc='lower left',
#                    borderpad=0.8,
#                    fontsize=legend_font,
#                    fancybox=True)
        
        
        
        

#         ax1.grid(b=None,which='major',axis='both')
        ylim(0., 1.)
#         xlim(0, 2.5)

        if plot_detection_fraction_velcut_both_save:
            print 'saving...'
            savefig('{0}/detection_fraction_min0_50_100_200_300_both.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()


##########################################################################################
##########################################################################################
##########################################################################################



#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


if __name__=="__main__":
    # do the work
    main()
    