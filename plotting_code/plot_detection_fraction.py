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
    for i in resampled_medians:
        med = np.median(i)
        resampled_medians.append(med)
    
    median_err = std(resampled_medians)
        

    return mean_err, median_err


    
def main():
    # plot detection fraction as a function of both impact parameter and likelihood - 
    # does not work right now
    plot_detection_fraction_both = True
    plot_detection_fraction_both_save = True
    
    # plot detection fraction as a function of impact parameter
    plot_detection_fraction_impact = False
    plot_detection_fraction_impact_save = False
    
    # plot detection fraction as a function of likelihood
    plot_detection_fraction_likelihood = False
    plot_detection_fraction_likelihood_save = False
    
    # plot detection fraction as a function of likelihood
    plot_detection_fraction_likelihood_inc = True
    plot_detection_fraction_likelihood_inc_save = True

    # plot detection fraction as a function of likelihood
    plot_detection_fraction_impact_inc = True
    plot_detection_fraction_impact_inc_save = True

    # plot detection fraction as a function of likelihood
    plot_detection_fraction_likelihood_lstar = True
    plot_detection_fraction_likelihood_lstar_save = True

    # plot detection fraction as a function of likelihood
    plot_detection_fraction_impact_lstar = True
    plot_detection_fraction_impact_lstar_save = True
    
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
        
#         detection_fraction_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction3.p'
        detection_fraction_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_lstarcut_{0}.p'.format(lstar_cut)

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
    
    # now impact parameter non-detection inclinations 
    dv400_imp1000_non_inc = np.array(full_dict['dv400_imp1000_non_inc'])
    dv400_imp750_non_inc = np.array(full_dict['dv400_imp750_non_inc'])
    dv400_imp500_non_inc = np.array(full_dict['dv400_imp500_non_inc'])
    dv400_imp400_non_inc = np.array(full_dict['dv400_imp400_non_inc'])
    dv400_imp300_non_inc = np.array(full_dict['dv400_imp300_non_inc'])
    dv400_imp200_non_inc = np.array(full_dict['dv400_imp200_non_inc'])
    dv400_imp100_non_inc = np.array(full_dict['dv400_imp100_non_inc'])
    dv400_imp50_non_inc = np.array(full_dict['dv400_imp50_non_inc'])
    dv400_imp25_non_inc = np.array(full_dict['dv400_imp25_non_inc'])


    # impact parameter detection lstars
    dv400_imp1000_det_lstar = np.array(full_dict['dv400_imp1000_det_lstar'])
    dv400_imp750_det_lstar = np.array(full_dict['dv400_imp750_det_lstar'])
    dv400_imp500_det_lstar = np.array(full_dict['dv400_imp500_det_lstar'])
    dv400_imp400_det_lstar = np.array(full_dict['dv400_imp400_det_lstar'])
    dv400_imp300_det_lstar = np.array(full_dict['dv400_imp300_det_lstar'])
    dv400_imp200_det_lstar = np.array(full_dict['dv400_imp200_det_lstar'])
    dv400_imp100_det_lstar = np.array(full_dict['dv400_imp100_det_lstar'])
    dv400_imp50_det_lstar = np.array(full_dict['dv400_imp50_det_lstar'])
    dv400_imp25_det_lstar = np.array(full_dict['dv400_imp25_det_lstar'])
    
    # now impact parameter non-detection lstars
    dv400_imp1000_non_lstar = np.array(full_dict['dv400_imp1000_non_lstar'])
    dv400_imp750_non_lstar = np.array(full_dict['dv400_imp750_non_lstar'])
    dv400_imp500_non_lstar = np.array(full_dict['dv400_imp500_non_lstar'])
    dv400_imp400_non_lstar = np.array(full_dict['dv400_imp400_non_lstar'])
    dv400_imp300_non_lstar = np.array(full_dict['dv400_imp300_non_lstar'])
    dv400_imp200_non_lstar = np.array(full_dict['dv400_imp200_non_lstar'])
    dv400_imp100_non_lstar = np.array(full_dict['dv400_imp100_non_lstar'])
    dv400_imp50_non_lstar = np.array(full_dict['dv400_imp50_non_lstar'])
    dv400_imp25_non_lstar = np.array(full_dict['dv400_imp25_non_lstar'])


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

    # now for likelihood detections inclinations
    dv400_l0001_det_inc = np.array(full_dict['dv400_l0001_det_inc'])
    dv400_l0005_det_inc = np.array(full_dict['dv400_l0005_det_inc'])
    dv400_l001_det_inc = np.array(full_dict['dv400_l001_det_inc'])
    dv400_l005_det_inc = np.array(full_dict['dv400_l005_det_inc'])
    dv400_l01_det_inc = np.array(full_dict['dv400_l01_det_inc'])
    dv400_l05_det_inc = np.array(full_dict['dv400_l05_det_inc'])
    dv400_l1_det_inc = np.array(full_dict['dv400_l1_det_inc'])
    dv400_l5_det_inc = np.array(full_dict['dv400_l5_det_inc'])
    dv400_l75_det_inc = np.array(full_dict['dv400_l75_det_inc'])

    # now for likelihood non-detections inclinations
    dv400_l0001_non_inc = np.array(full_dict['dv400_l0001_non_inc'])
    dv400_l0005_non_inc = np.array(full_dict['dv400_l0005_non_inc'])
    dv400_l001_non_inc = np.array(full_dict['dv400_l001_non_inc'])
    dv400_l005_non_inc = np.array(full_dict['dv400_l005_non_inc'])
    dv400_l01_non_inc = np.array(full_dict['dv400_l01_non_inc'])
    dv400_l05_non_inc = np.array(full_dict['dv400_l05_non_inc'])
    dv400_l1_non_inc = np.array(full_dict['dv400_l1_non_inc'])
    dv400_l5_non_inc = np.array(full_dict['dv400_l5_non_inc'])
    dv400_l75_non_inc = np.array(full_dict['dv400_l75_non_inc'])
    
    # now for likelihood detections lstars
    dv400_l0001_det_lstar = np.array(full_dict['dv400_l0001_det_lstar'])
    dv400_l0005_det_lstar = np.array(full_dict['dv400_l0005_det_lstar'])
    dv400_l001_det_lstar = np.array(full_dict['dv400_l001_det_lstar'])
    dv400_l005_det_lstar = np.array(full_dict['dv400_l005_det_lstar'])
    dv400_l01_det_lstar = np.array(full_dict['dv400_l01_det_lstar'])
    dv400_l05_det_lstar = np.array(full_dict['dv400_l05_det_lstar'])
    dv400_l1_det_lstar = np.array(full_dict['dv400_l1_det_lstar'])
    dv400_l5_det_lstar = np.array(full_dict['dv400_l5_det_lstar'])
    dv400_l75_det_lstar = np.array(full_dict['dv400_l75_det_lstar'])

    # now for likelihood non-detections lstars
    dv400_l0001_non_lstar = np.array(full_dict['dv400_l0001_non_lstar'])
    dv400_l0005_non_lstar = np.array(full_dict['dv400_l0005_non_lstar'])
    dv400_l001_non_lstar = np.array(full_dict['dv400_l001_non_lstar'])
    dv400_l005_non_lstar = np.array(full_dict['dv400_l005_non_lstar'])
    dv400_l01_non_lstar = np.array(full_dict['dv400_l01_non_lstar'])
    dv400_l05_non_lstar = np.array(full_dict['dv400_l05_non_lstar'])
    dv400_l1_non_lstar = np.array(full_dict['dv400_l1_non_lstar'])
    dv400_l5_non_lstar = np.array(full_dict['dv400_l5_non_lstar'])
    dv400_l75_non_lstar = np.array(full_dict['dv400_l75_non_lstar'])
    
    
    
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
    
    dv400_l0001_det_inc = dv400_l0001_det_inc.astype(np.float)
    dv400_l0005_det_inc = dv400_l0005_det_inc.astype(np.float)
    dv400_l001_det_inc = dv400_l001_det_inc.astype(np.float)
    dv400_l005_det_inc = dv400_l005_det_inc.astype(np.float)
    dv400_l01_det_inc = dv400_l01_det_inc.astype(np.float)
    dv400_l05_det_inc = dv400_l05_det_inc.astype(np.float)
    dv400_l1_det_inc = dv400_l1_det_inc.astype(np.float)
    dv400_l5_det_inc = dv400_l5_det_inc.astype(np.float)
    dv400_l75_det_inc = dv400_l75_det_inc.astype(np.float)
    
    dv400_l0001_non_inc = dv400_l0001_non_inc.astype(np.float)
    dv400_l0005_non_inc = dv400_l0005_non_inc.astype(np.float)
    dv400_l001_non_inc = dv400_l001_non_inc.astype(np.float)
    dv400_l005_non_inc = dv400_l005_non_inc.astype(np.float)
    dv400_l01_non_inc = dv400_l01_non_inc.astype(np.float)
    dv400_l05_non_inc = dv400_l05_non_inc.astype(np.float)
    dv400_l1_non_inc = dv400_l1_non_inc.astype(np.float)
    dv400_l5_non_inc = dv400_l5_non_inc.astype(np.float)
    dv400_l75_non_inc = dv400_l75_non_inc.astype(np.float)
    
    
    
    #######
    # cut out '-99' from all the data
    # likelihood version first - inclinations
    dv400_l0001_det_inc = list(filter(lambda x: x!= -99., dv400_l0001_det_inc))
    dv400_l0005_det_inc = list(filter(lambda x: x!= -99., dv400_l0005_det_inc))
    dv400_l001_det_inc = list(filter(lambda x: x!= -99., dv400_l001_det_inc))
    dv400_l005_det_inc = list(filter(lambda x: x!= -99., dv400_l005_det_inc))
    dv400_l01_det_inc = list(filter(lambda x: x!= -99., dv400_l01_det_inc))
    dv400_l05_det_inc = list(filter(lambda x: x!= -99., dv400_l05_det_inc))
    dv400_l1_det_inc = list(filter(lambda x: x!= -99., dv400_l1_det_inc))
    dv400_l5_det_inc = list(filter(lambda x: x!= -99., dv400_l5_det_inc))
    dv400_l75_det_inc = list(filter(lambda x: x!= -99., dv400_l75_det_inc))

    dv400_l0001_non_inc = list(filter(lambda x: x!= -99., dv400_l0001_non_inc))
    dv400_l0005_non_inc = list(filter(lambda x: x!= -99., dv400_l0005_non_inc))
    dv400_l001_non_inc = list(filter(lambda x: x!= -99., dv400_l001_non_inc))
    dv400_l005_non_inc = list(filter(lambda x: x!= -99., dv400_l005_non_inc))
    dv400_l01_non_inc = list(filter(lambda x: x!= -99., dv400_l01_non_inc))
    dv400_l05_non_inc = list(filter(lambda x: x!= -99., dv400_l05_non_inc))
    dv400_l1_non_inc = list(filter(lambda x: x!= -99., dv400_l1_non_inc))
    dv400_l5_non_inc = list(filter(lambda x: x!= -99., dv400_l5_non_inc))
    dv400_l75_non_inc = list(filter(lambda x: x!= -99., dv400_l75_non_inc))
    
    # likelihood version - lstars
    dv400_l0001_det_lstar = list(filter(lambda x: x!= -99., dv400_l0001_det_lstar))
    dv400_l0005_det_lstar = list(filter(lambda x: x!= -99., dv400_l0005_det_lstar))
    dv400_l001_det_lstar = list(filter(lambda x: x!= -99., dv400_l001_det_lstar))
    dv400_l005_det_lstar = list(filter(lambda x: x!= -99., dv400_l005_det_lstar))
    dv400_l01_det_lstar = list(filter(lambda x: x!= -99., dv400_l01_det_lstar))
    dv400_l05_det_lstar = list(filter(lambda x: x!= -99., dv400_l05_det_lstar))
    dv400_l1_det_lstar = list(filter(lambda x: x!= -99., dv400_l1_det_lstar))
    dv400_l5_det_lstar = list(filter(lambda x: x!= -99., dv400_l5_det_lstar))
    dv400_l75_det_lstar = list(filter(lambda x: x!= -99., dv400_l75_det_lstar))

    dv400_l0001_non_lstar = list(filter(lambda x: x!= -99., dv400_l0001_non_lstar))
    dv400_l0005_non_lstar = list(filter(lambda x: x!= -99., dv400_l0005_non_lstar))
    dv400_l001_non_lstar = list(filter(lambda x: x!= -99., dv400_l001_non_lstar))
    dv400_l005_non_lstar = list(filter(lambda x: x!= -99., dv400_l005_non_lstar))
    dv400_l01_non_lstar = list(filter(lambda x: x!= -99., dv400_l01_non_lstar))
    dv400_l05_non_lstar = list(filter(lambda x: x!= -99., dv400_l05_non_lstar))
    dv400_l1_non_lstar = list(filter(lambda x: x!= -99., dv400_l1_non_lstar))
    dv400_l5_non_lstar = list(filter(lambda x: x!= -99., dv400_l5_non_lstar))
    dv400_l75_non_lstar = list(filter(lambda x: x!= -99., dv400_l75_non_lstar))
    
    # now impact version - inclinations
    dv400_imp1000_det_inc = list(filter(lambda x: x!= -99., dv400_imp1000_det_inc))
    dv400_imp750_det_inc = list(filter(lambda x: x!= -99., dv400_imp750_det_inc))
    dv400_imp500_det_inc = list(filter(lambda x: x!= -99., dv400_imp500_det_inc))
    dv400_imp400_det_inc = list(filter(lambda x: x!= -99., dv400_imp400_det_inc))
    dv400_imp300_det_inc = list(filter(lambda x: x!= -99., dv400_imp300_det_inc))
    dv400_imp200_det_inc = list(filter(lambda x: x!= -99., dv400_imp200_det_inc))
    dv400_imp100_det_inc = list(filter(lambda x: x!= -99., dv400_imp100_det_inc))
    dv400_imp50_det_inc = list(filter(lambda x: x!= -99., dv400_imp50_det_inc))
    dv400_imp25_det_inc = list(filter(lambda x: x!= -99., dv400_imp25_det_inc))

    dv400_imp1000_non_inc = list(filter(lambda x: x!= -99., dv400_imp1000_non_inc))
    dv400_imp750_non_inc = list(filter(lambda x: x!= -99., dv400_imp750_non_inc))
    dv400_imp500_non_inc = list(filter(lambda x: x!= -99., dv400_imp500_non_inc))
    dv400_imp400_non_inc = list(filter(lambda x: x!= -99., dv400_imp400_non_inc))
    dv400_imp300_non_inc = list(filter(lambda x: x!= -99., dv400_imp300_non_inc))
    dv400_imp200_non_inc = list(filter(lambda x: x!= -99., dv400_imp200_non_inc))
    dv400_imp100_non_inc = list(filter(lambda x: x!= -99., dv400_imp100_non_inc))
    dv400_imp50_non_inc = list(filter(lambda x: x!= -99., dv400_imp50_non_inc))
    dv400_imp25_non_inc = list(filter(lambda x: x!= -99., dv400_imp25_non_inc))
    
    # now impact version - lstars
    dv400_imp1000_det_lstar = list(filter(lambda x: x!= -99., dv400_imp1000_det_lstar))
    dv400_imp750_det_lstar = list(filter(lambda x: x!= -99., dv400_imp750_det_lstar))
    dv400_imp500_det_lstar = list(filter(lambda x: x!= -99., dv400_imp500_det_lstar))
    dv400_imp400_det_lstar = list(filter(lambda x: x!= -99., dv400_imp400_det_lstar))
    dv400_imp300_det_lstar = list(filter(lambda x: x!= -99., dv400_imp300_det_lstar))
    dv400_imp200_det_lstar = list(filter(lambda x: x!= -99., dv400_imp200_det_lstar))
    dv400_imp100_det_lstar = list(filter(lambda x: x!= -99., dv400_imp100_det_lstar))
    dv400_imp50_det_lstar = list(filter(lambda x: x!= -99., dv400_imp50_det_lstar))
    dv400_imp25_det_lstar = list(filter(lambda x: x!= -99., dv400_imp25_det_lstar))

    dv400_imp1000_non_lstar = list(filter(lambda x: x!= -99., dv400_imp1000_non_lstar))
    dv400_imp750_non_lstar = list(filter(lambda x: x!= -99., dv400_imp750_non_lstar))
    dv400_imp500_non_lstar = list(filter(lambda x: x!= -99., dv400_imp500_non_lstar))
    dv400_imp400_non_lstar = list(filter(lambda x: x!= -99., dv400_imp400_non_lstar))
    dv400_imp300_non_lstar = list(filter(lambda x: x!= -99., dv400_imp300_non_lstar))
    dv400_imp200_non_lstar = list(filter(lambda x: x!= -99., dv400_imp200_non_lstar))
    dv400_imp100_non_lstar = list(filter(lambda x: x!= -99., dv400_imp100_non_lstar))
    dv400_imp50_non_lstar = list(filter(lambda x: x!= -99., dv400_imp50_non_lstar))
    dv400_imp25_non_lstar = list(filter(lambda x: x!= -99., dv400_imp25_non_lstar))
    
    
    
    print 'dv400_imp1000_det_inc: ',dv400_imp1000_det_inc


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
    
    print 'dv400_imp25_det_inc: ',dv400_imp25_det_inc
    print 'dv400_imp25_non_inc: ',dv400_imp25_non_inc
    print 'dv400_l75_det_inc: ',dv400_l75_det_inc
    print 'dv400_l75_non_inc: ',dv400_l75_non_inc
    print
    
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
    print 'Detection fraction for 0.0001 L mean inc: ', bmean(np.array(dv400_l0001_det_inc)), ', ',bmean(np.array(dv400_l0001_non_inc))
    print 'Detection fraction for 0.0005 L mean inc: ', bmean(np.array(dv400_l0005_det_inc)), ', ',bmean(np.array(dv400_l0005_non_inc))
    print 'Detection fraction for 0.001 L mean inc: ', bmean(np.array(dv400_l001_det_inc)), ', ',bmean(np.array(dv400_l001_non_inc))
    print 'Detection fraction for 0.00 L mean inc: ', bmean(np.array(dv400_l005_det_inc)), ', ',bmean(np.array(dv400_l005_non_inc))
    print 'Detection fraction for 0.01 L mean inc: ', bmean(np.array(dv400_l01_det_inc)), ', ',bmean(np.array(dv400_l01_non_inc))
    print 'Detection fraction for 0.05 L mean inc: ', bmean(np.array(dv400_l05_det_inc)), ', ',bmean(np.array(dv400_l05_non_inc))
    print 'Detection fraction for 0.1 L mean inc: ', bmean(np.array(dv400_l1_det_inc)), ', ',bmean(np.array(dv400_l1_non_inc))
    print 'Detection fraction for 0.5 L mean inc: ', bmean(np.array(dv400_l5_det_inc)), ', ',bmean(np.array(dv400_l5_non_inc))
    print 'Detection fraction for 0.75 L mean inc: ', bmean(np.array(dv400_l75_det_inc)), ', ',bmean(np.array(dv400_l75_non_inc))



#########################################################################################
#########################################################################################


#########################################################################################
#########################################################################################
    # bootstrap this shit
    
#     reps = 1000
#     xb = np.random.choice(x, (n, reps), replace=True)
#     yb = 1/np.arange(1, n+1)[:, None] * np.cumsum(xb, axis=0)
#     upper, lower = np.percentile(yb, [2.5, 97.5], axis=1)

    reps = 10000
    
    # detection fraction - likelihood
    dv400_l0001_det_meanerr, dv400_l0001_det_medianerr = return_bootstrap_errors(dv400_l0001_det, reps)
    dv400_l0005_det_meanerr, dv400_l0005_det_medianerr = return_bootstrap_errors(dv400_l0005_det, reps)
    dv400_l001_det_meanerr, dv400_l001_det_medianerr = return_bootstrap_errors(dv400_l001_det, reps)
    dv400_l005_det_meanerr, dv400_l005_det_medianerr = return_bootstrap_errors(dv400_l005_det, reps)
    dv400_l01_det_meanerr, dv400_l01_det_medianerr = return_bootstrap_errors(dv400_l01_det, reps)
    dv400_l05_det_meanerr, dv400_l05_det_medianerr = return_bootstrap_errors(dv400_l05_det, reps)
    dv400_l1_det_meanerr, dv400_l1_det_medianerr = return_bootstrap_errors(dv400_l1_det, reps)
    dv400_l5_det_meanerr, dv400_l5_det_medianerr = return_bootstrap_errors(dv400_l5_det, reps)
    dv400_l75_det_meanerr, dv400_l75_det_medianerr = return_bootstrap_errors(dv400_l75_det, reps)
    
    dv400_l0001_non_meanerr, dv400_l0001_non_medianerr = return_bootstrap_errors(dv400_l0001_non, reps)
    dv400_l0005_non_meanerr, dv400_l0005_non_medianerr = return_bootstrap_errors(dv400_l0005_non, reps)
    dv400_l001_non_meanerr, dv400_l001_non_medianerr = return_bootstrap_errors(dv400_l001_non, reps)
    dv400_l005_non_meanerr, dv400_l005_non_medianerr = return_bootstrap_errors(dv400_l005_non, reps)
    dv400_l01_non_meanerr, dv400_l01_non_medianerr = return_bootstrap_errors(dv400_l01_non, reps)
    dv400_l05_non_meanerr, dv400_l05_non_medianerr = return_bootstrap_errors(dv400_l05_non, reps)
    dv400_l1_non_meanerr, dv400_l1_non_medianerr = return_bootstrap_errors(dv400_l1_non, reps)
    dv400_l5_non_meanerr, dv400_l5_non_medianerr = return_bootstrap_errors(dv400_l5_non, reps)
    dv400_l75_non_meanerr, dv400_l75_non_medianerr = return_bootstrap_errors(dv400_l75_non, reps)
    

    # detection fraction - impact
    dv400_imp1000_det_meanerr, dv400_imp1000_det_medianerr = return_bootstrap_errors(dv400_imp1000_det_lstar, reps)
    dv400_imp750_det_meanerr, dv400_imp750_det_medianerr = return_bootstrap_errors(dv400_imp750_det_lstar, reps)
    dv400_imp500_det_meanerr, dv400_imp500_det_medianerr = return_bootstrap_errors(dv400_imp500_det_lstar, reps)
    dv400_imp400_det_meanerr, dv400_imp400_det_medianerr = return_bootstrap_errors(dv400_imp400_det_lstar, reps)
    dv400_imp300_det_meanerr, dv400_imp300_det_medianerr = return_bootstrap_errors(dv400_imp300_det_lstar, reps)
    dv400_imp200_det_meanerr, dv400_imp200_det_medianerr = return_bootstrap_errors(dv400_imp200_det_lstar, reps)
    dv400_imp100_det_meanerr, dv400_imp100_det_medianerr = return_bootstrap_errors(dv400_imp100_det_lstar, reps)
    dv400_imp50_det_meanerr, dv400_imp50_det_medianerr = return_bootstrap_errors(dv400_imp50_det_lstar, reps)
    dv400_imp25_det_meanerr, dv400_imp25_det_medianerr = return_bootstrap_errors(dv400_imp25_det_lstar, reps)
    
    dv400_imp1000_non_meanerr, dv400_imp1000_non_medianerr = return_bootstrap_errors(dv400_imp1000_non_lstar, reps)
    dv400_imp750_non_meanerr, dv400_imp750_non_medianerr = return_bootstrap_errors(dv400_imp750_non_lstar, reps)
    dv400_imp500_non_meanerr, dv400_imp500_non_medianerr = return_bootstrap_errors(dv400_imp500_non_lstar, reps)
    dv400_imp400_non_meanerr, dv400_imp400_non_medianerr = return_bootstrap_errors(dv400_imp400_non_lstar, reps)
    dv400_imp300_non_meanerr, dv400_imp300_non_medianerr = return_bootstrap_errors(dv400_imp300_non_lstar, reps)
    dv400_imp200_non_meanerr, dv400_imp200_non_medianerr = return_bootstrap_errors(dv400_imp200_non_lstar, reps)
    dv400_imp100_non_meanerr, dv400_imp100_non_medianerr = return_bootstrap_errors(dv400_imp100_non_lstar, reps)
    dv400_imp50_non_meanerr, dv400_imp50_non_medianerr = return_bootstrap_errors(dv400_imp50_non_lstar, reps)
    dv400_imp25_non_meanerr, dv400_imp25_non_medianerr = return_bootstrap_errors(dv400_imp25_non_lstar, reps)
    
    
    
    # likelihood - Lstars
    dv400_l0001_det_lstar_meanerr, dv400_l0001_det_lstar_medianerr = return_bootstrap_errors(dv400_l0001_det_lstar, reps)
    dv400_l0005_det_lstar_meanerr, dv400_l0005_det_lstar_medianerr = return_bootstrap_errors(dv400_l0005_det_lstar, reps)
    dv400_l001_det_lstar_meanerr, dv400_l001_det_lstar_medianerr = return_bootstrap_errors(dv400_l001_det_lstar, reps)
    dv400_l005_det_lstar_meanerr, dv400_l005_det_lstar_medianerr = return_bootstrap_errors(dv400_l005_det_lstar, reps)
    dv400_l01_det_lstar_meanerr, dv400_l01_det_lstar_medianerr = return_bootstrap_errors(dv400_l01_det_lstar, reps)
    dv400_l05_det_lstar_meanerr, dv400_l05_det_lstar_medianerr = return_bootstrap_errors(dv400_l05_det_lstar, reps)
    dv400_l1_det_lstar_meanerr, dv400_l1_det_lstar_medianerr = return_bootstrap_errors(dv400_l1_det_lstar, reps)
    dv400_l5_det_lstar_meanerr, dv400_l5_det_lstar_medianerr = return_bootstrap_errors(dv400_l5_det_lstar, reps)
    dv400_l75_det_lstar_meanerr, dv400_l75_det_lstar_medianerr = return_bootstrap_errors(dv400_l75_det_lstar, reps)

    dv400_l0001_non_lstar_meanerr, dv400_l0001_non_lstar_medianerr = return_bootstrap_errors(dv400_l0001_non_lstar, reps)
    dv400_l0005_non_lstar_meanerr, dv400_l0005_non_lstar_medianerr = return_bootstrap_errors(dv400_l0005_non_lstar, reps)
    dv400_l001_non_lstar_meanerr, dv400_l001_non_lstar_medianerr = return_bootstrap_errors(dv400_l001_non_lstar, reps)
    dv400_l005_non_lstar_meanerr, dv400_l005_non_lstar_medianerr = return_bootstrap_errors(dv400_l005_non_lstar, reps)
    dv400_l01_non_lstar_meanerr, dv400_l01_non_lstar_medianerr = return_bootstrap_errors(dv400_l01_non_lstar, reps)
    dv400_l05_non_lstar_meanerr, dv400_l05_non_lstar_medianerr = return_bootstrap_errors(dv400_l05_non_lstar, reps)
    dv400_l1_non_lstar_meanerr, dv400_l1_non_lstar_medianerr = return_bootstrap_errors(dv400_l1_non_lstar, reps)
    dv400_l5_non_lstar_meanerr, dv400_l5_non_lstar_medianerr = return_bootstrap_errors(dv400_l5_non_lstar, reps)
    dv400_l75_non_lstar_meanerr, dv400_l75_non_lstar_medianerr = return_bootstrap_errors(dv400_l75_non_lstar, reps)

    # likelihood - inclinations
    dv400_l0001_det_inc_meanerr, dv400_l0001_det_inc_medianerr = return_bootstrap_errors(dv400_l0001_det_inc, reps)
    dv400_l0005_det_inc_meanerr, dv400_l0005_det_inc_medianerr = return_bootstrap_errors(dv400_l0005_det_inc, reps)
    dv400_l001_det_inc_meanerr, dv400_l001_det_inc_medianerr = return_bootstrap_errors(dv400_l001_det_inc, reps)
    dv400_l005_det_inc_meanerr, dv400_l005_det_inc_medianerr = return_bootstrap_errors(dv400_l005_det_inc, reps)
    dv400_l01_det_inc_meanerr, dv400_l01_det_inc_medianerr = return_bootstrap_errors(dv400_l01_det_inc, reps)
    dv400_l05_det_inc_meanerr, dv400_l05_det_inc_medianerr = return_bootstrap_errors(dv400_l05_det_inc, reps)
    dv400_l1_det_inc_meanerr, dv400_l1_det_inc_medianerr = return_bootstrap_errors(dv400_l1_det_inc, reps)
    dv400_l5_det_inc_meanerr, dv400_l5_det_inc_medianerr = return_bootstrap_errors(dv400_l5_det_inc, reps)
    dv400_l75_det_inc_meanerr, dv400_l75_det_inc_medianerr = return_bootstrap_errors(dv400_l75_det_inc, reps)

    dv400_l0001_non_inc_meanerr, dv400_l0001_non_inc_medianerr = return_bootstrap_errors(dv400_l0001_non_inc, reps)
    dv400_l0005_non_inc_meanerr, dv400_l0005_non_inc_medianerr = return_bootstrap_errors(dv400_l0005_non_inc, reps)
    dv400_l001_non_inc_meanerr, dv400_l001_non_inc_medianerr = return_bootstrap_errors(dv400_l001_non_inc, reps)
    dv400_l005_non_inc_meanerr, dv400_l005_non_inc_medianerr = return_bootstrap_errors(dv400_l005_non_inc, reps)
    dv400_l01_non_inc_meanerr, dv400_l01_non_inc_medianerr = return_bootstrap_errors(dv400_l01_non_inc, reps)
    dv400_l05_non_inc_meanerr, dv400_l05_non_inc_medianerr = return_bootstrap_errors(dv400_l05_non_inc, reps)
    dv400_l1_non_inc_meanerr, dv400_l1_non_inc_medianerr = return_bootstrap_errors(dv400_l1_non_inc, reps)
    dv400_l5_non_inc_meanerr, dv400_l5_non_inc_medianerr = return_bootstrap_errors(dv400_l5_non_inc, reps)
    dv400_l75_non_inc_meanerr, dv400_l75_non_inc_medianerr = return_bootstrap_errors(dv400_l75_non_inc, reps)


    # impact - Lstars
    dv400_imp1000_det_lstar_meanerr, dv400_imp1000_det_lstar_medianerr = return_bootstrap_errors(dv400_imp1000_det_lstar, reps)
    dv400_imp750_det_lstar_meanerr, dv400_imp750_det_lstar_medianerr = return_bootstrap_errors(dv400_imp750_det_lstar, reps)
    dv400_imp500_det_lstar_meanerr, dv400_imp500_det_lstar_medianerr = return_bootstrap_errors(dv400_imp500_det_lstar, reps)
    dv400_imp400_det_lstar_meanerr, dv400_imp400_det_lstar_medianerr = return_bootstrap_errors(dv400_imp400_det_lstar, reps)
    dv400_imp300_det_lstar_meanerr, dv400_imp300_det_lstar_medianerr = return_bootstrap_errors(dv400_imp300_det_lstar, reps)
    dv400_imp200_det_lstar_meanerr, dv400_imp200_det_lstar_medianerr = return_bootstrap_errors(dv400_imp200_det_lstar, reps)
    dv400_imp100_det_lstar_meanerr, dv400_imp100_det_lstar_medianerr = return_bootstrap_errors(dv400_imp100_det_lstar, reps)
    dv400_imp50_det_lstar_meanerr, dv400_imp50_det_lstar_medianerr = return_bootstrap_errors(dv400_imp50_det_lstar, reps)
    dv400_imp25_det_lstar_meanerr, dv400_imp25_det_lstar_medianerr = return_bootstrap_errors(dv400_imp25_det_lstar, reps)

    dv400_imp1000_non_lstar_meanerr, dv400_imp1000_non_lstar_medianerr = return_bootstrap_errors(dv400_imp1000_non_lstar, reps)
    dv400_imp750_non_lstar_meanerr, dv400_imp750_non_lstar_medianerr = return_bootstrap_errors(dv400_imp750_non_lstar, reps)
    dv400_imp500_non_lstar_meanerr, dv400_imp500_non_lstar_medianerr = return_bootstrap_errors(dv400_imp500_non_lstar, reps)
    dv400_imp400_non_lstar_meanerr, dv400_imp400_non_lstar_medianerr = return_bootstrap_errors(dv400_imp400_non_lstar, reps)
    dv400_imp300_non_lstar_meanerr, dv400_imp300_non_lstar_medianerr = return_bootstrap_errors(dv400_imp300_non_lstar, reps)
    dv400_imp200_non_lstar_meanerr, dv400_imp200_non_lstar_medianerr = return_bootstrap_errors(dv400_imp200_non_lstar, reps)
    dv400_imp100_non_lstar_meanerr, dv400_imp100_non_lstar_medianerr = return_bootstrap_errors(dv400_imp100_non_lstar, reps)
    dv400_imp50_non_lstar_meanerr, dv400_imp50_non_lstar_medianerr = return_bootstrap_errors(dv400_imp50_non_lstar, reps)
    dv400_imp25_non_lstar_meanerr, dv400_imp25_non_lstar_medianerr = return_bootstrap_errors(dv400_imp25_non_lstar, reps)

    # impact - inclinations
    dv400_imp1000_det_inc_meanerr, dv400_imp1000_det_inc_medianerr = return_bootstrap_errors(dv400_imp1000_det_inc, reps)
    dv400_imp750_det_inc_meanerr, dv400_imp750_det_inc_medianerr = return_bootstrap_errors(dv400_imp750_det_inc, reps)
    dv400_imp500_det_inc_meanerr, dv400_imp500_det_inc_medianerr = return_bootstrap_errors(dv400_imp500_det_inc, reps)
    dv400_imp400_det_inc_meanerr, dv400_imp400_det_inc_medianerr = return_bootstrap_errors(dv400_imp400_det_inc, reps)
    dv400_imp300_det_inc_meanerr, dv400_imp300_det_inc_medianerr = return_bootstrap_errors(dv400_imp300_det_inc, reps)
    dv400_imp200_det_inc_meanerr, dv400_imp200_det_inc_medianerr = return_bootstrap_errors(dv400_imp200_det_inc, reps)
    dv400_imp100_det_inc_meanerr, dv400_imp100_det_inc_medianerr = return_bootstrap_errors(dv400_imp100_det_inc, reps)
    dv400_imp50_det_inc_meanerr, dv400_imp50_det_inc_medianerr = return_bootstrap_errors(dv400_imp50_det_inc, reps)
    dv400_imp25_det_inc_meanerr, dv400_imp25_det_inc_medianerr = return_bootstrap_errors(dv400_imp25_det_inc, reps)

    dv400_imp1000_non_inc_meanerr, dv400_imp1000_non_inc_medianerr = return_bootstrap_errors(dv400_imp1000_non_inc, reps)
    dv400_imp750_non_inc_meanerr, dv400_imp750_non_inc_medianerr = return_bootstrap_errors(dv400_imp750_non_inc, reps)
    dv400_imp500_non_inc_meanerr, dv400_imp500_non_inc_medianerr = return_bootstrap_errors(dv400_imp500_non_inc, reps)
    dv400_imp400_non_inc_meanerr, dv400_imp400_non_inc_medianerr = return_bootstrap_errors(dv400_imp400_non_inc, reps)
    dv400_imp300_non_inc_meanerr, dv400_imp300_non_inc_medianerr = return_bootstrap_errors(dv400_imp300_non_inc, reps)
    dv400_imp200_non_inc_meanerr, dv400_imp200_non_inc_medianerr = return_bootstrap_errors(dv400_imp200_non_inc, reps)
    dv400_imp100_non_inc_meanerr, dv400_imp100_non_inc_medianerr = return_bootstrap_errors(dv400_imp100_non_inc, reps)
    dv400_imp50_non_inc_meanerr, dv400_imp50_non_inc_medianerr = return_bootstrap_errors(dv400_imp50_non_inc, reps)
    dv400_imp25_non_inc_meanerr, dv400_imp25_non_inc_medianerr = return_bootstrap_errors(dv400_imp25_non_inc, reps)
    



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
            

        x = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 0.75]
        y_det_mean = [bmean(dv400_l0001_det_inc),
                    bmean(dv400_l0005_det_inc),
                    bmean(dv400_l001_det_inc),
                    bmean(dv400_l005_det_inc),
                    bmean(dv400_l01_det_inc),
                    bmean(dv400_l05_det_inc),
                    bmean(dv400_l1_det_inc),
                    bmean(dv400_l5_det_inc),
                    bmean(dv400_l75_det_inc)]
                    
        # bootstrap errors in MEAN detection inclination
        y_det_err_mean = [dv400_l0001_det_inc_meanerr,
                        dv400_l0005_det_inc_meanerr,
                        dv400_l001_det_inc_meanerr,
                        dv400_l005_det_inc_meanerr,
                        dv400_l01_det_inc_meanerr,
                        dv400_l05_det_inc_meanerr,
                        dv400_l1_det_inc_meanerr,
                        dv400_l5_det_inc_meanerr,
                        dv400_l75_det_inc_meanerr]
            
        
        y_non_mean = [bmean(dv400_l0001_non_inc),
                    bmean(dv400_l0005_non_inc),
                    bmean(dv400_l001_non_inc),
                    bmean(dv400_l005_non_inc), 
                    bmean(dv400_l01_non_inc),
                    bmean(dv400_l05_non_inc),
                    bmean(dv400_l1_non_inc),
                    bmean(dv400_l5_non_inc),
                    bmean(dv400_l75_non_inc)]
        
        # bootstrap errors in MEAN non-detection inclination
        y_non_err_mean = [dv400_l0001_non_inc_meanerr,
                        dv400_l0005_non_inc_meanerr,
                        dv400_l001_non_inc_meanerr,
                        dv400_l005_non_inc_meanerr,
                        dv400_l01_non_inc_meanerr,
                        dv400_l05_non_inc_meanerr,
                        dv400_l1_non_inc_meanerr,
                        dv400_l5_non_inc_meanerr,
                        dv400_l75_non_inc_meanerr]
                    
                
        y_det_median = [bmedian(dv400_l0001_det_inc),
                    bmedian(dv400_l0005_det_inc),
                    bmedian(dv400_l001_det_inc),
                    bmedian(dv400_l005_det_inc), 
                    bmedian(dv400_l01_det_inc),
                    bmedian(dv400_l05_det_inc),
                    bmedian(dv400_l1_det_inc),
                    bmedian(dv400_l5_det_inc),
                    bmedian(dv400_l75_det_inc)]
                    
        # bootstrap errors in MEDIAN detection inclination
        y_det_err_median = [dv400_l0001_det_inc_medianerr,
                        dv400_l0005_det_inc_medianerr,
                        dv400_l001_det_inc_medianerr,
                        dv400_l005_det_inc_medianerr,
                        dv400_l01_det_inc_medianerr,
                        dv400_l05_det_inc_medianerr,
                        dv400_l1_det_inc_medianerr,
                        dv400_l5_det_inc_medianerr,
                        dv400_l75_det_inc_medianerr]
                    
            
        y_non_median = [bmedian(dv400_l0001_non_inc),
                    bmedian(dv400_l0005_non_inc),
                    bmedian(dv400_l001_non_inc),
                    bmedian(dv400_l005_non_inc), 
                    bmedian(dv400_l01_non_inc),
                    bmedian(dv400_l05_non_inc),
                    bmedian(dv400_l1_non_inc),
                    bmedian(dv400_l5_non_inc),
                    bmedian(dv400_l75_non_inc)]

        # bootstrap errors in MEDIAN non-detection inclination
        y_det_err_median = [dv400_l0001_non_inc_medianerr,
                        dv400_l0005_non_inc_medianerr,
                        dv400_l001_non_inc_medianerr,
                        dv400_l005_non_inc_medianerr,
                        dv400_l01_non_inc_medianerr,
                        dv400_l05_non_inc_medianerr,
                        dv400_l1_non_inc_medianerr,
                        dv400_l5_non_inc_medianerr,
                        dv400_l75_non_inc_medianerr]
                    




        # MEAN inclination for detections with errors
        ax1.errorbar(x,
                    y_det_mean,
                    yerr=y_det_meanerr,
                    marker=symbol_det_mean,
                    c=color_det_mean,
                    ms=markerSize,
                    markeredgecolor='black',
                    lw = lw,
                    alpha=alpha_det_mean,
                    label=label_det_mean)



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
            savefig('{0}/detection_fraction_likelihood_inc_median_lstarcut{1}_2.pdf'.format(saveDirectory, lstar_cut),format='pdf',bbox_inches='tight')
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
            savefig('{0}/detection_fraction_impact_inc_median_lstarcut{1}_2.pdf'.format(saveDirectory, lstar_cut),format='pdf',bbox_inches='tight')
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
        

        alpha_likelihood = 0.9
        alpha_impact = 0.9
        markerSize = 12
        lw = 3.

        
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
        frac_imp100 = float(dv400_imp100_det) / float(dv400_imp100_det + dv400_imp100_non)
        frac_imp50  = float(dv400_imp50_det)  / float(dv400_imp50_det  + dv400_imp50_non)
        
        try:
            frac_imp25  = float(dv400_imp25_det)  / float(dv400_imp25_det  + dv400_imp25_non)
        except Exception,e:
            print 'Error: ',e
            frac_imp25 = 0

        impact_x = [25, 50, 100, 200, 300, 400, 500, 750, 1000]
        impact_y = [frac_imp25, frac_imp50, frac_imp100, frac_imp200, frac_imp300, frac_imp400, frac_imp500, frac_imp750, frac_imp1000]


        frac_l0001 = float(dv400_l0001_det) / float(dv400_l0001_det + dv400_l0001_non)
        frac_l0005 = float(dv400_l0005_det) / float(dv400_l0005_det + dv400_l0005_non)
        frac_l001 = float(dv400_l001_det) / float(dv400_l001_det + dv400_l001_non)
        frac_l005 = float(dv400_l005_det) / float(dv400_l005_det + dv400_l005_non)
        frac_l01  = float(dv400_l01_det)  / float(dv400_l01_det  + dv400_l01_non)
        frac_l05  = float(dv400_l05_det)  / float(dv400_l05_det  + dv400_l05_non)
        frac_l1   = float(dv400_l1_det)   / float(dv400_l1_det   + dv400_l1_non)
        frac_l5   = float(dv400_l5_det)   / float(dv400_l5_det   + dv400_l5_non)
        frac_l75   = float(dv400_l75_det)   / float(dv400_l75_det   + dv400_l75_non)

        likelihood_x = [0.75, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001]
        likelihood_y = [frac_l75, frac_l5, frac_l1, frac_l05, frac_l01, frac_l005, frac_l001, frac_l0005, frac_l0001]


        # impact detection fraction
        plot1 = ax1.plot(impact_x,
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
        plot2 = ax2.plot(likelihood_x,
                likelihood_y,
                marker=symbol_likelihood,
                c=color_likelihood,
                ms=markerSize,
                markeredgecolor='black',
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
        plots = plot1 + plot2
        labs = [p.get_label() for p in plots]
        ax1.legend(plots, labs, prop={'size':12},loc='lower left',fancybox=True)


#         legend(scatterpoints=1,prop={'size':12},loc='lower left',fancybox=True)
#         leg.get_frame().set_alpha(0.5)

#         ax1.grid(b=None,which='major',axis='both')
        ylim(0.4, 1.)
#         xlim(0, 2.5)

        if plot_detection_fraction_both_save:
            savefig('{0}/detection_fraction_both_lstarcut_{1}_new2.pdf'.format(saveDirectory, lstar_cut),format='pdf',bbox_inches='tight')
        else:
            show()


##########################################################################################
##########################################################################################
##########################################################################################




##########################################################################################
##########################################################################################
##########################################################################################

    if plot_detection_fraction_likelihood_lstar:
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
        
        print 'len(dv400_l5_det_lstar) : ',len(dv400_l1_det_lstar)
        print 'len(dv400_l5_det_lstar) : ',len(dv400_l5_det_lstar)
        print 'len(dv400_l75_det_lstar) : ',len(dv400_l75_det_lstar)
        print
        print
        print 'len(dv400_l5_non_lstar) : ',len(dv400_l1_non_lstar)
        print 'len(dv400_l5_non_lstar) : ',len(dv400_l5_non_lstar)
        print 'len(dv400_l75_non_lstar) : ',len(dv400_l75_non_lstar)
        print
        print
            

        x = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 0.75]
        if use_mean:
            y_det = [bmean(dv400_l0001_det_lstar),
                    bmean(dv400_l0005_det_lstar),
                    bmean(dv400_l001_det_lstar),
                    bmean(dv400_l005_det_lstar),
                    bmean(dv400_l01_det_lstar),
                    bmean(dv400_l05_det_lstar),
                    bmean(dv400_l1_det_lstar),
                    bmean(dv400_l5_det_lstar),
                    bmean(dv400_l75_det_lstar)]
                
                
            y_non = [bmean(dv400_l0001_non_lstar),
                    bmean(dv400_l0005_non_lstar),
                    bmean(dv400_l001_non_lstar),
                    bmean(dv400_l005_non_lstar),
                    bmean(dv400_l01_non_lstar),
                    bmean(dv400_l05_non_lstar),
                    bmean(dv400_l1_non_lstar),
                    bmean(dv400_l5_non_lstar),
                    bmean(dv400_l75_non_lstar)]
                    
        else:
            y_det = [bmedian(dv400_l0001_det_lstar),
                    bmedian(dv400_l0005_det_lstar),
                    bmedian(dv400_l001_det_lstar),
                    bmedian(dv400_l005_det_lstar),
                    bmedian(dv400_l01_det_lstar),
                    bmedian(dv400_l05_det_lstar),
                    bmedian(dv400_l1_det_lstar),
                    bmedian(dv400_l5_det_lstar),
                    bmedian(dv400_l75_det_lstar)]
                    
            y_det_err = [std(dv400_l0001_det_lstar),
                        std(dv400_l0005_det_lstar),
                        std(dv400_l001_det_lstar),
                        std(dv400_l005_det_lstar),
                        std(dv400_l01_det_lstar),
                        std(dv400_l05_det_lstar),
                        std(dv400_l1_det_lstar),
                        std(dv400_l5_det_lstar),
                        std(dv400_l75_det_lstar)]
                        
            y_non = [bmedian(dv400_l0001_non_lstar),
                    bmedian(dv400_l0005_non_lstar),
                    bmedian(dv400_l001_non_lstar),
                    bmedian(dv400_l005_non_lstar),
                    bmedian(dv400_l01_non_lstar),
                    bmedian(dv400_l05_non_lstar),
                    bmedian(dv400_l1_non_lstar),
                    bmedian(dv400_l5_non_lstar),
                    bmedian(dv400_l75_non_lstar)]

            y_non_err = [std(dv400_l0001_non_lstar),
                        std(dv400_l0005_non_lstar),
                        std(dv400_l001_non_lstar),
                        std(dv400_l005_non_lstar),
                        std(dv400_l01_non_lstar),
                        std(dv400_l05_non_lstar),
                        std(dv400_l1_non_lstar),
                        std(dv400_l5_non_lstar),
                        std(dv400_l75_non_lstar)]


        print 'y_det: ',y_det
        print
        print 'y_non: ',y_non
        print
        print

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

#         ax1.errorbar(x,
#                     y_det,
#                     yerr=y_det_err,
#                     marker=symbol_det,
#                     c=color_det,
#                     ms=markerSize,
#                     markeredgecolor='black',
#                     lw = lw,
#                     alpha=alpha_det,
#                     label=label_det)

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

#         ax1.errorbar(x,
#                     y_non,
#                     yerr=y_non_err,
#                     marker=symbol_non,
#                     c=color_non,
#                     ms=markerSize,
#                     markeredgecolor='black',
#                     lw = lw,
#                     alpha=alpha_non,
#                     label=label_non)


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
        
        ax1.set_ylabel(r'$\rm L^{\**}$')
        
        leg = ax1.legend(scatterpoints=1,prop={'size':12},loc='upper right',fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax1.grid(b=None,which='major',axis='both')
#         ylim(0., .)
        xlim(0.0001, 1.)
        ax1.invert_xaxis()


        if plot_detection_fraction_likelihood_lstar_save:
            savefig('{0}/detection_fraction_likelihood_lstar_median_lstarcut{1}_2.pdf'.format(saveDirectory, lstar_cut),format='pdf',bbox_inches='tight')
        else:
            show()

##########################################################################################
##########################################################################################


##########################################################################################
##########################################################################################
    
    if plot_detection_fraction_impact_lstar:
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
        
        print 'len(dv400_imp100_det_lstar) : ',len(dv400_imp100_det_lstar)
        print 'len(dv400_imp50_det_lstar) : ',len(dv400_imp50_det_lstar)
        print 'len(dv400_imp25_det_lstar) : ',len(dv400_imp25_det_lstar)
        print
        print
        print 'len(dv400_imp100_non_lstar) : ',len(dv400_imp100_non_lstar)
        print 'len(dv400_imp50_non_lstar) : ',len(dv400_imp50_non_lstar)
        print 'len(dv400_imp25_non_lstar) : ',len(dv400_imp25_non_lstar)
        print
        print 'bmean(dv400_imp1000_det_lstar): ',bmean(dv400_imp1000_det_lstar)
        print 'bmean(dv400_imp1000_non_lstar): ',bmean(dv400_imp1000_non_lstar)
        print
        print 'bmean(dv400_imp25_det_lstar): ',bmean(dv400_imp25_det_lstar)
        print 'bmean(dv400_imp25_non_lstar): ',bmean(dv400_imp25_non_lstar)
        print
        print
        print 'dv400_imp25_det_lstar: ',dv400_imp25_det_lstar
        print 'dv400_imp50_det_lstar: ',dv400_imp50_det_lstar
        print
        print 'dv400_imp25_non_lstar: ',dv400_imp25_non_lstar
        print 'dv400_imp50_non_lstar: ',dv400_imp50_non_lstar
        print

        x = [1000, 750, 500, 400, 300, 200, 100, 50, 25]
        if use_mean:
            y_det = [bmean(dv400_imp1000_det_lstar),
                    bmean(dv400_imp750_det_lstar),
                    bmean(dv400_imp500_det_lstar),
                    bmean(dv400_imp400_det_lstar), 
                    bmean(dv400_imp300_det_lstar),
                    bmean(dv400_imp200_det_lstar),
                    bmean(dv400_imp100_det_lstar),
                    bmean(dv400_imp50_det_lstar),
                    bmean(dv400_imp25_det_lstar)]
                
                
            y_non = [bmean(dv400_imp1000_non_lstar),
                    bmean(dv400_imp750_non_lstar),
                    bmean(dv400_imp500_non_lstar),
                    bmean(dv400_imp400_non_lstar), 
                    bmean(dv400_imp300_non_lstar),
                    bmean(dv400_imp200_non_lstar),
                    bmean(dv400_imp100_non_lstar),
                    bmean(dv400_imp50_non_lstar),
                    bmean(dv400_imp25_non_lstar)]
                    
        else:
            y_det = [bmedian(dv400_imp1000_det_lstar),
                    bmedian(dv400_imp750_det_lstar),
                    bmedian(dv400_imp500_det_lstar),
                    bmedian(dv400_imp400_det_lstar), 
                    bmedian(dv400_imp300_det_lstar),
                    bmedian(dv400_imp200_det_lstar),
                    bmedian(dv400_imp100_det_lstar),
                    bmedian(dv400_imp50_det_lstar),
                    bmedian(dv400_imp25_det_lstar)]
                
                
            y_non = [bmedian(dv400_imp1000_non_lstar),
                    bmedian(dv400_imp750_non_lstar),
                    bmedian(dv400_imp500_non_lstar),
                    bmedian(dv400_imp400_non_lstar), 
                    bmedian(dv400_imp300_non_lstar),
                    bmedian(dv400_imp200_non_lstar),
                    bmedian(dv400_imp100_non_lstar),
                    bmedian(dv400_imp50_non_lstar),
                    bmedian(dv400_imp25_non_lstar)]

        print 'y_det: ',y_det
        print
        print 'y_non: ',y_non
        print
        print

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
        majorLocator   = MultipleLocator(0.01)
        majorFormatter = FormatStrFormatter(r'$\rm %.2f$')
        minorLocator   = MultipleLocator(0.005)
        ax1.yaxis.set_major_locator(majorLocator)
        ax1.yaxis.set_major_formatter(majorFormatter)
        ax1.yaxis.set_minor_locator(minorLocator)
        
        ax1.set_ylabel(r'$\rm L^{\**}$')

        leg = ax1.legend(scatterpoints=1,prop={'size':12},loc='upper right',fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax1.grid(b=None,which='major',axis='both')
        ylim(0., 0.075)
        xlim(0, 1000)

        if plot_detection_fraction_impact_lstar_save:
            savefig('{0}/detection_fraction_impact_lstar_median_lstarcut{1}_2.pdf'.format(saveDirectory, lstar_cut),format='pdf',bbox_inches='tight')
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
    