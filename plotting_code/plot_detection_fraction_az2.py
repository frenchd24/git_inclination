#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_detection_fraction_az.py, v 1.0 07/02/18

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
    for i in x_resample:
        med = np.median(i)
        resampled_medians.append(med)
    
    median_err = std(resampled_medians)
        

    return mean_err, median_err
    


def main():   
    # plot detection fraction as a function of azlination
    plot_detection_fraction_likelihood_az = False
    plot_detection_fraction_likelihood_az_save = False
    
    # plot detection fraction as a function of median azlination -> This one!
    plot_detection_fraction_like_az2 = False
    plot_detection_fraction_like_az_save2 = False

    # plot detection fraction as a function of azlination
    plot_detection_fraction_impact_az = False
    plot_detection_fraction_impact_az_save = False

    # plot detection fraction as a function of azlination for both impact and likelihood
    plot_detection_fraction_both_az = False
    plot_detection_fraction_both_az_save = False
    
    # plot detection fraction as a function of azlination for both impact and likelihood
    plot_detection_fraction_vs_az_like = True
    plot_detection_fraction_vs_az_like_save = True
    
    # plot detection fraction as a function of azlination for both impact and likelihood
    plot_detection_fraction_vs_az_imp = False
    plot_detection_fraction_vs_az_imp_save = False
    
    # plot_number = 1 for just the isolated sample, =2 adds the associated, =3 adds two+
    # =4 adds groups with 2 or more members
    Lstar_min = 0.5
    
    bootstrap = False
    
    # only consider galaxies within this velocity limit
    cut = 'minEW200'
    
    # which lstar cut subset to use?
    lstar_cut  = 'include3'

    # some colors
    color_blue = '#436bad'  # french blue
    color_red = '#ec2d01'   # tomato red

    if getpass.getuser() == 'frenchd':

#         gtPickleFilename = '/Users/frenchd/Research/GT_update2/pickleGT_filteredAll.p'

        saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/plotting_code/detection_fraction_figs/'
        
#         detection_fraction_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction3.p'
#         detection_fraction_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_lstarcut_{0}.p'.format(lstar_cut)
        detection_fraction_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_minEW0_closestonly2.p'
        
#         detection_fraction_vcut_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_vcut2500.p'
#         detection_fraction_vcut_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_lstarcut01_vcut2500.p'
#         detection_fraction_vcut_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_dmin75_vcut2500_minEW50.p'
        detection_fraction_min50_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_minEW50_closestonly2.p'
        detection_fraction_min100_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_minEW100_closestonly.p'
        detection_fraction_min200_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_minEW200_closestonly2.p'
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

    # impact parameter detection azlinations
    dv400_imp1000_det_az = np.array(full_dict['dv400_imp1000_det_az'])
    dv400_imp750_det_az = np.array(full_dict['dv400_imp750_det_az'])
    dv400_imp500_det_az = np.array(full_dict['dv400_imp500_det_az'])
    dv400_imp400_det_az = np.array(full_dict['dv400_imp400_det_az'])
    dv400_imp300_det_az = np.array(full_dict['dv400_imp300_det_az'])
    dv400_imp200_det_az = np.array(full_dict['dv400_imp200_det_az'])
    dv400_imp100_det_az = np.array(full_dict['dv400_imp100_det_az'])
    dv400_imp50_det_az = np.array(full_dict['dv400_imp50_det_az'])
    dv400_imp25_det_az = np.array(full_dict['dv400_imp25_det_az'])
    
    # now impact parameter non-detection azlinations 
    dv400_imp1000_non_az = np.array(full_dict['dv400_imp1000_non_az'])
    dv400_imp750_non_az = np.array(full_dict['dv400_imp750_non_az'])
    dv400_imp500_non_az = np.array(full_dict['dv400_imp500_non_az'])
    dv400_imp400_non_az = np.array(full_dict['dv400_imp400_non_az'])
    dv400_imp300_non_az = np.array(full_dict['dv400_imp300_non_az'])
    dv400_imp200_non_az = np.array(full_dict['dv400_imp200_non_az'])
    dv400_imp100_non_az = np.array(full_dict['dv400_imp100_non_az'])
    dv400_imp50_non_az = np.array(full_dict['dv400_imp50_non_az'])
    dv400_imp25_non_az = np.array(full_dict['dv400_imp25_non_az'])

    # minEW = 50
    # impact parameter detection azlinations
    dv400_imp1000_min50_det_az = np.array(full_dict_min50['dv400_imp1000_det_az'])
    dv400_imp750_min50_det_az = np.array(full_dict_min50['dv400_imp750_det_az'])
    dv400_imp500_min50_det_az = np.array(full_dict_min50['dv400_imp500_det_az'])
    dv400_imp400_min50_det_az = np.array(full_dict_min50['dv400_imp400_det_az'])
    dv400_imp300_min50_det_az = np.array(full_dict_min50['dv400_imp300_det_az'])
    dv400_imp200_min50_det_az = np.array(full_dict_min50['dv400_imp200_det_az'])
    dv400_imp100_min50_det_az = np.array(full_dict_min50['dv400_imp100_det_az'])
    dv400_imp50_min50_det_az = np.array(full_dict_min50['dv400_imp50_det_az'])
    dv400_imp25_min50_det_az = np.array(full_dict_min50['dv400_imp25_det_az'])
    
    # now impact parameter non-detection azlinations 
    dv400_imp1000_min50_non_az = np.array(full_dict_min50['dv400_imp1000_non_az'])
    dv400_imp750_min50_non_az = np.array(full_dict_min50['dv400_imp750_non_az'])
    dv400_imp500_min50_non_az = np.array(full_dict_min50['dv400_imp500_non_az'])
    dv400_imp400_min50_non_az = np.array(full_dict_min50['dv400_imp400_non_az'])
    dv400_imp300_min50_non_az = np.array(full_dict_min50['dv400_imp300_non_az'])
    dv400_imp200_min50_non_az = np.array(full_dict_min50['dv400_imp200_non_az'])
    dv400_imp100_min50_non_az = np.array(full_dict_min50['dv400_imp100_non_az'])
    dv400_imp50_min50_non_az = np.array(full_dict_min50['dv400_imp50_non_az'])
    dv400_imp25_min50_non_az = np.array(full_dict_min50['dv400_imp25_non_az'])

    # minEW = 100
    # impact parameter detection azlinations
#     dv400_imp1000_min100_det_az = np.array(full_dict_min100['dv400_imp1000_det_az'])
#     dv400_imp750_min100_det_az = np.array(full_dict_min100['dv400_imp750_det_az'])
#     dv400_imp500_min100_det_az = np.array(full_dict_min100['dv400_imp500_det_az'])
#     dv400_imp400_min100_det_az = np.array(full_dict_min100['dv400_imp400_det_az'])
#     dv400_imp300_min100_det_az = np.array(full_dict_min100['dv400_imp300_det_az'])
#     dv400_imp200_min100_det_az = np.array(full_dict_min100['dv400_imp200_det_az'])
#     dv400_imp100_min100_det_az = np.array(full_dict_min100['dv400_imp100_det_az'])
#     dv400_imp50_min100_det_az = np.array(full_dict_min100['dv400_imp50_det_az'])
#     dv400_imp25_min100_det_az = np.array(full_dict_min100['dv400_imp25_det_az'])
#     
#     # now impact parameter non-detection azlinations 
#     dv400_imp1000_min100_non_az = np.array(full_dict_min100['dv400_imp1000_non_az'])
#     dv400_imp750_min100_non_az = np.array(full_dict_min100['dv400_imp750_non_az'])
#     dv400_imp500_min100_non_az = np.array(full_dict_min100['dv400_imp500_non_az'])
#     dv400_imp400_min100_non_az = np.array(full_dict_min100['dv400_imp400_non_az'])
#     dv400_imp300_min100_non_az = np.array(full_dict_min100['dv400_imp300_non_az'])
#     dv400_imp200_min100_non_az = np.array(full_dict_min100['dv400_imp200_non_az'])
#     dv400_imp100_min100_non_az = np.array(full_dict_min100['dv400_imp100_non_az'])
#     dv400_imp50_min100_non_az = np.array(full_dict_min100['dv400_imp50_non_az'])
#     dv400_imp25_min100_non_az = np.array(full_dict_min100['dv400_imp25_non_az'])


    # minEW = 200
    # impact parameter detection azlinations
    dv400_imp1000_min200_det_az = np.array(full_dict_min200['dv400_imp1000_det_az'])
    dv400_imp750_min200_det_az = np.array(full_dict_min200['dv400_imp750_det_az'])
    dv400_imp500_min200_det_az = np.array(full_dict_min200['dv400_imp500_det_az'])
    dv400_imp400_min200_det_az = np.array(full_dict_min200['dv400_imp400_det_az'])
    dv400_imp300_min200_det_az = np.array(full_dict_min200['dv400_imp300_det_az'])
    dv400_imp200_min200_det_az = np.array(full_dict_min200['dv400_imp200_det_az'])
    dv400_imp100_min200_det_az = np.array(full_dict_min200['dv400_imp100_det_az'])
    dv400_imp50_min200_det_az = np.array(full_dict_min200['dv400_imp50_det_az'])
    dv400_imp25_min200_det_az = np.array(full_dict_min200['dv400_imp25_det_az'])
    
    # now impact parameter non-detection azlinations 
    dv400_imp1000_min200_non_az = np.array(full_dict_min200['dv400_imp1000_non_az'])
    dv400_imp750_min200_non_az = np.array(full_dict_min200['dv400_imp750_non_az'])
    dv400_imp500_min200_non_az = np.array(full_dict_min200['dv400_imp500_non_az'])
    dv400_imp400_min200_non_az = np.array(full_dict_min200['dv400_imp400_non_az'])
    dv400_imp300_min200_non_az = np.array(full_dict_min200['dv400_imp300_non_az'])
    dv400_imp200_min200_non_az = np.array(full_dict_min200['dv400_imp200_non_az'])
    dv400_imp100_min200_non_az = np.array(full_dict_min200['dv400_imp100_non_az'])
    dv400_imp50_min200_non_az = np.array(full_dict_min200['dv400_imp50_non_az'])
    dv400_imp25_min200_non_az = np.array(full_dict_min200['dv400_imp25_non_az'])


    # minEW = 300
    # impact parameter detection azlinations
#     dv400_imp1000_min300_det_az = np.array(full_dict_min300['dv400_imp1000_det_az'])
#     dv400_imp750_min300_det_az = np.array(full_dict_min300['dv400_imp750_det_az'])
#     dv400_imp500_min300_det_az = np.array(full_dict_min300['dv400_imp500_det_az'])
#     dv400_imp400_min300_det_az = np.array(full_dict_min300['dv400_imp400_det_az'])
#     dv400_imp300_min300_det_az = np.array(full_dict_min300['dv400_imp300_det_az'])
#     dv400_imp200_min300_det_az = np.array(full_dict_min300['dv400_imp200_det_az'])
#     dv400_imp100_min300_det_az = np.array(full_dict_min300['dv400_imp100_det_az'])
#     dv400_imp50_min300_det_az = np.array(full_dict_min300['dv400_imp50_det_az'])
#     dv400_imp25_min300_det_az = np.array(full_dict_min300['dv400_imp25_det_az'])
#     
#     # now impact parameter non-detection azlinations 
#     dv400_imp1000_min300_non_az = np.array(full_dict_min300['dv400_imp1000_non_az'])
#     dv400_imp750_min300_non_az = np.array(full_dict_min300['dv400_imp750_non_az'])
#     dv400_imp500_min300_non_az = np.array(full_dict_min300['dv400_imp500_non_az'])
#     dv400_imp400_min300_non_az = np.array(full_dict_min300['dv400_imp400_non_az'])
#     dv400_imp300_min300_non_az = np.array(full_dict_min300['dv400_imp300_non_az'])
#     dv400_imp200_min300_non_az = np.array(full_dict_min300['dv400_imp200_non_az'])
#     dv400_imp100_min300_non_az = np.array(full_dict_min300['dv400_imp100_non_az'])
#     dv400_imp50_min300_non_az = np.array(full_dict_min300['dv400_imp50_non_az'])
#     dv400_imp25_min300_non_az = np.array(full_dict_min300['dv400_imp25_non_az'])

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


    # now for likelihood detections azlinations
    dv400_l0001_det_az = np.array(full_dict['dv400_l0001_det_az'])
    dv400_l0005_det_az = np.array(full_dict['dv400_l0005_det_az'])
    dv400_l001_det_az = np.array(full_dict['dv400_l001_det_az'])
    dv400_l005_det_az = np.array(full_dict['dv400_l005_det_az'])
    dv400_l01_det_az = np.array(full_dict['dv400_l01_det_az'])
    dv400_l05_det_az = np.array(full_dict['dv400_l05_det_az'])
    dv400_l1_det_az = np.array(full_dict['dv400_l1_det_az'])
    dv400_l5_det_az = np.array(full_dict['dv400_l5_det_az'])
    dv400_l75_det_az = np.array(full_dict['dv400_l75_det_az'])

    # now for likelihood non-detections azlinations
    dv400_l0001_non_az = np.array(full_dict['dv400_l0001_non_az'])
    dv400_l0005_non_az = np.array(full_dict['dv400_l0005_non_az'])
    dv400_l001_non_az = np.array(full_dict['dv400_l001_non_az'])
    dv400_l005_non_az = np.array(full_dict['dv400_l005_non_az'])
    dv400_l01_non_az = np.array(full_dict['dv400_l01_non_az'])
    dv400_l05_non_az = np.array(full_dict['dv400_l05_non_az'])
    dv400_l1_non_az = np.array(full_dict['dv400_l1_non_az'])
    dv400_l5_non_az = np.array(full_dict['dv400_l5_non_az'])
    dv400_l75_non_az = np.array(full_dict['dv400_l75_non_az'])
    
    # minEW = 50
    # now for likelihood detections azlinations
    dv400_l0001_min50_det_az = np.array(full_dict_min50['dv400_l0001_det_az'])
    dv400_l0005_min50_det_az = np.array(full_dict_min50['dv400_l0005_det_az'])
    dv400_l001_min50_det_az = np.array(full_dict_min50['dv400_l001_det_az'])
    dv400_l005_min50_det_az = np.array(full_dict_min50['dv400_l005_det_az'])
    dv400_l01_min50_det_az = np.array(full_dict_min50['dv400_l01_det_az'])
    dv400_l05_min50_det_az = np.array(full_dict_min50['dv400_l05_det_az'])
    dv400_l1_min50_det_az = np.array(full_dict_min50['dv400_l1_det_az'])
    dv400_l5_min50_det_az = np.array(full_dict_min50['dv400_l5_det_az'])
    dv400_l75_min50_det_az = np.array(full_dict_min50['dv400_l75_det_az'])

    # now for likelihood non-detections azlinations
    dv400_l0001_min50_non_az = np.array(full_dict_min50['dv400_l0001_non_az'])
    dv400_l0005_min50_non_az = np.array(full_dict_min50['dv400_l0005_non_az'])
    dv400_l001_min50_non_az = np.array(full_dict_min50['dv400_l001_non_az'])
    dv400_l005_min50_non_az = np.array(full_dict_min50['dv400_l005_non_az'])
    dv400_l01_min50_non_az = np.array(full_dict_min50['dv400_l01_non_az'])
    dv400_l05_min50_non_az = np.array(full_dict_min50['dv400_l05_non_az'])
    dv400_l1_min50_non_az = np.array(full_dict_min50['dv400_l1_non_az'])
    dv400_l5_min50_non_az = np.array(full_dict_min50['dv400_l5_non_az'])
    dv400_l75_min50_non_az = np.array(full_dict_min50['dv400_l75_non_az'])
    
    # minEW = 100
    # now for likelihood detections azlinations
#     dv400_l0001_min100_det_az = np.array(full_dict_min100['dv400_l0001_det_az'])
#     dv400_l0005_min100_det_az = np.array(full_dict_min100['dv400_l0005_det_az'])
#     dv400_l001_min100_det_az = np.array(full_dict_min100['dv400_l001_det_az'])
#     dv400_l005_min100_det_az = np.array(full_dict_min100['dv400_l005_det_az'])
#     dv400_l01_min100_det_az = np.array(full_dict_min100['dv400_l01_det_az'])
#     dv400_l05_min100_det_az = np.array(full_dict_min100['dv400_l05_det_az'])
#     dv400_l1_min100_det_az = np.array(full_dict_min100['dv400_l1_det_az'])
#     dv400_l5_min100_det_az = np.array(full_dict_min100['dv400_l5_det_az'])
#     dv400_l75_min100_det_az = np.array(full_dict_min100['dv400_l75_det_az'])
# 
#     # now for likelihood non-detections azlinations
#     dv400_l0001_min100_non_az = np.array(full_dict_min100['dv400_l0001_non_az'])
#     dv400_l0005_min100_non_az = np.array(full_dict_min100['dv400_l0005_non_az'])
#     dv400_l001_min100_non_az = np.array(full_dict_min100['dv400_l001_non_az'])
#     dv400_l005_min100_non_az = np.array(full_dict_min100['dv400_l005_non_az'])
#     dv400_l01_min100_non_az = np.array(full_dict_min100['dv400_l01_non_az'])
#     dv400_l05_min100_non_az = np.array(full_dict_min100['dv400_l05_non_az'])
#     dv400_l1_min100_non_az = np.array(full_dict_min100['dv400_l1_non_az'])
#     dv400_l5_min100_non_az = np.array(full_dict_min100['dv400_l5_non_az'])
#     dv400_l75_min100_non_az = np.array(full_dict_min100['dv400_l75_non_az'])

    # minEW = 200
    # now for likelihood detections azlinations
    dv400_l0001_min200_det_az = np.array(full_dict_min200['dv400_l0001_det_az'])
    dv400_l0005_min200_det_az = np.array(full_dict_min200['dv400_l0005_det_az'])
    dv400_l001_min200_det_az = np.array(full_dict_min200['dv400_l001_det_az'])
    dv400_l005_min200_det_az = np.array(full_dict_min200['dv400_l005_det_az'])
    dv400_l01_min200_det_az = np.array(full_dict_min200['dv400_l01_det_az'])
    dv400_l05_min200_det_az = np.array(full_dict_min200['dv400_l05_det_az'])
    dv400_l1_min200_det_az = np.array(full_dict_min200['dv400_l1_det_az'])
    dv400_l5_min200_det_az = np.array(full_dict_min200['dv400_l5_det_az'])
    dv400_l75_min200_det_az = np.array(full_dict_min200['dv400_l75_det_az'])

    # now for likelihood non-detections azlinations
    dv400_l0001_min200_non_az = np.array(full_dict_min200['dv400_l0001_non_az'])
    dv400_l0005_min200_non_az = np.array(full_dict_min200['dv400_l0005_non_az'])
    dv400_l001_min200_non_az = np.array(full_dict_min200['dv400_l001_non_az'])
    dv400_l005_min200_non_az = np.array(full_dict_min200['dv400_l005_non_az'])
    dv400_l01_min200_non_az = np.array(full_dict_min200['dv400_l01_non_az'])
    dv400_l05_min200_non_az = np.array(full_dict_min200['dv400_l05_non_az'])
    dv400_l1_min200_non_az = np.array(full_dict_min200['dv400_l1_non_az'])
    dv400_l5_min200_non_az = np.array(full_dict_min200['dv400_l5_non_az'])
    dv400_l75_min200_non_az = np.array(full_dict_min200['dv400_l75_non_az'])

    # minEW = 300
    # now for likelihood detections azlinations
#     dv400_l0001_min300_det_az = np.array(full_dict_min300['dv400_l0001_det_az'])
#     dv400_l0005_min300_det_az = np.array(full_dict_min300['dv400_l0005_det_az'])
#     dv400_l001_min300_det_az = np.array(full_dict_min300['dv400_l001_det_az'])
#     dv400_l005_min300_det_az = np.array(full_dict_min300['dv400_l005_det_az'])
#     dv400_l01_min300_det_az = np.array(full_dict_min300['dv400_l01_det_az'])
#     dv400_l05_min300_det_az = np.array(full_dict_min300['dv400_l05_det_az'])
#     dv400_l1_min300_det_az = np.array(full_dict_min300['dv400_l1_det_az'])
#     dv400_l5_min300_det_az = np.array(full_dict_min300['dv400_l5_det_az'])
#     dv400_l75_min300_det_az = np.array(full_dict_min300['dv400_l75_det_az'])
# 
#     # now for likelihood non-detections azlinations
#     dv400_l0001_min300_non_az = np.array(full_dict_min300['dv400_l0001_non_az'])
#     dv400_l0005_min300_non_az = np.array(full_dict_min300['dv400_l0005_non_az'])
#     dv400_l001_min300_non_az = np.array(full_dict_min300['dv400_l001_non_az'])
#     dv400_l005_min300_non_az = np.array(full_dict_min300['dv400_l005_non_az'])
#     dv400_l01_min300_non_az = np.array(full_dict_min300['dv400_l01_non_az'])
#     dv400_l05_min300_non_az = np.array(full_dict_min300['dv400_l05_non_az'])
#     dv400_l1_min300_non_az = np.array(full_dict_min300['dv400_l1_non_az'])
#     dv400_l5_min300_non_az = np.array(full_dict_min300['dv400_l5_non_az'])
#     dv400_l75_min300_non_az = np.array(full_dict_min300['dv400_l75_non_az'])
    
    # turn everything into floats
    dv400_imp1000_det_az = dv400_imp1000_det_az.astype(np.float)
    dv400_imp750_det_az = dv400_imp750_det_az.astype(np.float)
    dv400_imp500_det_az = dv400_imp500_det_az.astype(np.float)
    dv400_imp400_det_az = dv400_imp400_det_az.astype(np.float)
    dv400_imp300_det_az = dv400_imp300_det_az.astype(np.float)
    dv400_imp200_det_az = dv400_imp200_det_az.astype(np.float)
    dv400_imp100_det_az = dv400_imp100_det_az.astype(np.float)
    dv400_imp50_det_az = dv400_imp50_det_az.astype(np.float)
    dv400_imp25_det_az = dv400_imp25_det_az.astype(np.float)

    dv400_imp1000_non_az = dv400_imp1000_non_az.astype(np.float)
    dv400_imp750_non_az = dv400_imp750_non_az.astype(np.float)
    dv400_imp500_non_az = dv400_imp500_non_az.astype(np.float)
    dv400_imp400_non_az = dv400_imp400_non_az.astype(np.float)
    dv400_imp300_non_az = dv400_imp300_non_az.astype(np.float)
    dv400_imp200_non_az = dv400_imp200_non_az.astype(np.float)
    dv400_imp100_non_az = dv400_imp100_non_az.astype(np.float)
    dv400_imp50_non_az = dv400_imp50_non_az.astype(np.float)
    dv400_imp25_non_az = dv400_imp25_non_az.astype(np.float)
    
    dv400_l0001_det_az = dv400_l0001_det_az.astype(np.float)
    dv400_l0005_det_az = dv400_l0005_det_az.astype(np.float)
    dv400_l001_det_az = dv400_l001_det_az.astype(np.float)
    dv400_l005_det_az = dv400_l005_det_az.astype(np.float)
    dv400_l01_det_az = dv400_l01_det_az.astype(np.float)
    dv400_l05_det_az = dv400_l05_det_az.astype(np.float)
    dv400_l1_det_az = dv400_l1_det_az.astype(np.float)
    dv400_l5_det_az = dv400_l5_det_az.astype(np.float)
    dv400_l75_det_az = dv400_l75_det_az.astype(np.float)
    
    dv400_l0001_non_az = dv400_l0001_non_az.astype(np.float)
    dv400_l0005_non_az = dv400_l0005_non_az.astype(np.float)
    dv400_l001_non_az = dv400_l001_non_az.astype(np.float)
    dv400_l005_non_az = dv400_l005_non_az.astype(np.float)
    dv400_l01_non_az = dv400_l01_non_az.astype(np.float)
    dv400_l05_non_az = dv400_l05_non_az.astype(np.float)
    dv400_l1_non_az = dv400_l1_non_az.astype(np.float)
    dv400_l5_non_az = dv400_l5_non_az.astype(np.float)
    dv400_l75_non_az = dv400_l75_non_az.astype(np.float)
    
    # minEW = 50; turn everything in to floats
    dv400_imp1000_min50_det_az = dv400_imp1000_min50_det_az.astype(np.float)
    dv400_imp750_min50_det_az = dv400_imp750_min50_det_az.astype(np.float)
    dv400_imp500_min50_det_az = dv400_imp500_min50_det_az.astype(np.float)
    dv400_imp400_min50_det_az = dv400_imp400_min50_det_az.astype(np.float)
    dv400_imp300_min50_det_az = dv400_imp300_min50_det_az.astype(np.float)
    dv400_imp200_min50_det_az = dv400_imp200_min50_det_az.astype(np.float)
    dv400_imp100_min50_det_az = dv400_imp100_min50_det_az.astype(np.float)
    dv400_imp50_min50_det_az = dv400_imp50_min50_det_az.astype(np.float)
    dv400_imp25_min50_det_az = dv400_imp25_min50_det_az.astype(np.float)

    dv400_imp1000_min50_non_az = dv400_imp1000_min50_non_az.astype(np.float)
    dv400_imp750_min50_non_az = dv400_imp750_min50_non_az.astype(np.float)
    dv400_imp500_min50_non_az = dv400_imp500_min50_non_az.astype(np.float)
    dv400_imp400_min50_non_az = dv400_imp400_min50_non_az.astype(np.float)
    dv400_imp300_min50_non_az = dv400_imp300_min50_non_az.astype(np.float)
    dv400_imp200_min50_non_az = dv400_imp200_min50_non_az.astype(np.float)
    dv400_imp100_min50_non_az = dv400_imp100_min50_non_az.astype(np.float)
    dv400_imp50_min50_non_az = dv400_imp50_min50_non_az.astype(np.float)
    dv400_imp25_min50_non_az = dv400_imp25_min50_non_az.astype(np.float)
    
    dv400_l0001_min50_det_az = dv400_l0001_min50_det_az.astype(np.float)
    dv400_l0005_min50_det_az = dv400_l0005_min50_det_az.astype(np.float)
    dv400_l001_min50_det_az = dv400_l001_min50_det_az.astype(np.float)
    dv400_l005_min50_det_az = dv400_l005_min50_det_az.astype(np.float)
    dv400_l01_min50_det_az = dv400_l01_min50_det_az.astype(np.float)
    dv400_l05_min50_det_az = dv400_l05_min50_det_az.astype(np.float)
    dv400_l1_min50_det_az = dv400_l1_min50_det_az.astype(np.float)
    dv400_l5_min50_det_az = dv400_l5_min50_det_az.astype(np.float)
    dv400_l75_min50_det_az = dv400_l75_min50_det_az.astype(np.float)
    
    dv400_l0001_min50_non_az = dv400_l0001_min50_non_az.astype(np.float)
    dv400_l0005_min50_non_az = dv400_l0005_min50_non_az.astype(np.float)
    dv400_l001_min50_non_az = dv400_l001_min50_non_az.astype(np.float)
    dv400_l005_min50_non_az = dv400_l005_min50_non_az.astype(np.float)
    dv400_l01_min50_non_az = dv400_l01_min50_non_az.astype(np.float)
    dv400_l05_min50_non_az = dv400_l05_min50_non_az.astype(np.float)
    dv400_l1_min50_non_az = dv400_l1_min50_non_az.astype(np.float)
    dv400_l5_min50_non_az = dv400_l5_min50_non_az.astype(np.float)
    dv400_l75_min50_non_az = dv400_l75_min50_non_az.astype(np.float)
    
    # minEW = 100; turn everything in to floats
#     dv400_imp1000_min100_det_az = dv400_imp1000_min100_det_az.astype(np.float)
#     dv400_imp750_min100_det_az = dv400_imp750_min100_det_az.astype(np.float)
#     dv400_imp500_min100_det_az = dv400_imp500_min100_det_az.astype(np.float)
#     dv400_imp400_min100_det_az = dv400_imp400_min100_det_az.astype(np.float)
#     dv400_imp300_min100_det_az = dv400_imp300_min100_det_az.astype(np.float)
#     dv400_imp200_min100_det_az = dv400_imp200_min100_det_az.astype(np.float)
#     dv400_imp100_min100_det_az = dv400_imp100_min100_det_az.astype(np.float)
#     dv400_imp50_min100_det_az = dv400_imp50_min100_det_az.astype(np.float)
#     dv400_imp25_min100_det_az = dv400_imp25_min100_det_az.astype(np.float)
# 
#     dv400_imp1000_min100_non_az = dv400_imp1000_min100_non_az.astype(np.float)
#     dv400_imp750_min100_non_az = dv400_imp750_min100_non_az.astype(np.float)
#     dv400_imp500_min100_non_az = dv400_imp500_min100_non_az.astype(np.float)
#     dv400_imp400_min100_non_az = dv400_imp400_min100_non_az.astype(np.float)
#     dv400_imp300_min100_non_az = dv400_imp300_min100_non_az.astype(np.float)
#     dv400_imp200_min100_non_az = dv400_imp200_min100_non_az.astype(np.float)
#     dv400_imp100_min100_non_az = dv400_imp100_min100_non_az.astype(np.float)
#     dv400_imp50_min100_non_az = dv400_imp50_min100_non_az.astype(np.float)
#     dv400_imp25_min100_non_az = dv400_imp25_min100_non_az.astype(np.float)
#     
#     dv400_l0001_min100_det_az = dv400_l0001_min100_det_az.astype(np.float)
#     dv400_l0005_min100_det_az = dv400_l0005_min100_det_az.astype(np.float)
#     dv400_l001_min100_det_az = dv400_l001_min100_det_az.astype(np.float)
#     dv400_l005_min100_det_az = dv400_l005_min100_det_az.astype(np.float)
#     dv400_l01_min100_det_az = dv400_l01_min100_det_az.astype(np.float)
#     dv400_l05_min100_det_az = dv400_l05_min100_det_az.astype(np.float)
#     dv400_l1_min100_det_az = dv400_l1_min100_det_az.astype(np.float)
#     dv400_l5_min100_det_az = dv400_l5_min100_det_az.astype(np.float)
#     dv400_l75_min100_det_az = dv400_l75_min100_det_az.astype(np.float)
#     
#     dv400_l0001_min100_non_az = dv400_l0001_min100_non_az.astype(np.float)
#     dv400_l0005_min100_non_az = dv400_l0005_min100_non_az.astype(np.float)
#     dv400_l001_min100_non_az = dv400_l001_min100_non_az.astype(np.float)
#     dv400_l005_min100_non_az = dv400_l005_min100_non_az.astype(np.float)
#     dv400_l01_min100_non_az = dv400_l01_min100_non_az.astype(np.float)
#     dv400_l05_min100_non_az = dv400_l05_min100_non_az.astype(np.float)
#     dv400_l1_min100_non_az = dv400_l1_min100_non_az.astype(np.float)
#     dv400_l5_min100_non_az = dv400_l5_min100_non_az.astype(np.float)
#     dv400_l75_min100_non_az = dv400_l75_min100_non_az.astype(np.float)
    
    # minEW = 200; turn everything in to floats
    dv400_imp1000_min200_det_az = dv400_imp1000_min200_det_az.astype(np.float)
    dv400_imp750_min200_det_az = dv400_imp750_min200_det_az.astype(np.float)
    dv400_imp500_min200_det_az = dv400_imp500_min200_det_az.astype(np.float)
    dv400_imp400_min200_det_az = dv400_imp400_min200_det_az.astype(np.float)
    dv400_imp300_min200_det_az = dv400_imp300_min200_det_az.astype(np.float)
    dv400_imp200_min200_det_az = dv400_imp200_min200_det_az.astype(np.float)
    dv400_imp100_min200_det_az = dv400_imp100_min200_det_az.astype(np.float)
    dv400_imp50_min200_det_az = dv400_imp50_min200_det_az.astype(np.float)
    dv400_imp25_min200_det_az = dv400_imp25_min200_det_az.astype(np.float)

    dv400_imp1000_min200_non_az = dv400_imp1000_min200_non_az.astype(np.float)
    dv400_imp750_min200_non_az = dv400_imp750_min200_non_az.astype(np.float)
    dv400_imp500_min200_non_az = dv400_imp500_min200_non_az.astype(np.float)
    dv400_imp400_min200_non_az = dv400_imp400_min200_non_az.astype(np.float)
    dv400_imp300_min200_non_az = dv400_imp300_min200_non_az.astype(np.float)
    dv400_imp200_min200_non_az = dv400_imp200_min200_non_az.astype(np.float)
    dv400_imp100_min200_non_az = dv400_imp100_min200_non_az.astype(np.float)
    dv400_imp50_min200_non_az = dv400_imp50_min200_non_az.astype(np.float)
    dv400_imp25_min200_non_az = dv400_imp25_min200_non_az.astype(np.float)
    
    dv400_l0001_min200_det_az = dv400_l0001_min200_det_az.astype(np.float)
    dv400_l0005_min200_det_az = dv400_l0005_min200_det_az.astype(np.float)
    dv400_l001_min200_det_az = dv400_l001_min200_det_az.astype(np.float)
    dv400_l005_min200_det_az = dv400_l005_min200_det_az.astype(np.float)
    dv400_l01_min200_det_az = dv400_l01_min200_det_az.astype(np.float)
    dv400_l05_min200_det_az = dv400_l05_min200_det_az.astype(np.float)
    dv400_l1_min200_det_az = dv400_l1_min200_det_az.astype(np.float)
    dv400_l5_min200_det_az = dv400_l5_min200_det_az.astype(np.float)
    dv400_l75_min200_det_az = dv400_l75_min200_det_az.astype(np.float)
    
    dv400_l0001_min200_non_az = dv400_l0001_min200_non_az.astype(np.float)
    dv400_l0005_min200_non_az = dv400_l0005_min200_non_az.astype(np.float)
    dv400_l001_min200_non_az = dv400_l001_min200_non_az.astype(np.float)
    dv400_l005_min200_non_az = dv400_l005_min200_non_az.astype(np.float)
    dv400_l01_min200_non_az = dv400_l01_min200_non_az.astype(np.float)
    dv400_l05_min200_non_az = dv400_l05_min200_non_az.astype(np.float)
    dv400_l1_min200_non_az = dv400_l1_min200_non_az.astype(np.float)
    dv400_l5_min200_non_az = dv400_l5_min200_non_az.astype(np.float)
    dv400_l75_min200_non_az = dv400_l75_min200_non_az.astype(np.float)

    # minEW = 300; turn everything in to floats
#     dv400_imp1000_min300_det_az = dv400_imp1000_min300_det_az.astype(np.float)
#     dv400_imp750_min300_det_az = dv400_imp750_min300_det_az.astype(np.float)
#     dv400_imp500_min300_det_az = dv400_imp500_min300_det_az.astype(np.float)
#     dv400_imp400_min300_det_az = dv400_imp400_min300_det_az.astype(np.float)
#     dv400_imp300_min300_det_az = dv400_imp300_min300_det_az.astype(np.float)
#     dv400_imp200_min300_det_az = dv400_imp200_min300_det_az.astype(np.float)
#     dv400_imp100_min300_det_az = dv400_imp100_min300_det_az.astype(np.float)
#     dv400_imp50_min300_det_az = dv400_imp50_min300_det_az.astype(np.float)
#     dv400_imp25_min300_det_az = dv400_imp25_min300_det_az.astype(np.float)
# 
#     dv400_imp1000_min300_non_az = dv400_imp1000_min300_non_az.astype(np.float)
#     dv400_imp750_min300_non_az = dv400_imp750_min300_non_az.astype(np.float)
#     dv400_imp500_min300_non_az = dv400_imp500_min300_non_az.astype(np.float)
#     dv400_imp400_min300_non_az = dv400_imp400_min300_non_az.astype(np.float)
#     dv400_imp300_min300_non_az = dv400_imp300_min300_non_az.astype(np.float)
#     dv400_imp200_min300_non_az = dv400_imp200_min300_non_az.astype(np.float)
#     dv400_imp100_min300_non_az = dv400_imp100_min300_non_az.astype(np.float)
#     dv400_imp50_min300_non_az = dv400_imp50_min300_non_az.astype(np.float)
#     dv400_imp25_min300_non_az = dv400_imp25_min300_non_az.astype(np.float)
#     
#     dv400_l0001_min300_det_az = dv400_l0001_min300_det_az.astype(np.float)
#     dv400_l0005_min300_det_az = dv400_l0005_min300_det_az.astype(np.float)
#     dv400_l001_min300_det_az = dv400_l001_min300_det_az.astype(np.float)
#     dv400_l005_min300_det_az = dv400_l005_min300_det_az.astype(np.float)
#     dv400_l01_min300_det_az = dv400_l01_min300_det_az.astype(np.float)
#     dv400_l05_min300_det_az = dv400_l05_min300_det_az.astype(np.float)
#     dv400_l1_min300_det_az = dv400_l1_min300_det_az.astype(np.float)
#     dv400_l5_min300_det_az = dv400_l5_min300_det_az.astype(np.float)
#     dv400_l75_min300_det_az = dv400_l75_min300_det_az.astype(np.float)
#     
#     dv400_l0001_min300_non_az = dv400_l0001_min300_non_az.astype(np.float)
#     dv400_l0005_min300_non_az = dv400_l0005_min300_non_az.astype(np.float)
#     dv400_l001_min300_non_az = dv400_l001_min300_non_az.astype(np.float)
#     dv400_l005_min300_non_az = dv400_l005_min300_non_az.astype(np.float)
#     dv400_l01_min300_non_az = dv400_l01_min300_non_az.astype(np.float)
#     dv400_l05_min300_non_az = dv400_l05_min300_non_az.astype(np.float)
#     dv400_l1_min300_non_az = dv400_l1_min300_non_az.astype(np.float)
#     dv400_l5_min300_non_az = dv400_l5_min300_non_az.astype(np.float)
#     dv400_l75_min300_non_az = dv400_l75_min300_non_az.astype(np.float)


    #######
    # cut out '-99' from all the data
    # likelihood version first - azlinations
    dv400_l0001_det_az = list(filter(lambda x: x!= -99., dv400_l0001_det_az))
    dv400_l0005_det_az = list(filter(lambda x: x!= -99., dv400_l0005_det_az))
    dv400_l001_det_az = list(filter(lambda x: x!= -99., dv400_l001_det_az))
    dv400_l005_det_az = list(filter(lambda x: x!= -99., dv400_l005_det_az))
    dv400_l01_det_az = list(filter(lambda x: x!= -99., dv400_l01_det_az))
    dv400_l05_det_az = list(filter(lambda x: x!= -99., dv400_l05_det_az))
    dv400_l1_det_az = list(filter(lambda x: x!= -99., dv400_l1_det_az))
    dv400_l5_det_az = list(filter(lambda x: x!= -99., dv400_l5_det_az))
    dv400_l75_det_az = list(filter(lambda x: x!= -99., dv400_l75_det_az))

    dv400_l0001_non_az = list(filter(lambda x: x!= -99., dv400_l0001_non_az))
    dv400_l0005_non_az = list(filter(lambda x: x!= -99., dv400_l0005_non_az))
    dv400_l001_non_az = list(filter(lambda x: x!= -99., dv400_l001_non_az))
    dv400_l005_non_az = list(filter(lambda x: x!= -99., dv400_l005_non_az))
    dv400_l01_non_az = list(filter(lambda x: x!= -99., dv400_l01_non_az))
    dv400_l05_non_az = list(filter(lambda x: x!= -99., dv400_l05_non_az))
    dv400_l1_non_az = list(filter(lambda x: x!= -99., dv400_l1_non_az))
    dv400_l5_non_az = list(filter(lambda x: x!= -99., dv400_l5_non_az))
    dv400_l75_non_az = list(filter(lambda x: x!= -99., dv400_l75_non_az))
    
    # now impact version - azlinations
    dv400_imp1000_det_az = list(filter(lambda x: x!= -99., dv400_imp1000_det_az))
    dv400_imp750_det_az = list(filter(lambda x: x!= -99., dv400_imp750_det_az))
    dv400_imp500_det_az = list(filter(lambda x: x!= -99., dv400_imp500_det_az))
    dv400_imp400_det_az = list(filter(lambda x: x!= -99., dv400_imp400_det_az))
    dv400_imp300_det_az = list(filter(lambda x: x!= -99., dv400_imp300_det_az))
    dv400_imp200_det_az = list(filter(lambda x: x!= -99., dv400_imp200_det_az))
    dv400_imp100_det_az = list(filter(lambda x: x!= -99., dv400_imp100_det_az))
    dv400_imp50_det_az = list(filter(lambda x: x!= -99., dv400_imp50_det_az))
    dv400_imp25_det_az = list(filter(lambda x: x!= -99., dv400_imp25_det_az))

    dv400_imp1000_non_az = list(filter(lambda x: x!= -99., dv400_imp1000_non_az))
    dv400_imp750_non_az = list(filter(lambda x: x!= -99., dv400_imp750_non_az))
    dv400_imp500_non_az = list(filter(lambda x: x!= -99., dv400_imp500_non_az))
    dv400_imp400_non_az = list(filter(lambda x: x!= -99., dv400_imp400_non_az))
    dv400_imp300_non_az = list(filter(lambda x: x!= -99., dv400_imp300_non_az))
    dv400_imp200_non_az = list(filter(lambda x: x!= -99., dv400_imp200_non_az))
    dv400_imp100_non_az = list(filter(lambda x: x!= -99., dv400_imp100_non_az))
    dv400_imp50_non_az = list(filter(lambda x: x!= -99., dv400_imp50_non_az))
    dv400_imp25_non_az = list(filter(lambda x: x!= -99., dv400_imp25_non_az))
    
    
    # minEW = 50
    # likelihood version first - azlinations
    dv400_l0001_min50_det_az = list(filter(lambda x: x!= -99., dv400_l0001_min50_det_az))
    dv400_l0005_min50_det_az = list(filter(lambda x: x!= -99., dv400_l0005_min50_det_az))
    dv400_l001_min50_det_az = list(filter(lambda x: x!= -99., dv400_l001_min50_det_az))
    dv400_l005_min50_det_az = list(filter(lambda x: x!= -99., dv400_l005_min50_det_az))
    dv400_l01_min50_det_az = list(filter(lambda x: x!= -99., dv400_l01_min50_det_az))
    dv400_l05_min50_det_az = list(filter(lambda x: x!= -99., dv400_l05_min50_det_az))
    dv400_l1_min50_det_az = list(filter(lambda x: x!= -99., dv400_l1_min50_det_az))
    dv400_l5_min50_det_az = list(filter(lambda x: x!= -99., dv400_l5_min50_det_az))
    dv400_l75_min50_det_az = list(filter(lambda x: x!= -99., dv400_l75_min50_det_az))

    dv400_l0001_min50_non_az = list(filter(lambda x: x!= -99., dv400_l0001_min50_non_az))
    dv400_l0005_min50_non_az = list(filter(lambda x: x!= -99., dv400_l0005_min50_non_az))
    dv400_l001_min50_non_az = list(filter(lambda x: x!= -99., dv400_l001_min50_non_az))
    dv400_l005_min50_non_az = list(filter(lambda x: x!= -99., dv400_l005_min50_non_az))
    dv400_l01_min50_non_az = list(filter(lambda x: x!= -99., dv400_l01_min50_non_az))
    dv400_l05_min50_non_az = list(filter(lambda x: x!= -99., dv400_l05_min50_non_az))
    dv400_l1_min50_non_az = list(filter(lambda x: x!= -99., dv400_l1_min50_non_az))
    dv400_l5_min50_non_az = list(filter(lambda x: x!= -99., dv400_l5_min50_non_az))
    dv400_l75_min50_non_az = list(filter(lambda x: x!= -99., dv400_l75_min50_non_az))
    
    # now impact version - azlinations
    dv400_imp1000_min50_det_az = list(filter(lambda x: x!= -99., dv400_imp1000_min50_det_az))
    dv400_imp750_min50_det_az = list(filter(lambda x: x!= -99., dv400_imp750_min50_det_az))
    dv400_imp500_min50_det_az = list(filter(lambda x: x!= -99., dv400_imp500_min50_det_az))
    dv400_imp400_min50_det_az = list(filter(lambda x: x!= -99., dv400_imp400_min50_det_az))
    dv400_imp300_min50_det_az = list(filter(lambda x: x!= -99., dv400_imp300_min50_det_az))
    dv400_imp200_min50_det_az = list(filter(lambda x: x!= -99., dv400_imp200_min50_det_az))
    dv400_imp100_min50_det_az = list(filter(lambda x: x!= -99., dv400_imp100_min50_det_az))
    dv400_imp50_min50_det_az = list(filter(lambda x: x!= -99., dv400_imp50_min50_det_az))
    dv400_imp25_min50_det_az = list(filter(lambda x: x!= -99., dv400_imp25_min50_det_az))

    dv400_imp1000_min50_non_az = list(filter(lambda x: x!= -99., dv400_imp1000_min50_non_az))
    dv400_imp750_min50_non_az = list(filter(lambda x: x!= -99., dv400_imp750_min50_non_az))
    dv400_imp500_min50_non_az = list(filter(lambda x: x!= -99., dv400_imp500_min50_non_az))
    dv400_imp400_min50_non_az = list(filter(lambda x: x!= -99., dv400_imp400_min50_non_az))
    dv400_imp300_min50_non_az = list(filter(lambda x: x!= -99., dv400_imp300_min50_non_az))
    dv400_imp200_min50_non_az = list(filter(lambda x: x!= -99., dv400_imp200_min50_non_az))
    dv400_imp100_min50_non_az = list(filter(lambda x: x!= -99., dv400_imp100_min50_non_az))
    dv400_imp50_min50_non_az = list(filter(lambda x: x!= -99., dv400_imp50_min50_non_az))
    dv400_imp25_min50_non_az = list(filter(lambda x: x!= -99., dv400_imp25_min50_non_az))
    
    # minEW = 100
    # likelihood version first - azlinations
#     dv400_l0001_min100_det_az = list(filter(lambda x: x!= -99., dv400_l0001_min100_det_az))
#     dv400_l0005_min100_det_az = list(filter(lambda x: x!= -99., dv400_l0005_min100_det_az))
#     dv400_l001_min100_det_az = list(filter(lambda x: x!= -99., dv400_l001_min100_det_az))
#     dv400_l005_min100_det_az = list(filter(lambda x: x!= -99., dv400_l005_min100_det_az))
#     dv400_l01_min100_det_az = list(filter(lambda x: x!= -99., dv400_l01_min100_det_az))
#     dv400_l05_min100_det_az = list(filter(lambda x: x!= -99., dv400_l05_min100_det_az))
#     dv400_l1_min100_det_az = list(filter(lambda x: x!= -99., dv400_l1_min100_det_az))
#     dv400_l5_min100_det_az = list(filter(lambda x: x!= -99., dv400_l5_min100_det_az))
#     dv400_l75_min100_det_az = list(filter(lambda x: x!= -99., dv400_l75_min100_det_az))
# 
#     dv400_l0001_min100_non_az = list(filter(lambda x: x!= -99., dv400_l0001_min100_non_az))
#     dv400_l0005_min100_non_az = list(filter(lambda x: x!= -99., dv400_l0005_min100_non_az))
#     dv400_l001_min100_non_az = list(filter(lambda x: x!= -99., dv400_l001_min100_non_az))
#     dv400_l005_min100_non_az = list(filter(lambda x: x!= -99., dv400_l005_min100_non_az))
#     dv400_l01_min100_non_az = list(filter(lambda x: x!= -99., dv400_l01_min100_non_az))
#     dv400_l05_min100_non_az = list(filter(lambda x: x!= -99., dv400_l05_min100_non_az))
#     dv400_l1_min100_non_az = list(filter(lambda x: x!= -99., dv400_l1_min100_non_az))
#     dv400_l5_min100_non_az = list(filter(lambda x: x!= -99., dv400_l5_min100_non_az))
#     dv400_l75_min100_non_az = list(filter(lambda x: x!= -99., dv400_l75_min100_non_az))
#     
#     # now impact version - azlinations
#     dv400_imp1000_min100_det_az = list(filter(lambda x: x!= -99., dv400_imp1000_min100_det_az))
#     dv400_imp750_min100_det_az = list(filter(lambda x: x!= -99., dv400_imp750_min100_det_az))
#     dv400_imp500_min100_det_az = list(filter(lambda x: x!= -99., dv400_imp500_min100_det_az))
#     dv400_imp400_min100_det_az = list(filter(lambda x: x!= -99., dv400_imp400_min100_det_az))
#     dv400_imp300_min100_det_az = list(filter(lambda x: x!= -99., dv400_imp300_min100_det_az))
#     dv400_imp200_min100_det_az = list(filter(lambda x: x!= -99., dv400_imp200_min100_det_az))
#     dv400_imp100_min100_det_az = list(filter(lambda x: x!= -99., dv400_imp100_min100_det_az))
#     dv400_imp50_min100_det_az = list(filter(lambda x: x!= -99., dv400_imp50_min100_det_az))
#     dv400_imp25_min100_det_az = list(filter(lambda x: x!= -99., dv400_imp25_min100_det_az))
# 
#     dv400_imp1000_min100_non_az = list(filter(lambda x: x!= -99., dv400_imp1000_min100_non_az))
#     dv400_imp750_min100_non_az = list(filter(lambda x: x!= -99., dv400_imp750_min100_non_az))
#     dv400_imp500_min100_non_az = list(filter(lambda x: x!= -99., dv400_imp500_min100_non_az))
#     dv400_imp400_min100_non_az = list(filter(lambda x: x!= -99., dv400_imp400_min100_non_az))
#     dv400_imp300_min100_non_az = list(filter(lambda x: x!= -99., dv400_imp300_min100_non_az))
#     dv400_imp200_min100_non_az = list(filter(lambda x: x!= -99., dv400_imp200_min100_non_az))
#     dv400_imp100_min100_non_az = list(filter(lambda x: x!= -99., dv400_imp100_min100_non_az))
#     dv400_imp50_min100_non_az = list(filter(lambda x: x!= -99., dv400_imp50_min100_non_az))
#     dv400_imp25_min100_non_az = list(filter(lambda x: x!= -99., dv400_imp25_min100_non_az))
    
    # minEW = 200
    # likelihood version first - azlinations
    dv400_l0001_min200_det_az = list(filter(lambda x: x!= -99., dv400_l0001_min200_det_az))
    dv400_l0005_min200_det_az = list(filter(lambda x: x!= -99., dv400_l0005_min200_det_az))
    dv400_l001_min200_det_az = list(filter(lambda x: x!= -99., dv400_l001_min200_det_az))
    dv400_l005_min200_det_az = list(filter(lambda x: x!= -99., dv400_l005_min200_det_az))
    dv400_l01_min200_det_az = list(filter(lambda x: x!= -99., dv400_l01_min200_det_az))
    dv400_l05_min200_det_az = list(filter(lambda x: x!= -99., dv400_l05_min200_det_az))
    dv400_l1_min200_det_az = list(filter(lambda x: x!= -99., dv400_l1_min200_det_az))
    dv400_l5_min200_det_az = list(filter(lambda x: x!= -99., dv400_l5_min200_det_az))
    dv400_l75_min200_det_az = list(filter(lambda x: x!= -99., dv400_l75_min200_det_az))

    dv400_l0001_min200_non_az = list(filter(lambda x: x!= -99., dv400_l0001_min200_non_az))
    dv400_l0005_min200_non_az = list(filter(lambda x: x!= -99., dv400_l0005_min200_non_az))
    dv400_l001_min200_non_az = list(filter(lambda x: x!= -99., dv400_l001_min200_non_az))
    dv400_l005_min200_non_az = list(filter(lambda x: x!= -99., dv400_l005_min200_non_az))
    dv400_l01_min200_non_az = list(filter(lambda x: x!= -99., dv400_l01_min200_non_az))
    dv400_l05_min200_non_az = list(filter(lambda x: x!= -99., dv400_l05_min200_non_az))
    dv400_l1_min200_non_az = list(filter(lambda x: x!= -99., dv400_l1_min200_non_az))
    dv400_l5_min200_non_az = list(filter(lambda x: x!= -99., dv400_l5_min200_non_az))
    dv400_l75_min200_non_az = list(filter(lambda x: x!= -99., dv400_l75_min200_non_az))
    
    # now impact version - azlinations
    dv400_imp1000_min200_det_az = list(filter(lambda x: x!= -99., dv400_imp1000_min200_det_az))
    dv400_imp750_min200_det_az = list(filter(lambda x: x!= -99., dv400_imp750_min200_det_az))
    dv400_imp500_min200_det_az = list(filter(lambda x: x!= -99., dv400_imp500_min200_det_az))
    dv400_imp400_min200_det_az = list(filter(lambda x: x!= -99., dv400_imp400_min200_det_az))
    dv400_imp300_min200_det_az = list(filter(lambda x: x!= -99., dv400_imp300_min200_det_az))
    dv400_imp200_min200_det_az = list(filter(lambda x: x!= -99., dv400_imp200_min200_det_az))
    dv400_imp100_min200_det_az = list(filter(lambda x: x!= -99., dv400_imp100_min200_det_az))
    dv400_imp50_min200_det_az = list(filter(lambda x: x!= -99., dv400_imp50_min200_det_az))
    dv400_imp25_min200_det_az = list(filter(lambda x: x!= -99., dv400_imp25_min200_det_az))

    dv400_imp1000_min200_non_az = list(filter(lambda x: x!= -99., dv400_imp1000_min200_non_az))
    dv400_imp750_min200_non_az = list(filter(lambda x: x!= -99., dv400_imp750_min200_non_az))
    dv400_imp500_min200_non_az = list(filter(lambda x: x!= -99., dv400_imp500_min200_non_az))
    dv400_imp400_min200_non_az = list(filter(lambda x: x!= -99., dv400_imp400_min200_non_az))
    dv400_imp300_min200_non_az = list(filter(lambda x: x!= -99., dv400_imp300_min200_non_az))
    dv400_imp200_min200_non_az = list(filter(lambda x: x!= -99., dv400_imp200_min200_non_az))
    dv400_imp100_min200_non_az = list(filter(lambda x: x!= -99., dv400_imp100_min200_non_az))
    dv400_imp50_min200_non_az = list(filter(lambda x: x!= -99., dv400_imp50_min200_non_az))
    dv400_imp25_min200_non_az = list(filter(lambda x: x!= -99., dv400_imp25_min200_non_az))
    
    # minEW = 300
    # likelihood version first - azlinations
#     dv400_l0001_min300_det_az = list(filter(lambda x: x!= -99., dv400_l0001_min300_det_az))
#     dv400_l0005_min300_det_az = list(filter(lambda x: x!= -99., dv400_l0005_min300_det_az))
#     dv400_l001_min300_det_az = list(filter(lambda x: x!= -99., dv400_l001_min300_det_az))
#     dv400_l005_min300_det_az = list(filter(lambda x: x!= -99., dv400_l005_min300_det_az))
#     dv400_l01_min300_det_az = list(filter(lambda x: x!= -99., dv400_l01_min300_det_az))
#     dv400_l05_min300_det_az = list(filter(lambda x: x!= -99., dv400_l05_min300_det_az))
#     dv400_l1_min300_det_az = list(filter(lambda x: x!= -99., dv400_l1_min300_det_az))
#     dv400_l5_min300_det_az = list(filter(lambda x: x!= -99., dv400_l5_min300_det_az))
#     dv400_l75_min300_det_az = list(filter(lambda x: x!= -99., dv400_l75_min300_det_az))
# 
#     dv400_l0001_min300_non_az = list(filter(lambda x: x!= -99., dv400_l0001_min300_non_az))
#     dv400_l0005_min300_non_az = list(filter(lambda x: x!= -99., dv400_l0005_min300_non_az))
#     dv400_l001_min300_non_az = list(filter(lambda x: x!= -99., dv400_l001_min300_non_az))
#     dv400_l005_min300_non_az = list(filter(lambda x: x!= -99., dv400_l005_min300_non_az))
#     dv400_l01_min300_non_az = list(filter(lambda x: x!= -99., dv400_l01_min300_non_az))
#     dv400_l05_min300_non_az = list(filter(lambda x: x!= -99., dv400_l05_min300_non_az))
#     dv400_l1_min300_non_az = list(filter(lambda x: x!= -99., dv400_l1_min300_non_az))
#     dv400_l5_min300_non_az = list(filter(lambda x: x!= -99., dv400_l5_min300_non_az))
#     dv400_l75_min300_non_az = list(filter(lambda x: x!= -99., dv400_l75_min300_non_az))
#    
#     # now impact version - azlinations
#     dv400_imp1000_min300_det_az = list(filter(lambda x: x!= -99., dv400_imp1000_min300_det_az))
#     dv400_imp750_min300_det_az = list(filter(lambda x: x!= -99., dv400_imp750_min300_det_az))
#     dv400_imp500_min300_det_az = list(filter(lambda x: x!= -99., dv400_imp500_min300_det_az))
#     dv400_imp400_min300_det_az = list(filter(lambda x: x!= -99., dv400_imp400_min300_det_az))
#     dv400_imp300_min300_det_az = list(filter(lambda x: x!= -99., dv400_imp300_min300_det_az))
#     dv400_imp200_min300_det_az = list(filter(lambda x: x!= -99., dv400_imp200_min300_det_az))
#     dv400_imp100_min300_det_az = list(filter(lambda x: x!= -99., dv400_imp100_min300_det_az))
#     dv400_imp50_min300_det_az = list(filter(lambda x: x!= -99., dv400_imp50_min300_det_az))
#     dv400_imp25_min300_det_az = list(filter(lambda x: x!= -99., dv400_imp25_min300_det_az))
# 
#     dv400_imp1000_min300_non_az = list(filter(lambda x: x!= -99., dv400_imp1000_min300_non_az))
#     dv400_imp750_min300_non_az = list(filter(lambda x: x!= -99., dv400_imp750_min300_non_az))
#     dv400_imp500_min300_non_az = list(filter(lambda x: x!= -99., dv400_imp500_min300_non_az))
#     dv400_imp400_min300_non_az = list(filter(lambda x: x!= -99., dv400_imp400_min300_non_az))
#     dv400_imp300_min300_non_az = list(filter(lambda x: x!= -99., dv400_imp300_min300_non_az))
#     dv400_imp200_min300_non_az = list(filter(lambda x: x!= -99., dv400_imp200_min300_non_az))
#     dv400_imp100_min300_non_az = list(filter(lambda x: x!= -99., dv400_imp100_min300_non_az))
#     dv400_imp50_min300_non_az = list(filter(lambda x: x!= -99., dv400_imp50_min300_non_az))
#     dv400_imp25_min300_non_az = list(filter(lambda x: x!= -99., dv400_imp25_min300_non_az))
    
    
    print 'dv400_imp25_det_az: ',dv400_imp25_det_az
    print 'dv400_imp25_non_az: ',dv400_imp25_non_az
    print 'dv400_l75_det_az: ',dv400_l75_det_az
    print 'dv400_l75_non_az: ',dv400_l75_non_az
    print
    
    print
    print
    print 'Detection fraction for 1000 kpc mean az: ', bmean(np.array(dv400_imp1000_det_az)), ', ',bmean(np.array(dv400_imp1000_non_az))
    print 'Detection fraction for 750 kpc mean az: ', bmean(np.array(dv400_imp750_det_az)), ', ',bmean(np.array(dv400_imp750_non_az))
    print 'Detection fraction for 500 kpc mean az: ', bmean(np.array(dv400_imp500_det_az)), ', ',bmean(np.array(dv400_imp500_non_az))
    print 'Detection fraction for 400 kpc mean az: ', bmean(np.array(dv400_imp400_det_az)), ', ',bmean(np.array(dv400_imp400_non_az))
    print 'Detection fraction for 300 kpc mean az: ', bmean(np.array(dv400_imp300_det_az)), ', ',bmean(np.array(dv400_imp300_non_az))
    print 'Detection fraction for 200 kpc mean az: ', bmean(np.array(dv400_imp200_det_az)), ', ',bmean(np.array(dv400_imp200_non_az))
    print 'Detection fraction for 100 kpc mean az: ', bmean(np.array(dv400_imp100_det_az)), ', ',bmean(np.array(dv400_imp100_non_az))
    print 'Detection fraction for 50 kpc mean az: ', bmean(np.array(dv400_imp50_det_az)), ', ',bmean(np.array(dv400_imp50_non_az))
    print 'Detection fraction for 25 kpc mean az: ', bmean(np.array(dv400_imp25_det_az)), ', ',bmean(np.array(dv400_imp25_non_az))
    print
    print
    print
    print
    print 'Detection fraction for 0.0001 L mean az: ', bmean(np.array(dv400_l0001_det_az)), ', ',bmean(np.array(dv400_l0001_non_az))
    print 'Detection fraction for 0.0005 L mean az: ', bmean(np.array(dv400_l0005_det_az)), ', ',bmean(np.array(dv400_l0005_non_az))
    print 'Detection fraction for 0.001 L mean az: ', bmean(np.array(dv400_l001_det_az)), ', ',bmean(np.array(dv400_l001_non_az))
    print 'Detection fraction for 0.00 L mean az: ', bmean(np.array(dv400_l005_det_az)), ', ',bmean(np.array(dv400_l005_non_az))
    print 'Detection fraction for 0.01 L mean az: ', bmean(np.array(dv400_l01_det_az)), ', ',bmean(np.array(dv400_l01_non_az))
    print 'Detection fraction for 0.05 L mean az: ', bmean(np.array(dv400_l05_det_az)), ', ',bmean(np.array(dv400_l05_non_az))
    print 'Detection fraction for 0.1 L mean az: ', bmean(np.array(dv400_l1_det_az)), ', ',bmean(np.array(dv400_l1_non_az))
    print 'Detection fraction for 0.5 L mean az: ', bmean(np.array(dv400_l5_det_az)), ', ',bmean(np.array(dv400_l5_non_az))
    print 'Detection fraction for 0.75 L mean az: ', bmean(np.array(dv400_l75_det_az)), ', ',bmean(np.array(dv400_l75_non_az))



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
    
    print
    print "Starting bootstrap..."
    
    if bootstrap:
    
        # likelihood - azlinations
        dv400_l0001_det_az_meanerr, dv400_l0001_det_az_medianerr = return_bootstrap_errors(dv400_l0001_det_az, reps)
        dv400_l0005_det_az_meanerr, dv400_l0005_det_az_medianerr = return_bootstrap_errors(dv400_l0005_det_az, reps)
        dv400_l001_det_az_meanerr, dv400_l001_det_az_medianerr = return_bootstrap_errors(dv400_l001_det_az, reps)
        dv400_l005_det_az_meanerr, dv400_l005_det_az_medianerr = return_bootstrap_errors(dv400_l005_det_az, reps)
        dv400_l01_det_az_meanerr, dv400_l01_det_az_medianerr = return_bootstrap_errors(dv400_l01_det_az, reps)
        dv400_l05_det_az_meanerr, dv400_l05_det_az_medianerr = return_bootstrap_errors(dv400_l05_det_az, reps)
        dv400_l1_det_az_meanerr, dv400_l1_det_az_medianerr = return_bootstrap_errors(dv400_l1_det_az, reps)
        dv400_l5_det_az_meanerr, dv400_l5_det_az_medianerr = return_bootstrap_errors(dv400_l5_det_az, reps)
        dv400_l75_det_az_meanerr, dv400_l75_det_az_medianerr = return_bootstrap_errors(dv400_l75_det_az, reps)

        dv400_l0001_non_az_meanerr, dv400_l0001_non_az_medianerr = return_bootstrap_errors(dv400_l0001_non_az, reps)
        dv400_l0005_non_az_meanerr, dv400_l0005_non_az_medianerr = return_bootstrap_errors(dv400_l0005_non_az, reps)
        dv400_l001_non_az_meanerr, dv400_l001_non_az_medianerr = return_bootstrap_errors(dv400_l001_non_az, reps)
        dv400_l005_non_az_meanerr, dv400_l005_non_az_medianerr = return_bootstrap_errors(dv400_l005_non_az, reps)
        dv400_l01_non_az_meanerr, dv400_l01_non_az_medianerr = return_bootstrap_errors(dv400_l01_non_az, reps)
        dv400_l05_non_az_meanerr, dv400_l05_non_az_medianerr = return_bootstrap_errors(dv400_l05_non_az, reps)
        dv400_l1_non_az_meanerr, dv400_l1_non_az_medianerr = return_bootstrap_errors(dv400_l1_non_az, reps)
        dv400_l5_non_az_meanerr, dv400_l5_non_az_medianerr = return_bootstrap_errors(dv400_l5_non_az, reps)
        dv400_l75_non_az_meanerr, dv400_l75_non_az_medianerr = return_bootstrap_errors(dv400_l75_non_az, reps)

        print
        print 'Finished likelihood - azlinations bootstraping...'
    
        # impact - azlinations
        dv400_imp1000_det_az_meanerr, dv400_imp1000_det_az_medianerr = return_bootstrap_errors(dv400_imp1000_det_az, reps)
        dv400_imp750_det_az_meanerr, dv400_imp750_det_az_medianerr = return_bootstrap_errors(dv400_imp750_det_az, reps)
        dv400_imp500_det_az_meanerr, dv400_imp500_det_az_medianerr = return_bootstrap_errors(dv400_imp500_det_az, reps)
        dv400_imp400_det_az_meanerr, dv400_imp400_det_az_medianerr = return_bootstrap_errors(dv400_imp400_det_az, reps)
        dv400_imp300_det_az_meanerr, dv400_imp300_det_az_medianerr = return_bootstrap_errors(dv400_imp300_det_az, reps)
        dv400_imp200_det_az_meanerr, dv400_imp200_det_az_medianerr = return_bootstrap_errors(dv400_imp200_det_az, reps)
        dv400_imp100_det_az_meanerr, dv400_imp100_det_az_medianerr = return_bootstrap_errors(dv400_imp100_det_az, reps)
        dv400_imp50_det_az_meanerr, dv400_imp50_det_az_medianerr = return_bootstrap_errors(dv400_imp50_det_az, reps)
        dv400_imp25_det_az_meanerr, dv400_imp25_det_az_medianerr = return_bootstrap_errors(dv400_imp25_det_az, reps)

        dv400_imp1000_non_az_meanerr, dv400_imp1000_non_az_medianerr = return_bootstrap_errors(dv400_imp1000_non_az, reps)
        dv400_imp750_non_az_meanerr, dv400_imp750_non_az_medianerr = return_bootstrap_errors(dv400_imp750_non_az, reps)
        dv400_imp500_non_az_meanerr, dv400_imp500_non_az_medianerr = return_bootstrap_errors(dv400_imp500_non_az, reps)
        dv400_imp400_non_az_meanerr, dv400_imp400_non_az_medianerr = return_bootstrap_errors(dv400_imp400_non_az, reps)
        dv400_imp300_non_az_meanerr, dv400_imp300_non_az_medianerr = return_bootstrap_errors(dv400_imp300_non_az, reps)
        dv400_imp200_non_az_meanerr, dv400_imp200_non_az_medianerr = return_bootstrap_errors(dv400_imp200_non_az, reps)
        dv400_imp100_non_az_meanerr, dv400_imp100_non_az_medianerr = return_bootstrap_errors(dv400_imp100_non_az, reps)
        dv400_imp50_non_az_meanerr, dv400_imp50_non_az_medianerr = return_bootstrap_errors(dv400_imp50_non_az, reps)
        dv400_imp25_non_az_meanerr, dv400_imp25_non_az_medianerr = return_bootstrap_errors(dv400_imp25_non_az, reps)
    
        print
        print 'Finished impact - azlinations bootstraping...'
        print 'dv400_imp300_non_az_meanerr, dv400_imp300_non_az_medianerr = ', dv400_imp300_non_az_meanerr, dv400_imp300_non_az_medianerr


        # minEW = 200
        # likelihood - azlinations
        dv400_l0001_min200_det_az_meanerr, dv400_l0001_min200_det_az_medianerr = return_bootstrap_errors(dv400_l0001_min200_det_az, reps)
        dv400_l0005_min200_det_az_meanerr, dv400_l0005_min200_det_az_medianerr = return_bootstrap_errors(dv400_l0005_min200_det_az, reps)
        dv400_l001_min200_det_az_meanerr, dv400_l001_min200_det_az_medianerr = return_bootstrap_errors(dv400_l001_min200_det_az, reps)
        dv400_l005_min200_det_az_meanerr, dv400_l005_min200_det_az_medianerr = return_bootstrap_errors(dv400_l005_min200_det_az, reps)
        dv400_l01_min200_det_az_meanerr, dv400_l01_min200_det_az_medianerr = return_bootstrap_errors(dv400_l01_min200_det_az, reps)
        dv400_l05_min200_det_az_meanerr, dv400_l05_min200_det_az_medianerr = return_bootstrap_errors(dv400_l05_min200_det_az, reps)
        dv400_l1_min200_det_az_meanerr, dv400_l1_min200_det_az_medianerr = return_bootstrap_errors(dv400_l1_min200_det_az, reps)
        dv400_l5_min200_det_az_meanerr, dv400_l5_min200_det_az_medianerr = return_bootstrap_errors(dv400_l5_min200_det_az, reps)
        dv400_l75_min200_det_az_meanerr, dv400_l75_min200_det_az_medianerr = return_bootstrap_errors(dv400_l75_min200_det_az, reps)

        dv400_l0001_min200_non_az_meanerr, dv400_l0001_min200_non_az_medianerr = return_bootstrap_errors(dv400_l0001_min200_non_az, reps)
        dv400_l0005_min200_non_az_meanerr, dv400_l0005_min200_non_az_medianerr = return_bootstrap_errors(dv400_l0005_min200_non_az, reps)
        dv400_l001_min200_non_az_meanerr, dv400_l001_min200_non_az_medianerr = return_bootstrap_errors(dv400_l001_min200_non_az, reps)
        dv400_l005_min200_non_az_meanerr, dv400_l005_min200_non_az_medianerr = return_bootstrap_errors(dv400_l005_min200_non_az, reps)
        dv400_l01_min200_non_az_meanerr, dv400_l01_min200_non_az_medianerr = return_bootstrap_errors(dv400_l01_min200_non_az, reps)
        dv400_l05_min200_non_az_meanerr, dv400_l05_min200_non_az_medianerr = return_bootstrap_errors(dv400_l05_min200_non_az, reps)
        dv400_l1_min200_non_az_meanerr, dv400_l1_min200_non_az_medianerr = return_bootstrap_errors(dv400_l1_min200_non_az, reps)
        dv400_l5_min200_non_az_meanerr, dv400_l5_min200_non_az_medianerr = return_bootstrap_errors(dv400_l5_min200_non_az, reps)
        dv400_l75_min200_non_az_meanerr, dv400_l75_min200_non_az_medianerr = return_bootstrap_errors(dv400_l75_min200_non_az, reps)

        print
        print 'Finished minEW = 200 likelihood - azlinations bootstraping...'
    
        # impact - azlinations
        dv400_imp1000_min200_det_az_meanerr, dv400_imp1000_min200_det_az_medianerr = return_bootstrap_errors(dv400_imp1000_min200_det_az, reps)
        dv400_imp750_min200_det_az_meanerr, dv400_imp750_min200_det_az_medianerr = return_bootstrap_errors(dv400_imp750_min200_det_az, reps)
        dv400_imp500_min200_det_az_meanerr, dv400_imp500_min200_det_az_medianerr = return_bootstrap_errors(dv400_imp500_min200_det_az, reps)
        dv400_imp400_min200_det_az_meanerr, dv400_imp400_min200_det_az_medianerr = return_bootstrap_errors(dv400_imp400_min200_det_az, reps)
        dv400_imp300_min200_det_az_meanerr, dv400_imp300_min200_det_az_medianerr = return_bootstrap_errors(dv400_imp300_min200_det_az, reps)
        dv400_imp200_min200_det_az_meanerr, dv400_imp200_min200_det_az_medianerr = return_bootstrap_errors(dv400_imp200_min200_det_az, reps)
        dv400_imp100_min200_det_az_meanerr, dv400_imp100_min200_det_az_medianerr = return_bootstrap_errors(dv400_imp100_min200_det_az, reps)
        dv400_imp50_min200_det_az_meanerr, dv400_imp50_min200_det_az_medianerr = return_bootstrap_errors(dv400_imp50_min200_det_az, reps)
        dv400_imp25_min200_det_az_meanerr, dv400_imp25_min200_det_az_medianerr = return_bootstrap_errors(dv400_imp25_min200_det_az, reps)

        dv400_imp1000_min200_non_az_meanerr, dv400_imp1000_min200_non_az_medianerr = return_bootstrap_errors(dv400_imp1000_min200_non_az, reps)
        dv400_imp750_min200_non_az_meanerr, dv400_imp750_min200_non_az_medianerr = return_bootstrap_errors(dv400_imp750_min200_non_az, reps)
        dv400_imp500_min200_non_az_meanerr, dv400_imp500_min200_non_az_medianerr = return_bootstrap_errors(dv400_imp500_min200_non_az, reps)
        dv400_imp400_min200_non_az_meanerr, dv400_imp400_min200_non_az_medianerr = return_bootstrap_errors(dv400_imp400_min200_non_az, reps)
        dv400_imp300_min200_non_az_meanerr, dv400_imp300_min200_non_az_medianerr = return_bootstrap_errors(dv400_imp300_min200_non_az, reps)
        dv400_imp200_min200_non_az_meanerr, dv400_imp200_min200_non_az_medianerr = return_bootstrap_errors(dv400_imp200_min200_non_az, reps)
        dv400_imp100_min200_non_az_meanerr, dv400_imp100_min200_non_az_medianerr = return_bootstrap_errors(dv400_imp100_min200_non_az, reps)
        dv400_imp50_min200_non_az_meanerr, dv400_imp50_min200_non_az_medianerr = return_bootstrap_errors(dv400_imp50_min200_non_az, reps)
        dv400_imp25_min200_non_az_meanerr, dv400_imp25_min200_non_az_medianerr = return_bootstrap_errors(dv400_imp25_min200_non_az, reps)
    
        print
        print 'Finished minEW = 200 impact - azlinations bootstraping...'
        print 'dv400_imp300_min200_non_az_meanerr, dv400_imp300_min200_non_az_medianerr = ', dv400_imp300_min200_non_az_meanerr, dv400_imp300_min200_non_az_medianerr




##########################################################################################
##########################################################################################
    
    if plot_detection_fraction_likelihood_az:
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
        

        alpha_det_mean = 0.9
        alpha_non_mean = 0.9
        alpha_det_median = 0.6
        alpha_non_median = 0.6
        markerSize = 10
        lw = 2.
        
        binSize = 100
        bins = arange(0, 600, binSize)

        
        label_det_mean = r'$\rm Detections - Mean~az.$'
        label_non_mean = r'$\rm Non-Detections - Mean ~az.$'
        
        label_det_median = r'$\rm Detections - Median ~az.$'
        label_non_median = r'$\rm Non-Detections - Median ~az.$'

        symbol_det_mean = 'D'
        symbol_non_mean = 'X'
        symbol_det_median = 'D'
        symbol_non_median = 'X'
        
        color_det_mean = color_blue
        color_non_mean = color_red
        color_det_median = color_purple2
        color_non_median = color_orange

        ls_det_mean = 'solid'
        ls_non_mean = 'dashed'
        ls_det_median = 'solid'
        ls_non_median = 'dashed'

        maxEW = 15000.

##########################################################################################
        # do the plotting
        
        print 'len(dv400_l5_det_az) : ',len(dv400_l1_det_az)
        print 'len(dv400_l5_det_az) : ',len(dv400_l5_det_az)
        print 'len(dv400_l75_det_az) : ',len(dv400_l75_det_az)
        print
        print
        print 'len(dv400_l5_non_az) : ',len(dv400_l1_non_az)
        print 'len(dv400_l5_non_az) : ',len(dv400_l5_non_az)
        print 'len(dv400_l75_non_az) : ',len(dv400_l75_non_az)
        print            
            

        x = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 0.75]
        y_det_mean = [bmean(dv400_l0001_det_az),
                    bmean(dv400_l0005_det_az),
                    bmean(dv400_l001_det_az),
                    bmean(dv400_l005_det_az),
                    bmean(dv400_l01_det_az),
                    bmean(dv400_l05_det_az),
                    bmean(dv400_l1_det_az),
                    bmean(dv400_l5_det_az),
                    bmean(dv400_l75_det_az)]
                    
        # bootstrap errors in MEAN detection azlination
        y_det_err_mean = [dv400_l0001_det_az_meanerr,
                        dv400_l0005_det_az_meanerr,
                        dv400_l001_det_az_meanerr,
                        dv400_l005_det_az_meanerr,
                        dv400_l01_det_az_meanerr,
                        dv400_l05_det_az_meanerr,
                        dv400_l1_det_az_meanerr,
                        dv400_l5_det_az_meanerr,
                        dv400_l75_det_az_meanerr]
            
        
        y_non_mean = [bmean(dv400_l0001_non_az),
                    bmean(dv400_l0005_non_az),
                    bmean(dv400_l001_non_az),
                    bmean(dv400_l005_non_az), 
                    bmean(dv400_l01_non_az),
                    bmean(dv400_l05_non_az),
                    bmean(dv400_l1_non_az),
                    bmean(dv400_l5_non_az),
                    bmean(dv400_l75_non_az)]
        
        # bootstrap errors in MEAN non-detection azlination
        y_non_err_mean = [dv400_l0001_non_az_meanerr,
                        dv400_l0005_non_az_meanerr,
                        dv400_l001_non_az_meanerr,
                        dv400_l005_non_az_meanerr,
                        dv400_l01_non_az_meanerr,
                        dv400_l05_non_az_meanerr,
                        dv400_l1_non_az_meanerr,
                        dv400_l5_non_az_meanerr,
                        dv400_l75_non_az_meanerr]
                    
                
        y_det_median = [bmedian(dv400_l0001_det_az),
                    bmedian(dv400_l0005_det_az),
                    bmedian(dv400_l001_det_az),
                    bmedian(dv400_l005_det_az), 
                    bmedian(dv400_l01_det_az),
                    bmedian(dv400_l05_det_az),
                    bmedian(dv400_l1_det_az),
                    bmedian(dv400_l5_det_az),
                    bmedian(dv400_l75_det_az)]
                    
        # bootstrap errors in MEDIAN detection azlination
        y_det_err_median = [dv400_l0001_det_az_medianerr,
                        dv400_l0005_det_az_medianerr,
                        dv400_l001_det_az_medianerr,
                        dv400_l005_det_az_medianerr,
                        dv400_l01_det_az_medianerr,
                        dv400_l05_det_az_medianerr,
                        dv400_l1_det_az_medianerr,
                        dv400_l5_det_az_medianerr,
                        dv400_l75_det_az_medianerr]
                    
            
        y_non_median = [bmedian(dv400_l0001_non_az),
                    bmedian(dv400_l0005_non_az),
                    bmedian(dv400_l001_non_az),
                    bmedian(dv400_l005_non_az), 
                    bmedian(dv400_l01_non_az),
                    bmedian(dv400_l05_non_az),
                    bmedian(dv400_l1_non_az),
                    bmedian(dv400_l5_non_az),
                    bmedian(dv400_l75_non_az)]

        # bootstrap errors in MEDIAN non-detection azlination
        y_non_err_median = [dv400_l0001_non_az_medianerr,
                        dv400_l0005_non_az_medianerr,
                        dv400_l001_non_az_medianerr,
                        dv400_l005_non_az_medianerr,
                        dv400_l01_non_az_medianerr,
                        dv400_l05_non_az_medianerr,
                        dv400_l1_non_az_medianerr,
                        dv400_l5_non_az_medianerr,
                        dv400_l75_non_az_medianerr]
                    

        ###########

        # MEAN azlination for detections with errors
        ax1.errorbar(x,
                    y_det_mean,
                    yerr=y_det_err_mean,
                    marker=symbol_det_mean,
                    c=color_det_mean,
                    ms=markerSize,
                    markeredgecolor='black',
                    ls = ls_det_mean,
                    lw = lw,
                    alpha=alpha_det_mean,
                    label=label_det_mean)

        # MEDIAN azlination for detections with errors
        ax1.errorbar(x,
                    y_det_median,
                    yerr=y_det_err_median,
                    marker=symbol_det_median,
                    c=color_det_median,
                    ms=markerSize,
                    markeredgecolor='black',
                    ls = ls_det_median,
                    lw = lw,
                    alpha=alpha_det_median,
                    label=label_det_median)


        # MEAN azlination for non-detections with errors
        ax1.errorbar(x,
                    y_non_mean,
                    yerr=y_non_err_mean,
                    marker=symbol_non_mean,
                    c=color_non_mean,
                    ms=markerSize,
                    markeredgecolor='black',
                    ls = ls_non_mean,
                    lw = lw,
                    alpha=alpha_non_mean,
                    label=label_non_mean)

        # MEDIAN azlination for non-detections with errors
        ax1.errorbar(x,
                    y_non_median,
                    yerr=y_non_err_median,
                    marker=symbol_non_median,
                    c=color_non_median,
                    ms=markerSize,
                    markeredgecolor='black',
                    ls = ls_non_median,
                    lw = lw,
                    alpha=alpha_non_median,
                    label=label_non_median)


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
        
        ax1.set_ylabel(r'$\rm Azimuth~[deg]$')
        
        leg = ax1.legend(scatterpoints=1,prop={'size':12},loc='lower left',fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax1.grid(b=None,which='major',axis='both')
        ylim(0., 90.)
        xlim(0.001, 1.)
        ax1.invert_xaxis()


        if plot_detection_fraction_likelihood_az_save:
            savefig('{0}/detection_fraction_likelihood_az_errors_lstarcut{1}.pdf'.format(saveDirectory, lstar_cut),format='pdf',bbox_inches='tight')
        else:
            show()

##########################################################################################
##########################################################################################


##########################################################################################
##########################################################################################
    
    if plot_detection_fraction_impact_az:
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
        

        alpha_det_mean = 0.9
        alpha_non_mean = 0.9
        alpha_det_median = 0.6
        alpha_non_median = 0.6
        markerSize = 10
        lw = 2.
        
        binSize = 100
        bins = arange(0, 600, binSize)

        
        label_det_mean = r'$\rm Detections - Mean ~az.$'
        label_non_mean = r'$\rm Non-Detections - Mean ~az.$'
        
        label_det_median = r'$\rm Detections - Median ~az.$'
        label_non_median = r'$\rm Non-Detections - Median ~az.$'

        symbol_det_mean = 'D'
        symbol_non_mean = 'X'
        symbol_det_median = 'D'
        symbol_non_median = 'X'
        
        color_det_mean = color_blue
        color_non_mean = color_red
        color_det_median = color_purple2
        color_non_median = color_orange

        ls_det_mean = 'solid'
        ls_non_mean = 'dashed'
        ls_det_median = 'solid'
        ls_non_median = 'dashed'

        maxEW = 15000.

##########################################################################################
        # do the plotting
        
        print 'len(dv400_imp100_det_az) : ',len(dv400_imp100_det_az)
        print 'len(dv400_imp50_det_az) : ',len(dv400_imp50_det_az)
        print 'len(dv400_imp25_det_az) : ',len(dv400_imp25_det_az)
        print
        print
        print 'len(dv400_imp100_non_az) : ',len(dv400_imp100_non_az)
        print 'len(dv400_imp50_non_az) : ',len(dv400_imp50_non_az)
        print 'len(dv400_imp25_non_az) : ',len(dv400_imp25_non_az)
        print


        x = [1000, 750, 500, 400, 300, 200, 100, 50, 25]
        y_det_mean = [bmean(dv400_imp1000_det_az),
                    bmean(dv400_imp750_det_az),
                    bmean(dv400_imp500_det_az),
                    bmean(dv400_imp400_det_az),
                    bmean(dv400_imp300_det_az),
                    bmean(dv400_imp200_det_az),
                    bmean(dv400_imp100_det_az),
                    bmean(dv400_imp50_det_az),
                    bmean(dv400_imp25_det_az)]


        y_det_err_mean = [bmean(dv400_imp1000_det_az_meanerr),
                        bmean(dv400_imp750_det_az_meanerr),
                        bmean(dv400_imp500_det_az_meanerr),
                        bmean(dv400_imp400_det_az_meanerr),
                        bmean(dv400_imp300_det_az_meanerr),
                        bmean(dv400_imp200_det_az_meanerr),
                        bmean(dv400_imp100_det_az_meanerr),
                        bmean(dv400_imp50_det_az_meanerr),
                        bmean(dv400_imp25_det_az_meanerr)]


        y_non_mean = [bmean(dv400_imp1000_non_az),
                    bmean(dv400_imp750_non_az),
                    bmean(dv400_imp500_non_az),
                    bmean(dv400_imp400_non_az),
                    bmean(dv400_imp300_non_az),
                    bmean(dv400_imp200_non_az),
                    bmean(dv400_imp100_non_az),
                    bmean(dv400_imp50_non_az),
                    bmean(dv400_imp25_non_az)]
                
        y_non_err_mean = [bmean(dv400_imp1000_non_az_meanerr),
                        bmean(dv400_imp750_non_az_meanerr),
                        bmean(dv400_imp500_non_az_meanerr),
                        bmean(dv400_imp400_non_az_meanerr),
                        bmean(dv400_imp300_non_az_meanerr),
                        bmean(dv400_imp200_non_az_meanerr),
                        bmean(dv400_imp100_non_az_meanerr),
                        bmean(dv400_imp50_non_az_meanerr),
                        bmean(dv400_imp25_non_az_meanerr)]
                
                
                
        y_det_median = [bmedian(dv400_imp1000_det_az),
                        bmedian(dv400_imp750_det_az),
                        bmedian(dv400_imp500_det_az),
                        bmedian(dv400_imp400_det_az),
                        bmedian(dv400_imp300_det_az),
                        bmedian(dv400_imp200_det_az),
                        bmedian(dv400_imp100_det_az),
                        bmedian(dv400_imp50_det_az),
                        bmedian(dv400_imp25_det_az)]
            
        y_det_err_median = [bmedian(dv400_imp1000_det_az_medianerr),
                            bmedian(dv400_imp750_det_az_medianerr),
                            bmedian(dv400_imp500_det_az_medianerr),
                            bmedian(dv400_imp400_det_az_medianerr),
                            bmedian(dv400_imp300_det_az_medianerr),
                            bmedian(dv400_imp200_det_az_medianerr),
                            bmedian(dv400_imp100_det_az_medianerr),
                            bmedian(dv400_imp50_det_az_medianerr),
                            bmedian(dv400_imp25_det_az_medianerr)]

                
        y_non_median = [bmedian(dv400_imp1000_non_az),
                        bmedian(dv400_imp750_non_az),
                        bmedian(dv400_imp500_non_az),
                        bmedian(dv400_imp400_non_az), 
                        bmedian(dv400_imp300_non_az),
                        bmedian(dv400_imp200_non_az),
                        bmedian(dv400_imp100_non_az),
                        bmedian(dv400_imp50_non_az),
                        bmedian(dv400_imp25_non_az)]

        y_non_err_median = [bmedian(dv400_imp1000_non_az_medianerr),
                            bmedian(dv400_imp750_non_az_medianerr),
                            bmedian(dv400_imp500_non_az_medianerr),
                            bmedian(dv400_imp400_non_az_medianerr), 
                            bmedian(dv400_imp300_non_az_medianerr),
                            bmedian(dv400_imp200_non_az_medianerr),
                            bmedian(dv400_imp100_non_az_medianerr),
                            bmedian(dv400_imp50_non_az_medianerr),
                            bmedian(dv400_imp25_non_az_medianerr)]
                

        ###########

        # MEAN azlination for detections with errors
        ax1.errorbar(x,
                    y_det_mean,
                    yerr=y_det_err_mean,
                    marker=symbol_det_mean,
                    c=color_det_mean,
                    ms=markerSize,
                    markeredgecolor='black',
                    ls = ls_det_mean,
                    lw = lw,
                    alpha=alpha_det_mean,
                    label=label_det_mean)

        # MEDIAN azlination for detections with errors
        ax1.errorbar(x,
                    y_det_median,
                    yerr=y_det_err_median,
                    marker=symbol_det_median,
                    c=color_det_median,
                    ms=markerSize,
                    markeredgecolor='black',
                    ls = ls_det_median,
                    lw = lw,
                    alpha=alpha_det_median,
                    label=label_det_median)


        # MEAN azlination for non-detections with errors
        ax1.errorbar(x,
                    y_non_mean,
                    yerr=y_non_err_mean,
                    marker=symbol_non_mean,
                    c=color_non_mean,
                    ms=markerSize,
                    markeredgecolor='black',
                    ls = ls_non_mean,
                    lw = lw,
                    alpha=alpha_non_mean,
                    label=label_non_mean)

        # MEDIAN azlination for non-detections with errors
        ax1.errorbar(x,
                    y_non_median,
                    yerr=y_non_err_median,
                    marker=symbol_non_median,
                    c=color_non_median,
                    ms=markerSize,
                    markeredgecolor='black',
                    ls = ls_non_median,
                    lw = lw,
                    alpha=alpha_non_median,
                    label=label_non_median)
                

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
        
        ax1.set_ylabel(r'$\rm Azimuth~[deg]$')

        leg = ax1.legend(scatterpoints=1,prop={'size':12},loc='lower left',fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax1.grid(b=None,which='major',axis='both')
        ylim(0., 90.)
        xlim(0, 1000)

        if plot_detection_fraction_impact_az_save:
            savefig('{0}/detection_fraction_impact_az_errors_lstarcut{1}.pdf'.format(saveDirectory, lstar_cut),format='pdf',bbox_inches='tight')
        else:
            show()

##########################################################################################
##########################################################################################



##########################################################################################
##########################################################################################
    
    if plot_detection_fraction_like_az2:
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
        
        markerSize = 10
        lw = 2.
        
        binSize = 100
        bins = arange(0, 600, binSize)

        label_like_det = r'$\rm Detections$'
        label_like_non = r'$\rm Non-detections$'

        label_like_min200_det = r'$\rm Detections;~EW_{min} = 200$'
        label_like_min200_non = r'$\rm Non-detections;~EW_{min} = 200$'

        label_min200 = r'$\rm EW_{min} = 200; Median ~Az.$'


        ls_det_mean = 'solid'
        ls_non_mean = 'dashed'

        alpha_like_det = 0.99
        alpha_like_min200_det = 0.5
        alpha_like_non = 0.99
        alpha_like_min200_non = 0.5

        legend_size = 10
        legend_font = 10

        markerSize_like = 12
        markerSize_like_min200 = 6

        lw_like = 2.0
        lw_like_min200 = 0.8


        symbol_like_det = 'D'
        symbol_like_min200_det = 'D'
        symbol_like_non = 'X'
        symbol_like_min200_non = 'X'
        

        color_like_det = color_blue
        color_like_non = color_red
        color_like_min200_det = color_blue
        color_like_min200_non = color_red

        ls_like_det = 'solid'
        ls_like_non = 'dashed'
        ls_like_min200_det = 'solid'
        ls_like_min200_non = 'dashed'

##########################################################################################
        # do the plotting
        
        x_like = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 0.75]

########
        # now likelihoods
        y_like_det_median = [bmedian(dv400_l0001_det_az),
                                    bmedian(dv400_l0005_det_az),
                                    bmedian(dv400_l001_det_az),
                                    bmedian(dv400_l005_det_az), 
                                    bmedian(dv400_l01_det_az),
                                    bmedian(dv400_l05_det_az),
                                    bmedian(dv400_l1_det_az),
                                    bmedian(dv400_l5_det_az),
                                    bmedian(dv400_l75_det_az)]
                    
        # bootstrap errors in MEDIAN detection azlination
        y_like_det_err_median = [dv400_l0001_det_az_medianerr,
                                        dv400_l0005_det_az_medianerr,
                                        dv400_l001_det_az_medianerr,
                                        dv400_l005_det_az_medianerr,
                                        dv400_l01_det_az_medianerr,
                                        dv400_l05_det_az_medianerr,
                                        dv400_l1_det_az_medianerr,
                                        dv400_l5_det_az_medianerr,
                                        dv400_l75_det_az_medianerr]
                    
            
        y_like_non_median = [bmedian(dv400_l0001_non_az),
                                    bmedian(dv400_l0005_non_az),
                                    bmedian(dv400_l001_non_az),
                                    bmedian(dv400_l005_non_az), 
                                    bmedian(dv400_l01_non_az),
                                    bmedian(dv400_l05_non_az),
                                    bmedian(dv400_l1_non_az),
                                    bmedian(dv400_l5_non_az),
                                    bmedian(dv400_l75_non_az)]

        # bootstrap errors in MEDIAN non-detection azlination
        y_like_non_err_median = [dv400_l0001_non_az_medianerr,
                                        dv400_l0005_non_az_medianerr,
                                        dv400_l001_non_az_medianerr,
                                        dv400_l005_non_az_medianerr,
                                        dv400_l01_non_az_medianerr,
                                        dv400_l05_non_az_medianerr,
                                        dv400_l1_non_az_medianerr,
                                        dv400_l5_non_az_medianerr,
                                        dv400_l75_non_az_medianerr]
        
        
        # minEW = 200; likelihoods
        y_like_min200_det_median = [bmedian(dv400_l0001_min200_det_az),
                                    bmedian(dv400_l0005_min200_det_az),
                                    bmedian(dv400_l001_min200_det_az),
                                    bmedian(dv400_l005_min200_det_az), 
                                    bmedian(dv400_l01_min200_det_az),
                                    bmedian(dv400_l05_min200_det_az),
                                    bmedian(dv400_l1_min200_det_az),
                                    bmedian(dv400_l5_min200_det_az),
                                    bmedian(dv400_l75_min200_det_az)]
                    
        # bootstrap errors in MEDIAN detection azlination
        y_like_min200_det_err_median = [dv400_l0001_min200_det_az_medianerr,
                                        dv400_l0005_min200_det_az_medianerr,
                                        dv400_l001_min200_det_az_medianerr,
                                        dv400_l005_min200_det_az_medianerr,
                                        dv400_l01_min200_det_az_medianerr,
                                        dv400_l05_min200_det_az_medianerr,
                                        dv400_l1_min200_det_az_medianerr,
                                        dv400_l5_min200_det_az_medianerr,
                                        dv400_l75_min200_det_az_medianerr]
                    
            
        y_like_min200_non_median = [bmedian(dv400_l0001_min200_non_az),
                                    bmedian(dv400_l0005_min200_non_az),
                                    bmedian(dv400_l001_min200_non_az),
                                    bmedian(dv400_l005_min200_non_az), 
                                    bmedian(dv400_l01_min200_non_az),
                                    bmedian(dv400_l05_min200_non_az),
                                    bmedian(dv400_l1_min200_non_az),
                                    bmedian(dv400_l5_min200_non_az),
                                    bmedian(dv400_l75_min200_non_az)]

        # bootstrap errors in MEDIAN non-detection azlination
        y_like_min200_non_err_median = [dv400_l0001_min200_non_az_medianerr,
                                        dv400_l0005_min200_non_az_medianerr,
                                        dv400_l001_min200_non_az_medianerr,
                                        dv400_l005_min200_non_az_medianerr,
                                        dv400_l01_min200_non_az_medianerr,
                                        dv400_l05_min200_non_az_medianerr,
                                        dv400_l1_min200_non_az_medianerr,
                                        dv400_l5_min200_non_az_medianerr,
                                        dv400_l75_min200_non_az_medianerr]


        
        # MEDIAN azlination for detections with errors
        ax1.errorbar(x_like,
                    y_like_det_median,
                    yerr=y_like_det_err_median,
                    marker=symbol_like_det,
                    c=color_like_det,
                    ms=markerSize_like,
                    markeredgecolor='black',
                    ls = ls_like_det,
                    lw = lw_like,
                    alpha=alpha_like_det,
                    label=label_like_det)
        
        
        # MEDIAN azlination for non-detections with errors
        ax1.errorbar(x_like,
                    y_like_non_median,
                    yerr=y_like_non_err_median,
                    marker=symbol_like_non,
                    c=color_like_non,
                    ms=markerSize_like,
                    markeredgecolor='black',
                    ls = ls_like_non,
                    lw = lw_like,
                    alpha=alpha_like_non,
                    label=label_like_non)
                    
        #
        # minEW = 200
        # MEDIAN azlination for detections with errors
        ax1.errorbar(x_like,
                    y_like_min200_det_median,
                    yerr=y_like_min200_det_err_median,
                    marker=symbol_like_min200_det,
                    c=color_like_min200_det,
                    ms=markerSize_like_min200,
                    markeredgecolor='black',
                    ls = ls_like_min200_det,
                    lw = lw_like_min200,
                    alpha=alpha_like_min200_det,
                    label=label_like_min200_det)
        
        
        # MEDIAN azlination for non-detections with errors
        ax1.errorbar(x_like,
                    y_like_min200_non_median,
                    yerr=y_like_min200_non_err_median,
                    marker=symbol_like_min200_non,
                    c=color_like_min200_non,
                    ms=markerSize_like_min200,
                    markeredgecolor='black',
                    ls = ls_like_min200_non,
                    lw = lw_like_min200,
                    alpha=alpha_like_min200_non,
                    label=label_like_min200_non)
                    
                
        ax1.set_xlabel(r'$\rm \mathcal{L}$')
        ax1.set_xscale("log")
        ax1.invert_xaxis()

        # y-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax1.yaxis.set_major_locator(majorLocator)
        ax1.yaxis.set_major_formatter(majorFormatter)
        ax1.yaxis.set_minor_locator(minorLocator)
        
        ax1.set_ylabel(r'$\rm Azimuth~[deg]$')

#         leg = ax1.legend(scatterpoints=1,prop={'size':12},loc='lower left',fancybox=True)
#         leg.get_frame().set_alpha(0.5)



        try:
            import matplotlib.patches as mpatches
            import matplotlib.lines as mlines
                              
            likelihood_det = mlines.Line2D([], [], color=color_like_det, marker=symbol_like_det,lw=lw_like,
                                    markeredgecolor='black', markersize=legend_size, label=label_like_det)

                              
            likelihood_non = mlines.Line2D([], [], color=color_like_non, marker=symbol_like_non,lw=lw_like,
                                    markeredgecolor='black', markersize=legend_size, label=label_like_non)
                                    
            likelihood_min200_det = mlines.Line2D([], [], color=color_like_min200_det, marker=symbol_like_min200_det,
                                    lw=lw_like_min200, markeredgecolor='black', markersize=markerSize_like_min200, 
                                    alpha = alpha_like_min200_det,label=label_like_min200_det)

                              
            likelihood_min200_non = mlines.Line2D([], [], color=color_like_min200_non, marker=symbol_like_min200_non,
                                    lw=lw_like_min200, markeredgecolor='black', markersize=markerSize_like_min200, 
                                    alpha = alpha_like_min200_det, label=label_like_min200_non)
                                    
                                

    #         min200 = mlines.Line2D([], [], color=color_impact_min200, marker=None,lw=lw_min200,
    #                                   markersize=legend_size, markeredgecolor='black', label=label_min200)
                                  
            plt.legend(handles=[likelihood_det, likelihood_non, likelihood_min200_det, likelihood_min200_non],loc='lower right',
                                    borderpad=0.8, fontsize=legend_font, fancybox=True)
        except Exception, e:
            print
            print "Legend issue: ",e
            print
        

        ax1.grid(b=None,which='major',axis='both')
        ylim(10., 90.)
#         xlim(0, 1000)

        if plot_detection_fraction_like_az_save2:
            savefig('{0}/detection_fraction_like_az_errors_minEW0_200.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()




##########################################################################################
##########################################################################################
    
    if plot_detection_fraction_both_az:
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
        

        alpha_det_mean = 0.9
        alpha_non_mean = 0.9
        alpha_det_median = 0.6
        alpha_non_median = 0.6
        markerSize = 10
        lw = 2.
        
        binSize = 100
        bins = arange(0, 600, binSize)

        label_like_det = r'$\rm \mathcal{L} - Detections; Median ~Az.$'
        label_imp_det = r'$\rm \rho - Detections; Median ~Az.$'
        label_like_non = r'$\rm \mathcal{L} - Non-detections; Median ~Az.$'
        label_imp_non = r'$\rm \rho - Non-detection; Median ~Az.$'

        label_imp_min200_det = r'$\rm EW_{min} = 200; Median ~Az.$'
        label_imp_min200_non = r'$\rm EW_{min} = 200; Median ~Az.$'
        label_like_min200_det = r'$\rm EW_{min} = 200; Median ~Az.$'
        label_like_min200_non = r'$\rm EW_{min} = 200; Median ~Az.$'

        label_min200 = r'$\rm EW_{min} = 200; Median ~Az.$'
        
        symbol_det_median = 'D'
        symbol_non_median = 'X'
        
        color_det_mean = color_blue
        color_non_mean = color_red
        color_det_median = color_purple2
        color_non_median = color_orange

        ls_det_mean = 'solid'
        ls_non_mean = 'dashed'

        alpha_like_det = 0.99
        alpha_like_min200_det = 0.5
        alpha_like_non = 0.99
        alpha_like_min200_non = 0.5

        alpha_imp_det = 0.99
        alpha_imp_min200_det = 0.5
        alpha_imp_non = 0.99
        alpha_imp_min200_non = 0.5

        legend_size = 10
        legend_font = 10

        markerSize_like = 12
        markerSize_imp = 12
        markerSize_like_min200 = 6
        markerSize_imp_min200 = 6

        lw_like = 2.0
        lw_imp = 2.0
        lw_like_min200 = 0.8
        lw_imp_min200 = 0.8


        symbol_like_det = 'D'
        symbol_like_min200_det = 'D'
        symbol_like_non = 'X'
        symbol_like_min200_non = 'X'
        
        symbol_imp_det = 'o'
        symbol_imp_min200_det = 'o'
        symbol_imp_non = 'P'
        symbol_imp_min200_non = 'P'

        color_like_det = color_blue
        color_like_non = color_red
        color_like_min200_det = color_blue
        color_like_min200_non = color_red

        color_imp_det = color_purple2
        color_imp_non = color_orange
        color_imp_min200_det = color_purple2
        color_imp_min200_non = color_orange

        ls_like_det = 'solid'
        ls_like_non = 'dashed'
        ls_like_min200_det = 'solid'
        ls_like_min200_non = 'dashed'

        ls_imp_det = 'solid'
        ls_imp_non = 'dashed'
        ls_imp_min200_det = 'solid'
        ls_imp_min200_non = 'dashed'


##########################################################################################
        # do the plotting
        
        print 'len(dv400_imp100_det_az) : ',len(dv400_imp100_det_az)
        print 'len(dv400_imp50_det_az) : ',len(dv400_imp50_det_az)
        print 'len(dv400_imp25_det_az) : ',len(dv400_imp25_det_az)
        print
        print
        print 'len(dv400_imp100_non_az) : ',len(dv400_imp100_non_az)
        print 'len(dv400_imp50_non_az) : ',len(dv400_imp50_non_az)
        print 'len(dv400_imp25_non_az) : ',len(dv400_imp25_non_az)
        print


        x_imp = [1000, 750, 500, 400, 300, 200, 100, 50, 25]
        x_like = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 0.75]
        
        y_imp_det_median = [bmedian(dv400_imp1000_det_az),
                        bmedian(dv400_imp750_det_az),
                        bmedian(dv400_imp500_det_az),
                        bmedian(dv400_imp400_det_az),
                        bmedian(dv400_imp300_det_az),
                        bmedian(dv400_imp200_det_az),
                        bmedian(dv400_imp100_det_az),
                        bmedian(dv400_imp50_det_az),
                        bmedian(dv400_imp25_det_az)]
            
        y_imp_det_err_median = [bmedian(dv400_imp1000_det_az_medianerr),
                            bmedian(dv400_imp750_det_az_medianerr),
                            bmedian(dv400_imp500_det_az_medianerr),
                            bmedian(dv400_imp400_det_az_medianerr),
                            bmedian(dv400_imp300_det_az_medianerr),
                            bmedian(dv400_imp200_det_az_medianerr),
                            bmedian(dv400_imp100_det_az_medianerr),
                            bmedian(dv400_imp50_det_az_medianerr),
                            bmedian(dv400_imp25_det_az_medianerr)]

                
        y_imp_non_median = [bmedian(dv400_imp1000_non_az),
                        bmedian(dv400_imp750_non_az),
                        bmedian(dv400_imp500_non_az),
                        bmedian(dv400_imp400_non_az), 
                        bmedian(dv400_imp300_non_az),
                        bmedian(dv400_imp200_non_az),
                        bmedian(dv400_imp100_non_az),
                        bmedian(dv400_imp50_non_az),
                        bmedian(dv400_imp25_non_az)]

        y_imp_non_err_median = [bmedian(dv400_imp1000_non_az_medianerr),
                            bmedian(dv400_imp750_non_az_medianerr),
                            bmedian(dv400_imp500_non_az_medianerr),
                            bmedian(dv400_imp400_non_az_medianerr), 
                            bmedian(dv400_imp300_non_az_medianerr),
                            bmedian(dv400_imp200_non_az_medianerr),
                            bmedian(dv400_imp100_non_az_medianerr),
                            bmedian(dv400_imp50_non_az_medianerr),
                            bmedian(dv400_imp25_non_az_medianerr)]
                
###############
        y_imp_min200_det_median = [bmedian(dv400_imp1000_min200_det_az),
                                    bmedian(dv400_imp750_min200_det_az),
                                    bmedian(dv400_imp500_min200_det_az),
                                    bmedian(dv400_imp400_min200_det_az),
                                    bmedian(dv400_imp300_min200_det_az),
                                    bmedian(dv400_imp200_min200_det_az),
                                    bmedian(dv400_imp100_min200_det_az),
                                    bmedian(dv400_imp50_min200_det_az),
                                    bmedian(dv400_imp25_min200_det_az)]
            
        y_imp_min200_det_err_median = [bmedian(dv400_imp1000_min200_det_az_medianerr),
                                        bmedian(dv400_imp750_min200_det_az_medianerr),
                                        bmedian(dv400_imp500_min200_det_az_medianerr),
                                        bmedian(dv400_imp400_min200_det_az_medianerr),
                                        bmedian(dv400_imp300_min200_det_az_medianerr),
                                        bmedian(dv400_imp200_min200_det_az_medianerr),
                                        bmedian(dv400_imp100_min200_det_az_medianerr),
                                        bmedian(dv400_imp50_min200_det_az_medianerr),
                                        bmedian(dv400_imp25_min200_det_az_medianerr)]

                
        y_imp_min200_non_median = [bmedian(dv400_imp1000_min200_non_az),
                                    bmedian(dv400_imp750_min200_non_az),
                                    bmedian(dv400_imp500_min200_non_az),
                                    bmedian(dv400_imp400_min200_non_az), 
                                    bmedian(dv400_imp300_min200_non_az),
                                    bmedian(dv400_imp200_min200_non_az),
                                    bmedian(dv400_imp100_min200_non_az),
                                    bmedian(dv400_imp50_min200_non_az),
                                    bmedian(dv400_imp25_min200_non_az)]

        y_imp_min200_non_err_median = [bmedian(dv400_imp1000_min200_non_az_medianerr),
                                        bmedian(dv400_imp750_min200_non_az_medianerr),
                                        bmedian(dv400_imp500_min200_non_az_medianerr),
                                        bmedian(dv400_imp400_min200_non_az_medianerr), 
                                        bmedian(dv400_imp300_min200_non_az_medianerr),
                                        bmedian(dv400_imp200_min200_non_az_medianerr),
                                        bmedian(dv400_imp100_min200_non_az_medianerr),
                                        bmedian(dv400_imp50_min200_non_az_medianerr),
                                        bmedian(dv400_imp25_min200_non_az_medianerr)]


########
        # now likelihoods
        y_like_det_median = [bmedian(dv400_l0001_det_az),
                                    bmedian(dv400_l0005_det_az),
                                    bmedian(dv400_l001_det_az),
                                    bmedian(dv400_l005_det_az), 
                                    bmedian(dv400_l01_det_az),
                                    bmedian(dv400_l05_det_az),
                                    bmedian(dv400_l1_det_az),
                                    bmedian(dv400_l5_det_az),
                                    bmedian(dv400_l75_det_az)]
                    
        # bootstrap errors in MEDIAN detection azlination
        y_like_det_err_median = [dv400_l0001_det_az_medianerr,
                                        dv400_l0005_det_az_medianerr,
                                        dv400_l001_det_az_medianerr,
                                        dv400_l005_det_az_medianerr,
                                        dv400_l01_det_az_medianerr,
                                        dv400_l05_det_az_medianerr,
                                        dv400_l1_det_az_medianerr,
                                        dv400_l5_det_az_medianerr,
                                        dv400_l75_det_az_medianerr]
                    
            
        y_like_non_median = [bmedian(dv400_l0001_non_az),
                                    bmedian(dv400_l0005_non_az),
                                    bmedian(dv400_l001_non_az),
                                    bmedian(dv400_l005_non_az), 
                                    bmedian(dv400_l01_non_az),
                                    bmedian(dv400_l05_non_az),
                                    bmedian(dv400_l1_non_az),
                                    bmedian(dv400_l5_non_az),
                                    bmedian(dv400_l75_non_az)]

        # bootstrap errors in MEDIAN non-detection azlination
        y_like_non_err_median = [dv400_l0001_non_az_medianerr,
                                        dv400_l0005_non_az_medianerr,
                                        dv400_l001_non_az_medianerr,
                                        dv400_l005_non_az_medianerr,
                                        dv400_l01_non_az_medianerr,
                                        dv400_l05_non_az_medianerr,
                                        dv400_l1_non_az_medianerr,
                                        dv400_l5_non_az_medianerr,
                                        dv400_l75_non_az_medianerr]
        
        
        # minEW = 200; likelihoods
        y_like_min200_det_median = [bmedian(dv400_l0001_min200_det_az),
                                    bmedian(dv400_l0005_min200_det_az),
                                    bmedian(dv400_l001_min200_det_az),
                                    bmedian(dv400_l005_min200_det_az), 
                                    bmedian(dv400_l01_min200_det_az),
                                    bmedian(dv400_l05_min200_det_az),
                                    bmedian(dv400_l1_min200_det_az),
                                    bmedian(dv400_l5_min200_det_az),
                                    bmedian(dv400_l75_min200_det_az)]
                    
        # bootstrap errors in MEDIAN detection azlination
        y_like_min200_det_err_median = [dv400_l0001_min200_det_az_medianerr,
                                        dv400_l0005_min200_det_az_medianerr,
                                        dv400_l001_min200_det_az_medianerr,
                                        dv400_l005_min200_det_az_medianerr,
                                        dv400_l01_min200_det_az_medianerr,
                                        dv400_l05_min200_det_az_medianerr,
                                        dv400_l1_min200_det_az_medianerr,
                                        dv400_l5_min200_det_az_medianerr,
                                        dv400_l75_min200_det_az_medianerr]
                    
            
        y_like_min200_non_median = [bmedian(dv400_l0001_min200_non_az),
                                    bmedian(dv400_l0005_min200_non_az),
                                    bmedian(dv400_l001_min200_non_az),
                                    bmedian(dv400_l005_min200_non_az), 
                                    bmedian(dv400_l01_min200_non_az),
                                    bmedian(dv400_l05_min200_non_az),
                                    bmedian(dv400_l1_min200_non_az),
                                    bmedian(dv400_l5_min200_non_az),
                                    bmedian(dv400_l75_min200_non_az)]

        # bootstrap errors in MEDIAN non-detection azlination
        y_like_min200_non_err_median = [dv400_l0001_min200_non_az_medianerr,
                                        dv400_l0005_min200_non_az_medianerr,
                                        dv400_l001_min200_non_az_medianerr,
                                        dv400_l005_min200_non_az_medianerr,
                                        dv400_l01_min200_non_az_medianerr,
                                        dv400_l05_min200_non_az_medianerr,
                                        dv400_l1_min200_non_az_medianerr,
                                        dv400_l5_min200_non_az_medianerr,
                                        dv400_l75_min200_non_az_medianerr]


        ##########
        # impact axis
        # MEDIAN azlination for detections with errors
        ax1.errorbar(x_imp,
                    y_imp_det_median,
                    yerr=y_imp_det_err_median,
                    marker=symbol_imp_det,
                    c=color_imp_det,
                    ms=markerSize_imp,
                    markeredgecolor='black',
                    ls = ls_imp_det,
                    lw = lw_imp,
                    alpha=alpha_imp_det,
                    label=label_imp_det)

        # MEDIAN azlination for non-detections with errors
        ax1.errorbar(x_imp,
                    y_imp_non_median,
                    yerr=y_imp_non_err_median,
                    marker=symbol_imp_non,
                    c=color_imp_non,
                    ms=markerSize_imp,
                    markeredgecolor='black',
                    ls = ls_imp_non,
                    lw = lw_imp,
                    alpha=alpha_imp_non,
                    label=label_imp_non)

        #
        # now minEW = 200
        # MEDIAN azlination for detections with errors
        ax1.errorbar(x_imp,
                    y_imp_min200_det_median,
                    yerr=y_imp_min200_det_err_median,
                    marker=symbol_imp_min200_det,
                    c=color_imp_min200_det,
                    ms=markerSize_imp_min200,
                    markeredgecolor='black',
                    ls = ls_imp_min200_det,
                    lw = lw_imp_min200,
                    alpha=alpha_imp_min200_det,
                    label=label_imp_min200_det)

        # MEDIAN azlination for non-detections with errors
        ax1.errorbar(x_imp,
                    y_imp_min200_non_median,
                    yerr=y_imp_min200_non_err_median,
                    marker=symbol_imp_min200_non,
                    c=color_imp_min200_non,
                    ms=markerSize_imp_min200,
                    markeredgecolor='black',
                    ls = ls_imp_min200_non,
                    lw = lw_imp_min200,
                    alpha=alpha_imp_min200_non,
                    label=label_imp_min200_non)


        ax1.set_xlabel(r'$\rm \rho~[kpc]$')
        xlim(0, 1000)

        # x-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)
        ax1.xaxis.set_major_locator(majorLocator)
        ax1.xaxis.set_major_formatter(majorFormatter)
        ax1.xaxis.set_minor_locator(minorLocator)

        # now do the likelihood axis
        # share a y-axis, have different top and bottom x-scales
        ax2 = ax1.twiny()
        
        # MEDIAN azlination for detections with errors
        ax2.errorbar(x_like,
                    y_like_det_median,
                    yerr=y_like_det_err_median,
                    marker=symbol_like_det,
                    c=color_like_det,
                    ms=markerSize_like,
                    markeredgecolor='black',
                    ls = ls_like_det,
                    lw = lw_like,
                    alpha=alpha_like_det,
                    label=label_like_det)
        
        
        # MEDIAN azlination for non-detections with errors
        ax2.errorbar(x_like,
                    y_like_non_median,
                    yerr=y_like_non_err_median,
                    marker=symbol_like_non,
                    c=color_like_non,
                    ms=markerSize_like,
                    markeredgecolor='black',
                    ls = ls_like_non,
                    lw = lw_like,
                    alpha=alpha_like_non,
                    label=label_like_non)
                    
        #
        # minEW = 200
        # MEDIAN azlination for detections with errors
        ax2.errorbar(x_like,
                    y_like_min200_det_median,
                    yerr=y_like_min200_det_err_median,
                    marker=symbol_like_min200_det,
                    c=color_like_min200_det,
                    ms=markerSize_like_min200,
                    markeredgecolor='black',
                    ls = ls_like_min200_det,
                    lw = lw_like_min200,
                    alpha=alpha_like_min200_det,
                    label=label_like_min200_det)
        
        
        # MEDIAN azlination for non-detections with errors
        ax2.errorbar(x_like,
                    y_like_min200_non_median,
                    yerr=y_like_min200_non_err_median,
                    marker=symbol_like_min200_non,
                    c=color_like_min200_non,
                    ms=markerSize_like_min200,
                    markeredgecolor='black',
                    ls = ls_like_min200_non,
                    lw = lw_like_min200,
                    alpha=alpha_like_min200_non,
                    label=label_like_min200_non)
                    
                
        ax2.set_xlabel(r'$\rm \mathcal{L}$')
        ax2.set_xscale("log")
        ax2.invert_xaxis()

        # y-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax1.yaxis.set_major_locator(majorLocator)
        ax1.yaxis.set_major_formatter(majorFormatter)
        ax1.yaxis.set_minor_locator(minorLocator)
        
        ax1.set_ylabel(r'$\rm Azimuth~[deg]$')

#         leg = ax1.legend(scatterpoints=1,prop={'size':12},loc='lower left',fancybox=True)
#         leg.get_frame().set_alpha(0.5)



        try:
            import matplotlib.patches as mpatches
            import matplotlib.lines as mlines
                              
            impact_det = mlines.Line2D([], [], color=color_imp_det, marker=symbol_imp_det,lw=lw_imp,
                                      markersize=legend_size, markeredgecolor='black', label=label_imp_det)
                              
            likelihood_det = mlines.Line2D([], [], color=color_like_det, marker=symbol_like_det,lw=lw_like,
                                    markeredgecolor='black', markersize=legend_size, label=label_like_det)
                                
            impact_non = mlines.Line2D([], [], color=color_imp_non, marker=symbol_imp_non,lw=lw_imp,
                                      markersize=legend_size, markeredgecolor='black', label=label_imp_non)
                              
            likelihood_non = mlines.Line2D([], [], color=color_like_non, marker=symbol_like_non,lw=lw_like,
                                    markeredgecolor='black', markersize=legend_size, label=label_like_det)
                                

    #         min200 = mlines.Line2D([], [], color=color_impact_min200, marker=None,lw=lw_min200,
    #                                   markersize=legend_size, markeredgecolor='black', label=label_min200)
                                  
            plt.legend(handles=[impact_det, likelihood_det, impact_non, likelihood_non],loc='lower right',
                                    borderpad=0.8, fontsize=legend_font, fancybox=True)
        except Exception, e:
            print
            print "Legend issue: ",e
            print
        

        ax1.grid(b=None,which='major',axis='both')
        ylim(10., 90.)
#         xlim(0, 1000)

        if plot_detection_fraction_both_az_save:
            savefig('{0}/detection_fraction_both_az_errors_minEW0_200.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()




##########################################################################################
##########################################################################################
    
    if plot_detection_fraction_vs_az_like:
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
        

        label_l1_outside = r'$\rm \mathcal{L} < 0.1; \rho/R_{vir} > 1.5$'
        label_l1_inside= r'$\rm \mathcal{L} \ge 0.1; \rho/R_{vir} \leq 1.5$'

        alpha_l1_outside = 0.99
        alpha_l1_inside = 0.99


        legend_size = 12
        legend_font = 12

        markerSize_l1_outside = 12
        markerSize_l1_inside = 12

        lw_l1_outside = 1.5
        lw_l1_inside = 1.5

        symbol_l1_outside = 'P'
        symbol_l1_inside = 'D'

        color_l1_outside = color_orange
        color_l1_inside = color_blue

        ls_l1_outside = 'dashed'
        ls_l1_inside = 'solid'


##########################################################################################
        # do the plotting
        
        print 'len(dv400_imp100_det_az) : ',len(dv400_imp100_det_az)
        print 'len(dv400_imp50_det_az) : ',len(dv400_imp50_det_az)
        print 'len(dv400_imp25_det_az) : ',len(dv400_imp25_det_az)
        print
        print
        print 'len(dv400_imp100_non_az) : ',len(dv400_imp100_non_az)
        print 'len(dv400_imp50_non_az) : ',len(dv400_imp50_non_az)
        print 'len(dv400_imp25_non_az) : ',len(dv400_imp25_non_az)
        print


        x = np.arange(7.5, 97.5, 15)
        
        # L = 0.01 ~ 2.15 rho/rvir
        # L = 0.05 ~ 1.73 rho/rvir
        # L = 0.1 ~ 1.5 rho/rvir
        # L = 0.5 ~ 0.83 rho/rvir
        # L = 0.75 ~ 0.54 rho/vir

#######
        l001_det_15 = []
        l001_det_30 = []
        l001_det_45 = []
        l001_det_60 = []
        l001_det_75 = []
        l001_det_90 = []

        l001_non_15 = []
        l001_non_30 = []
        l001_non_45 = []
        l001_non_60 = []
        l001_non_75 = []
        l001_non_90 = []

        l005_det_15 = []
        l005_det_30 = []
        l005_det_45 = []
        l005_det_60 = []
        l005_det_75 = []
        l005_det_90 = []

        l005_non_15 = []
        l005_non_30 = []
        l005_non_45 = []
        l005_non_60 = []
        l005_non_75 = []
        l005_non_90 = []

        l01_det_15 = []
        l01_det_30 = []
        l01_det_45 = []
        l01_det_60 = []
        l01_det_75 = []
        l01_det_90 = []

        l01_non_15 = []
        l01_non_30 = []
        l01_non_45 = []
        l01_non_60 = []
        l01_non_75 = []
        l01_non_90 = []
        
        l05_det_15 = []
        l05_det_30 = []
        l05_det_45 = []
        l05_det_60 = []
        l05_det_75 = []
        l05_det_90 = []

        l05_non_15 = []
        l05_non_30 = []
        l05_non_45 = []
        l05_non_60 = []
        l05_non_75 = []
        l05_non_90 = []
        
        l1_det_15 = []
        l1_det_30 = []
        l1_det_45 = []
        l1_det_60 = []
        l1_det_75 = []
        l1_det_90 = []

        l1_non_15 = []
        l1_non_30 = []
        l1_non_45 = []
        l1_non_60 = []
        l1_non_75 = []
        l1_non_90 = []
        
        l5_det_15 = []
        l5_det_30 = []
        l5_det_45 = []
        l5_det_60 = []
        l5_det_75 = []
        l5_det_90 = []

        l5_non_15 = []
        l5_non_30 = []
        l5_non_45 = []
        l5_non_60 = []
        l5_non_75 = []
        l5_non_90 = []
        
        l75_det_15 = []
        l75_det_30 = []
        l75_det_45 = []
        l75_det_60 = []
        l75_det_75 = []
        l75_det_90 = []

        l75_non_15 = []
        l75_non_30 = []
        l75_non_45 = []
        l75_non_60 = []
        l75_non_75 = []
        l75_non_90 = []
        



        for az in dv400_l001_det_az:
            az = float(az)
            if az <= 15:
                l001_det_15.append(az)

            elif az > 15 and az <= 30:
                l001_det_30.append(az)

            elif az > 30 and az <= 45:
                l001_det_45.append(az)

            elif az > 45 and az <= 60:
                l001_det_60.append(az)

            elif az > 60 and az <= 75:
                l001_det_75.append(az)

            elif az > 75 and az <= 90:
                l001_det_90.append(az)
            else:
                print 'wtf? l001_det_az = ',az
            
        for az in dv400_l001_non_az:
            az = float(az)
            if az <= 15:
                l001_non_15.append(az)

            elif az > 15 and az <= 30:
                l001_non_30.append(az)

            elif az > 30 and az <= 45:
                l001_non_45.append(az)

            elif az > 45 and az <= 60:
                l001_non_60.append(az)

            elif az > 60 and az <= 75:
                l001_non_75.append(az)

            elif az > 75 and az <= 90:
                l001_non_90.append(az)
            else:
                print 'wtf? l001_non_az = ',az
        
        
        for az in dv400_l005_det_az:
            az = float(az)
            if az <= 15:
                l005_det_15.append(az)

            elif az > 15 and az <= 30:
                l005_det_30.append(az)

            elif az > 30 and az <= 45:
                l005_det_45.append(az)

            elif az > 45 and az <= 60:
                l005_det_60.append(az)

            elif az > 60 and az <= 75:
                l005_det_75.append(az)

            elif az > 75 and az <= 90:
                l005_det_90.append(az)
            else:
                print 'wtf? l005_det_az = ',az
            
        for az in dv400_l005_non_az:
            az = float(az)
            if az <= 15:
                l005_non_15.append(az)

            elif az > 15 and az <= 30:
                l005_non_30.append(az)

            elif az > 30 and az <= 45:
                l005_non_45.append(az)

            elif az > 45 and az <= 60:
                l005_non_60.append(az)

            elif az > 60 and az <= 75:
                l005_non_75.append(az)

            elif az > 75 and az <= 90:
                l005_non_90.append(az)
            else:
                print 'wtf? l005_non_az = ',az
        
        
        for az in dv400_l01_det_az:
            az = float(az)
            if az <= 15:
                l01_det_15.append(az)

            elif az > 15 and az <= 30:
                l01_det_30.append(az)

            elif az > 30 and az <= 45:
                l01_det_45.append(az)

            elif az > 45 and az <= 60:
                l01_det_60.append(az)

            elif az > 60 and az <= 75:
                l01_det_75.append(az)

            elif az > 75 and az <= 90:
                l01_det_90.append(az)
            else:
                print 'wtf? l01_det_az = ',az
            
        for az in dv400_l01_non_az:
            az = float(az)
            if az <= 15:
                l01_non_15.append(az)

            elif az > 15 and az <= 30:
                l01_non_30.append(az)

            elif az > 30 and az <= 45:
                l01_non_45.append(az)

            elif az > 45 and az <= 60:
                l01_non_60.append(az)

            elif az > 60 and az <= 75:
                l01_non_75.append(az)

            elif az > 75 and az <= 90:
                l01_non_90.append(az)
            else:
                print 'wtf? l01_non_az = ',az

        for az in dv400_l05_det_az:
            az = float(az)
            if az <= 15:
                l05_det_15.append(az)

            elif az > 15 and az <= 30:
                l05_det_30.append(az)

            elif az > 30 and az <= 45:
                l05_det_45.append(az)

            elif az > 45 and az <= 60:
                l05_det_60.append(az)

            elif az > 60 and az <= 75:
                l05_det_75.append(az)

            elif az > 75 and az <= 90:
                l05_det_90.append(az)
            else:
                print 'wtf? l05_det_az = ',az


        for az in dv400_l05_non_az:
            az = float(az)
            if az <= 15:
                l05_non_15.append(az)

            elif az > 15 and az <= 30:
                l05_non_30.append(az)

            elif az > 30 and az <= 45:
                l05_non_45.append(az)

            elif az > 45 and az <= 60:
                l05_non_60.append(az)

            elif az > 60 and az <= 75:
                l05_non_75.append(az)

            elif az > 75 and az <= 90:
                l05_non_90.append(az)
            else:
                print 'wtf? l05_non_az = ',az


        for az in dv400_l1_det_az:
            az = float(az)
            if az <= 15:
                l1_det_15.append(az)

            elif az > 15 and az <= 30:
                l1_det_30.append(az)

            elif az > 30 and az <= 45:
                l1_det_45.append(az)

            elif az > 45 and az <= 60:
                l1_det_60.append(az)

            elif az > 60 and az <= 75:
                l1_det_75.append(az)

            elif az > 75 and az <= 90:
                l1_det_90.append(az)
            else:
                print 'wtf? l1_det_az = ',az

        for az in dv400_l1_non_az:
            az = float(az)
            if az <= 15:
                l1_non_15.append(az)

            elif az > 15 and az <= 30:
                l1_non_30.append(az)

            elif az > 30 and az <= 45:
                l1_non_45.append(az)

            elif az > 45 and az <= 60:
                l1_non_60.append(az)

            elif az > 60 and az <= 75:
                l1_non_75.append(az)

            elif az > 75 and az <= 90:
                l1_non_90.append(az)
            else:
                print 'wtf? l1_non_az = ',az

        for az in dv400_l5_det_az:
            az = float(az)
            if az <= 15:
                l5_det_15.append(az)

            elif az > 15 and az <= 30:
                l5_det_30.append(az)

            elif az > 30 and az <= 45:
                l5_det_45.append(az)

            elif az > 45 and az <= 60:
                l5_det_60.append(az)

            elif az > 60 and az <= 75:
                l5_det_75.append(az)

            elif az > 75 and az <= 90:
                l5_det_90.append(az)
            else:
                print 'wtf? l5_det_az = ',az

        for az in dv400_l5_non_az:
            az = float(az)
            if az <= 15:
                l5_non_15.append(az)

            elif az > 15 and az <= 30:
                l5_non_30.append(az)

            elif az > 30 and az <= 45:
                l5_non_45.append(az)

            elif az > 45 and az <= 60:
                l5_non_60.append(az)

            elif az > 60 and az <= 75:
                l5_non_75.append(az)

            elif az > 75 and az <= 90:
                l5_non_90.append(az)
            else:
                print 'wtf? l5_non_az = ',az


        for az in dv400_l75_det_az:
            az = float(az)
            if az <= 15:
                l75_det_15.append(az)

            elif az > 15 and az <= 30:
                l75_det_30.append(az)

            elif az > 30 and az <= 45:
                l75_det_45.append(az)

            elif az > 45 and az <= 60:
                l75_det_60.append(az)

            elif az > 60 and az <= 75:
                l75_det_75.append(az)

            elif az > 75 and az <= 90:
                l75_det_90.append(az)
            else:
                print 'wtf? l75_det_az = ',az

        for az in dv400_l75_non_az:
            az = float(az)
            if az <= 15:
                l75_non_15.append(az)

            elif az > 15 and az <= 30:
                l75_non_30.append(az)

            elif az > 30 and az <= 45:
                l75_non_45.append(az)

            elif az > 45 and az <= 60:
                l75_non_60.append(az)

            elif az > 60 and az <= 75:
                l75_non_75.append(az)

            elif az > 75 and az <= 90:
                l75_non_90.append(az)
            else:
                print 'wtf? l75_non_az = ',az


        y_l001_det = [len(l001_det_15),
                     len(l001_det_30),
                     len(l001_det_45),
                     len(l001_det_60),
                     len(l001_det_75),
                     len(l001_det_90)]
                     
        y_l001_non = [len(l001_non_15),
                     len(l001_non_30),
                     len(l001_non_45),
                     len(l001_non_60),
                     len(l001_non_75),
                     len(l001_non_90)]

        y_l005_det = [len(l005_det_15),
                     len(l005_det_30),
                     len(l005_det_45),
                     len(l005_det_60),
                     len(l005_det_75),
                     len(l005_det_90)]
                     
        y_l005_non = [len(l005_non_15),
                     len(l005_non_30),
                     len(l005_non_45),
                     len(l005_non_60),
                     len(l005_non_75),
                     len(l005_non_90)]

        y_l01_det = [len(l01_det_15),
                     len(l01_det_30),
                     len(l01_det_45),
                     len(l01_det_60),
                     len(l01_det_75),
                     len(l01_det_90)]
                     
        y_l01_non = [len(l01_non_15),
                     len(l01_non_30),
                     len(l01_non_45),
                     len(l01_non_60),
                     len(l01_non_75),
                     len(l01_non_90)]
                     
        y_l05_det = [len(l05_det_15),
                     len(l05_det_30),
                     len(l05_det_45),
                     len(l05_det_60),
                     len(l05_det_75),
                     len(l05_det_90)]
                     
        y_l05_non = [len(l05_non_15),
                     len(l05_non_30),
                     len(l05_non_45),
                     len(l05_non_60),
                     len(l05_non_75),
                     len(l05_non_90)]
                     
        y_l1_det = [len(l1_det_15),
                     len(l1_det_30),
                     len(l1_det_45),
                     len(l1_det_60),
                     len(l1_det_75),
                     len(l1_det_90)]
                     
        y_l1_non = [len(l1_non_15),
                     len(l1_non_30),
                     len(l1_non_45),
                     len(l1_non_60),
                     len(l1_non_75),
                     len(l1_non_90)]
                     
        y_l5_det = [len(l5_det_15),
                     len(l5_det_30),
                     len(l5_det_45),
                     len(l5_det_60),
                     len(l5_det_75),
                     len(l5_det_90)]
                     
        y_l5_non = [len(l5_non_15),
                     len(l5_non_30),
                     len(l5_non_45),
                     len(l5_non_60),
                     len(l5_non_75),
                     len(l5_non_90)]
        
        y_l75_det = [len(l75_det_15),
                     len(l75_det_30),
                     len(l75_det_45),
                     len(l75_det_60),
                     len(l75_det_75),
                     len(l75_det_90)]
                     
        y_l75_non = [len(l75_non_15),
                     len(l75_non_30),
                     len(l75_non_45),
                     len(l75_non_60),
                     len(l75_non_75),
                     len(l75_non_90)]
        
        print 'y_l01_det: ',y_l01_det
        print 'y_l01_non: ',y_l01_non
        print 'y_l05_det: ',y_l05_det
        print 'y_l05_non: ',y_l05_non
        print 'y_l1_det: ',y_l1_det
        print 'y_l1_non: ',y_l1_non
        print 'y_l5_det: ',y_l5_det
        print 'y_l5_non: ',y_l5_non
        print 'y_l75_det: ',y_l75_det
        print 'y_l75_non: ',y_l75_non
        print
        
        # make them all floats
        
        y_l001_det = np.array(y_l001_det).astype(np.float)
        y_l005_det = np.array(y_l005_det).astype(np.float)
        y_l01_det = np.array(y_l01_det).astype(np.float)
        y_l05_det = np.array(y_l05_det).astype(np.float)
        y_l1_det = np.array(y_l1_det).astype(np.float)
        y_l5_det = np.array(y_l5_det).astype(np.float)
        y_l75_det = np.array(y_l75_det).astype(np.float)

        y_l001_non = np.array(y_l001_non).astype(np.float)
        y_l005_non = np.array(y_l005_non).astype(np.float)
        y_l01_non = np.array(y_l01_non).astype(np.float)
        y_l05_non = np.array(y_l05_non).astype(np.float)
        y_l1_non = np.array(y_l1_non).astype(np.float)
        y_l5_non = np.array(y_l5_non).astype(np.float)
        y_l75_non = np.array(y_l75_non).astype(np.float)
        
        
        # combine them into just 2 sets
        y_l1_det_inside = y_l75_det + y_l5_det + y_l1_det
        y_l1_det_outside = y_l05_det + y_l01_det + y_l005_det +y_l001_det
        
        y_l1_non_inside = y_l75_non + y_l5_non + y_l1_non
        y_l1_non_outside = y_l05_non + y_l01_non + y_l005_non + y_l001_non
        
        print 'y_l1_det_inside: ',y_l1_det_inside
        print 'y_l1_det_outside: ',y_l1_det_outside
        print 'y_l1_non_inside: ',y_l1_non_inside
        print 'y_l1_non_outside: ',y_l1_non_outside
        print
        print
        
        
        def detection_fraction_error(det, non):
            first_term = (non * det**2)/(non + det)**4
        
            second_term = det * (-det/(non + det)**2 + 1/(non + det))**2
            
            err = np.sqrt(first_term + second_term)
            
            return err
        
        detection_fraction_l1_outside = np.array(y_l1_det_outside) / np.array(np.array(y_l1_det_outside) + np.array(y_l1_non_outside))
        detection_fraction_l1_inside = np.array(y_l1_det_inside) / np.array(np.array(y_l1_det_inside) + np.array(y_l1_non_inside))


        detection_fraction_l1_outside_err = [detection_fraction_error(det, non) for det, non in zip(y_l1_det_outside, y_l1_non_outside)]
        detection_fraction_l1_inside_err = [detection_fraction_error(det, non) for det, non in zip(y_l1_det_inside, y_l1_non_inside)]

        print
        print 'detection_fraction_l1_outside_err: ',detection_fraction_l1_outside_err
        print 'detection_fraction_l1_inside_err: ',detection_fraction_l1_inside_err
        print
        print


        ##########
        print 'x: ',x
        print 'detection_fraction_l1_outside: ',detection_fraction_l1_outside
        print 'detection_fraction_l1_inside: ',detection_fraction_l1_inside
        ax1.errorbar(x,
                    detection_fraction_l1_outside,
                    yerr=detection_fraction_l1_outside_err,
                    marker=symbol_l1_outside,
                    capsize=3,
                    capthick=1,
                    c=color_l1_outside,
                    ms=markerSize_l1_outside,
                    markeredgecolor='black',
                    ls = ls_l1_outside,
                    lw = lw_l1_outside,
                    alpha=alpha_l1_outside,
                    label=label_l1_outside)

        ax1.errorbar(x,
                    detection_fraction_l1_inside,
                    yerr=detection_fraction_l1_inside_err,
                    marker=symbol_l1_inside,
                    c=color_l1_inside,
                    ms=markerSize_l1_inside,
                    markeredgecolor='black',
                    ls = ls_l1_inside,
                    lw = lw_l1_inside,
                    alpha=alpha_l1_inside,
                    label=label_l1_inside)


        ax1.set_xlabel(r'$\rm Azimuth~[deg]$')
        ax1.set_ylabel(r'$\rm Detection~Fraction$')

        xlim(0, 90)

        # x-axis
        majorLocator   = MultipleLocator(15)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax1.xaxis.set_major_locator(majorLocator)
        ax1.xaxis.set_major_formatter(majorFormatter)
        ax1.xaxis.set_minor_locator(minorLocator)


        # y-axis
        majorLocator   = MultipleLocator(0.1)
        majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
        minorLocator   = MultipleLocator(0.05)
        ax1.yaxis.set_major_locator(majorLocator)
        ax1.yaxis.set_major_formatter(majorFormatter)
        ax1.yaxis.set_minor_locator(minorLocator)

        ylim(0.4,1.0)

        try:
            import matplotlib.patches as mpatches
            import matplotlib.lines as mlines
                              
            l1_outside = mlines.Line2D([], [], color=color_l1_outside, marker=symbol_l1_outside,lw=lw_l1_outside,
                                      markersize=legend_size, markeredgecolor='black', label=label_l1_outside)
                                
            l1_inside = mlines.Line2D([], [], color=color_l1_inside, marker=symbol_l1_inside, lw=lw_l1_inside,
                                      markersize=legend_size, markeredgecolor='black', label=label_l1_inside)
                                  
            plt.legend(handles=[l1_inside, l1_outside],loc='lower right',
                                    borderpad=0.8, fontsize=legend_font, fancybox=True)
        except Exception, e:
            print
            print "Legend issue: ",e
            print

        ax1.grid(b=None,which='major',axis='both')


        if plot_detection_fraction_vs_az_like_save:
            savefig('{0}/detection_fraction_vs_az_like_sep15Rvir.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()




##########################################################################################
##########################################################################################
    
    if plot_detection_fraction_vs_az_imp:
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
        
        # L = 0.01 ~ 2.15 rho/rvir
        # L = 0.05 ~ 1.73 rho/rvir
        # L = 0.1 ~ 1.5 rho/rvir
        # L = 0.5 ~ 0.83 rho/rvir
        # L = 0.75 ~ 0.54 rho/vir

        label_imp1000 = r'$\rm \rho = 1000$'
        label_imp500 = r'$\rm \rho = 500$'
        label_imp300 = r'$\rm \rho = 300$'
        label_imp200 = r'$\rm \rho = 200$'
        label_imp100 = r'$\rm \rho = 100$'

#         label_l01 = '.01'
#         label_l05 = '.05'
#         label_l1 = '.1'
#         label_l5 = '.5'
#         label_l75 = '.75'

        alpha_imp1000 = 0.9
        alpha_imp500 = 0.9
        alpha_imp300 = 0.9
        alpha_imp200 = 0.9
        alpha_imp100 = 0.9

        legend_size = 10
        legend_font = 10

        markerSize_imp1000 = 10
        markerSize_imp500 = 10
        markerSize_imp300 = 10
        markerSize_imp200 = 10
        markerSize_imp100 = 10

        lw_imp1000 = 1.3
        lw_imp500 = 1.3
        lw_imp300 = 1.3
        lw_imp200 = 1.3
        lw_imp100 = 1.3

        symbol_imp1000 = 'X'
        symbol_imp500 = 'P'
        symbol_imp300 = 's'
        symbol_imp200 = 'd'
        symbol_imp100 = 'D'

        color_imp1000 = color_red
        color_imp500 = color_yellow
        color_imp300 = color_green
        color_imp200 = color_purple2
        color_imp100 = color_blue

        ls_imp1000 = 'solid'
        ls_imp500 = 'solid'
        ls_imp300 = 'solid'
        ls_imp200 = 'solid'
        ls_imp100 = 'solid'


##########################################################################################
        # do the plotting
        
        print 'len(dv400_imp100_det_az) : ',len(dv400_imp100_det_az)
        print 'len(dv400_imp50_det_az) : ',len(dv400_imp50_det_az)
        print 'len(dv400_imp25_det_az) : ',len(dv400_imp25_det_az)
        print
        print
        print 'len(dv400_imp100_non_az) : ',len(dv400_imp100_non_az)
        print 'len(dv400_imp50_non_az) : ',len(dv400_imp50_non_az)
        print 'len(dv400_imp25_non_az) : ',len(dv400_imp25_non_az)
        print


        x = np.arange(7.5, 97.5, 15)
        
        # L = 0.01 ~ 2.15 rho/rvir
        # L = 0.05 ~ 1.73 rho/rvir
        # L = 0.1 ~ 1.5 rho/rvir
        # L = 0.5 ~ 0.83 rho/rvir
        # L = 0.75 ~ 0.54 rho/vir

#######
        imp1000_det_15 = []
        imp1000_det_30 = []
        imp1000_det_45 = []
        imp1000_det_60 = []
        imp1000_det_75 = []
        imp1000_det_90 = []

        imp1000_non_15 = []
        imp1000_non_30 = []
        imp1000_non_45 = []
        imp1000_non_60 = []
        imp1000_non_75 = []
        imp1000_non_90 = []
        
        imp500_det_15 = []
        imp500_det_30 = []
        imp500_det_45 = []
        imp500_det_60 = []
        imp500_det_75 = []
        imp500_det_90 = []

        imp500_non_15 = []
        imp500_non_30 = []
        imp500_non_45 = []
        imp500_non_60 = []
        imp500_non_75 = []
        imp500_non_90 = []
        
        imp300_det_15 = []
        imp300_det_30 = []
        imp300_det_45 = []
        imp300_det_60 = []
        imp300_det_75 = []
        imp300_det_90 = []

        imp300_non_15 = []
        imp300_non_30 = []
        imp300_non_45 = []
        imp300_non_60 = []
        imp300_non_75 = []
        imp300_non_90 = []
        
        imp200_det_15 = []
        imp200_det_30 = []
        imp200_det_45 = []
        imp200_det_60 = []
        imp200_det_75 = []
        imp200_det_90 = []

        imp200_non_15 = []
        imp200_non_30 = []
        imp200_non_45 = []
        imp200_non_60 = []
        imp200_non_75 = []
        imp200_non_90 = []
        
        imp100_det_15 = []
        imp100_det_30 = []
        imp100_det_45 = []
        imp100_det_60 = []
        imp100_det_75 = []
        imp100_det_90 = []

        imp100_non_15 = []
        imp100_non_30 = []
        imp100_non_45 = []
        imp100_non_60 = []
        imp100_non_75 = []
        imp100_non_90 = []
        
        
        for az in dv400_imp1000_det_az:
            az = float(az)
            if az <= 15:
                imp1000_det_15.append(az)

            elif az > 15 and az <= 30:
                imp1000_det_30.append(az)

            elif az > 30 and az <= 45:
                imp1000_det_45.append(az)

            elif az > 45 and az <= 60:
                imp1000_det_60.append(az)

            elif az > 60 and az <= 75:
                imp1000_det_75.append(az)

            elif az > 75 and az <= 90:
                imp1000_det_90.append(az)
            else:
                print 'wtf? imp1000_det_az = ',az
            
        for az in dv400_imp1000_non_az:
            az = float(az)
            if az <= 15:
                imp1000_non_15.append(az)

            elif az > 15 and az <= 30:
                imp1000_non_30.append(az)

            elif az > 30 and az <= 45:
                imp1000_non_45.append(az)

            elif az > 45 and az <= 60:
                imp1000_non_60.append(az)

            elif az > 60 and az <= 75:
                imp1000_non_75.append(az)

            elif az > 75 and az <= 90:
                imp1000_non_90.append(az)
            else:
                print 'wtf? imp1000_non_az = ',az

        for az in dv400_imp500_det_az:
            az = float(az)
            if az <= 15:
                imp500_det_15.append(az)

            elif az > 15 and az <= 30:
                imp500_det_30.append(az)

            elif az > 30 and az <= 45:
                imp500_det_45.append(az)

            elif az > 45 and az <= 60:
                imp500_det_60.append(az)

            elif az > 60 and az <= 75:
                imp500_det_75.append(az)

            elif az > 75 and az <= 90:
                imp500_det_90.append(az)
            else:
                print 'wtf? imp500_det_az = ',az


        for az in dv400_imp500_non_az:
            az = float(az)
            if az <= 15:
                imp500_non_15.append(az)

            elif az > 15 and az <= 30:
                imp500_non_30.append(az)

            elif az > 30 and az <= 45:
                imp500_non_45.append(az)

            elif az > 45 and az <= 60:
                imp500_non_60.append(az)

            elif az > 60 and az <= 75:
                imp500_non_75.append(az)

            elif az > 75 and az <= 90:
                imp500_non_90.append(az)
            else:
                print 'wtf? imp500_non_az = ',az


        for az in dv400_imp300_det_az:
            az = float(az)
            if az <= 15:
                imp300_det_15.append(az)

            elif az > 15 and az <= 30:
                imp300_det_30.append(az)

            elif az > 30 and az <= 45:
                imp300_det_45.append(az)

            elif az > 45 and az <= 60:
                imp300_det_60.append(az)

            elif az > 60 and az <= 75:
                imp300_det_75.append(az)

            elif az > 75 and az <= 90:
                imp300_det_90.append(az)
            else:
                print 'wtf? imp300_det_az = ',az

        for az in dv400_imp300_non_az:
            az = float(az)
            if az <= 15:
                imp300_non_15.append(az)

            elif az > 15 and az <= 30:
                imp300_non_30.append(az)

            elif az > 30 and az <= 45:
                imp300_non_45.append(az)

            elif az > 45 and az <= 60:
                imp300_non_60.append(az)

            elif az > 60 and az <= 75:
                imp300_non_75.append(az)

            elif az > 75 and az <= 90:
                imp300_non_90.append(az)
            else:
                print 'wtf? imp300_non_az = ',az

        for az in dv400_imp200_det_az:
            az = float(az)
            if az <= 15:
                imp200_det_15.append(az)

            elif az > 15 and az <= 30:
                imp200_det_30.append(az)

            elif az > 30 and az <= 45:
                imp200_det_45.append(az)

            elif az > 45 and az <= 60:
                imp200_det_60.append(az)

            elif az > 60 and az <= 75:
                imp200_det_75.append(az)

            elif az > 75 and az <= 90:
                imp200_det_90.append(az)
            else:
                print 'wtf? imp200_det_az = ',az

        for az in dv400_imp200_non_az:
            az = float(az)
            if az <= 15:
                imp200_non_15.append(az)

            elif az > 15 and az <= 30:
                imp200_non_30.append(az)

            elif az > 30 and az <= 45:
                imp200_non_45.append(az)

            elif az > 45 and az <= 60:
                imp200_non_60.append(az)

            elif az > 60 and az <= 75:
                imp200_non_75.append(az)

            elif az > 75 and az <= 90:
                imp200_non_90.append(az)
            else:
                print 'wtf? imp200_non_az = ',az


        for az in dv400_imp100_det_az:
            az = float(az)
            if az <= 15:
                imp100_det_15.append(az)

            elif az > 15 and az <= 30:
                imp100_det_30.append(az)

            elif az > 30 and az <= 45:
                imp100_det_45.append(az)

            elif az > 45 and az <= 60:
                imp100_det_60.append(az)

            elif az > 60 and az <= 75:
                imp100_det_75.append(az)

            elif az > 75 and az <= 90:
                imp100_det_90.append(az)
            else:
                print 'wtf? imp100_det_az = ',az

        for az in dv400_imp100_non_az:
            az = float(az)
            if az <= 15:
                imp100_non_15.append(az)

            elif az > 15 and az <= 30:
                imp100_non_30.append(az)

            elif az > 30 and az <= 45:
                imp100_non_45.append(az)

            elif az > 45 and az <= 60:
                imp100_non_60.append(az)

            elif az > 60 and az <= 75:
                imp100_non_75.append(az)

            elif az > 75 and az <= 90:
                imp100_non_90.append(az)
            else:
                print 'wtf? imp100_non_az = ',az



        y_imp1000_det = [len(imp1000_det_15),
                     len(imp1000_det_30),
                     len(imp1000_det_45),
                     len(imp1000_det_60),
                     len(imp1000_det_75),
                     len(imp1000_det_90)]
                     
        y_imp1000_non = [len(imp1000_non_15),
                     len(imp1000_non_30),
                     len(imp1000_non_45),
                     len(imp1000_non_60),
                     len(imp1000_non_75),
                     len(imp1000_non_90)]
                     
        y_imp500_det = [len(imp500_det_15),
                     len(imp500_det_30),
                     len(imp500_det_45),
                     len(imp500_det_60),
                     len(imp500_det_75),
                     len(imp500_det_90)]
                     
        y_imp500_non = [len(imp500_non_15),
                     len(imp500_non_30),
                     len(imp500_non_45),
                     len(imp500_non_60),
                     len(imp500_non_75),
                     len(imp500_non_90)]
                     
        y_imp300_det = [len(imp300_det_15),
                     len(imp300_det_30),
                     len(imp300_det_45),
                     len(imp300_det_60),
                     len(imp300_det_75),
                     len(imp300_det_90)]
                     
        y_imp300_non = [len(imp300_non_15),
                     len(imp300_non_30),
                     len(imp300_non_45),
                     len(imp300_non_60),
                     len(imp300_non_75),
                     len(imp300_non_90)]
                     
        y_imp200_det = [len(imp200_det_15),
                     len(imp200_det_30),
                     len(imp200_det_45),
                     len(imp200_det_60),
                     len(imp200_det_75),
                     len(imp200_det_90)]
                     
        y_imp200_non = [len(imp200_non_15),
                     len(imp200_non_30),
                     len(imp200_non_45),
                     len(imp200_non_60),
                     len(imp200_non_75),
                     len(imp200_non_90)]
        
        y_imp100_det = [len(imp100_det_15),
                     len(imp100_det_30),
                     len(imp100_det_45),
                     len(imp100_det_60),
                     len(imp100_det_75),
                     len(imp100_det_90)]
                     
        y_imp100_non = [len(imp100_non_15),
                     len(imp100_non_30),
                     len(imp100_non_45),
                     len(imp100_non_60),
                     len(imp100_non_75),
                     len(imp100_non_90)]
        
        print 'y_imp1000_det: ',y_imp1000_det
        print 'y_imp1000_non: ',y_imp1000_non
        print 'y_imp500_det: ',y_imp500_det
        print 'y_imp500_non: ',y_imp500_non
        print 'y_imp300_det: ',y_imp300_det
        print 'y_imp300_non: ',y_imp300_non
        print 'y_imp200_det: ',y_imp200_det
        print 'y_imp200_non: ',y_imp200_non
        print 'y_imp100_det: ',y_imp100_det
        print 'y_imp100_non: ',y_imp100_non
        print
        
        # make them all floats
        y_imp1000_det = np.array(y_imp1000_det).astype(np.float)
        y_imp500_det = np.array(y_imp500_det).astype(np.float)
        y_imp300_det = np.array(y_imp300_det).astype(np.float)
        y_imp200_det = np.array(y_imp200_det).astype(np.float)
        y_imp100_det = np.array(y_imp100_det).astype(np.float)

        y_imp1000_non = np.array(y_imp1000_non).astype(np.float)
        y_imp500_non = np.array(y_imp500_non).astype(np.float)
        y_imp300_non = np.array(y_imp300_non).astype(np.float)
        y_imp200_non = np.array(y_imp200_non).astype(np.float)
        y_imp100_non = np.array(y_imp100_non).astype(np.float)
        
        detection_fraction_imp1000 = np.array(y_imp1000_det) / np.array(np.array(y_imp1000_det) + np.array(y_imp1000_non))
        detection_fraction_imp500 = np.array(y_imp500_det) / np.array(np.array(y_imp500_det) + np.array(y_imp500_non))
        detection_fraction_imp300 = np.array(y_imp300_det) / np.array(np.array(y_imp300_det) + np.array(y_imp300_non))
        detection_fraction_imp200 = np.array(y_imp200_det) / np.array(np.array(y_imp200_det) + np.array(y_imp200_non))
        detection_fraction_imp100 = np.array(y_imp100_det) / np.array(np.array(y_imp100_det) + np.array(y_imp100_non))

        detection_fraction_imp1000_err = np.sqrt(y_imp1000_det) / np.array(np.array(y_imp1000_det) + np.array(y_imp1000_non))
        detection_fraction_imp500_err = np.sqrt(y_imp500_det) / np.array(np.array(y_imp500_det) + np.array(y_imp500_non))
        detection_fraction_imp300_err = np.sqrt(y_imp300_det) / np.array(np.array(y_imp300_det) + np.array(y_imp300_non))
        detection_fraction_imp200_err = np.sqrt(y_imp200_det) / np.array(np.array(y_imp200_det) + np.array(y_imp200_non))
        detection_fraction_imp100_err = np.sqrt(y_imp100_det) / np.array(np.array(y_imp100_det) + np.array(y_imp100_non))



        ##########
        print 'x: ',x
        print 'detection_fraction_imp1000: ',detection_fraction_imp1000
        print 'detection_fraction_imp1000_err: ',detection_fraction_imp1000_err
        ax1.errorbar(x,
                    detection_fraction_imp1000,
#                     yerr=detection_fraction_imp1000_err,
                    marker=symbol_imp1000,
                    c=color_imp1000,
                    ms=markerSize_imp1000,
                    markeredgecolor='black',
                    ls = ls_imp1000,
                    lw = lw_imp1000,
                    alpha=alpha_imp1000,
                    label=label_imp1000)

        ax1.errorbar(x,
                    detection_fraction_imp500,
#                     yerr=detection_fraction_imp500_err,
                    marker=symbol_imp500,
                    c=color_imp500,
                    ms=markerSize_imp500,
                    markeredgecolor='black',
                    ls = ls_imp500,
                    lw = lw_imp500,
                    alpha=alpha_imp500,
                    label=label_imp500)

        ax1.errorbar(x,
                    detection_fraction_imp300,
#                     yerr=detection_fraction_imp300_err,
                    marker=symbol_imp300,
                    c=color_imp300,
                    ms=markerSize_imp300,
                    markeredgecolor='black',
                    ls = ls_imp300,
                    lw = lw_imp300,
                    alpha=alpha_imp300,
                    label=label_imp300)

        ax1.errorbar(x,
                    detection_fraction_imp200,
#                     yerr=detection_fraction_imp200_err,
                    marker=symbol_imp200,
                    c=color_imp200,
                    ms=markerSize_imp200,
                    markeredgecolor='black',
                    ls = ls_imp200,
                    lw = lw_imp200,
                    alpha=alpha_imp200,
                    label=label_imp200)
                    
        ax1.errorbar(x,
                    detection_fraction_imp100,
#                     yerr=detection_fraction_imp100_err,
                    marker=symbol_imp100,
                    c=color_imp100,
                    ms=markerSize_imp100,
                    markeredgecolor='black',
                    ls = ls_imp100,
                    lw = lw_imp100,
                    alpha=alpha_imp100,
                    label=label_imp100)


        ax1.set_xlabel(r'$\rm Azimuth~[deg]$')
        ax1.set_ylabel(r'$\rm Detection~Fraction$')

        xlim(0, 90)

        # x-axis
        majorLocator   = MultipleLocator(15)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax1.xaxis.set_major_locator(majorLocator)
        ax1.xaxis.set_major_formatter(majorFormatter)
        ax1.xaxis.set_minor_locator(minorLocator)


        # y-axis
        majorLocator   = MultipleLocator(0.1)
        majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
        minorLocator   = MultipleLocator(0.05)
        ax1.yaxis.set_major_locator(majorLocator)
        ax1.yaxis.set_major_formatter(majorFormatter)
        ax1.yaxis.set_minor_locator(minorLocator)

        ylim(0,1)

        try:
            import matplotlib.patches as mpatches
            import matplotlib.lines as mlines
                              
            imp1000 = mlines.Line2D([], [], color=color_imp1000, marker=symbol_imp1000,lw=lw_imp1000,
                                      markersize=legend_size, markeredgecolor='black', label=label_imp1000)
                              
            imp500 = mlines.Line2D([], [], color=color_imp500, marker=symbol_imp500,lw=lw_imp500,
                                    markeredgecolor='black', markersize=legend_size, label=label_imp500)
                                
            imp300 = mlines.Line2D([], [], color=color_imp300, marker=symbol_imp300,lw=lw_imp300,
                                      markersize=legend_size, markeredgecolor='black', label=label_imp300)
                              
            imp200 = mlines.Line2D([], [], color=color_imp200, marker=symbol_imp200,lw=lw_imp200,
                                    markeredgecolor='black', markersize=legend_size, label=label_imp200)
                                
            imp100 = mlines.Line2D([], [], color=color_imp100, marker=symbol_imp100,lw=lw_imp100,
                                    markeredgecolor='black', markersize=legend_size, label=label_imp100)
                                  
            plt.legend(handles=[imp1000, imp500, imp300, imp200, imp100],loc='lower right',
                                    borderpad=0.8, fontsize=legend_font, fancybox=True)
        except Exception, e:
            print
            print "Legend issue: ",e
            print

        ax1.grid(b=None,which='major',axis='both')


        if plot_detection_fraction_vs_az_imp_save:
            savefig('{0}/detection_fraction_vs_az_imp_errors.pdf'.format(saveDirectory),format='pdf',bbox_inches='tight')
        else:
            show()





#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


if __name__=="__main__":
    # do the work
    main()
    