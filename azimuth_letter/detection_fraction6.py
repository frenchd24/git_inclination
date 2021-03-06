#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)


$Id: detection_fraction5.py v 5.0 11/07/2018

Determine the detection fraction. Correlate the sightlines with galaxies, and figure out 
for what level of likelihood an absorber is always found, etc.


v5: Update for azimuth letter (11/07/18)

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
import correlateSingle11 as correlateSingle


# from astropy.io.votable import parse,tree

# from vo.table import parse
# import vo.tree


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

from matplotlib import rc


###########################################################################
# here are the dictionaries to be pickled. each one will the lists which follow
all = {}
isolated = {}
L_isolated = {}
L_associated_isolated = {}
L_associated = {}
L_nonassociated = {}
L_two = {}
L_three_plus = {}
L_group = {}
L_summed = {}



def calculate_likelihood(impact, R_vir, vel_dif):
    # likelihood with virial radius
    likelihood = math.exp(-(impact/R_vir)**2) * math.exp(-(vel_dif/200.)**2)
    
    return likelihood
    
    

def calculate_cusLikelihood(impact, MajDiam, vel_dif):
    # try this "sphere of influence" value
    m15 = MajDiam**1.5

    # likelihood with m15 instead
    likelihoodm15 = math.exp(-(impact/m15)**2) * math.exp(-(vel_dif/200.)**2)   

    return likelihoodm15
    
    
    
    
def bmean(a):
    # compute the mean of 'a' but just return -99.99 if the array is empty

    try:
        return np.mean(a)
    except Exception, e:
        print 'Error: ',e
        return -99.99
    
    
def is_null(x):
    # return True if x is any of -99, '-99', -99., '-99.', , -99.99, '-99.99', 'x', or blank
    return_val = False
    
    x = str(x).strip()
    if isNumber(x):
        if x == -99 or x == -99.99 or str(x) == '-99' or str(x) =='-99.99' or x == -99. or x == '-99.':
            return_val = True
    else:
        return_val = True
    
    return return_val

    
def is_closest_impact(Name, Vhel, Lya_v, impact, other_list, max_deltav):
    # check if there are other closer galaxies within 400 km/s of this one
    closest = True
    for other in other_list:
        other_impact, other_vhel, other_name, other_rvir = other
        if abs(float(other_vhel) - float(Lya_v)) <= max_deltav:
            if float(other_impact) < float(impact):
                if other_name != Name:
                    closest = False
    
    return closest


def is_closest_likelihood(Name, Vhel, Lya_v, l_used, other_list, max_deltav):
    # check if there are higher likelihood galaxies than this one
    closest = True
    for other in other_list:
        other_impact, other_vhel, other_name, other_rvir = other
        
        # velocity difference between absorption and other galaxy
        other_vel_dif = float(Lya_v) - other_vhel
        
        # calculate likelihood for other galaxy
        other_l_used = calculate_likelihood(other_impact, other_rvir, other_vel_dif)

        if other_l_used > l_used:
            if other_name != Name:
                closest = False
    
    return closest

    
def is_closest_impact_rvir(Name, Vhel, Lya_v, l_used, other_list, max_deltav):
    # check if there are higher likelihood galaxies than this one
    closest = True
    for other in other_list:
        other_impact, other_vhel, other_name, other_rvir = other
        
        # velocity difference between absorption and other galaxy
        other_vel_dif = float(Lya_v) - other_vhel
        
        # calculate likelihood for other galaxy
        other_l_used = other_impact/other_rvir

        if abs(other_vel_dif) <= max_deltav:
            if other_l_used < l_used:
                if other_name != Name:
                    closest = False
    
    return closest

    
def main():
    # correlation options
    maxSep = 1000.
    agnSeparation = 4000.
    minVcorr = 450.
    minSize = 0.
    max_deltav = 400.
#     rigor = 5
    
    # set Lstar min and max values. Set Lstar_max = False to disregard
    Lstar_min = 0
    Lstar_max = False # set to False to not require an Lstar measurement
    
    # maximum velocity cut
    vcut = False
    
#     min_EW = 300
    min_EW = 200

    # minimum diameter. False to bypass, a number > 0 to impose a diameter measurement needed
#     d_min = 0.01 # must have a diameter > 0.01
#     d_min = False # does not need a diameter measurement at all
    d_min = 5.0
    
    
    # sort based on likelihood cus instead of the regular one?
    use_likelihood_cus = False
    
    # double l if impact <= 1 R_vir?
    double_l_within_rvir = False
    
    # print a ton of shit out for testing if True
    verbose = False
    
    # cutoff - will only return this many targers - use for testing together with verbose. 
    cutoff = 1000000.
    
    
    if getpass.getuser() == 'frenchd':
#         filename = '/Users/frenchd/Research/inclination/git_inclination/maps/LG_correlation_combined5_14_edit.csv'
#         pickleFilename = '/Users/frenchd/Research/inclination/git_inclination/picklePilot_plusSALT_14.p'

#         targetlist_filename = '/Users/frenchd/Research/correlation/TARGETLIST_10_17_17_TOTAL.csv'
        targetlist_filename = '/Users/frenchd/Research/correlation/TARGETLIST_06_06_18_TOTAL.csv'
        gtPickleFilename = '/Users/frenchd/Research/GT_update2/pickleGT_filteredAll.p'

#         filename = '/Users/frenchd/Research/inclination/git_inclination/correlatedTargetList_5_29_18_measurements.csv'
#         filename = '/Users/frenchd/Research/inclination/git_inclination/targets/correlatedTargetList_5_29_18_measurements_copy_short.csv'
#         filename = '/Users/frenchd/Research/inclination/git_inclination/targets/correlatedTargetList_5_29_18_measurements_copy.csv'
        filename = '/Users/frenchd/Research/inclination/git_inclination/targets/correlatedTargetList_5_29_18_measurements_copy.csv'

        # pickle files
#         detection_fraction_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_lstarcut_all-1.p'
#         detection_fraction_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_dmin75_vcut2500_minEW50.p'
#         detection_fraction_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_minEW300_closestonly2.p'
#         detection_fraction_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_minEW50.p'
        detection_fraction_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_minEW200_dmin5_Lstar0-False.p'

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    
    theFile = open(filename,'rU')
    reader = csv.DictReader(theFile)
    
    # open the pickle files
    detection_fraction = open(detection_fraction_filename,'wt')

    
##########################################################################################
    # grab the target coordinates
    targetFile = open(targetlist_filename,'rU')
    targetReader = csv.DictReader(targetFile)
    
    target_coords = {}
    for t in targetReader:
        name = t['targetName']
        ra = t['degreesRA']
        dec = t['degreesDec']
        
        if not target_coords.has_key(name):
            target_coords[name] = {'RA':ra,'Dec':dec}
    
    targetFile.close()
    
##########################################################################################
##########################################################################################    
    
    # loop through the measurement list and make a dictionary containing a list of all the
    # lines and their velocity, EW, Na and b-parameter for each target (the dict key)
    
    target_dict = {}
    
    # everything lists
    # detections
    dv400_all_det_imps = []
    dv400_all_det_incs = []
    dv400_all_det_azs = []
    dv400_all_det_rvirs = []
    dv400_all_det_lstars = []
    dv400_all_det_mtypes = []
    dv400_all_det_majdiams = []
    dv400_all_det_bmags = []
    dv400_all_det_group_nums = []
    dv400_all_det_group_mems = []
    dv400_all_det_group_dists = []
    dv400_all_det_likelihoods = []

    # non-detections
    dv400_all_non_imps = []
    dv400_all_non_incs = []
    dv400_all_non_azs = []
    dv400_all_non_rvirs = []
    dv400_all_non_lstars = []
    dv400_all_non_mtypes = []
    dv400_all_non_majdiams = []
    dv400_all_non_bmags = []
    dv400_all_non_group_nums = []
    dv400_all_non_group_mems = []
    dv400_all_non_group_dists = []
    dv400_all_non_likelihoods = []
    
    
    
    # lists of \Delta v for galaxies within 500, 400, 300, 200, 100, 50, 25 kpc of a target
    dv_1000 = []
    dv_750 = []
    dv_500 = []
    dv_400 = []
    dv_300 = []
    dv_200 = []
    dv_100 = []
    dv_50 = []
    dv_25 = []
    
    # count of detections within each impact parameter window
    dv400_imp1000_det = 0
    dv400_imp750_det = 0
    dv400_imp500_det = 0
    dv400_imp400_det = 0
    dv400_imp300_det = 0
    dv400_imp200_det = 0
    dv400_imp100_det = 0
    dv400_imp50_det = 0
    dv400_imp25_det = 0
    
    
    # count of non-detections within each impact parameter window
    dv400_imp1000_non = 0
    dv400_imp750_non = 0
    dv400_imp500_non = 0
    dv400_imp400_non = 0
    dv400_imp300_non = 0
    dv400_imp200_non = 0
    dv400_imp100_non = 0
    dv400_imp50_non = 0
    dv400_imp25_non = 0
    

    
    # inclinations of galaxies near detections within each impact parameter window
    dv400_imp1000_det_inc = []
    dv400_imp750_det_inc = []
    dv400_imp500_det_inc = []
    dv400_imp400_det_inc = []
    dv400_imp300_det_inc = []
    dv400_imp200_det_inc = []
    dv400_imp100_det_inc = []
    dv400_imp50_det_inc = []
    dv400_imp25_det_inc = []
    
    # inclinations of galaxies near non-detections within each impact parameter window
    dv400_imp1000_non_inc = []
    dv400_imp750_non_inc = []
    dv400_imp500_non_inc = []
    dv400_imp400_non_inc = []
    dv400_imp300_non_inc = []
    dv400_imp200_non_inc = []
    dv400_imp100_non_inc = []
    dv400_imp50_non_inc = []
    dv400_imp25_non_inc = []
    
    
    # Lstars of galaxies near detections within each impact parameter window
    dv400_imp1000_det_lstar = []
    dv400_imp750_det_lstar = []
    dv400_imp500_det_lstar = []
    dv400_imp400_det_lstar = []
    dv400_imp300_det_lstar = []
    dv400_imp200_det_lstar = []
    dv400_imp100_det_lstar = []
    dv400_imp50_det_lstar = []
    dv400_imp25_det_lstar = []
    
    # Lstars of galaxies near non-detections within each impact parameter window
    dv400_imp1000_non_lstar = []
    dv400_imp750_non_lstar = []
    dv400_imp500_non_lstar = []
    dv400_imp400_non_lstar = []
    dv400_imp300_non_lstar = []
    dv400_imp200_non_lstar = []
    dv400_imp100_non_lstar = []
    dv400_imp50_non_lstar = []
    dv400_imp25_non_lstar = []
    
    
    # az of galaxies near detections within each impact parameter window
    dv400_imp1000_det_az = []
    dv400_imp750_det_az = []
    dv400_imp500_det_az = []
    dv400_imp400_det_az = []
    dv400_imp300_det_az = []
    dv400_imp200_det_az = []
    dv400_imp100_det_az = []
    dv400_imp50_det_az = []
    dv400_imp25_det_az = []
    
    # az of galaxies near non-detections within each impact parameter window
    dv400_imp1000_non_az = []
    dv400_imp750_non_az = []
    dv400_imp500_non_az = []
    dv400_imp400_non_az = []
    dv400_imp300_non_az = []
    dv400_imp200_non_az = []
    dv400_imp100_non_az = []
    dv400_imp50_non_az = []
    dv400_imp25_non_az = []

##########################################################################################
    # now for likelihood thresholds
    # lists of \Delta v for galaxies within 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 0.75 L of a target
    dv_l0001 = []
    dv_l0005 = []
    dv_l001 = []
    dv_l005 = []
    dv_l01 = []
    dv_l05 = []
    dv_l1 = []
    dv_l5 = []
    dv_l75 = []
    
    # count of detections within each likelihood window
    dv400_l0001_det = 0
    dv400_l0005_det = 0
    dv400_l001_det = 0
    dv400_l005_det = 0
    dv400_l01_det = 0
    dv400_l05_det = 0
    dv400_l1_det = 0
    dv400_l5_det = 0
    dv400_l75_det = 0
    
    # count of non-detections within each likelihood window
    dv400_l0001_non = 0
    dv400_l0005_non = 0
    dv400_l001_non = 0
    dv400_l005_non = 0
    dv400_l01_non = 0
    dv400_l05_non = 0
    dv400_l1_non = 0
    dv400_l5_non = 0
    dv400_l75_non = 0
    
    
    # inclinations of galaxies near detections within each likelihood window
    dv400_l0001_det_inc = []
    dv400_l0005_det_inc = []
    dv400_l001_det_inc = []
    dv400_l005_det_inc = []
    dv400_l01_det_inc = []
    dv400_l05_det_inc = []
    dv400_l1_det_inc = []
    dv400_l5_det_inc = []
    dv400_l75_det_inc = []
    
    # inclinations of galaxies near non-detections within each likelihood window
    dv400_l0001_non_inc = []
    dv400_l0005_non_inc = []
    dv400_l001_non_inc = []
    dv400_l005_non_inc = []
    dv400_l01_non_inc = []
    dv400_l05_non_inc = []
    dv400_l1_non_inc = []
    dv400_l5_non_inc = []
    dv400_l75_non_inc = []
    
    
    # lstars of galaxies near detections within each likelihood window
    dv400_l0001_det_lstar = []
    dv400_l0005_det_lstar = []
    dv400_l001_det_lstar = []
    dv400_l005_det_lstar = []
    dv400_l01_det_lstar = []
    dv400_l05_det_lstar = []
    dv400_l1_det_lstar = []
    dv400_l5_det_lstar = []
    dv400_l75_det_lstar = []
    
    # lstars of galaxies near non-detections within each likelihood window
    dv400_l0001_non_lstar = []
    dv400_l0005_non_lstar = []
    dv400_l001_non_lstar = []
    dv400_l005_non_lstar = []
    dv400_l01_non_lstar = []
    dv400_l05_non_lstar = []
    dv400_l1_non_lstar = []
    dv400_l5_non_lstar = []
    dv400_l75_non_lstar = []
    
    
    # az of galaxies near detections within each likelihood window
    dv400_l0001_det_az = []
    dv400_l0005_det_az = []
    dv400_l001_det_az = []
    dv400_l005_det_az = []
    dv400_l01_det_az = []
    dv400_l05_det_az = []
    dv400_l1_det_az = []
    dv400_l5_det_az = []
    dv400_l75_det_az = []
    
    # az of galaxies near non-detections within each likelihood window
    dv400_l0001_non_az = []
    dv400_l0005_non_az = []
    dv400_l001_non_az = []
    dv400_l005_non_az = []
    dv400_l01_non_az = []
    dv400_l05_non_az = []
    dv400_l1_non_az = []
    dv400_l5_non_az = []
    dv400_l75_non_az = []
    
    
##########################################################################################
    # now for impact / R_vir
    # lists of \Delta v for galaxies within 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5 imp/Rvir of a target
    dv_rvir025 = []
    dv_rvir05 = []
    dv_rvir075 = []
    dv_rvir1 = []
    dv_rvir15 = []
    dv_rvir2 = []
    dv_rvir25 = []
    
    # count of detections within each imp/Rvir window
    dv400_rvir025_det = 0
    dv400_rvir05_det = 0
    dv400_rvir075_det = 0
    dv400_rvir1_det = 0
    dv400_rvir15_det = 0
    dv400_rvir2_det = 0
    dv400_rvir25_det = 0
    
    # count of non-detections within each impact/Rvir window
    dv400_rvir025_non = 0
    dv400_rvir05_non = 0
    dv400_rvir075_non = 0
    dv400_rvir1_non = 0
    dv400_rvir15_non = 0
    dv400_rvir2_non = 0
    dv400_rvir25_non = 0


    # no galaxy detections
    impact_no_galaxies = []
    likelihood_no_galaxies = []
    impact_rvir_no_galaxies = []

    print
    print 'Starting loop!'
    print
    
    total = 1551
    counter = 0
    stopCount = 50000
    for i in reader:
        counter +=1
        target = i['target']

        identified = i['identified']
        Lya_v = i['Lya_v']
        e_Lya_v = i['e_v_center']
        Lya_W = i['Lya_W']
        e_Lya_W = i['e_Lya_W']
        Na = i['Na']
        e_Na = i['e_Na']
        b = i['b']
        e_b = i['e_b']
        W = i['W']
        e_W = i['e_W']
        z_target = i['z']
        
        if not isNumber(Lya_v):
            Lya_v = 'x'

        if not isNumber(Lya_W):
            Lya_W = 'x'
            
        if not isNumber(Na):
            Na = 'x'
            
        if not isNumber(b):
            b = 'x'


        if isNumber(Lya_W):
            if float(Lya_W) < min_EW:
                Lya_v = 'x'
                Lya_W = 'x'
                Na = 'x'
                b = 'x'

        
        proceed = False
        # 'lya' are properly identified lines, 'maybe' are for possible lines, not yet
        # identified, and 'yes' are for identified targets with no Lya in the spectra
        if bfind(identified.lower(), 'lya') or bfind(identified.lower(), 'maybe') or bfind(identified.lower(), 'yes'):
            proceed = True

        if counter >= stopCount:
            proceed = False
            
                
        if proceed:
            try:
                RA_target = target_coords[target]['RA']
                Dec_target = target_coords[target]['Dec']
            except Exception, e:
                proceed = False
                print
                print "Target not found: ",target
            
            
            if target_dict.has_key(target):
            
                info = target_dict[target]
                
                # add this line velocity to the list
                Lya_vs = info['Lya_vs']
                Lya_vs.append(Lya_v)
                info['Lya_vs'] = Lya_vs
                
                # add line EW
                Lya_Ws = info['Lya_Ws']
                Lya_Ws.append(Lya_W)
                info['Lya_Ws'] = Lya_Ws
                
                # add line Na 
                Nas = info['Nas']
                Nas.append(Na)
                info['Nas'] = Nas
                
                # add line b 
                bs = info['bs']
                bs.append(b)
                info['bs'] = bs


            else:
                info = {}
                info['Lya_vs'] = [Lya_v]
                info['Lya_Ws'] = [Lya_W]
                info['Nas'] = [Na]
                info['bs'] = [b]
                info['z_target'] = z_target
                info['RA_target'] = RA_target
                info['Dec_target'] = Dec_target
                
                target_dict[target] = info

    print
    print 'Finished compiling target stuff'
    print
    print 'len(target_dict) : ',len(target_dict)
    print
#     print 'target_dict: ',target_dict
    print
    print 'Starting correlation...'
    print
    
    counter = 0
    for target in target_dict:
        counter +=1
        sys.stdout.write("\r Percent Complete: {0} / {1}".format(counter, total))
        sys.stdout.flush()
    
        info = target_dict[target]
        Lya_vs = info['Lya_vs']
        Lya_vs_like = info['Lya_vs']
        Lya_vs_rvir = info['Lya_vs']


        # do the correlation: returns a dictionary where the keys are all the correlated
        # galaxies and the values are a dictionary of all the galaxy info
#         correlation_original = correlateSingle.correlateTarget(target, maxSep, agnSeparation, minVcorr, minSize, slow=False, searchAll=False, returnDict=True)
        
        if counter < cutoff:
            correlation = correlateSingle.correlateTarget(target, maxSep, agnSeparation, minVcorr, minSize, slow=False, searchAll=False, returnDict=True)
        else:
            break

        # only look at the closest galaxy
#         impact_sorted_list = []
#         for c in correlation_original:
#             impact = correlation_original[c]['impact']
#             MajDiam = correlation_original[c]['MajDiam']
#             Vhel = correlation_original[c]['Vhel']
# 
#             if isNumber(MajDiam):
#                 if float(MajDiam) >= d_min:
#                     impact_sorted_list.append([impact, correlation_original[c]])
#        
#         impact_sorted_list.sort()
#        
#         c2 = impact_sorted_list[0][1]
#        
#         correlation = {}
#         correlation[c2['Name']] = c2

        other_list = []
        
        if verbose:
            print
            print 'target = {0}  ----'.format(target)
            
        # make a list of all galaxies in this correlation, decide which galaxies pass the 
        # given criteria for inclusion
        
        # the list of galaxy names which pass the test
        passing_galaxies = []
        
        for c in correlation:
            passing = True
            target = correlation[c]['target']
            Name = correlation[c]['Name']
            Vhel = correlation[c]['Vhel']
            impact = correlation[c]['impact']
            R_vir = correlation[c]['R_vir']
            Lstar_med = correlation[c]['Lstar_med']
            Lstar_sdss = correlation[c]['Lstar_sdss']
            MajDiam = correlation[c]['MajDiam']

            if str(Lstar_med) == '-99.99':
                Lstar = Lstar_sdss
            else:
                Lstar = Lstar_med
                
            # restrict galaxies based on Lstar min and max?
            if Lstar_max:
                if not is_null(Lstar):
                    Lstar = float(Lstar)
                    if float(Lstar) >=0:
                        if float(Lstar) < Lstar_min or float(Lstar) > Lstar_max:
                            passing = False
                        
                    else:
                        passing = False
                else:
                    passing = False
            
            # restrict galaxies based on velocity?
            if vcut:
                if Vhel > vcut:
                    passing = False
            
            # restrict based on size?
            if d_min:
                if isNumber(MajDiam):
                    if float(MajDiam) < d_min:
                        passing = False
                else:
                    passing = False

            if passing:
                other_list.append([float(impact), float(Vhel), Name, float(R_vir)])
                passing_galaxies.append(Name)
                
                if verbose:
                    print 'Name = {0}'.format(Name)
                    print 'impact = {0} , vhel = {1}'.format(impact, Vhel)
                    print
        
        
        # now do the real work
        for c in correlation:
#             proceed = True
            
            # unpack everything
            target = correlation[c]['target']
            z_target = correlation[c]['z_target']
            RA_target = correlation[c]['RA_target']
            Dec_target = correlation[c]['Dec_target']
            Name = correlation[c]['Name']
            RAdeg = correlation[c]['RAdeg']
            DEdeg = correlation[c]['DEdeg']
            impact = correlation[c]['impact']
            azimuth = correlation[c]['azimuth']
            PA = correlation[c]['PA']
            inc = correlation[c]['inc']
            adjustedInc = correlation[c]['adjustedInc']

            R_vir = correlation[c]['R_vir']
            MajDiam = correlation[c]['MajDiam']
            MType = correlation[c]['MType']
            Vhel = correlation[c]['Vhel']
            vcorr = correlation[c]['vcorr']
            bestDist = correlation[c]['bestDist']
            e_bestDist = correlation[c]['e_bestDist']
            group_num = correlation[c]['group_num']
            group_mem = correlation[c]['group_mem']
            group_dist = correlation[c]['group_dist']
            Lstar_med = correlation[c]['Lstar_med']
            e_Lstar_med = correlation[c]['e_Lstar_med']
            Lstar_sdss = correlation[c]['Lstar_sdss']
            e_Lstar_sdss = correlation[c]['e_Lstar_sdss']
            Bmag = correlation[c]['Bmag']
            Bmag_sdss = correlation[c]['Bmag_sdss']

            impact = float(impact)
            R_vir = float(R_vir)
            Vhel = float(Vhel)
            adjustedInc = float(adjustedInc)
            azimuth = float(azimuth)
            
            if str(Lstar_med) == '-99.99':
                Lstar = Lstar_sdss
            else:
                Lstar = Lstar_med
            
            if is_null(Bmag):
                Bmag_used = Bmag_sdss
            else:
                Bmag_used = Bmag
                
            Lstar = float(Lstar)

            
#             restrict galaxies based on Lstar min and max?
#             if Lstar_max:
#                 if not is_null(Lstar):
#                     Lstar = float(Lstar)
#                     if float(Lstar) >=0:
#                         if float(Lstar) < Lstar_min or float(Lstar) > Lstar_max:
#                             proceed = False
#                         
#                     else:
#                         proceed = False
#                 else:
#                     proceed = False
#             
#             restrict galaxies based on velocity?
#             if vcut:
#                 if Vhel > vcut:
#                     proceed = False
#             
#             restrict based on size?
#             if d_min:
#                 if isNumber(MajDiam):
#                     if float(MajDiam) < d_min:
#                         proceed = False
#                 else:
#                     proceed = False

            # did this galaxy pass the previous tests?
            if Name in passing_galaxies:
                proceed = True
            else:
                proceed = False

            if verbose:
                print 'Now starting {0}'.format(Name)
                print 'proceed = {0}'.format(proceed)
                print

            if proceed:
                # grab everything
                if float(impact) <= maxSep:
                    detection = False
                    for Lya_v in Lya_vs:
                        if isNumber(Lya_v):
                            dv = float(Lya_v) - Vhel
                        
                            if abs(dv) <= max_deltav:
                                detection = True
                            
                            if detection:
                                dv400_all_det_imps.append(float(impact))
                                dv400_all_det_incs.append(float(adjustedInc))
                                dv400_all_det_azs.append(float(azimuth))
                                dv400_all_det_rvirs.append(float(R_vir))
                                dv400_all_det_lstars.append(float(Lstar))
                                dv400_all_det_mtypes.append(MType)
                                dv400_all_det_majdiams.append(MajDiam)
                                dv400_all_det_bmags.append(Bmag_used)
                                dv400_all_det_group_nums.append(group_num)
                                dv400_all_det_group_mems.append(group_mem)
                                dv400_all_det_group_dists.append(group_dist)
                            else:
                                dv400_all_non_imps.append(float(impact))
                                dv400_all_non_incs.append(float(adjustedInc))
                                dv400_all_non_azs.append(float(azimuth))
                                dv400_all_non_rvirs.append(float(R_vir))
                                dv400_all_non_lstars.append(float(Lstar))
                                dv400_all_non_mtypes.append(MType)
                                dv400_all_non_majdiams.append(MajDiam)
                                dv400_all_non_bmags.append(Bmag_used)
                                dv400_all_non_group_nums.append(group_num)
                                dv400_all_non_group_mems.append(group_mem)
                                dv400_all_non_group_dists.append(group_dist)
                            

                # if there's a galaxy within 25 kpc, search for a line 
                if float(impact) <= 25:
                    detection = False
                    for Lya_v in Lya_vs:
                        if isNumber(Lya_v):
                            dv = float(Lya_v) - Vhel
                            dv_25.append(dv)
            
                            if abs(dv) <= max_deltav:
                                if is_closest_impact(Name, Vhel, Lya_v, impact, other_list, max_deltav):
                                    detection = True
    #                                 Lya_vs.pop(Lya_vs.index(Lya_v))

                                    if verbose:
                                        print 'Detection <= 25kpc - {0}'.format(Name)
                                        print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)

                    if detection:
                        dv400_imp25_det +=1
                        if not is_null(adjustedInc):
                            if float(adjustedInc) < 0:
                                print 'WTF: inc = {0} , Name = {1}'.format(adjustedInc, Name)
                                
                            dv400_imp25_det_inc.append(adjustedInc)
                        
                        if not is_null(Lstar):
                            dv400_imp25_det_lstar.append(Lstar)
                            
                        if not is_null(azimuth):
                            dv400_imp25_det_az.append(azimuth)
                    else:
                        dv400_imp25_non +=1
                        if not is_null(adjustedInc):
                            if float(adjustedInc) < 0:
                                print 'WTF: inc = {0} , Name = {1}'.format(adjustedInc, Name)
                            dv400_imp25_non_inc.append(adjustedInc)
                        
                        if not is_null(Lstar):
                            dv400_imp25_non_lstar.append(Lstar)
                            
                        if not is_null(azimuth):
                            dv400_imp25_non_az.append(azimuth)

                        print 'non-detection within 25: {0} - {1}'.format(target, Name)
            
            
                # if there's a galaxy within 50 kpc, search for a line 
                elif float(impact) <= 50:
                    detection = False
                    for Lya_v in Lya_vs:
                        if isNumber(Lya_v):
                            dv = float(Lya_v) - Vhel
                            dv_50.append(dv)
            
                            if abs(dv) <= max_deltav:
                                if is_closest_impact(Name, Vhel, Lya_v, impact, other_list, max_deltav):
                                    detection = True
    #                                 Lya_vs.pop(Lya_vs.index(Lya_v))

                                    if verbose:
                                        print 'Detection <= 50kpc - {0}'.format(Name)
                                        print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)

                    if detection:
                        dv400_imp50_det +=1
                        if not is_null(adjustedInc):
                            dv400_imp50_det_inc.append(adjustedInc)
                        
                        if not is_null(Lstar):
                            dv400_imp50_det_lstar.append(Lstar)
                            
                        if not is_null(azimuth):
                            dv400_imp50_det_az.append(azimuth)
                    else:
                        dv400_imp50_non +=1
                        if not is_null(adjustedInc):
                            dv400_imp50_non_inc.append(adjustedInc)
                        
                        if not is_null(Lstar):
                            dv400_imp50_non_lstar.append(Lstar)
                            
                        if not is_null(azimuth):
                            dv400_imp50_non_az.append(azimuth)


                # if there's a galaxy within 100 kpc, search for a line 
                elif float(impact) <= 100:
                    detection = False
                    for Lya_v in Lya_vs:
                        if isNumber(Lya_v):
                            dv = float(Lya_v) - Vhel
                            dv_100.append(dv)
                        
                            if abs(dv) <= max_deltav:
                                if is_closest_impact(Name, Vhel, Lya_v, impact, other_list, max_deltav):
                                    detection = True
    #                                 Lya_vs.pop(Lya_vs.index(Lya_v))

                                    if verbose:
                                        print 'Detection <= 100kpc - {0}'.format(Name)
                                        print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)

                    if detection:
                        dv400_imp100_det +=1
                        if not is_null(adjustedInc):
                            dv400_imp100_det_inc.append(adjustedInc)
                        
                        if not is_null(Lstar): 
                            dv400_imp100_det_lstar.append(Lstar)
                            
                        if not is_null(azimuth): 
                            dv400_imp100_det_az.append(azimuth)
                    else:
                        dv400_imp100_non +=1
                        if not is_null(adjustedInc):
                            dv400_imp100_non_inc.append(adjustedInc)
                        
                        if not is_null(Lstar):
                            dv400_imp100_non_lstar.append(Lstar)
                            
                        if not is_null(azimuth): 
                            dv400_imp100_non_az.append(azimuth)
                            

                # if there's a galaxy within 200 kpc, search for a line 
                elif float(impact) <= 200:
                    detection = False
                    for Lya_v in Lya_vs:
                        if isNumber(Lya_v):
                            dv = float(Lya_v) - Vhel
                            dv_200.append(dv)

                            if abs(dv) <= max_deltav:
                                if is_closest_impact(Name, Vhel, Lya_v, impact, other_list, max_deltav):
                                    detection = True
    #                                 Lya_vs.pop(Lya_vs.index(Lya_v))

                                    if verbose:
                                        print 'Detection <= 200kpc - {0}'.format(Name)
                                        print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)

                    if detection:
                        dv400_imp200_det +=1
                        if not is_null(adjustedInc):
                            dv400_imp200_det_inc.append(adjustedInc)
                        
                        if not is_null(Lstar):
                            dv400_imp200_det_lstar.append(Lstar)
                            
                        if not is_null(azimuth):
                            dv400_imp200_det_az.append(azimuth)
                    else:
                        dv400_imp200_non +=1
                        if not is_null(adjustedInc):
                            dv400_imp200_non_inc.append(adjustedInc)
                            
                        if not is_null(Lstar):
                            dv400_imp200_non_lstar.append(Lstar)
                            
                        if not is_null(azimuth):
                            dv400_imp200_non_az.append(azimuth)


                # if there's a galaxy within 300 kpc, search for a line 
                elif float(impact) <= 300:
                    detection = False
                    for Lya_v in Lya_vs:
                        if isNumber(Lya_v):
                            dv = float(Lya_v) - Vhel
                            dv_300.append(dv)
                        
                            if abs(dv) <= max_deltav:
                                if is_closest_impact(Name, Vhel, Lya_v, impact, other_list, max_deltav):
                                    detection = True
    #                                 Lya_vs.pop(Lya_vs.index(Lya_v))

                                    if verbose:
                                        print 'Detection <= 300kpc - {0}'.format(Name)
                                        print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)
                                    
                    if detection:
                        dv400_imp300_det +=1
                        if not is_null(adjustedInc):
                            dv400_imp300_det_inc.append(adjustedInc)
                        
                        if not is_null(Lstar):
                            dv400_imp300_det_lstar.append(Lstar)
                            
                        if not is_null(azimuth):
                            dv400_imp300_det_az.append(azimuth)
                    else:
                        dv400_imp300_non +=1
                        if not is_null(adjustedInc):
                            dv400_imp300_non_inc.append(adjustedInc)
                        
                        if not is_null(Lstar):
                            dv400_imp300_non_lstar.append(Lstar)
                            
                        if not is_null(azimuth):
                            dv400_imp300_non_az.append(azimuth)


                # if there's a galaxy within 400 kpc, search for a line
                elif float(impact) <= 400:
                    detection = False
                    for Lya_v in Lya_vs:
                        if isNumber(Lya_v):
                            dv = float(Lya_v) - Vhel
                            dv_400.append(dv)
            
                            if abs(dv) <= max_deltav:
                                if is_closest_impact(Name, Vhel, Lya_v, impact, other_list, max_deltav):
                                    detection = True
    #                                 Lya_vs.pop(Lya_vs.index(Lya_v))

                                    if verbose:
                                        print 'Detection <= 400kpc - {0}'.format(Name)
                                        print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)

                    if detection:
                        dv400_imp400_det +=1
                        if not is_null(adjustedInc):
                            dv400_imp400_det_inc.append(adjustedInc)
                        
                        if not is_null(Lstar):
                            dv400_imp400_det_lstar.append(Lstar)
                            
                        if not is_null(azimuth):
                            dv400_imp400_det_az.append(azimuth)
                    else:
                        dv400_imp400_non +=1
                        if not is_null(adjustedInc):
                            dv400_imp400_non_inc.append(adjustedInc)
                            
                        if not is_null(Lstar):
                            dv400_imp400_non_lstar.append(Lstar)
                            
                        if not is_null(azimuth):
                            dv400_imp400_non_az.append(azimuth)


                # if there's a galaxy within 500 kpc, search for a line 
                elif float(impact) <= 500:
                    detection = False
                    for Lya_v in Lya_vs:
                        if isNumber(Lya_v):
                            dv = float(Lya_v) - Vhel
                            dv_500.append(dv)
                        
                            if abs(dv) <= max_deltav:
                                if is_closest_impact(Name, Vhel, Lya_v, impact, other_list, max_deltav):
                                    detection = True
    #                                 Lya_vs.pop(Lya_vs.index(Lya_v))

                                    if verbose:
                                        print 'Detection <= 500kpc - {0}'.format(Name)
                                        print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)

                    if detection:
                        dv400_imp500_det +=1
                        if not is_null(adjustedInc):
                            dv400_imp500_det_inc.append(adjustedInc)
                        
                        if not is_null(Lstar):
                            dv400_imp500_det_lstar.append(Lstar)
                            
                        if not is_null(azimuth):
                            dv400_imp500_det_az.append(azimuth)
                    else:
                        dv400_imp500_non +=1
                        if not is_null(adjustedInc):
                            dv400_imp500_non_inc.append(adjustedInc)
                         
                        if not is_null(Lstar):
                            dv400_imp500_non_lstar.append(Lstar)
                            
                        if not is_null(azimuth):
                            dv400_imp500_non_az.append(azimuth)


                # if there's a galaxy within 750 kpc, search for a line 
                elif float(impact) <= 750:
                    detection = False
                    for Lya_v in Lya_vs:
                        if isNumber(Lya_v):
                            dv = float(Lya_v) - Vhel
                            dv_750.append(dv)
                        
                            if abs(dv) <= max_deltav:
                                if is_closest_impact(Name, Vhel, Lya_v, impact, other_list, max_deltav):
                                    detection = True
    #                                 Lya_vs.pop(Lya_vs.index(Lya_v))

                                    if verbose:
                                        print 'Detection <= 750kpc - {0}'.format(Name)
                                        print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)

                    if detection:
                        dv400_imp750_det +=1
                        if not is_null(adjustedInc):
                            dv400_imp750_det_inc.append(adjustedInc)
                            
                        if not is_null(Lstar):
                            dv400_imp750_det_lstar.append(Lstar)
                            
                        if not is_null(azimuth):
                            dv400_imp750_det_az.append(azimuth)
                    else:
                        dv400_imp750_non +=1
                        if not is_null(adjustedInc):
                            dv400_imp750_non_inc.append(adjustedInc)
                        
                        if not is_null(Lstar):
                            dv400_imp750_non_lstar.append(Lstar)
                            
                        if not is_null(azimuth):
                            dv400_imp750_non_az.append(azimuth)
                    
                    
                # if there's a galaxy within 1000 kpc, search for a line 
                elif float(impact) <= 1000:
                    detection = False
                    for Lya_v in Lya_vs:
                        if isNumber(Lya_v):
                            dv = float(Lya_v) - Vhel
                            dv_1000.append(dv)
                        
                            if abs(dv) <= max_deltav:
                                if is_closest_impact(Name, Vhel, Lya_v, impact, other_list, max_deltav):
                                    detection = True
    #                                 Lya_vs.pop(Lya_vs.index(Lya_v))

                                    if verbose:
                                        print 'Detection <= 1000kpc - {0}'.format(Name)
                                        print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)

                    if detection:
                        dv400_imp1000_det +=1
                        if not is_null(adjustedInc):
                            dv400_imp1000_det_inc.append(adjustedInc)
                            
                        if not is_null(Lstar):
                            dv400_imp1000_det_lstar.append(Lstar)

                        if not is_null(azimuth):
                            dv400_imp1000_det_az.append(azimuth)
                    else:
                        dv400_imp1000_non +=1
                        if not is_null(adjustedInc):
                            dv400_imp1000_non_inc.append(adjustedInc)
                            
                        if not is_null(Lstar):
                            dv400_imp1000_non_lstar.append(Lstar)

                        if not is_null(azimuth):
                            dv400_imp1000_non_az.append(azimuth)

                else:
                    for v in Lya_vs:
                        impact_no_galaxies.append(v)
            
    ##########################################################################################            
                # now do it for likelihood
            
            if not is_null(MajDiam):
                # "sphere of influence" value for likelihood_custom
                MajDiam = float(MajDiam)
                impact = float(impact)
                R_vir = float(R_vir)
            
                cus = MajDiam**1.5

                # first for the virial radius
                likelihood = math.exp(-(impact/R_vir)**2)
        
                # then for the second 'virial like' m15 radius
                likelihood_cus = math.exp(-(impact/cus)**2)

                # sort based on which likelihood metric?
                if use_likelihood_cus:
                    l_used = likelihood_cus
                else:
                    l_used = likelihood
            
                # multiply the likelihood by two if within 1 R_vir?
                if double_l_within_rvir:
                    if impact <= R_vir:
                        l_used = l_used * 2

                
                if proceed:
                    # grab everything
                    if float(impact) <= maxSep:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                                        
                                if abs(dv) <= max_deltav:
                                    detection = True
                            
                                if detection:
                                    dv400_all_det_likelihoods.append(float(l_used))
                                else:
                                    dv400_all_non_likelihoods.append(float(l_used))
                            

                    # if the likelihood is greater than 0.75, see if there's a corresponding
                    if l_used >= 0.75:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_l75.append(dv)
                        
                                if abs(dv) <= max_deltav:
                                    if is_closest_likelihood(Name, Vhel, Lya_v, l_used, other_list, max_deltav):
                                        detection = True
    #                                     Lya_vs_like.pop(Lya_vs_like.index(Lya_v))

                                        if verbose:
                                            print 'Detection <= 0.75L - {0}'.format(Name)
                                            print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)
                                        
                        if detection:
                            dv400_l75_det +=1
                            if not is_null(adjustedInc):
                                dv400_l75_det_inc.append(adjustedInc)
                            
                            if not is_null(Lstar):
                                dv400_l75_det_lstar.append(Lstar)
                                
                            if not is_null(azimuth):
                                dv400_l75_det_az.append(azimuth)
                        else:
                            dv400_l75_non +=1
                            if not is_null(adjustedInc):
                                dv400_l75_non_inc.append(adjustedInc)
                            
                            if not is_null(Lstar):
                                dv400_l75_non_lstar.append(Lstar)
                                
                            if not is_null(azimuth):
                                dv400_l75_non_az.append(azimuth)


                    # if the likelihood is greater than 0.5, see if there's a corresponding
                    elif l_used >= 0.5:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_l5.append(dv)
                        
                                if abs(dv) <= max_deltav:
                                    if is_closest_likelihood(Name, Vhel, Lya_v, l_used, other_list, max_deltav):
                                        detection = True
    #                                     Lya_vs_like.pop(Lya_vs_like.index(Lya_v))
                                        if verbose:
                                            print 'Detection <= 0.5L - {0}'.format(Name)
                                            print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)
                        if detection:
                            dv400_l5_det +=1
                            if not is_null(adjustedInc):
                                dv400_l5_det_inc.append(adjustedInc)
                            
                            if not is_null(Lstar):
                                dv400_l5_det_lstar.append(Lstar)
                                
                            if not is_null(azimuth):
                                dv400_l5_det_az.append(azimuth)
                        else:
                            dv400_l5_non +=1
                            if not is_null(adjustedInc):
                                dv400_l5_non_inc.append(adjustedInc)
                            
                            if not is_null(Lstar):
                                dv400_l5_non_lstar.append(Lstar)
                                
                            if not is_null(azimuth):
                                dv400_l5_non_az.append(azimuth)


                    # if the likelihood is greater than 0.1, see if there's a corresponding absorber
                    elif l_used >= 0.1:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_l1.append(dv)
                        
                                if abs(dv) <= max_deltav:
                                    if is_closest_likelihood(Name, Vhel, Lya_v, l_used, other_list, max_deltav):
                                        detection = True
    #                                     Lya_vs_like.pop(Lya_vs_like.index(Lya_v))
                                        if verbose:
                                            print 'Detection <= 0.1L - {0}'.format(Name)
                                            print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)
                                        
                        if detection:
                            dv400_l1_det +=1
                            if not is_null(adjustedInc):
                                dv400_l1_det_inc.append(adjustedInc)
                            
                            if not is_null(Lstar):
                                dv400_l1_det_lstar.append(Lstar)
                                
                            if not is_null(azimuth):
                                dv400_l1_det_az.append(azimuth)
                        else:
                            dv400_l1_non +=1
                            if not is_null(adjustedInc):
                                dv400_l1_non_inc.append(adjustedInc)
                            
                            if not is_null(Lstar):
                                dv400_l1_non_lstar.append(Lstar)
                                
                            if not is_null(azimuth):
                                dv400_l1_non_az.append(azimuth)
            
                        
                    # if the likelihood is greater than 0.05, see if there's a corresponding absorber
                    elif l_used >= 0.05:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_l05.append(dv)
                        
                                if abs(dv) <= max_deltav:
                                    if is_closest_likelihood(Name, Vhel, Lya_v, l_used, other_list, max_deltav):
                                        detection = True
    #                                     Lya_vs_like.pop(Lya_vs_like.index(Lya_v))
                                        if verbose:
                                            print 'Detection <= 0.05L - {0}'.format(Name)
                                            print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)
                                        
                        if detection:
                            dv400_l05_det +=1
                            if not is_null(adjustedInc):
                                dv400_l05_det_inc.append(adjustedInc)
                            
                            if not is_null(Lstar):
                                dv400_l05_det_lstar.append(Lstar)
                                
                            if not is_null(azimuth):
                                dv400_l05_det_az.append(azimuth)
                        else:
                            dv400_l05_non +=1
                            if not is_null(adjustedInc):
                                dv400_l05_non_inc.append(adjustedInc)
                            
                            if not is_null(Lstar):
                                dv400_l05_non_lstar.append(Lstar)
                                
                            if not is_null(azimuth):
                                dv400_l05_non_az.append(azimuth)


                    # if the likelihood is greater than 0.01, see if there's a corresponding absorber
                    elif l_used >= 0.01:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_l01.append(dv)
                        
                                if abs(dv) <= max_deltav:
                                    if is_closest_likelihood(Name, Vhel, Lya_v, l_used, other_list, max_deltav):
                                        detection = True
    #                                     Lya_vs_like.pop(Lya_vs_like.index(Lya_v))
                                        if verbose:
                                            print 'Detection <= 0.01L - {0}'.format(Name)
                                            print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)
                                        
                        if detection:
                            dv400_l01_det +=1
                            if not is_null(adjustedInc):
                                dv400_l01_det_inc.append(adjustedInc)
                            
                            if not is_null(Lstar):
                                dv400_l01_det_lstar.append(Lstar)
                                
                            if not is_null(azimuth):
                                dv400_l01_det_az.append(azimuth)
                        else:
                            dv400_l01_non +=1
                            if not is_null(adjustedInc):
                                dv400_l01_non_inc.append(adjustedInc)
                            
                            if not is_null(Lstar):
                                dv400_l01_non_lstar.append(Lstar)
                                
                            if not is_null(azimuth):
                                dv400_l01_non_az.append(azimuth)


                    # if the likelihood is greater than 0.005, see if there's a corresponding absorber
                    elif l_used >= 0.005:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_l005.append(dv)
                        
                                if abs(dv) <= max_deltav:
                                    if is_closest_likelihood(Name, Vhel, Lya_v, l_used, other_list, max_deltav):
                                        detection = True
    #                                     Lya_vs_like.pop(Lya_vs_like.index(Lya_v))
                                        if verbose:
                                            print 'Detection <= 0.005L - {0}'.format(Name)
                                            print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)
                                        
                        if detection:
                            dv400_l005_det +=1
                            if not is_null(adjustedInc):
                                dv400_l005_det_inc.append(adjustedInc)
                            
                            if not is_null(Lstar): 
                                dv400_l005_det_lstar.append(Lstar)
                                
                            if not is_null(azimuth): 
                                dv400_l005_det_az.append(azimuth)
                        else:
                            dv400_l005_non +=1
                            if not is_null(adjustedInc):
                                dv400_l005_non_inc.append(adjustedInc)
                            
                            if not is_null(Lstar):
                                dv400_l005_non_lstar.append(Lstar)
                                
                            if not is_null(azimuth): 
                                dv400_l005_non_az.append(azimuth)


                    # if the likelihood is greater than 0.001, see if there's a corresponding absorber
                    elif l_used >= 0.001:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_l001.append(dv)
                        
                                if abs(dv) <= max_deltav:
                                    if is_closest_likelihood(Name, Vhel, Lya_v, l_used, other_list, max_deltav):
                                        detection = True
    #                                     Lya_vs_like.pop(Lya_vs_like.index(Lya_v))
                                        if verbose:
                                            print 'Detection <= 0.001L - {0}'.format(Name)
                                            print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)
                                        
                        if detection:
                            dv400_l001_det +=1
                            if not is_null(adjustedInc):
                                dv400_l001_det_inc.append(adjustedInc)
                            
                            if not is_null(Lstar):
                                dv400_l001_det_lstar.append(Lstar)
                                
                            if not is_null(azimuth):
                                dv400_l001_det_az.append(azimuth)
                        else:
                            dv400_l001_non +=1
                            if not is_null(adjustedInc):
                                dv400_l001_non_inc.append(adjustedInc)
                            
                            if not is_null(Lstar):
                                dv400_l001_non_lstar.append(Lstar)
                                
                            if not is_null(azimuth):
                                dv400_l001_non_az.append(azimuth)


                    # if the likelihood is greater than 0.0005, see if there's a corresponding absorber
                    elif l_used >= 0.0005:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_l0005.append(dv)
                        
                                if abs(dv) <= max_deltav:
                                    if is_closest_likelihood(Name, Vhel, Lya_v, l_used, other_list, max_deltav):
                                        detection = True
    #                                     Lya_vs_like.pop(Lya_vs_like.index(Lya_v))
                                        if verbose:
                                            print 'Detection <= 0.0005L - {0}'.format(Name)
                                            print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)
                                        
                        if detection:
                            dv400_l0005_det +=1
                            if not is_null(adjustedInc):
                                dv400_l0005_det_inc.append(adjustedInc)
                            
                            if not is_null(Lstar):
                                dv400_l0005_det_lstar.append(Lstar)
                                
                            if not is_null(azimuth):
                                dv400_l0005_det_az.append(azimuth)
                        else:
                            dv400_l0005_non +=1
                            if not is_null(adjustedInc):
                                dv400_l0005_non_inc.append(adjustedInc)
                            
                            if not is_null(Lstar):
                                dv400_l0005_non_lstar.append(Lstar)
                                
                            if not is_null(azimuth):
                                dv400_l0005_non_az.append(azimuth)


                    # if the likelihood is greater than 0.0001, see if there's a corresponding absorber
                    elif l_used >= 0.0001:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_l0001.append(dv)
                        
                                if abs(dv) <= max_deltav:
                                    if is_closest_likelihood(Name, Vhel, Lya_v, l_used, other_list, max_deltav):
                                        detection = True
    #                                     Lya_vs_like.pop(Lya_vs_like.index(Lya_v))
                                        if verbose:
                                            print 'Detection <= 0.0001L - {0}'.format(Name)
                                            print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)

                        if detection:
                            dv400_l0001_det +=1
                            if not is_null(adjustedInc):
                                dv400_l0001_det_inc.append(adjustedInc)
                            
                            if not is_null(Lstar): 
                                dv400_l0001_det_lstar.append(Lstar)
                                
                            if not is_null(azimuth): 
                                dv400_l0001_det_az.append(azimuth)
                        else:
                            dv400_l0001_non +=1
                            if not is_null(adjustedInc):
                                dv400_l0001_non_inc.append(adjustedInc)
                            
                            if not is_null(Lstar):
                                dv400_l0001_non_lstar.append(Lstar)
                                
                            if not is_null(azimuth): 
                                dv400_l0001_non_az.append(azimuth)
                                
                    else:
                        for v in Lya_vs_like:
                            likelihood_no_galaxies.append(v)

                        
    ##########################################################################################            
                # now do it for imp/R_vir
                #
                # Smaller l_used here should mean higher detection fraction
                
                l_used = impact/R_vir
                                    
                if proceed:
                
                    # if the imp/Rvir is less than 0.25, see if there's a corresponding line
                    if l_used <= 0.25:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_rvir025.append(dv)
                        
                                if abs(dv) <= max_deltav:
                                    if is_closest_impact_rvir(Name, Vhel, Lya_v, l_used, other_list, max_deltav):
                                        detection = True
    #                                     Lya_vs_rvir.pop(Lya_vs_rvir.index(Lya_v))
                                        if verbose:
                                            print 'Detection <= 0.25imp/rvir - {0}'.format(Name)
                                            print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)

                        if detection:
                            dv400_rvir025_det +=1
                        else:
                            dv400_rvir025_non +=1
            
            
                    # if the imp/Rvir is less than 0.5, see if there's a corresponding line
                    elif l_used <= 0.5:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_rvir05.append(dv)
                        
                                if abs(dv) <= max_deltav:
                                    if is_closest_impact_rvir(Name, Vhel, Lya_v, l_used, other_list, max_deltav):
                                        detection = True
    #                                     Lya_vs_rvir.pop(Lya_vs_rvir.index(Lya_v))
                                        if verbose:
                                            print 'Detection <= 0.5imp/rvir - {0}'.format(Name)
                                            print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)

                        if detection:
                            dv400_rvir05_det +=1
                        else:
                            dv400_rvir05_non +=1


                    # if the imp/Rvir is less than 0.75, see if there's a corresponding line
                    elif l_used <= 0.75:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_rvir075.append(dv)
                        
                                if abs(dv) <= max_deltav:
                                    if is_closest_impact_rvir(Name, Vhel, Lya_v, l_used, other_list, max_deltav):
                                        detection = True
    #                                     Lya_vs_rvir.pop(Lya_vs_rvir.index(Lya_v))
                                        if verbose:
                                            print 'Detection <= 0.75imp/rvir - {0}'.format(Name)
                                            print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)
                                        
                        if detection:
                            dv400_rvir075_det +=1
                        else:
                            dv400_rvir075_non +=1


                    # if the imp/Rvir is greater than 1.0, see if there's a corresponding line
                    elif l_used <= 1.0:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_rvir1.append(dv)
                        
                                if abs(dv) <= max_deltav:
                                    if is_closest_impact_rvir(Name, Vhel, Lya_v, l_used, other_list, max_deltav):
                                        detection = True
    #                                     Lya_vs_rvir.pop(Lya_vs_rvir.index(Lya_v))
                                        if verbose:
                                            print 'Detection <= 1.0imp/rvir - {0}'.format(Name)
                                            print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)
                                        
                        if detection:
                            dv400_rvir1_det +=1
                        else:
                            dv400_rvir1_non +=1


                    # if the imp/Rvir is greater than 1.5, see if there's a corresponding line
                    elif l_used <= 1.5:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_rvir15.append(dv)
                        
                                if abs(dv) <= max_deltav:
                                    if is_closest_impact_rvir(Name, Vhel, Lya_v, l_used, other_list, max_deltav):
                                        detection = True
    #                                     Lya_vs_rvir.pop(Lya_vs_rvir.index(Lya_v))
                                        if verbose:
                                            print 'Detection <= 1.5imp/rvir - {0}'.format(Name)
                                            print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)

                        if detection:
                            dv400_rvir15_det +=1
                        else:
                            dv400_rvir15_non +=1


                    # if the imp/Rvir is greater than 2.0, see if there's a corresponding line
                    elif l_used <= 2.0:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_rvir2.append(dv)
                        
                                if abs(dv) <= max_deltav:
                                    if is_closest_impact_rvir(Name, Vhel, Lya_v, l_used, other_list, max_deltav):
                                        detection = True
    #                                     Lya_vs_rvir.pop(Lya_vs_rvir.index(Lya_v))
                                        if verbose:
                                            print 'Detection <= 2.0imp/rvir - {0}'.format(Name)
                                            print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)

                        if detection:
                            dv400_rvir2_det +=1
                        else:
                            dv400_rvir2_non +=1


                    # if the imp/Rvir is greater than 2.5, see if there's a corresponding line
                    elif l_used <= 2.5:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_rvir25.append(dv)
                        
                                if abs(dv) <= max_deltav:
                                    if is_closest_impact_rvir(Name, Vhel, Lya_v, l_used, other_list, max_deltav):
                                        detection = True
    #                                     Lya_vs_rvir.pop(Lya_vs_rvir.index(Lya_v))
                                        if verbose:
                                            print 'Detection <= 2.5imp/rvir - {0}'.format(Name)
                                            print 'Lya_v = {0}, dv = {1}'.format(Lya_v, dv)

                        if detection:
                            dv400_rvir25_det +=1
                        else:
                            dv400_rvir25_non +=1
                    
                    else:
                        for v in Lya_vs_rvir:
                            impact_rvir_no_galaxies.append(v)

    try:
        print 'finished! '
        print
        print 'Detection fraction for 1000 kpc: ', float(dv400_imp1000_det) / (dv400_imp1000_det + dv400_imp1000_non)
        print 'Detection fraction for 750 kpc: ', float(dv400_imp750_det) / (dv400_imp750_det + dv400_imp750_non)
        print 'Detection fraction for 500 kpc: ', float(dv400_imp500_det) / (dv400_imp500_det + dv400_imp500_non)
        print 'Detection fraction for 400 kpc: ', float(dv400_imp400_det) / (dv400_imp400_det + dv400_imp400_non)
        print 'Detection fraction for 300 kpc: ', float(dv400_imp300_det) / (dv400_imp300_det + dv400_imp300_non)
        print 'Detection fraction for 200 kpc: ', float(dv400_imp200_det) / (dv400_imp200_det + dv400_imp200_non)
        print 'Detection fraction for 100 kpc: ', float(dv400_imp100_det) / (dv400_imp100_det + dv400_imp100_non)
        print 'Detection fraction for 50 kpc: ', float(dv400_imp50_det)  / (dv400_imp50_det  + dv400_imp50_non)
        print 'Detection fraction for 25 kpc: ', float(dv400_imp25_det)  / (dv400_imp25_det  + dv400_imp25_non)
        print
        print
        print 'Mean inc for 1000 kpc detections vs non: ',bmean(dv400_imp1000_det_inc), ', ',bmean(dv400_imp1000_non_inc)
        print 'Mean inc for 750 kpc detections vs non: ',bmean(dv400_imp750_det_inc), ', ',bmean(dv400_imp750_non_inc)
        print 'Mean inc for 500 kpc detections vs non: ',bmean(dv400_imp500_det_inc), ', ',bmean(dv400_imp500_non_inc)
        print 'Mean inc for 400 kpc detections vs non: ',bmean(dv400_imp400_det_inc), ', ',bmean(dv400_imp400_non_inc)
        print 'Mean inc for 300 kpc detections vs non: ',bmean(dv400_imp300_det_inc), ', ',bmean(dv400_imp300_non_inc)
        print 'Mean inc for 200 kpc detections vs non: ',bmean(dv400_imp200_det_inc), ', ',bmean(dv400_imp200_non_inc)
        print 'Mean inc for 100 kpc detections vs non: ',bmean(dv400_imp100_det_inc), ', ',bmean(dv400_imp100_non_inc)
        print 'Mean inc for 50 kpc detections vs non: ',bmean(dv400_imp50_det_inc), ', ',bmean(dv400_imp50_non_inc)
        print 'Mean inc for 25 kpc detections vs non: ',bmean(dv400_imp25_det_inc), ', ',bmean(dv400_imp25_non_inc)
        print
        print
        print
        print
        print
        print 'Detection fraction for 0.0001 L: ', float(dv400_l0001_det) / (dv400_l0001_det   + dv400_l0001_non)
        print 'Detection fraction for 0.0005 L: ', float(dv400_l0005_det) / (dv400_l0005_det   + dv400_l0005_non)
        print 'Detection fraction for 0.001 L: ', float(dv400_l001_det) / (dv400_l001_det   + dv400_l001_non)
        print 'Detection fraction for 0.005 L: ', float(dv400_l005_det) / (dv400_l005_det   + dv400_l005_non)
        print 'Detection fraction for 0.01 L: ', float(dv400_l01_det)   / (dv400_l01_det    + dv400_l01_non)
        print 'Detection fraction for 0.05 L: ', float(dv400_l05_det)   / (dv400_l05_det    + dv400_l05_non)
        print 'Detection fraction for 0.1 L: ', float(dv400_l1_det)     / (dv400_l1_det     + dv400_l1_non)
        print 'Detection fraction for 0.5 L: ', float(dv400_l5_det)     / (dv400_l5_det     + dv400_l5_non)
        print 'Detection fraction for 0.75 L: ', float(dv400_l75_det)   / (dv400_l75_det    + dv400_l75_non)
        print
        print
        print
        print 'Mean inc for 0.0001 L detections vs non: ',bmean(dv400_l0001_det_inc), ', ',bmean(dv400_l0001_non_inc)
        print 'Mean inc for 0.0005 L detections vs non: ',bmean(dv400_l0005_det_inc), ', ',bmean(dv400_l0005_non_inc)
        print 'Mean inc for 0.001 L detections vs non: ',bmean(dv400_l001_det_inc), ', ',bmean(dv400_l001_non_inc)
        print 'Mean inc for 0.005 L detections vs non: ',bmean(dv400_l005_det_inc), ', ',bmean(dv400_l005_non_inc)
        print 'Mean inc for 0.01 L detections vs non: ',bmean(dv400_l01_det_inc), ', ',bmean(dv400_l01_non_inc)
        print 'Mean inc for 0.05 L detections vs non: ',bmean(dv400_l05_det_inc), ', ',bmean(dv400_l05_non_inc)
        print 'Mean inc for 0.1 L detections vs non: ',bmean(dv400_l1_det_inc), ', ',bmean(dv400_l1_non_inc)
        print 'Mean inc for 0.5 L detections vs non: ',bmean(dv400_l5_det_inc), ', ',bmean(dv400_l5_non_inc)
        print 'Mean inc for 0.75 L detections vs non: ',bmean(dv400_l75_det_inc), ', ',bmean(dv400_l75_non_inc)
        print
        print
        print '------- Lstar stuff -------'
        print 'dv400_l0001_det_lstar: ', dv400_l0001_det_lstar
        print 'dv400_l0001_non_lstar: ',dv400_l0001_non_lstar
        print
        print 'dv400_l0005_det_lstar: ',dv400_l0005_det_lstar
        print 'dv400_l0005_non_lstar: ',dv400_l0005_non_lstar
        print
        print 'dv400_l001_det_lstar: ',dv400_l001_det_lstar
        print 'dv400_l001_non_lstar: ',dv400_l001_non_lstar
        print
        print 'dv400_l005_det_lstar: ',dv400_l005_det_lstar
        print 'dv400_l005_non_lstar: ',dv400_l005_non_lstar
        print
        print 'dv400_l01_det_lstar: ',dv400_l01_det_lstar
        print 'dv400_l01_non_lstar: ',dv400_l01_non_lstar
        print
        print 'dv400_l05_det_lstar: ',dv400_l05_det_lstar
        print 'dv400_l05_non_lstar: ',dv400_l05_non_lstar
        print
        print 'dv400_l1_det_lstar: ',dv400_l1_det_lstar
        print 'dv400_l1_non_lstar: ',dv400_l1_non_lstar
        print
        print 'dv400_l5_det_lstar: ',dv400_l5_det_lstar
        print 'dv400_l5_non_lstar: ',dv400_l5_non_lstar
        print
        print 'dv400_l75_det_lstar: ',dv400_l75_det_lstar
        print 'dv400_l75_non_lstar: ',dv400_l75_non_lstar
        print
        print
        print 'Mean lstar for 0.0001 L detections vs non: ',bmean(dv400_l0001_det_lstar), ', ',bmean(dv400_l0001_non_lstar)
        print 'Mean lstar for 0.0005 L detections vs non: ',bmean(dv400_l0005_det_lstar), ', ',bmean(dv400_l0005_non_lstar)
        print 'Mean lstar for 0.001 L detections vs non: ',bmean(dv400_l001_det_lstar), ', ',bmean(dv400_l001_non_lstar)
        print 'Mean lstar for 0.005 L detections vs non: ',bmean(dv400_l005_det_lstar), ', ',bmean(dv400_l005_non_lstar)
        print 'Mean lstar for 0.01 L detections vs non: ',bmean(dv400_l01_det_lstar), ', ',bmean(dv400_l01_non_lstar)
        print 'Mean lstar for 0.05 L detections vs non: ',bmean(dv400_l05_det_lstar), ', ',bmean(dv400_l05_non_lstar)
        print 'Mean lstar for 0.1 L detections vs non: ',bmean(dv400_l1_det_lstar), ', ',bmean(dv400_l1_non_lstar)
        print 'Mean lstar for 0.5 L detections vs non: ',bmean(dv400_l5_det_lstar), ', ',bmean(dv400_l5_non_lstar)
        print 'Mean lstar for 0.75 L detections vs non: ',bmean(dv400_l75_det_lstar), ', ',bmean(dv400_l75_non_lstar)
        print
        print
        print
        print 'Detection fraction for 0.25 imp/vir: ', float(dv400_rvir025_det) / (dv400_rvir025_det    + dv400_rvir025_non)
        print 'Detection fraction for 0.5 imp/vir: ', float(dv400_rvir05_det)   / (dv400_rvir05_det     + dv400_rvir05_non)
        print 'Detection fraction for 0.75 imp/vir: ', float(dv400_rvir075_det) / (dv400_rvir075_det    + dv400_rvir075_non)
        print 'Detection fraction for 1.0 imp/vir: ', float(dv400_rvir1_det)    / (dv400_rvir1_det      + dv400_rvir1_non)
        print 'Detection fraction for 1.5 imp/vir: ', float(dv400_rvir15_det)   / (dv400_rvir15_det     + dv400_rvir15_non)
        print 'Detection fraction for 2.0 imp/vir: ', float(dv400_rvir2_det)    / (dv400_rvir2_det      + dv400_rvir2_non)
        print 'Detection fraction for 2.5 imp/vir: ', float(dv400_rvir25_det)   / (dv400_rvir25_det     + dv400_rvir25_non)
        print
        print
        print 'dv400_imp500_det: ',dv400_imp500_det
        print 'dv400_imp500_non: ',dv400_imp500_non
        print
        print 'dv400_imp400_det: ',dv400_imp400_det
        print 'dv400_imp400_non: ',dv400_imp400_non
        print
        print
        print
        print 'dv400_imp1000_det, dv400_imp1000_non = {0}'.format(dv400_imp1000_det, dv400_imp1000_non)
        print 'dv400_imp750_det, dv400_imp750_non = {0}'.format(dv400_imp750_det, dv400_imp750_non)
        print 'dv400_imp500_det, dv400_imp500_non = {0}'.format(dv400_imp500_det, dv400_imp500_non)
        print 'dv400_imp400_det, dv400_imp400_non = {0}'.format(dv400_imp400_det, dv400_imp400_non)
        print 'dv400_imp300_det, dv400_imp300_non = {0}'.format(dv400_imp300_det, dv400_imp300_non)
        print 'dv400_imp200_det, dv400_imp200_non = {0}'.format(dv400_imp200_det, dv400_imp200_non)
        print 'dv400_imp100_det, dv400_imp100_non = {0}'.format(dv400_imp100_det, dv400_imp100_non)
        print 'dv400_imp50_det, dv400_imp50_non = {0}'.format(dv400_imp50_det, dv400_imp50_non)
        print 'dv400_imp25_det, dv400_imp25_non = {0}'.format(dv400_imp25_det, dv400_imp25_non)
        print
        print
        print 'dv400_l0001_det, dv400_l0001_non = {0}'.format(dv400_l0001_det, dv400_l0001_non)
        print 'dv400_l0005_det, dv400_l0005_non: {0}'.format(dv400_l0005_det, dv400_l0005_non)
        print 'dv400_l001_det, dv400_l001_non: {0} '.format(dv400_l001_det, dv400_l001_non)
        print 'dv400_l005_det, dv400_l005_non: {0} '.format(dv400_l005_det, dv400_l005_non)
        print 'dv400_l01_det, dv400_l01_non: {0} '.format(dv400_l01_det, dv400_l01_non)
        print 'dv400_l05_det, dv400_l05_non: {0} '.format(dv400_l05_det, dv400_l05_non)
        print 'dv400_l1_det, dv400_l1_non: {0} '.format(dv400_l1_det, dv400_l1_non)
        print 'dv400_l5_det, dv400_l5_non: {0} '.format(dv400_l5_det, dv400_l5_non)
        print 'dv400_l75_det, dv400_l75_non: {0} '.format(dv400_l75_det, dv400_l75_non)
        print
        print
        
#         print 'no galaxies: '
#         print 'impact: ', impact_no_galaxies
#         print
#         print
#         print 'Likelihood: ',likelihood_no_galaxies
#         print
#         print
#         print 'impact/R_Vir: ',impact_rvir_no_galaxies
#         print
    
    except Exception, e:
        print 'Error: ',e
        print
        
    
    full_dict = {}
    
    # everything lists
    # detections
    full_dict['dv400_all_det_imps'] = dv400_all_det_imps
    full_dict['dv400_all_det_incs'] = dv400_all_det_incs
    full_dict['dv400_all_det_azs'] = dv400_all_det_azs
    full_dict['dv400_all_det_rvirs'] = dv400_all_det_rvirs
    full_dict['dv400_all_det_lstars'] = dv400_all_det_lstars
    full_dict['dv400_all_det_mtypes'] = dv400_all_det_mtypes
    full_dict['dv400_all_det_majdiams'] = dv400_all_det_majdiams
    full_dict['dv400_all_det_bmags'] = dv400_all_det_bmags
    full_dict['dv400_all_det_group_nums'] = dv400_all_det_group_nums
    full_dict['dv400_all_det_group_mems'] = dv400_all_det_group_mems
    full_dict['dv400_all_det_group_dists'] = dv400_all_det_group_dists
    full_dict['dv400_all_det_likelihoods'] = dv400_all_det_likelihoods

    # non-detections
    full_dict['dv400_all_non_imps'] = dv400_all_non_imps
    full_dict['dv400_all_non_incs'] = dv400_all_non_incs
    full_dict['dv400_all_non_azs'] = dv400_all_non_azs
    full_dict['dv400_all_non_rvirs'] = dv400_all_non_rvirs
    full_dict['dv400_all_non_lstars'] = dv400_all_non_lstars
    full_dict['dv400_all_non_mtypes'] = dv400_all_non_mtypes
    full_dict['dv400_all_non_majdiams'] = dv400_all_non_majdiams
    full_dict['dv400_all_non_bmags'] = dv400_all_non_bmags
    full_dict['dv400_all_non_group_nums'] = dv400_all_non_group_nums
    full_dict['dv400_all_non_group_mems'] = dv400_all_non_group_mems
    full_dict['dv400_all_non_group_dists'] = dv400_all_non_group_dists
    full_dict['dv400_all_non_likelihoods'] = dv400_all_non_likelihoods
    
    
    # \Delta v lists for impact parameter
    full_dict['dv_1000'] = dv_1000
    full_dict['dv_750'] = dv_750
    full_dict['dv_500'] = dv_500
    full_dict['dv_400'] = dv_400
    full_dict['dv_300'] = dv_300
    full_dict['dv_200'] = dv_200
    full_dict['dv_100'] = dv_100
    full_dict['dv_50'] = dv_50
    full_dict['dv_25'] = dv_25

    # now impact parameter detection counts
    full_dict['dv400_imp1000_det'] = dv400_imp1000_det
    full_dict['dv400_imp750_det'] = dv400_imp750_det
    full_dict['dv400_imp500_det'] = dv400_imp500_det
    full_dict['dv400_imp400_det'] = dv400_imp400_det
    full_dict['dv400_imp300_det'] = dv400_imp300_det
    full_dict['dv400_imp200_det'] = dv400_imp200_det
    full_dict['dv400_imp100_det'] = dv400_imp100_det
    full_dict['dv400_imp50_det'] = dv400_imp50_det
    full_dict['dv400_imp25_det'] = dv400_imp25_det
    
    # now impact parameter non-detection counts
    full_dict['dv400_imp1000_non'] = dv400_imp1000_non
    full_dict['dv400_imp750_non'] = dv400_imp750_non
    full_dict['dv400_imp500_non'] = dv400_imp500_non
    full_dict['dv400_imp400_non'] = dv400_imp400_non
    full_dict['dv400_imp300_non'] = dv400_imp300_non
    full_dict['dv400_imp200_non'] = dv400_imp200_non
    full_dict['dv400_imp100_non'] = dv400_imp100_non
    full_dict['dv400_imp50_non'] = dv400_imp50_non
    full_dict['dv400_imp25_non'] = dv400_imp25_non

    # now impact parameter detection inc 
    full_dict['dv400_imp1000_det_inc'] = dv400_imp1000_det_inc
    full_dict['dv400_imp750_det_inc'] = dv400_imp750_det_inc
    full_dict['dv400_imp500_det_inc'] = dv400_imp500_det_inc
    full_dict['dv400_imp400_det_inc'] = dv400_imp400_det_inc
    full_dict['dv400_imp300_det_inc'] = dv400_imp300_det_inc
    full_dict['dv400_imp200_det_inc'] = dv400_imp200_det_inc
    full_dict['dv400_imp100_det_inc'] = dv400_imp100_det_inc
    full_dict['dv400_imp50_det_inc'] = dv400_imp50_det_inc
    full_dict['dv400_imp25_det_inc'] = dv400_imp25_det_inc
    
    # now impact parameter non-detection inc
    full_dict['dv400_imp1000_non_inc'] = dv400_imp1000_non_inc
    full_dict['dv400_imp750_non_inc'] = dv400_imp750_non_inc
    full_dict['dv400_imp500_non_inc'] = dv400_imp500_non_inc
    full_dict['dv400_imp400_non_inc'] = dv400_imp400_non_inc
    full_dict['dv400_imp300_non_inc'] = dv400_imp300_non_inc
    full_dict['dv400_imp200_non_inc'] = dv400_imp200_non_inc
    full_dict['dv400_imp100_non_inc'] = dv400_imp100_non_inc
    full_dict['dv400_imp50_non_inc'] = dv400_imp50_non_inc
    full_dict['dv400_imp25_non_inc'] = dv400_imp25_non_inc

    # now impact parameter detection Lstar 
    full_dict['dv400_imp1000_det_lstar'] = dv400_imp1000_det_lstar
    full_dict['dv400_imp750_det_lstar'] = dv400_imp750_det_lstar
    full_dict['dv400_imp500_det_lstar'] = dv400_imp500_det_lstar
    full_dict['dv400_imp400_det_lstar'] = dv400_imp400_det_lstar
    full_dict['dv400_imp300_det_lstar'] = dv400_imp300_det_lstar
    full_dict['dv400_imp200_det_lstar'] = dv400_imp200_det_lstar
    full_dict['dv400_imp100_det_lstar'] = dv400_imp100_det_lstar
    full_dict['dv400_imp50_det_lstar'] = dv400_imp50_det_lstar
    full_dict['dv400_imp25_det_lstar'] = dv400_imp25_det_lstar
    
    # now impact parameter non-detection Lstar
    full_dict['dv400_imp1000_non_lstar'] = dv400_imp1000_non_lstar
    full_dict['dv400_imp750_non_lstar'] = dv400_imp750_non_lstar
    full_dict['dv400_imp500_non_lstar'] = dv400_imp500_non_lstar
    full_dict['dv400_imp400_non_lstar'] = dv400_imp400_non_lstar
    full_dict['dv400_imp300_non_lstar'] = dv400_imp300_non_lstar
    full_dict['dv400_imp200_non_lstar'] = dv400_imp200_non_lstar
    full_dict['dv400_imp100_non_lstar'] = dv400_imp100_non_lstar
    full_dict['dv400_imp50_non_lstar'] = dv400_imp50_non_lstar
    full_dict['dv400_imp25_non_lstar'] = dv400_imp25_non_lstar


    # now impact parameter detection az 
    full_dict['dv400_imp1000_det_az'] = dv400_imp1000_det_az
    full_dict['dv400_imp750_det_az'] = dv400_imp750_det_az
    full_dict['dv400_imp500_det_az'] = dv400_imp500_det_az
    full_dict['dv400_imp400_det_az'] = dv400_imp400_det_az
    full_dict['dv400_imp300_det_az'] = dv400_imp300_det_az
    full_dict['dv400_imp200_det_az'] = dv400_imp200_det_az
    full_dict['dv400_imp100_det_az'] = dv400_imp100_det_az
    full_dict['dv400_imp50_det_az'] = dv400_imp50_det_az
    full_dict['dv400_imp25_det_az'] = dv400_imp25_det_az
    
    # now impact parameter non-detection az
    full_dict['dv400_imp1000_non_az'] = dv400_imp1000_non_az
    full_dict['dv400_imp750_non_az'] = dv400_imp750_non_az
    full_dict['dv400_imp500_non_az'] = dv400_imp500_non_az
    full_dict['dv400_imp400_non_az'] = dv400_imp400_non_az
    full_dict['dv400_imp300_non_az'] = dv400_imp300_non_az
    full_dict['dv400_imp200_non_az'] = dv400_imp200_non_az
    full_dict['dv400_imp100_non_az'] = dv400_imp100_non_az
    full_dict['dv400_imp50_non_az'] = dv400_imp50_non_az
    full_dict['dv400_imp25_non_az'] = dv400_imp25_non_az




    # now for likelihood thresholds
    full_dict['dv_l0001'] = dv_l0001
    full_dict['dv_l0005'] = dv_l0005
    full_dict['dv_l001'] = dv_l001
    full_dict['dv_l005'] = dv_l005
    full_dict['dv_l01'] = dv_l01
    full_dict['dv_l05'] = dv_l05
    full_dict['dv_l1'] = dv_l1
    full_dict['dv_l5'] = dv_l5
    full_dict['dv_l75'] = dv_l75
    
    # now for likelihood detections
    full_dict['dv400_l0001_det'] = dv400_l0001_det
    full_dict['dv400_l0005_det'] = dv400_l0005_det
    full_dict['dv400_l001_det'] = dv400_l001_det
    full_dict['dv400_l005_det'] = dv400_l005_det
    full_dict['dv400_l01_det'] = dv400_l01_det
    full_dict['dv400_l05_det'] = dv400_l05_det
    full_dict['dv400_l1_det'] = dv400_l1_det
    full_dict['dv400_l5_det'] = dv400_l5_det
    full_dict['dv400_l75_det'] = dv400_l75_det
    
    # now for likelihood non-detections
    full_dict['dv400_l0001_non'] = dv400_l0001_non
    full_dict['dv400_l0005_non'] = dv400_l0005_non
    full_dict['dv400_l001_non'] = dv400_l001_non
    full_dict['dv400_l005_non'] = dv400_l005_non
    full_dict['dv400_l01_non'] = dv400_l01_non
    full_dict['dv400_l05_non'] = dv400_l05_non
    full_dict['dv400_l1_non'] = dv400_l1_non
    full_dict['dv400_l5_non'] = dv400_l5_non
    full_dict['dv400_l75_non'] = dv400_l75_non
    
    # now for likelihood detection incs
    full_dict['dv400_l0001_det_inc'] = dv400_l0001_det_inc
    full_dict['dv400_l0005_det_inc'] = dv400_l0005_det_inc
    full_dict['dv400_l001_det_inc'] = dv400_l001_det_inc
    full_dict['dv400_l005_det_inc'] = dv400_l005_det_inc
    full_dict['dv400_l01_det_inc'] = dv400_l01_det_inc
    full_dict['dv400_l05_det_inc'] = dv400_l05_det_inc
    full_dict['dv400_l1_det_inc'] = dv400_l1_det_inc
    full_dict['dv400_l5_det_inc'] = dv400_l5_det_inc
    full_dict['dv400_l75_det_inc'] = dv400_l75_det_inc
    
    # now for likelihood non-detection incs
    full_dict['dv400_l0001_non_inc'] = dv400_l0001_non_inc
    full_dict['dv400_l0005_non_inc'] = dv400_l0005_non_inc
    full_dict['dv400_l001_non_inc'] = dv400_l001_non_inc
    full_dict['dv400_l005_non_inc'] = dv400_l005_non_inc
    full_dict['dv400_l01_non_inc'] = dv400_l01_non_inc
    full_dict['dv400_l05_non_inc'] = dv400_l05_non_inc
    full_dict['dv400_l1_non_inc'] = dv400_l1_non_inc
    full_dict['dv400_l5_non_inc'] = dv400_l5_non_inc
    full_dict['dv400_l75_non_inc'] = dv400_l75_non_inc
    
    # now for likelihood detection Lstars
    full_dict['dv400_l0001_det_lstar'] = dv400_l0001_det_lstar
    full_dict['dv400_l0005_det_lstar'] = dv400_l0005_det_lstar
    full_dict['dv400_l001_det_lstar'] = dv400_l001_det_lstar
    full_dict['dv400_l005_det_lstar'] = dv400_l005_det_lstar
    full_dict['dv400_l01_det_lstar'] = dv400_l01_det_lstar
    full_dict['dv400_l05_det_lstar'] = dv400_l05_det_lstar
    full_dict['dv400_l1_det_lstar'] = dv400_l1_det_lstar
    full_dict['dv400_l5_det_lstar'] = dv400_l5_det_lstar
    full_dict['dv400_l75_det_lstar'] = dv400_l75_det_lstar
    
    # now for likelihood non-detection Lstars
    full_dict['dv400_l0001_non_lstar'] = dv400_l0001_non_lstar
    full_dict['dv400_l0005_non_lstar'] = dv400_l0005_non_lstar
    full_dict['dv400_l001_non_lstar'] = dv400_l001_non_lstar
    full_dict['dv400_l005_non_lstar'] = dv400_l005_non_lstar
    full_dict['dv400_l01_non_lstar'] = dv400_l01_non_lstar
    full_dict['dv400_l05_non_lstar'] = dv400_l05_non_lstar
    full_dict['dv400_l1_non_lstar'] = dv400_l1_non_lstar
    full_dict['dv400_l5_non_lstar'] = dv400_l5_non_lstar
    full_dict['dv400_l75_non_lstar'] = dv400_l75_non_lstar
    
    # now for likelihood detection az
    full_dict['dv400_l0001_det_az'] = dv400_l0001_det_az
    full_dict['dv400_l0005_det_az'] = dv400_l0005_det_az
    full_dict['dv400_l001_det_az'] = dv400_l001_det_az
    full_dict['dv400_l005_det_az'] = dv400_l005_det_az
    full_dict['dv400_l01_det_az'] = dv400_l01_det_az
    full_dict['dv400_l05_det_az'] = dv400_l05_det_az
    full_dict['dv400_l1_det_az'] = dv400_l1_det_az
    full_dict['dv400_l5_det_az'] = dv400_l5_det_az
    full_dict['dv400_l75_det_az'] = dv400_l75_det_az
    
    # now for likelihood non-detection az
    full_dict['dv400_l0001_non_az'] = dv400_l0001_non_az
    full_dict['dv400_l0005_non_az'] = dv400_l0005_non_az
    full_dict['dv400_l001_non_az'] = dv400_l001_non_az
    full_dict['dv400_l005_non_az'] = dv400_l005_non_az
    full_dict['dv400_l01_non_az'] = dv400_l01_non_az
    full_dict['dv400_l05_non_az'] = dv400_l05_non_az
    full_dict['dv400_l1_non_az'] = dv400_l1_non_az
    full_dict['dv400_l5_non_az'] = dv400_l5_non_az
    full_dict['dv400_l75_non_az'] = dv400_l75_non_az
    
    
    
    
    # now for imp/rvir thresholds
    full_dict['dv_rvir025'] = dv_rvir025
    full_dict['dv_rvir05'] = dv_rvir05
    full_dict['dv_rvir075'] = dv_rvir075
    full_dict['dv_rvir1'] = dv_rvir1
    full_dict['dv_rvir15'] = dv_rvir15
    full_dict['dv_rvir2'] = dv_rvir2
    full_dict['dv_rvir25'] = dv_rvir25

    # now for imp/rvir detections
    full_dict['dv400_rvir025_det'] = dv400_rvir025_det
    full_dict['dv400_rvir05_det']  = dv400_rvir05_det
    full_dict['dv400_rvir075_det'] = dv400_rvir075_det
    full_dict['dv400_rvir1_det'] = dv400_rvir1_det
    full_dict['dv400_rvir15_det'] = dv400_rvir15_det
    full_dict['dv400_rvir2_det'] = dv400_rvir2_det
    full_dict['dv400_rvir25_det'] = dv400_rvir25_det

    # now for imp/rvir non-detections
    full_dict['dv400_rvir025_non'] = dv400_rvir025_non
    full_dict['dv400_rvir05_non']  = dv400_rvir05_non
    full_dict['dv400_rvir075_non'] = dv400_rvir075_non
    full_dict['dv400_rvir1_non'] = dv400_rvir1_non
    full_dict['dv400_rvir15_non'] = dv400_rvir15_non
    full_dict['dv400_rvir2_non'] = dv400_rvir2_non
    full_dict['dv400_rvir25_non'] = dv400_rvir25_non
    
    
    # no galaxy detections
    full_dict['impact_no_galaxies'] = impact_no_galaxies
    full_dict['likelihood_no_galaxies'] = likelihood_no_galaxies
    full_dict['impact_rvir_no_galaxies'] = impact_rvir_no_galaxies
    
##########################################################################################
##########################################################################################
##########################################################################################
    
    pickle.dump(full_dict, detection_fraction)
    detection_fraction.close()

    theFile.close()
    print
    print 'Done!'
    print
    print 'Filename is: ',detection_fraction_filename
    print
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    