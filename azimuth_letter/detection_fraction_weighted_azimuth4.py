#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)


$Id: detection_fraction_weighted_azimuth.py v 1.0 02/01/2019

Create a detection fraction vs azimuth plot where the mean azimuth in each bin is weighted
by the impact parameter of each galaxy


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


def return_bootstrap_errors(x, reps):
    # takes in x as an array of values, calculates bootstrap mean and median errors based 
    # on 'reps' number of iterations
    
    x = np.array(x)
    lenx = len(x)

    if lenx != 0:
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
        
    else:
        mean_err, median_err = -99, -99
        
    return mean_err, median_err




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

    # --- correlation options
    maxSep = 1000.
    agnSeparation = 4000.
    minVcorr = 450.
    minSize = 0.
    max_deltav = 400.
#     rigor = 5
    
    # --- set Lstar min and max values. Set Lstar_max = False to disregard
    Lstar_min = 0
    Lstar_max = False # set to False to not require an Lstar measurement
    
    # --- maximum velocity cut
    vcut = False
    
#     min_EW = 300
    min_EW = 0

    # --- minimum diameter. False to bypass, a number > 0 to impose a diameter measurement needed
#     d_min = 0.01 # must have a diameter > 0.01
#     d_min = False # does not need a diameter measurement at all
#     d_min = 0.01
    d_min = 0.10
    d_min_label = '01'

    # --- maximum velocity difference between galaxy and absorber
    vel_dif_limit = 400
    
    # --- where is the rvir_ratio limit for inside vs outside?
    rvir_limit = 1.5
    
    # --- sort based on likelihood cus instead of the regular one?
    use_likelihood_cus = False
    
    # --- double l if impact <= 1 R_vir?
    double_l_within_rvir = False
    
    # --- print a ton of shit out for testing if True
    verbose = False
    
    # --- cutoff - will only return this many targers - use for testing together with verbose. 
    cutoff = 1000000.
    
    # target_sightline can be either a single target name or 'all'
#     target_sightline = 'MRK876'
    target_sightline = 'all'
#     target_sightline = 'US2816'
#     target_sightline = 'Zw535.012'


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
        detection_fraction_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_az_minEW0_dminFalse_Lstar0-False.p'

        saveDirectory = '/Users/frenchd/Research/test/distributions/'
        
    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    
    theFile = open(filename,'rU')
    reader = csv.DictReader(theFile)
    
    # --- open the pickle files
#     detection_fraction = open(detection_fraction_filename,'wt')

    
##########################################################################################
    # --- grab the target coordinates
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
    
    # --- loop through the measurement list and make a dictionary containing a list of all the
    # --- lines and their velocity, EW, Na and b-parameter for each target (the dict key)
    
    target_dict = {}
    
    passing_galaxy_dict = {}
    
    zero_50_az = []
    fifty_100_az = []
    one_200_az = []
    two_300_az = []
    three_400_az = []
    gr_400_az = []

    zero_50_inc = []
    fifty_100_inc = []
    one_200_inc = []
    two_300_inc = []
    three_400_inc = []
    gr_400_inc = []

    zero_50_impact = []
    fifty_100_impact = []
    one_200_impact = []
    two_300_impact = []
    three_400_impact = []
    gr_400_impact = []
    
        
    # --- these are for the ones with an absorber nearby
    az_0_15 = []
    az_15_30 = []
    az_30_45 = []
    az_45_60 = []
    az_60_75 = []
    az_75_90 = []
    
    az_0_15_impact = []
    az_15_30_impact = []
    az_30_45_impact = []
    az_45_60_impact = []
    az_60_75_impact = []
    az_75_90_impact = []
    
    az_0_15_w = []
    az_15_30_w = []
    az_30_45_w = []
    az_45_60_w = []
    az_60_75_w = []
    az_75_90_w = []
    
    az_0_15_v = []
    az_15_30_v = []
    az_30_45_v = []
    az_45_60_v = []
    az_60_75_v = []
    az_75_90_v = []
    
    Names = []
    Vhels = []
    impacts = []
    R_virs = []
    Lstars = []
    MajDiams = []
    azimuths =  []
    Ws = []
    likelihoods = []    
    
    
    # --- these are all the rest
    all_az_0_15 = []
    all_az_15_30 = []
    all_az_30_45 = []
    all_az_45_60 = []
    all_az_60_75 = []
    all_az_75_90 = []
    
    all_az_0_15_impact = []
    all_az_15_30_impact = []
    all_az_30_45_impact = []
    all_az_45_60_impact = []
    all_az_60_75_impact = []
    all_az_75_90_impact = []
    
    all_az_0_15_R_vir = []
    all_az_15_30_R_vir = []
    all_az_30_45_R_vir = []
    all_az_45_60_R_vir = []
    all_az_60_75_R_vir = []
    all_az_75_90_R_vir = []
    
    all_az_0_15_v = []
    all_az_15_30_v = []
    all_az_30_45_v = []
    all_az_45_60_v = []
    all_az_60_75_v = []
    all_az_75_90_v = []
    
    all_names = []
    all_Vhels = []
    all_azimuths = []
    all_incs = []
    all_impacts = []
    all_R_virs = []
    all_Lstars = []
    all_MajDiams = []
    
    # --- detection lists
    inside_det_names = []
    inside_det_Vhels = []
    inside_det_azimuths = []
    inside_det_incs = []
    inside_det_impacts = []
    inside_det_R_virs = []
    inside_det_Lstars = []
    inside_det_MajDiams = []
    
    outside_det_names = []
    outside_det_Vhels = []
    outside_det_azimuths = []
    outside_det_incs = []
    outside_det_impacts = []
    outside_det_R_virs = []
    outside_det_Lstars = []
    outside_det_MajDiams = []
    
    # --- non-detection lists
    inside_non_names = []
    inside_non_Vhels = []
    inside_non_azimuths = []
    inside_non_incs = []
    inside_non_impacts = []
    inside_non_R_virs = []
    inside_non_Lstars = []
    inside_non_MajDiams = []
    
    outside_non_names = []
    outside_non_Vhels = []
    outside_non_azimuths = []
    outside_non_incs = []
    outside_non_impacts = []
    outside_non_R_virs = []
    outside_non_Lstars = []
    outside_non_MajDiams = []
    
    
    absorber_dict = {}
    
    # --- these are the dictionaries that will hold lists of highest likelihood and 
    # --- lowest impact parameter results
    L_dict = {}
    impact_dict = {}
    ratio_dict = {}

    impact_dict_non = {}
    ratio_dict_non = {}
    

    L_dict['Names'] = []
    L_dict['Vhels'] = []
    L_dict['impacts'] = []
    L_dict['R_virs'] = []
    L_dict['Lstars'] = []
    L_dict['MajDiams'] = []
    L_dict['azimuths'] = []
    L_dict['likelihoods'] = []
    L_dict['Ws'] = []
    L_dict['vs'] = []
    
    # for detections
    impact_dict['Names'] = []
    impact_dict['Vhels'] = []
    impact_dict['impacts'] = []
    impact_dict['R_virs'] = []
    impact_dict['Lstars'] = []
    impact_dict['MajDiams'] = []
    impact_dict['azimuths'] = []
    impact_dict['likelihoods'] = []
    impact_dict['Ws'] = []
    impact_dict['vs'] = []
    
    ratio_dict['Names'] = []
    ratio_dict['Vhels'] = []
    ratio_dict['impacts'] = []
    ratio_dict['R_virs'] = []
    ratio_dict['Lstars'] = []
    ratio_dict['MajDiams'] = []
    ratio_dict['azimuths'] = []
    ratio_dict['likelihoods'] = []
    ratio_dict['Ws'] = []
    ratio_dict['vs'] = []


    print
    print 'Starting loop!'
    print
    
    total = 1551
    counter = 0
    stopCount = 500000
    found_target = False
    for i in reader:
        counter +=1
        target = i['target']

        if target_sightline == 'all' or target == target_sightline:
                
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
            else:
                Lya_v = float(Lya_v)

            if not isNumber(Lya_W):
                Lya_W = 'x'
            else:
                Lya_W = float(Lya_W)
            
            if not isNumber(Na):
                Na = 'x'
            else:
                Na = float(Na)
            
            if not isNumber(b):
                b = 'x'
            else:
                b = float(b)


            if isNumber(Lya_W):
                if float(Lya_W) < min_EW:
                    Lya_v = 'x'
                    Lya_W = 'x'
                    Na = 'x'
                    b = 'x'

        
            proceed = False
            # --- 'lya' are properly identified lines, 'maybe' are for possible lines, not yet
            # --- identified, and 'yes' are for identified targets with no Lya in the spectra
            if bfind(identified.lower(), 'lya') or bfind(identified.lower(), 'maybe') or bfind(identified.lower(), 'yes'):
                proceed = True

            if counter >= stopCount:
                proceed = False
            
                
            if Lya_v == 'x' or Lya_W == 'x':
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
                
                    # --- add this line velocity to the list
                    Lya_vs = info['Lya_vs']
                    Lya_vs.append(Lya_v)
                    info['Lya_vs'] = Lya_vs
                
                    # --- add line EW
                    Lya_Ws = info['Lya_Ws']
                    Lya_Ws.append(Lya_W)
                    info['Lya_Ws'] = Lya_Ws
                
                    # --- add line Na 
                    Nas = info['Nas']
                    Nas.append(Na)
                    info['Nas'] = Nas
                
                    # --- add line b 
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
    
    
    counter = 0
    # --- now do the galaxy correlation based on the targets and absorbers collected above
    for target in target_dict:
        counter +=1
        sys.stdout.write("\r Percent Complete: {0} / {1}".format(counter, total))
        sys.stdout.flush()
        print 'on target = ',target
    
        # --- the collected info for this target
        info = target_dict[target]
        Lya_vs = info['Lya_vs']
        Lya_Ws = info['Lya_Ws']
        Nas = info['Nas']
        bs = info['bs']
        z_target = info['z_target']
        RA_target = info['RA_target']
        Dec_target = info['Dec_target']


        # --- do the correlation: returns a dictionary where the keys are all the 
        # --- correlated galaxies and the values are a dictionary of all the galaxy info        
        if counter < cutoff:
            correlation = correlateSingle.correlateTarget(target, maxSep, agnSeparation, minVcorr, minSize, slow=False, searchAll=False, returnDict=True)
        else:
            break

        other_list = []
        
        if verbose:
            print
            print 'target = {0}  ----'.format(target)
            
        # --- make a list of all galaxies in this correlation, decide which galaxies pass the 
        # --- given criteria for inclusion
        
        # --- the list of galaxy names which pass the test
        passing_galaxies = []
        L_list = []
        impact_list = []
        ratio_list = []

    
        # --- each 'c' here is the name of a galaxy returned by the correlation
        for c in correlation:
            print 'c: ',c
            passing = True
            target = correlation[c]['target']
            z_target = correlation[c]['z_target']
            RA_target = correlation[c]['RA_target']
            Dec_target = correlation[c]['Dec_target']
            
            Name = correlation[c]['Name']
            RAdeg = correlation[c]['RAdeg']
            DEdeg = correlation[c]['DEdeg']
            Vhel = correlation[c]['Vhel']
            vcorr = correlation[c]['vcorr']

            impact = correlation[c]['impact']
            PA = correlation[c]['PA']
            azimuth = correlation[c]['azimuth']
            inc = correlation[c]['inc']
            adjustedInc = correlation[c]['adjustedInc']
            MajDiam = correlation[c]['MajDiam']
            R_vir = correlation[c]['R_vir']
            MType = correlation[c]['MType']
            bestDist = correlation[c]['bestDist']
            e_bestDist = correlation[c]['e_bestDist']
            
#             group_num = correlation[c]['group_num']
#             group_mem = correlation[c]['group_mem']
#             group_dist = correlation[c]['group_dist']
            
            Lstar_med = correlation[c]['Lstar_med']
            e_Lstar_med = correlation[c]['e_Lstar_med']
            Lstar_sdss = correlation[c]['Lstar_sdss']
            e_Lstar_sdss = correlation[c]['e_Lstar_sdss']
            Bmag = correlation[c]['Bmag']
            Bmag_sdss = correlation[c]['Bmag_sdss']


            if str(Lstar_med) == '-99.99':
                Lstar = Lstar_sdss
            else:
                Lstar = Lstar_med
                
            # --- restrict galaxies based on Lstar min and max?
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
            
            # --- restrict galaxies based on velocity?
            if vcut:
                if Vhel > vcut:
                    passing = False
            
            # --- restrict based on size?
            if d_min:
                if isNumber(MajDiam):
                    if float(MajDiam) < d_min:
                        passing = False
                else:
                    passing = False

            if passing:
                other_list.append([float(impact), float(Vhel), Name, float(R_vir)])
                passing_galaxies.append(Name)

                # --- convert to floats
                Vhel = float(Vhel)
                impact = float(impact)
                R_vir = float(R_vir)
                Lstar = float(Lstar)
                MajDiam = float(MajDiam)
                azimuth = float(azimuth)
                inc = float(inc)

                # --- add needed values to the corresponding lists
                all_names.append(Name)
                all_Vhels.append(Vhel)
                all_azimuths.append(azimuth)
                all_incs.append(inc)
                all_impacts.append(impact)
                all_R_virs.append(R_vir)
                all_Lstars.append(Lstar)
                all_MajDiams.append(MajDiam)

                
                # --- get the distribution of az, vhel, r_vir, impact for all galaxies
                # --- within 1 Mpc of target sightline
                if azimuth <= 15:
                    all_az_0_15.append(azimuth)
                    all_az_0_15_v.append(Vhel)
                    all_az_0_15_R_vir.append(R_vir)
                    all_az_0_15_impact.append(impact)

                elif azimuth > 15 and azimuth <= 30:
                    all_az_15_30.append(azimuth)
                    all_az_15_30_R_vir.append(R_vir)
                    all_az_15_30_v.append(Vhel)
                    all_az_15_30_impact.append(impact)

                elif azimuth > 30 and azimuth <= 45:
                    all_az_30_45.append(azimuth)
                    all_az_30_45_R_vir.append(R_vir)
                    all_az_30_45_v.append(Vhel)
                    all_az_30_45_impact.append(impact)

                elif azimuth > 45 and azimuth <= 60:
                    all_az_45_60.append(azimuth)
                    all_az_45_60_R_vir.append(R_vir)
                    all_az_45_60_v.append(Vhel)
                    all_az_45_60_impact.append(impact)

                elif azimuth > 60 and azimuth <= 75:
                    all_az_60_75.append(azimuth)
                    all_az_60_75_R_vir.append(R_vir)
                    all_az_60_75_v.append(Vhel)
                    all_az_60_75_impact.append(impact)

                elif azimuth > 75 and azimuth <= 90:
                    all_az_45_60.append(azimuth)
                    all_az_45_60_R_vir.append(R_vir)
                    all_az_45_60_v.append(Vhel)
                    all_az_45_60_impact.append(impact)

                else:
                    print "Error, azimuth={0} does not fit into any bins".format(azimuth)


                # --- inside or outside of rvir_limit? add to corresponding lists 
                rvir_ratio = impact/R_vir
                
                inside_detection = False
                outside_detection = False
                
                # --- loop through the absorber list for this target
                for w, v in zip(Lya_Ws, Lya_vs):
                
                    # --- velocity difference
                    vel_dif = float(Vhel) - float(v)
                    
                    # --- velocity restriction
                    if abs(vel_dif) <= vel_dif_limit:

                        if float(w) >= min_EW:
                            likelihood = calculate_likelihood(impact, R_vir, vel_dif)
                            rvir_ratio = float(impact)/float(R_vir)

                            # --- add this absorber-galaxy system to the all lists
                            Names.append(Name)
                            Vhels.append(Vhel)
                            impacts.append(impact)
                            R_virs.append(R_vir)
                            Lstars.append(Lstar)
                            MajDiams.append(MajDiam)
                            azimuths.append(azimuth)
                            likelihoods.append(likelihood)
                            Ws.append(w)
                            

                            system_dict =  {'Name':Name,
                                            'Vhel':Vhel,
                                            'impact':impact,
                                            'R_vir':R_vir,
                                            'Lstar':Lstar,
                                            'MajDiam':MajDiam,
                                            'azimuth':azimuth,
                                            'likelihood':likelihood,
                                            'W':w,
                                            'v':v}
                            
                            L_list.append([likelihood,system_dict])
                            impact_list.append([impact,system_dict])
                            ratio_list.append([rvir_ratio,system_dict])

                            # --- azimuth bins
                            if azimuth <= 15:
                                az_0_15.append(azimuth)
                                az_0_15_w.append(w)
                                az_0_15_v.append(v)
                                az_0_15_impact.append(impact)

                            elif azimuth > 15 and azimuth <= 30:
                                az_15_30.append(azimuth)
                                az_15_30_w.append(w)
                                az_15_30_v.append(v)
                                az_15_30_impact.append(impact)

                            elif azimuth > 30 and azimuth <= 45:
                                az_30_45.append(azimuth)
                                az_30_45_w.append(w)
                                az_30_45_v.append(v)
                                az_30_45_impact.append(impact)

                            elif azimuth > 45 and azimuth <= 60:
                                az_45_60.append(azimuth)
                                az_45_60_w.append(w)
                                az_45_60_v.append(v)
                                az_45_60_impact.append(impact)

                            elif azimuth > 60 and azimuth <= 75:
                                az_60_75.append(azimuth)
                                az_60_75_w.append(w)
                                az_60_75_v.append(v)
                                az_60_75_impact.append(impact)

                            elif azimuth > 75 and azimuth <= 90:
                                az_45_60.append(azimuth)
                                az_45_60_w.append(w)
                                az_45_60_v.append(v)
                                az_45_60_impact.append(impact)

                            else:
                                print "Error, azimuth={0} does not fit into any bins".format(azimuth)


                            #  --- EW bins
                            if w >= 400:
                                gr_400_az.append(azimuth)
                                gr_400_inc.append(inc)
                                gr_400_impact.append(impact)
                            
                            elif w < 400 and  w >= 300:
                                three_400_az.append(azimuth)
                                three_400_inc.append(inc)
                                three_400_impact.append(impact)
                            
                            elif w < 300 and  w >= 200:
                                two_300_az.append(azimuth)
                                two_300_inc.append(inc)
                                two_300_impact.append(impact)
                            
                            elif w < 200 and w >= 100:
                                one_200_az.append(azimuth)
                                one_200_inc.append(inc)
                                one_200_impact.append(impact)
                            
                            elif w < 100 and w >= 50:
                                fifty_100_az.append(azimuth)
                                fifty_100_inc.append(inc)
                                fifty_100_impact.append(impact)
                        
                            else:
                                zero_50_az.append(azimuth)
                                zero_50_inc.append(inc)
                                zero_50_impact.append(impact)
                                
                            if rvir_ratio <= rvir_limit:
                                inside_detection = True
                            else:
                                outside_detection = True

                
                if inside_detection:
                    inside_det_names.append(Name)
                    inside_det_vhels.append(Vhel)
                    inside_det_impacts.append(impact)
                    inside_det_r_virs.append(R_vir)
                    inside_det_lstars.append(Lstar)
                    inside_det_majDiams.append(MajDiam)
                    inside_det_azimuths.append(azimuth)
                else:
                    inside_non_names.append(Name)
                    inside_non_vhels.append(Vhel)
                    inside_non_impacts.append(impact)
                    inside_non_r_virs.append(R_vir)
                    inside_non_lstars.append(Lstar)
                    inside_non_majDiams.append(MajDiam)
                    inside_non_azimuths.append(azimuth)
                
                if outside_detection:
                    outside_det_names.append(Name)
                    outside_det_vhels.append(Vhel)
                    outside_det_impacts.append(impact)
                    outside_det_r_virs.append(R_vir)
                    outside_det_lstars.append(Lstar)
                    outside_det_majDiams.append(MajDiam)
                    outside_det_azimuths.append(azimuth)
                else:
                    outside_non_names.append(Name)
                    outside_non_vhels.append(Vhel)
                    outside_non_impacts.append(impact)
                    outside_non_r_virs.append(R_vir)
                    outside_non_lstars.append(Lstar)
                    outside_non_majDiams.append(MajDiam)
                    outside_non_azimuths.append(azimuth)

                
                
                if verbose:
                    print 'Name = {0}'.format(Name)
                    print 'impact = {0} , vhel = {1}'.format(impact, Vhel)
                    print
        
        # --- now sort the impact and likelihood lists
        L_list.sort(reverse=True)
        impact_list.sort(reverse=False)
        ratio_list.sort(reverse=False)

        
        # --- append the highest likelihood results to the overall result dictionary
        if len(L_list) >0:
            print "L_list: ",L_list
            print
            print "L_list[0]: ",L_list[0]
            print
            L_dict['Names'].append(L_list[0][1]['Name'])
            L_dict['Vhels'].append(L_list[0][1]['Vhel'])
            L_dict['impacts'].append(L_list[0][1]['impact'])
            L_dict['R_virs'].append(L_list[0][1]['R_vir'])
            L_dict['Lstars'].append(L_list[0][1]['Lstar'])
            L_dict['MajDiams'].append(L_list[0][1]['MajDiam'])
            L_dict['azimuths'].append(L_list[0][1]['azimuth'])
            L_dict['likelihoods'].append(L_list[0][1]['likelihood'])
            L_dict['Ws'].append(L_list[0][1]['W'])
            L_dict['vs'].append(L_list[0][1]['v'])

        # --- append the lowest impact results to the overall result dictionary
        if len(impact_list) > 0:
            impact_dict['Names'].append(impact_list[0][1]['Name'])
            impact_dict['Vhels'].append(impact_list[0][1]['Vhel'])
            impact_dict['impacts'].append(impact_list[0][1]['impact'])
            impact_dict['R_virs'].append(impact_list[0][1]['R_vir'])
            impact_dict['Lstars'].append(impact_list[0][1]['Lstar'])
            impact_dict['MajDiams'].append(impact_list[0][1]['MajDiam'])
            impact_dict['azimuths'].append(impact_list[0][1]['azimuth'])
            impact_dict['likelihoods'].append(impact_list[0][1]['likelihood'])
            impact_dict['Ws'].append(impact_list[0][1]['W'])
            impact_dict['vs'].append(impact_list[0][1]['v'])
        
        
         # --- append the lowest rvir_ratio results to the dictionary if the next closest
         # --- is twice as far
        if len(ratio_list) >1:
            try:
                if 2*ratio_list[0] < ratio_list[1]:
                    ratio_dict['Names'].append(ratio_list[0][1]['Name'])
                    ratio_dict['Vhels'].append(ratio_list[0][1]['Vhel'])
                    ratio_dict['impacts'].append(ratio_list[0][1]['impact'])
                    ratio_dict['R_virs'].append(ratio_list[0][1]['R_vir'])
                    ratio_dict['Lstars'].append(ratio_list[0][1]['Lstar'])
                    ratio_dict['MajDiams'].append(ratio_list[0][1]['MajDiam'])
                    ratio_dict['azimuths'].append(ratio_list[0][1]['azimuth'])
                    ratio_dict['likelihoods'].append(ratio_list[0][1]['likelihood'])
                    ratio_dict['Ws'].append(ratio_list[0][1]['W'])
                    ratio_dict['vs'].append(ratio_list[0][1]['v'])
            except Exception, e:
                print 'empty'

        
    ################
    # --- finished with all loops at this point
    
    # --- plot the impact parameter histograms first


##########################################################################################
##########################################################################################

    # --- plot the azimuth distributions -> means and impact parameter weighted means
    plot_weighted_azimuth_dist = False
    
    if plot_weighted_azimuth_dist:
        # --- start creating figure
        fig = figure(figsize=(7.7,5.7))
        ax1 = fig.add_subplot(211)

        # --- colors for plotting
        color_purple = '#7570b3'
        color_purple2 = '#984ea3'
        color_purple3 = '#7570b3'
        color_green = '#1b9e77'
        color_orange = '#d95f02'
        color_pink = '#e7298a'
        color_lime = '#66a61e'
        color_yellow = '#e6ab02'
        color_brown = '#a6761d'
        color_coal = '#666666'
        color_blue = '#436bad'  # french blue
        color_red = '#ec2d01'   # tomato red

        galaxy_label = r'$\rm Galaxy$'
        absorber_label = r'$\rm Ly\alpha$'
    
        az_0_15_lw = 1.8
        az_15_30_lw = 1.8
        az_30_45_lw = 1.8
        az_45_60_lw = 1.8
        az_60_75_lw = 1.8
        az_75_90_lw = 1.8

        az_0_15_alpha = 0.8
        az_15_30_alpha = 0.8
        az_30_45_alpha = 0.8
        az_45_60_alpha = 0.8
        az_60_75_alpha = 0.8
        az_75_90_alpha = 0.8

        markerSize = 35

        ##########
    
        bins = arange(0, 100, 10)
        color_az_means = color_blue
        symbol_az_means = 'D'
        ls_az_means = '-'
        lw_az_means = 1.8
        alpha_az_means = 0.8
        label_az_means = r'$\rm Mean~Azimuth$'
    
        #
    
        color_az_weighted = color_red
        symbol_az_weighted = 'P'
        ls_az_weighted = '-'
        lw_az_weighted = 1.8
        alpha_az_weighted = 0.8
    #     label_az_weighted = r'$\rm R_{vir}/\rho ~Weighted~Azimuth$'
    #     label_az_weighted = r'$\rm \rho~Weighted~Azimuth$'
        label_az_weighted = r'$\rm \mathcal{L} ~Weighted~Azimuth$'


        #####
    
    #     weights = 1/np.array(all_impact)
    
    #     weights = np.array(all_R_vir)/np.array(all_impact)

        weights = np.array(likelihoods)
    
        ax1.hist(np.array(azimuths),
                bins = 6,
                range = (0, 90),
                weights = weights,
                histtype = 'bar',
                color = color_blue,
                label = label_az_weighted)
                      
        
        ax1.set_ylabel(r'$\rm Number$')
        ax1.set_xlabel(r'$\rm Azimuth~[deg]$')

        # x-axis
        majorLocator   = MultipleLocator(15)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax1.xaxis.set_major_locator(majorLocator)
        ax1.xaxis.set_major_formatter(majorFormatter)
        ax1.xaxis.set_minor_locator(minorLocator)

    #     ax1.grid(b=None, which='major', axis='both', alpha=0.5)

        ax1.legend(scatterpoints=1,prop={'size':12},loc='upper left',fancybox=True)

        title(target_sightline)

    
        ax2 = fig.add_subplot(212)

        ###########################################################################
            
        ax2.hist(np.array(azimuths),
                bins = 18,
                range = (0, 90),
                weights = weights,
                histtype = 'bar',
                color = color_green,
                label = label_az_weighted)
        
        ax2.set_ylabel(r'$\rm Number$')
        ax2.set_xlabel(r'$\rm Azimuth~[deg]$')

        # x-axis
        majorLocator   = MultipleLocator(15)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax2.xaxis.set_major_locator(majorLocator)
        ax2.xaxis.set_major_formatter(majorFormatter)
        ax2.xaxis.set_minor_locator(minorLocator)

        ax2.legend(scatterpoints=1,prop={'size':12},loc='upper left',fancybox=True)
        
        plot_name = '{0}_L_weighted_azimuth_hist_dmin_{1}_minEW_{2}.pdf'.format(target_sightline, d_min_label, min_EW)
        savefig('{0}/{1}'.format(saveDirectory, plot_name),format='pdf',bbox_inches='tight')


##########################################################################################
##########################################################################################
   
    # plot the azimuth within/outside 1 Rvir distributions
    plot_detection_raction_az_rvir_sep = True
    
    if plot_detection_raction_az_rvir_sep:
        # start creating figure
        fig = figure(figsize=(7.7,5.7))
        ax1 = fig.add_subplot(111)

        # colors for plotting
        color_purple = '#7570b3'
        color_purple2 = '#984ea3'
        color_purple3 = '#7570b3'
        color_green = '#1b9e77'
        color_orange = '#d95f02'
        color_pink = '#e7298a'
        color_lime = '#66a61e'
        color_yellow = '#e6ab02'
        color_brown = '#a6761d'
        color_coal = '#666666'
        color_blue = '#436bad'  # french blue
        color_red = '#ec2d01'   # tomato red

        galaxy_label = r'$\rm Galaxy$'
        absorber_label = r'$\rm Ly\alpha$'
        
        #######
        inside_marker = 'D'
        inside_color = color_blue
        inside_ls = '-'
        inside_lw = 1.8
        inside_alpha = 0.8
        inside_label = r'$\rm R_{vir}/\rho \leq 1.5~Detection~Fraction$'
    
        outside_marker = 'P'
        outside_color = color_red
        outside_ls = '-'
        outside_lw = 1.8
        outside_alpha = 0.8
        outside_label = r'$\rm R_{vir}/\rho > 1.5~Detection~Fraction$'
        #######


        ####
        # lists for plotting
        # --- detections
        inside_rvir_W = []
        inside_rvir_az = []
        inside_rvir_W_15 = []
        inside_rvir_az_15 = []
        inside_rvir_W_30 = []
        inside_rvir_az_30 = []
        inside_rvir_W_45 = []
        inside_rvir_az_45 = []
        inside_rvir_W_60 = []
        inside_rvir_az_60 = []
        inside_rvir_W_75 = []
        inside_rvir_az_75 = []
        inside_rvir_W_90 = []
        inside_rvir_az_90 = []
    
    
        outside_rvir_W = []
        outside_rvir_az = []
        outside_rvir_W_15 = []
        outside_rvir_az_15 = []
        outside_rvir_W_30 = []
        outside_rvir_az_30 = []
        outside_rvir_W_45 = []
        outside_rvir_az_45 = []
        outside_rvir_W_60 = []
        outside_rvir_az_60 = []
        outside_rvir_W_75 = []
        outside_rvir_az_75 = []
        outside_rvir_W_90 = []
        outside_rvir_az_90 = []

        # --- all galaxies
        inside_rvir_az_ag = []
        inside_rvir_az_15_ag = []
        inside_rvir_az_30_ag = []
        inside_rvir_az_45_ag = []
        inside_rvir_az_60_ag = []
        inside_rvir_az_75_ag = []
        inside_rvir_az_90_ag = []
    
        outside_rvir_az_ag = []
        outside_rvir_az_15_ag = []
        outside_rvir_az_30_ag = []
        outside_rvir_az_45_ag = []
        outside_rvir_az_60_ag = []
        outside_rvir_az_75_ag = []
        outside_rvir_az_90_ag = []
        
        R_virs = ratio_dict['R_virs']
        impacts = ratio_dict['impacts']
        azimuths = ratio_dict['azimuths']
        Ws = ratio_dict['Ws']
        
        # error function for detection fraction
        def detection_fraction_error(det, non):
            first_term = (non * det**2)/(non + det)**4
        
            second_term = det * (-det/(non + det)**2 + 1/(non + det))**2
            
            err = np.sqrt(first_term + second_term)
            
            return err
        
        
        ####
        # now separate inside vs outside R_vir
                        
        for rvir, impact, w, az in zip(R_virs, impacts, Ws, azimuths):
            if (impact / rvir) <= 2.0:
                inside_rvir_W.append(w)
                inside_rvir_az.append(az)
                
                if az <= 15:
                    inside_rvir_W_15.append(w)
                    inside_rvir_az_15.append(az)
                    
                if az <= 30 and az >15:
                    inside_rvir_W_30.append(w)
                    inside_rvir_az_30.append(az)
                
                if az <= 45 and az >30:
                    inside_rvir_W_45.append(w)
                    inside_rvir_az_45.append(az)
                    
                if az <= 60 and az >45:
                    inside_rvir_W_60.append(w)
                    inside_rvir_az_60.append(az)
                
                if az <= 75 and az >60:
                    inside_rvir_W_75.append(w)
                    inside_rvir_az_75.append(az)
                
                if az <= 90 and az >75:
                    inside_rvir_W_90.append(w)
                    inside_rvir_az_90.append(az)
                
            else:
                outside_rvir_W.append(w)
                outside_rvir_az.append(az)

                if az <= 15:
                    outside_rvir_W_15.append(w)
                    outside_rvir_az_15.append(az)
                    
                if az <= 30 and az >15:
                    outside_rvir_W_30.append(w)
                    outside_rvir_az_30.append(az)
                
                if az <= 45 and az >30:
                    outside_rvir_W_45.append(w)
                    outside_rvir_az_45.append(az)
                    
                if az <= 60 and az >45:
                    outside_rvir_W_60.append(w)
                    outside_rvir_az_60.append(az)
                
                if az <= 75 and az >60:
                    outside_rvir_W_75.append(w)
                    outside_rvir_az_75.append(az)
                
                if az <= 90 and az >75:
                    outside_rvir_W_90.append(w)
                    outside_rvir_az_90.append(az)
        
        #######
        # --- this is just for passing galaxies regardless of detections
        
        print 'R_virs: ',R_virs
        print
        print 'impacts: ',impacts
        print 
        print 'azimuths: ',azimuths
        
        for rvir, impact, az in zip(all_R_virs, all_impacts, all_azimuths):
            if (impact / rvir) <= 2.0:
                inside_rvir_az_ag.append(az)
                
                if az <= 15:
                    inside_rvir_az_15_ag.append(az)
                    
                if az <= 30 and az >15:
                    inside_rvir_az_30_ag.append(az)
                
                if az <= 45 and az >30:
                    inside_rvir_az_45_ag.append(az)
                    
                if az <= 60 and az >45:
                    inside_rvir_az_60_ag.append(az)
                
                if az <= 75 and az >60:
                    inside_rvir_az_75_ag.append(az)
                
                if az <= 90 and az >75:
                    inside_rvir_az_90_ag.append(az)
                
            else:
                outside_rvir_az_ag.append(az)

                if az <= 15:
                    outside_rvir_az_15_ag.append(az)
                    
                if az <= 30 and az >15:
                    outside_rvir_az_30_ag.append(az)
                
                if az <= 45 and az >30:
                    outside_rvir_az_45_ag.append(az)
                    
                if az <= 60 and az >45:
                    outside_rvir_az_60_ag.append(az)
                
                if az <= 75 and az >60:
                    outside_rvir_az_75_ag.append(az)
                
                if az <= 90 and az >75:
                    outside_rvir_az_90_ag.append(az)
    
        #####
        
        # compute detection fractions and mean azimuths for each
        # 0 - 15 deg
        try:
            print 'len(inside_rvir_az_15 vs ag): ',len(inside_rvir_az_15), float(len(inside_rvir_az_15_ag))
            inside_det_frac_15 = len(inside_rvir_az_15) / float(len(inside_rvir_az_15_ag))
            inside_det_frac_15_err = detection_fraction_error(det, non)
        except Exception, e:
            inside_det_frac_15 = 0
            print 'Error = {0}'.format(e)
            print 'len(inside_rvir_az_15) = {0}'.format(len(inside_rvir_az_15))
            print 'len(inside_rvir_az_15_ag) = {0}'.format(len(inside_rvir_az_15_ag))
        inside_mean_EW_15 = mean(inside_rvir_W_15)
        
        try:
            print 'len(outside_rvir_az_15 vs ag): ',len(outside_rvir_az_15), float(len(outside_rvir_az_15_ag))
            outside_det_frac_15 = len(outside_rvir_az_15) / float(len(outside_rvir_az_15_ag))
        except Exception, e:
            outside_det_frac_15 = 0
            print 'Error = {0}'.format(e)
            print 'len(outside_rvir_az_15) = {0}'.format(len(outside_rvir_az_15))
            print 'len(outside_rvir_az_15_ag) = {0}'.format(len(outside_rvir_az_15_ag))
        outside_mean_EW_15 = mean(outside_rvir_W_15)
        
        # 15 - 30 deg
        try:
            print 'len(inside_rvir_az_30 vs ag): ',len(inside_rvir_az_30), float(len(inside_rvir_az_30_ag))
            inside_det_frac_30 = len(inside_rvir_az_30) / float(len(inside_rvir_az_30_ag))
        except Exception, e:
            inside_det_frac_30 = 0
            print 'Error = {0}'.format(e)
            print 'len(inside_rvir_az_30) = {0}'.format(len(inside_rvir_az_30))
            print 'len(inside_rvir_az_30_ag) = {0}'.format(len(inside_rvir_az_30_ag))
        inside_mean_EW_30 = mean(inside_rvir_W_30)
        
        try:
            print 'len(outside_rvir_az_30 vs ag): ',len(outside_rvir_az_30), float(len(outside_rvir_az_30_ag))
            outside_det_frac_30 = len(outside_rvir_az_30) / float(len(outside_rvir_az_30_ag))
        except Exception, e:
            outside_det_frac_30 = 0
            print 'Error = {0}'.format(e)
            print 'len(outside_rvir_az_30) = {0}'.format(len(outside_rvir_az_30))
            print 'len(outside_rvir_az_30_ag) = {0}'.format(len(outside_rvir_az_30_ag))
        outside_mean_EW_30 = mean(outside_rvir_W_30)
        
        # 30 - 45 deg
        try:
            print 'len(inside_rvir_az_45 vs ag): ',len(inside_rvir_az_45), float(len(inside_rvir_az_45_ag))
            inside_det_frac_45 = len(inside_rvir_az_45) / float(len(inside_rvir_az_45_ag))
        except Exception, e:
            inside_det_frac_45 = 0
            print 'Error = {0}'.format(e)
            print 'len(inside_rvir_az_45) = {0}'.format(len(inside_rvir_az_45))
            print 'len(inside_rvir_az_45_ag) = {0}'.format(len(inside_rvir_az_45_ag))
        inside_mean_EW_45 = mean(inside_rvir_W_45)

        try:
            print 'len(outside_rvir_az_45 vs ag): ',len(outside_rvir_az_45), float(len(outside_rvir_az_45_ag))
            outside_det_frac_45 = len(outside_rvir_az_45) / float(len(outside_rvir_az_45_ag))
        except Exception, e:
            outside_det_frac_45 = 0
            print 'Error = {0}'.format(e)
            print 'len(outside_rvir_az_45) = {0}'.format(len(outside_rvir_az_45))
            print 'len(outside_rvir_az_45_ag) = {0}'.format(len(outside_rvir_az_45_ag))
        outside_mean_EW_45 = mean(outside_rvir_W_45)
        
        # 45 - 60 deg
        try:
            print 'len(inside_rvir_az_60 vs ag): ',len(inside_rvir_az_60), float(len(inside_rvir_az_60_ag))
            inside_det_frac_60 = len(inside_rvir_az_60) / float(len(inside_rvir_az_60_ag))
        except Exception, e:
            inside_det_frac_60 = 0
            print 'Error = {0}'.format(e)
            print 'len(inside_rvir_az_60) = {0}'.format(len(inside_rvir_az_60))
            print 'len(inside_rvir_az_60_ag) = {0}'.format(len(inside_rvir_az_60_ag))
        inside_mean_EW_60 = mean(inside_rvir_W_60)
        
        try:
            print 'len(outside_rvir_az_60 vs ag): ',len(outside_rvir_az_60), float(len(outside_rvir_az_60_ag))
            outside_det_frac_60 = len(outside_rvir_az_60) / float(len(outside_rvir_az_60_ag))
        except Exception, e:
            outside_det_frac_60 = 0
            print 'Error = {0}'.format(e)
            print 'len(outside_rvir_az_60) = {0}'.format(len(outside_rvir_az_60))
            print 'len(outside_rvir_az_60_ag) = {0}'.format(len(outside_rvir_az_60_ag))
        outside_mean_EW_60 = mean(outside_rvir_W_60)

        # 60 - 75 deg
        try:
            print 'len(inside_rvir_az_75 vs ag): ',len(inside_rvir_az_75), float(len(inside_rvir_az_75_ag))
            inside_det_frac_75 = len(inside_rvir_az_75) / float(len(inside_rvir_az_75_ag))
        except Exception, e:
            inside_det_frac_75 = 0
            print 'Error = {0}'.format(e)
            print 'len(inside_rvir_az_75) = {0}'.format(len(inside_rvir_az_75))
            print 'len(inside_rvir_az_75_ag) = {0}'.format(len(inside_rvir_az_75_ag))
        inside_mean_EW_75 = mean(inside_rvir_W_75)

        try:
            print 'len(outside_rvir_az_75 vs ag): ',len(outside_rvir_az_75), float(len(outside_rvir_az_75_ag))
            outside_det_frac_75 = len(outside_rvir_az_75) / float(len(outside_rvir_az_75_ag))
        except Exception, e:
            outside_det_frac_75 = 0
            print 'Error = {0}'.format(e)
            print 'len(outside_rvir_az_75) = {0}'.format(len(outside_rvir_az_75))
            print 'len(outside_rvir_az_75_ag) = {0}'.format(len(outside_rvir_az_75_ag))
        outside_mean_EW_75 = mean(outside_rvir_W_75)
        
        # 75 - 90 deg
        try:
            print 'len(inside_rvir_az_90 vs ag): ',len(inside_rvir_az_90), float(len(inside_rvir_az_90_ag))
            inside_det_frac_90 = len(inside_rvir_az_90) / float(len(inside_rvir_az_90_ag))
        except Exception, e:
            inside_det_frac_90 = 0
            print 'Error = {0}'.format(e)
            print 'len(inside_rvir_az_90) = {0}'.format(len(inside_rvir_az_90))
            print 'len(inside_rvir_az_90_ag) = {0}'.format(len(inside_rvir_az_90_ag))
        inside_mean_EW_90 = mean(inside_rvir_W_90)

        try:
            print 'len(outside_rvir_az_90 vs ag): ',len(outside_rvir_az_90), float(len(outside_rvir_az_90_ag))
            outside_det_frac_90 = len(outside_rvir_az_90) / float(len(outside_rvir_az_90_ag))
        except Exception, e:
            outside_det_frac_90 = 0
            print 'Error = {0}'.format(e)
            print 'len(outside_rvir_az_90) = {0}'.format(len(outside_rvir_az_90))
            print 'len(outside_rvir_az_90_ag) = {0}'.format(len(outside_rvir_az_90_ag))
        outside_mean_EW_90 = mean(outside_rvir_W_90)


        ##########
        # --- errors
        reps = 10000
        inside_rvir_az_15_err, inside_rvir_az_15_medianerr = return_bootstrap_errors(inside_rvir_az_15, reps)
        inside_rvir_az_30_err, inside_rvir_az_30_medianerr = return_bootstrap_errors(inside_rvir_az_30, reps)
        inside_rvir_az_45_err, inside_rvir_az_45_medianerr = return_bootstrap_errors(inside_rvir_az_45, reps)
        inside_rvir_az_60_err, inside_rvir_az_60_medianerr = return_bootstrap_errors(inside_rvir_az_60, reps)
        inside_rvir_az_75_err, inside_rvir_az_75_medianerr = return_bootstrap_errors(inside_rvir_az_75, reps)
        inside_rvir_az_90_err, inside_rvir_az_90_medianerr = return_bootstrap_errors(inside_rvir_az_90, reps)

        outside_rvir_az_15_err, outside_rvir_az_15_medianerr = return_bootstrap_errors(outside_rvir_az_15, reps)
        outside_rvir_az_30_err, outside_rvir_az_30_medianerr = return_bootstrap_errors(outside_rvir_az_30, reps)
        outside_rvir_az_45_err, outside_rvir_az_45_medianerr = return_bootstrap_errors(outside_rvir_az_45, reps)
        outside_rvir_az_60_err, outside_rvir_az_60_medianerr = return_bootstrap_errors(outside_rvir_az_60, reps)
        outside_rvir_az_75_err, outside_rvir_az_75_medianerr = return_bootstrap_errors(outside_rvir_az_75, reps)
        outside_rvir_az_90_err, outside_rvir_az_90_medianerr = return_bootstrap_errors(outside_rvir_az_90, reps)

        print "inside_rvir_az_15_err: ",inside_rvir_az_15_err
        print 'inside_rvir_az_30_err: ',inside_rvir_az_30_err
        print 'inside_rvir_az_45_err: ',inside_rvir_az_45_err
        print 'inside_rvir_az_60_err: ',inside_rvir_az_60_err
        print 'inside_rvir_az_75_err: ',inside_rvir_az_75_err
        print 'inside_rvir_az_90_err: ',inside_rvir_az_90_err
        print


        # now make arrays with this values
        inside_det_frac = [inside_det_frac_15, inside_det_frac_30, inside_det_frac_45, 
                            inside_det_frac_60, inside_det_frac_75, inside_det_frac_90]
                            
        outside_det_frac = [outside_det_frac_15, outside_det_frac_30, outside_det_frac_45,
                            outside_det_frac_60, outside_det_frac_75, outside_det_frac_90]
                            
        inside_mean_EW = [inside_mean_EW_15, inside_mean_EW_30, inside_mean_EW_45, 
                        inside_mean_EW_60, inside_mean_EW_75, inside_mean_EW_90]
                        
        outside_mean_EW = [outside_mean_EW_15, outside_mean_EW_30, outside_mean_EW_45,
                            outside_mean_EW_60, outside_mean_EW_75, outside_mean_EW_90]
        
        x_axis = np.array([15/2., 15+15/2., 30+15/2., 45+15/2., 60+15/2., 75+15/2.])

        detection_fraction_l1_outside_err = [detection_fraction_error(det, non) for det, non in zip(y_l1_det_outside, y_l1_non_outside)]
        detection_fraction_l1_inside_err = [detection_fraction_error(det, non) for det, non in zip(y_l1_det_inside, y_l1_non_inside)]
        
        inside_markersize = []
        inside_mean_EW = nan_to_num(inside_mean_EW)
        minSize = 50
        maxSize = 200
        largestEW = max(inside_mean_EW)
        smallestEW = min(inside_mean_EW)
        for w in inside_mean_EW:
            newSize = ((float(w) - smallestEW)/(largestEW - smallestEW)) * (maxSize - minSize) + minSize
            if not isNumber(newSize) or newSize == 0:
                newSize = minSize
            inside_markersize.append(newSize)

        outside_markersize = []
        outside_mean_EW = nan_to_num(outside_mean_EW)
        minSize = 50
        maxSize = 200
        largestEW = max(outside_mean_EW)
        smallestEW = min(outside_mean_EW)
        

        for w in outside_mean_EW:
            newSize = ((float(w) - smallestEW)/(largestEW - smallestEW)) * (maxSize - minSize) + minSize
            outside_markersize.append(newSize)
        
        inside_markersize = nan_to_num(inside_markersize)
        outside_markersize = nan_to_num(outside_markersize)

#         print 'inside_markersize: ',inside_markersize
#         print
#         print 'outside_markersize: ',outside_markersize
#         print
#         print 'inside_det_frac: ',inside_det_frac
#         print
#         print 'outside_det_frac: ',outside_det_frac
        
        
        # inside R_vir
#         plot(x_axis,
#             inside_det_frac,
#             marker = inside_marker,
#             c = inside_color,
#             ms = 1,
#             markeredgecolor = 'black',
#             ls = inside_ls,
#             lw = inside_lw,
#             alpha = inside_alpha,
#             label = inside_label)
            
        inside_errs = [inside_rvir_az_15_err, 
                        inside_rvir_az_30_err,
                        inside_rvir_az_45_err,
                        inside_rvir_az_60_err,
                        inside_rvir_az_75_err,
                        inside_rvir_az_90_err]

        outside_errs = [outside_rvir_az_15_err, 
                        outside_rvir_az_30_err,
                        outside_rvir_az_45_err,
                        outside_rvir_az_60_err,
                        outside_rvir_az_75_err,
                        outside_rvir_az_90_err]
            
        errorbar(x_axis,
            inside_det_frac,
            yerr = inside_errs,
            marker = inside_marker,
            c = inside_color,
            ms = 1,
            markeredgecolor = 'black',
            ls = inside_ls,
            lw = inside_lw,
            alpha = inside_alpha,
            label = inside_label)
            
        scatter(x_axis,
            inside_det_frac,
            marker = inside_marker,
            c = inside_color,
            s = inside_markersize,
            edgecolor = 'black',
            lw = inside_lw,
            alpha = inside_alpha)
        

        # outside R_vir
#         plot(x_axis,
#             outside_det_frac,
#             marker = outside_marker,
#             c = outside_color,
#             ms = 1,
#             markeredgecolor = 'black',
#             ls = outside_ls,
#             lw = outside_lw,
#             alpha = outside_alpha,
#             label = outside_label)
            
        errorbar(x_axis,
            outside_det_frac,
            yerr = outside_errs,
            marker = outside_marker,
            c = outside_color,
            ms = 1,
            markeredgecolor = 'black',
            ls = outside_ls,
            lw = outside_lw,
            alpha = outside_alpha,
            label = outside_label)
            

        scatter(x_axis,
            outside_det_frac,
            marker = outside_marker,
            c = outside_color,
            s = outside_markersize,
            edgecolor = 'black',
            lw = outside_lw,
            alpha = outside_alpha)
        
        ax1.set_ylabel(r'$\rm Detection~Fraction$')
        ax1.set_xlabel(r'$\rm Azimuth ~ [Deg]$')

        # x-axis
        majorLocator   = MultipleLocator(15)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
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
    
    #         leg = ax1.legend(scatterpoints=1,prop={'size':12},loc='lower left',fancybox=True)
    #         leg.get_frame().set_alpha(0.5)

        ax1.grid(b=None, which='major', axis='both', alpha=0.5)
        ax1.legend()
    
        ylim(0., 1)
        xlim(0, 90)
        title(target_sightline)

        plot_name = '{0}_azimuth_DF_ratio2_sep2_dmin_{1}_minEW_{2}_2.pdf'.format(target_sightline, d_min_label, min_EW)
        savefig('{0}/{1}'.format(saveDirectory, plot_name),format='pdf',bbox_inches='tight')


##########################################################################################
##########################################################################################


##########################################################################################
##########################################################################################
##########################################################################################
    
#     detection_fraction.close()
    theFile.close()
    print
    print 'Done!'
    print
        
###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    