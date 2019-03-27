#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)


$Id: detection_fraction_azimuth.py v 1.0 01/23/2019


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
    min_EW = 0

    # minimum diameter. False to bypass, a number > 0 to impose a diameter measurement needed
#     d_min = 0.01 # must have a diameter > 0.01
#     d_min = False # does not need a diameter measurement at all
    d_min = 0.01
    
    
    # sort based on likelihood cus instead of the regular one?
    use_likelihood_cus = False
    
    # double l if impact <= 1 R_vir?
    double_l_within_rvir = False
    
    # print a ton of shit out for testing if True
    verbose = False
    
    # cutoff - will only return this many targers - use for testing together with verbose. 
    cutoff = 1000000.
    
    # target_sightline can be either a single target name or 'all'
#     target_sightline = 'MRK876'
#     target_sightline = 'all'
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
    
    # open the pickle files
#     detection_fraction = open(detection_fraction_filename,'wt')

    
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
    
    passing_galaxy_dict = {}
    
    sub_50_az = []
    gr_50_az = []
    gr_100_az = []
    gr_200_az = []
    gr_300_az = []
    gr_400_az = []
    gr_500_az = []

    sub_50_inc = []
    gr_50_inc = []
    gr_100_inc = []
    gr_200_inc = []
    gr_300_inc = []
    gr_400_inc = []
    gr_500_inc = []    

    sub_50_impact = []
    gr_50_impact = []
    gr_100_impact = []
    gr_200_impact = []
    gr_300_impact = []
    gr_400_impact = []
    gr_500_impact = []
    
    all_impact = []
    all_az = []
    all_inc = []
    all_W = []
    all_R_vir = []
    
    
    print
    print 'Starting loop!'
    print
    
    total = 1551
    counter = 0
    stopCount = 2000
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
            # 'lya' are properly identified lines, 'maybe' are for possible lines, not yet
            # identified, and 'yes' are for identified targets with no Lya in the spectra
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

    
    counter = 0
    # now do the galaxy correlation based on the targets and absorbers collected above
    for target in target_dict:
        counter +=1
        sys.stdout.write("\r Percent Complete: {0} / {1}".format(counter, total))
        sys.stdout.flush()
        print 'on target = ',target
    
        # the collected info for this target
        info = target_dict[target]
        Lya_vs = info['Lya_vs']
        Lya_Ws = info['Lya_Ws']
        Nas = info['Nas']
        bs = info['bs']
        z_target = info['z_target']
        RA_target = info['RA_target']
        Dec_target = info['Dec_target']


        # do the correlation: returns a dictionary where the keys are all the correlated
        # galaxies and the values are a dictionary of all the galaxy info
#         correlation_original = correlateSingle.correlateTarget(target, maxSep, agnSeparation, minVcorr, minSize, slow=False, searchAll=False, returnDict=True)
        
        if counter < cutoff:
            correlation = correlateSingle.correlateTarget(target, maxSep, agnSeparation, minVcorr, minSize, slow=False, searchAll=False, returnDict=True)
        else:
            break

        other_list = []
        
        if verbose:
            print
            print 'target = {0}  ----'.format(target)
            
        # make a list of all galaxies in this correlation, decide which galaxies pass the 
        # given criteria for inclusion
        
        # the list of galaxy names which pass the test
        passing_galaxies = []
        
        Names = []
        Vhels = []
        impacts = []
        R_virs = []
        Lstars = []
        MajDiams = []
        azimuths =  []

        # each 'c' here is the name of a galaxy returned by the correlation
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
            
            group_num = correlation[c]['group_num']
            group_mem = correlation[c]['group_mem']
            group_dist = correlation[c]['group_dist']
            
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
                
                # convert to floats
                Vhel = float(Vhel)
                impact = float(impact)
                R_vir = float(R_vir)
                Lstar = float(Lstar)
                MajDiam = float(MajDiam)
                azimuth = float(azimuth)
                inc = float(inc)
                
                # add needed values to the corresponding lists
                Names.append(Name)
                Vhels.append(Vhel)
                impacts.append(impact)
                R_virs.append(R_vir)
                Lstars.append(Lstar)
                MajDiams.append(MajDiam)
                azimuths.append(azimuth)
                
#                 print
#                 print 'Lya_Ws: ',Lya_Ws
#                 print
#                 print 'Lya_vs: ',Lya_vs
                
                for w, v in zip(Lya_Ws, Lya_vs):
                    if abs(Vhel - v) <= 400:
                    
                        all_az.append(azimuth)
                        all_inc.append(inc)
                        all_impact.append(impact)
                        all_W.append(w)
                        all_R_vir.append(R_vir)

                        if w >= 500:
                            gr_500_az.append(azimuth)
                            gr_500_inc.append(inc)
                            gr_500_impact.append(impact)
                            
                        elif w >= 400:
                            gr_400_az.append(azimuth)
                            gr_400_inc.append(inc)
                            gr_400_impact.append(impact)
                            
                        elif w >= 300:
                            gr_300_az.append(azimuth)
                            gr_300_inc.append(inc)
                            gr_300_impact.append(impact)
                            
                        elif w >= 200:
                            gr_200_az.append(azimuth)
                            gr_200_inc.append(inc)
                            gr_200_impact.append(impact)
                            
                        elif w >= 100:
                            gr_100_az.append(azimuth)
                            gr_100_inc.append(inc)
                            gr_100_impact.append(impact)
                            
                        elif w >= 50:
                            gr_50_az.append(azimuth)
                            gr_50_inc.append(inc)
                            gr_50_impact.append(impact)
                        
                        else:
                            sub_50_az.append(azimuth)
                            sub_50_inc.append(inc)
                            sub_50_impact.append(impact)
                            
#                             print "wtf? w = {0}, v = {1}, target = {2}".format(w, v, target)
                
                
#                 passing_galaxy_dict[target] = {
                
                if verbose:
                    print 'Name = {0}'.format(Name)
                    print 'impact = {0} , vhel = {1}'.format(impact, Vhel)
                    print
        
    ################
    # finished with all loops at this point
    
    # plot the impact parameter histograms first
    
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
    
    lw_sub_50_impact = 1.8
    lw_gr_50_impact = 1.8
    lw_gr_100_impact = 1.8
    lw_gr_200_impact = 1.8
    lw_gr_300_impact = 1.8
    lw_gr_400_impact = 1.8
    lw_gr_500_impact = 1.8

    alpha_sub_50_impact = 0.8
    alpha_gr_50_impact = 0.8
    alpha_gr_100_impact = 0.8
    alpha_gr_200_impact = 0.8
    alpha_gr_300_impact = 0.8
    alpha_gr_400_impact = 0.8
    alpha_gr_500_impact = 0.8

    ##########
    
    bins = arange(0, 1000, 10)

    #####
    # sub_50_impact
    hist(sub_50_impact,
        bins=bins,
        histtype='step',
        cumulative=True,
        normed=1,
        color=color_purple3,
        lw=lw_sub_50_impact,
        alpha=alpha_sub_50_impact,
        label=r'$\rm sub\_50\_impact$')
    
    # gr_50_impact
    hist(gr_50_impact,
        bins=bins,
        histtype='step',
        cumulative=True,
        normed=1,
        color='black',
        lw=lw_gr_50_impact,
        alpha=alpha_gr_50_impact,
        label=r'$\rm gr\_50\_impact$')

    # gr_100_impact
    hist(gr_100_impact,
        bins=bins,
        histtype='step',
        cumulative=True,
        normed=1,
        color=color_green,
        lw=lw_gr_100_impact,
        alpha=alpha_gr_100_impact,
        label=r'$\rm gr\_100\_impact$')
    
    # gr_200_impact
    hist(gr_200_impact,
        bins=bins,
        histtype='step',
        cumulative=True,
        normed=1,
        color=color_orange,
        lw=lw_gr_200_impact,
        alpha=alpha_gr_200_impact,
        label=r'$\rm gr\_200\_impact$')
        
    # gr_300_impact
    hist(gr_300_impact,
        bins=bins,
        histtype='step',
        cumulative=True,
        normed=1,
        color=color_blue,
        lw=lw_gr_300_impact,
        alpha=alpha_gr_300_impact,
        label=r'$\rm gr\_300\_impact$')

    # gr_500_impact
    hist(gr_500_impact,
        bins=bins,
        histtype='step',
        cumulative=True,
        normed=1,
        color=color_red,
        lw=lw_gr_500_impact,
        alpha=alpha_gr_500_impact,
        label=r'$\rm gr\_500\_impact$')
        
    ax1.set_ylabel(r'$\rm Number$')
    ax1.set_xlabel(r'$\rm \rho ~ [kpc]$')

    # x-axis
    majorLocator   = MultipleLocator(100)
    majorFormatter = FormatStrFormatter(r'$\rm %d$')
    minorLocator   = MultipleLocator(50)
    ax1.xaxis.set_major_locator(majorLocator)
    ax1.xaxis.set_major_formatter(majorFormatter)
    ax1.xaxis.set_minor_locator(minorLocator)

    y_scale = int(len(gr_50_impact)/10.)
    
    if y_scale/2 <1:
        y_scale = 2
    
    print 'y_scale : ',y_scale

    # y-axis
#     majorLocator   = MultipleLocator(0.1)
#     majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
#     minorLocator   = MultipleLocator(0.05)
#     ax1.yaxis.set_major_locator(majorLocator)
#     ax1.yaxis.set_major_formatter(majorFormatter)
#     ax1.yaxis.set_minor_locator(minorLocator)
    
#         leg = ax1.legend(scatterpoints=1,prop={'size':12},loc='lower left',fancybox=True)
#         leg.get_frame().set_alpha(0.5)

    ax1.grid(b=None, which='major', axis='both', alpha=0.5)
    ax1.legend()
    
#     ylim(0., 1000)
#     xlim(0, 10000)
    title(target)
    
    savefig('{0}/{1}_impact_hists_CDF.pdf'.format(saveDirectory, target_sightline),format='pdf',bbox_inches='tight')

##########################################################################################
##########################################################################################

    # plot the azimuth distributions
    
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
    
    lw_sub_50_az = 1.8
    lw_gr_50_az = 1.8
    lw_gr_100_az = 1.8
    lw_gr_200_az = 1.8
    lw_gr_300_az = 1.8
    lw_gr_400_az = 1.8
    lw_gr_500_az = 1.8

    alpha_sub_50_az = 0.8
    alpha_gr_50_az = 0.8
    alpha_gr_100_az = 0.8
    alpha_gr_200_az = 0.8
    alpha_gr_300_az = 0.8
    alpha_gr_400_az = 0.8
    alpha_gr_500_az = 0.8

    ##########
    
    bins = arange(0, 100, 10)

    #####
    # sub_50_az
    hist(sub_50_az,
        bins=bins,
        histtype='step',
        normed=1,
        color=color_purple3,
        lw=lw_sub_50_az,
        alpha=alpha_sub_50_az,
        label=r'$\rm sub\_50\_az$')
    
    # gr_50_az
    hist(gr_50_az,
        bins=bins,
        histtype='step',
        normed=1,
        color='black',
        lw=lw_gr_50_az,
        alpha=alpha_gr_50_az,
        label=r'$\rm gr\_50\_az$')

    # gr_100_az
    hist(gr_100_az,
        bins=bins,
        histtype='step',
        normed=1,
        color=color_green,
        lw=lw_gr_100_az,
        alpha=alpha_gr_100_az,
        label=r'$\rm gr\_100\_az$')
    
    # gr_200_az
    hist(gr_200_az,
        bins=bins,
        histtype='step',
        normed=1,
        color=color_orange,
        lw=lw_gr_200_az,
        alpha=alpha_gr_200_az,
        label=r'$\rm gr\_200\_az$')
        
    # gr_300_az
    hist(gr_300_az,
        bins=bins,
        histtype='step',
        normed=1,
        color=color_blue,
        lw=lw_gr_300_az,
        alpha=alpha_gr_300_az,
        label=r'$\rm gr\_300\_az$')

    # gr_500_az
    hist(gr_500_az,
        bins=bins,
        histtype='step',
        normed=1,
        color=color_red,
        lw=lw_gr_500_az,
        alpha=alpha_gr_500_az,
        label=r'$\rm gr\_500\_az$')
        
    ax1.set_ylabel(r'$\rm Number$')
    ax1.set_xlabel(r'$\rm \rho ~ [kpc]$')

    # x-axis
    majorLocator   = MultipleLocator(10)
    majorFormatter = FormatStrFormatter(r'$\rm %d$')
    minorLocator   = MultipleLocator(5)
    ax1.xaxis.set_major_locator(majorLocator)
    ax1.xaxis.set_major_formatter(majorFormatter)
    ax1.xaxis.set_minor_locator(minorLocator)

    y_scale = int(len(gr_50_az)/10.)
    
    if y_scale/2 <1:
        y_scale = 2
    
    print 'y_scale : ',y_scale

    # y-axis
#     majorLocator   = MultipleLocator(0.1)
#     majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
#     minorLocator   = MultipleLocator(0.05)
#     ax1.yaxis.set_major_locator(majorLocator)
#     ax1.yaxis.set_major_formatter(majorFormatter)
#     ax1.yaxis.set_minor_locator(minorLocator)
    
#         leg = ax1.legend(scatterpoints=1,prop={'size':12},loc='lower left',fancybox=True)
#         leg.get_frame().set_alpha(0.5)

    ax1.grid(b=None, which='major', axis='both', alpha=0.5)
    ax1.legend()
    
#     ylim(0., 1000)
#     xlim(0, 10000)
    title(target)
    
    savefig('{0}/{1}_azimuth_hists.pdf'.format(saveDirectory, target_sightline),format='pdf',bbox_inches='tight')

##########################################################################################
##########################################################################################
   
    # plot the azimuth within/outside 1 Rvir distributions
    
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
    
    lw_inside_rvir = 1.8
    lw_outside_rvir = 1.8

    alpha_inside_rvir = 0.8
    alpha_outside_rvir = 0.8

    ####
    # lists for plotting
    inside_rvir_W = []
    inside_rvir_az = []
    
    outside_rvir_W = []
    outside_rvir_az = []
    
    ####
    # now separate inside vs outside R_vir
                        
    for rvir, impact, w, az in zip(all_R_vir, all_impact, all_W, all_az):
#         print 'rvir: ',rvir
#         print 'impact: ',impact
#         print '(rvir / impact)  : ',(rvir / impact)
#         print
        if (impact / rvir) <= 1.5:
            inside_rvir_W.append(w)
            inside_rvir_az.append(az)
        else:
            outside_rvir_W.append(w)
            outside_rvir_az.append(az)
            
    ##########
    
    
    
    bins = arange(0, 100, 10)

    #####
    # inside R_vir
    hist(inside_rvir_az,
        bins=bins,
        histtype='step',
        normed=1,
        color=color_blue,
        lw=lw_inside_rvir,
        alpha=alpha_inside_rvir,
        label=r'$\rm Az ~ (R_{{vir}} / \rho) \leq 1$')
    
    # outside R_vir
    hist(outside_rvir_az,
        bins=bins,
        histtype='step',
        normed=1,
        color=color_orange,
        lw=lw_outside_rvir,
        alpha=alpha_outside_rvir,
        label=r'$\rm Az ~ (R_{{vir}} / \rho) > 1$')
        
    ax1.set_ylabel(r'$\rm Number$')
    ax1.set_xlabel(r'$\rm Azimuth ~ [Deg]$')

    # x-axis
    majorLocator   = MultipleLocator(10)
    majorFormatter = FormatStrFormatter(r'$\rm %d$')
    minorLocator   = MultipleLocator(5)
    ax1.xaxis.set_major_locator(majorLocator)
    ax1.xaxis.set_major_formatter(majorFormatter)
    ax1.xaxis.set_minor_locator(minorLocator)

    y_scale = int(len(gr_50_az)/10.)
    
    if y_scale/2 <1:
        y_scale = 2
    
    print 'y_scale : ',y_scale

    # y-axis
#     majorLocator   = MultipleLocator(0.1)
#     majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
#     minorLocator   = MultipleLocator(0.05)
#     ax1.yaxis.set_major_locator(majorLocator)
#     ax1.yaxis.set_major_formatter(majorFormatter)
#     ax1.yaxis.set_minor_locator(minorLocator)
    
#         leg = ax1.legend(scatterpoints=1,prop={'size':12},loc='lower left',fancybox=True)
#         leg.get_frame().set_alpha(0.5)

    ax1.grid(b=None, which='major', axis='both', alpha=0.5)
    ax1.legend()
    
#     ylim(0., 1000)
#     xlim(0, 10000)
    title(target)
    
    savefig('{0}/{1}_az_rvir_hists_1r5.pdf'.format(saveDirectory, target_sightline),format='pdf',bbox_inches='tight')

##########################################################################################
##########################################################################################

    # plot the inclination distributions
    
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
    
    lw_sub_50_inc = 1.8
    lw_gr_50_inc = 1.8
    lw_gr_100_inc = 1.8
    lw_gr_200_inc = 1.8
    lw_gr_300_inc = 1.8
    lw_gr_400_inc = 1.8
    lw_gr_500_inc = 1.8

    alpha_sub_50_inc = 0.8
    alpha_gr_50_inc = 0.8
    alpha_gr_100_inc = 0.8
    alpha_gr_200_inc = 0.8
    alpha_gr_300_inc = 0.8
    alpha_gr_400_inc = 0.8
    alpha_gr_500_inc = 0.8

    ##########
    
    bins = arange(0, 100, 10)

    #####
    # sub_50_inc
    hist(sub_50_inc,
        bins=bins,
        histtype='step',
        normed=1,
        color=color_purple3,
        lw=lw_sub_50_inc,
        alpha=alpha_sub_50_inc,
        label=r'$\rm sub\_50\_inc$')
    
    # gr_50_inc
    hist(gr_50_inc,
        bins=bins,
        histtype='step',
        normed=1,
        color='black',
        lw=lw_gr_50_inc,
        alpha=alpha_gr_50_inc,
        label=r'$\rm gr\_50\_inc$')

    # gr_100_inc
    hist(gr_100_inc,
        bins=bins,
        histtype='step',
        normed=1,
        color=color_green,
        lw=lw_gr_100_inc,
        alpha=alpha_gr_100_inc,
        label=r'$\rm gr\_100\_inc$')
    
    # gr_200_inc
    hist(gr_200_inc,
        bins=bins,
        histtype='step',
        normed=1,
        color=color_orange,
        lw=lw_gr_200_inc,
        alpha=alpha_gr_200_inc,
        label=r'$\rm gr\_200\_inc$')
        
    # gr_300_inc
    hist(gr_300_inc,
        bins=bins,
        histtype='step',
        normed=1,
        color=color_blue,
        lw=lw_gr_300_inc,
        alpha=alpha_gr_300_inc,
        label=r'$\rm gr\_300\_inc$')

    # gr_500_inc
    hist(gr_500_inc,
        bins=bins,
        histtype='step',
        normed=1,
        color=color_red,
        lw=lw_gr_500_inc,
        alpha=alpha_gr_500_inc,
        label=r'$\rm gr\_500\_inc$')
        
    ax1.set_ylabel(r'$\rm Number$')
    ax1.set_xlabel(r'$\rm \rho ~ [kpc]$')

    # x-axis
    majorLocator   = MultipleLocator(10)
    majorFormatter = FormatStrFormatter(r'$\rm %d$')
    minorLocator   = MultipleLocator(5)
    ax1.xaxis.set_major_locator(majorLocator)
    ax1.xaxis.set_major_formatter(majorFormatter)
    ax1.xaxis.set_minor_locator(minorLocator)

    y_scale = int(len(gr_50_inc)/10.)
    
    if y_scale/2 <1:
        y_scale = 2
    
    print 'y_scale : ',y_scale

    # y-axis
#     majorLocator   = MultipleLocator(0.1)
#     majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
#     minorLocator   = MultipleLocator(0.05)
#     ax1.yaxis.set_major_locator(majorLocator)
#     ax1.yaxis.set_major_formatter(majorFormatter)
#     ax1.yaxis.set_minor_locator(minorLocator)
    
#         leg = ax1.legend(scatterpoints=1,prop={'size':12},loc='lower left',fancybox=True)
#         leg.get_frame().set_alpha(0.5)

    ax1.grid(b=None, which='major', axis='both', alpha=0.5)
    ax1.legend()
    
#     ylim(0., 1000)
#     xlim(0, 10000)
    title(target)
    
    savefig('{0}/{1}_inc_hists.pdf'.format(saveDirectory, target_sightline),format='pdf',bbox_inches='tight')



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
    