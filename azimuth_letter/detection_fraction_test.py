#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)


$Id: detection_fraction5.py v 5.0 11/07/2018

Determine the detection fraction. Correlate the sightlines with galaxies, and figure out 
for what level of likelihood an absorber is always found, etc.


v5: Update for azimuth letter (11/07/18)

v7: Update to include more information (galaxy, target name) (12/28/18)

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
    target_sightline = 'US2816'

    
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

        saveDirectory = '/Users/frenchd/Research/test/sightline_maps/'
        
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
    
    passing_galaxy_dict = {}

    print
    print 'Starting loop!'
    print
    
    total = 1551
    counter = 0
    stopCount = 50000
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
        
        # each 'c' here is the name of a galaxy returned by the correlation
        for c in correlation:
            print 'c: ',c
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
                
                Names.append(Name)
                Vhels.append(Vhel)
                impacts.append(impact)
                R_virs.append(R_vir)
                Lstars.append(Lstar)
                MajDiams.append(MajDiam)
                
#                 passing_galaxy_dict[target] = {
                
                if verbose:
                    print 'Name = {0}'.format(Name)
                    print 'impact = {0} , vhel = {1}'.format(impact, Vhel)
                    print
        
        
        # start creating figure
        fig = figure(figsize=(7.7,5.7))
        ax1 = fig.add_subplot(111)
        
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

        ##########
        # plot the target sightline
        target_x = [0, 10000]
        target_y = [0, 0]
        p1 = plot(target_x,
                target_y,
                lw=2,
                c='black')
                
        ##########
        # plot the absorber locations
        
        # scale marker size to EW
        print 'Lya_Ws: ',Lya_Ws
        print
        Lya_Ws = np.array(Lya_Ws)
        
        absorber_markerSizes = np.interp(Lya_Ws, (Lya_Ws.min(), Lya_Ws.max()), (200, 1000))
        
        print 'absorber_markerSizes: ',absorber_markerSizes
        print 'type(absorber_markerSizes) : ', type(absorber_markerSizes)
        print
                
        # put all the absorbers along the 0 y-axis
        absorber_y = np.zeros(len(Lya_vs))
        
        # plot the absorbers
        p2 = scatter(Lya_vs,
                    absorber_y,
                    marker="|",
                    s = absorber_markerSizes,
                    c = color_red)
#                     label = absorber_label)

        print 'Lya_vs: ',Lya_vs
        print
        print 'absorber_y: ',absorber_y
        print
        
        # now loop through all the galaxies and plot them
        print 'Starting galaxy plotting...'
        print
        for Name, Vhel, impact, R_vir, Lstar, MajDiam in zip(Names, Vhels, impacts, R_virs, Lstars, MajDiams):

            galaxy_markerSize = float(R_vir)
            Vhel = float(Vhel)
            
            print 'galaxy_markerSize: ',galaxy_markerSize
            print 'type(galaxy_markerSize) : ', type(galaxy_markerSize)
            print
    
            # plot impact parameter as a function of velocity for galaxies
            p1 = scatter(Vhel,
                        impact,
                        c = color_blue,
                        s = galaxy_markerSize,
                        marker = 'o')
            
            print 'Vhel: ', Vhel
            print 'impact: ',impact
            print
        
        ax1.set_ylabel(r'$\rm \rho~[kpc]$')
        ax1.set_xlabel(r'$\rm Vhel~[km s^{-1}]$')

        # x-axis
        majorLocator   = MultipleLocator(1000)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(500)
        ax1.xaxis.set_major_locator(majorLocator)
        ax1.xaxis.set_major_formatter(majorFormatter)
        ax1.xaxis.set_minor_locator(minorLocator)

        # y-axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)
        ax1.yaxis.set_major_locator(majorLocator)
        ax1.yaxis.set_major_formatter(majorFormatter)
        ax1.yaxis.set_minor_locator(minorLocator)
        
#         leg = ax1.legend(scatterpoints=1,prop={'size':12},loc='lower left',fancybox=True)
#         leg.get_frame().set_alpha(0.5)

        ax1.grid(b=None,which='major',axis='both')
        ylim(0., 1000)
        xlim(0, 10000)
        title(target)
        
        savefig('{0}/{1}_sightline_map.pdf'.format(saveDirectory, target),format='pdf',bbox_inches='tight')
            

##########################################################################################
##########################################################################################
##########################################################################################
    
    detection_fraction.close()
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
    