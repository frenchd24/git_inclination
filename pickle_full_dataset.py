#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)


$Id: pickle_full_dataset.py v 2.0 06/01/2018

Create a pickle file with already correlated results


Output bins:
1. Isolated (no galaxies within 500 kpc, 400 km/s)
2. No galaxies with L >= 0.01
3. 1 galaxy with L >= 0.01, isolated (no galaxies within 500 kpc, 400 km/s)
4. 1 galaxy with L >= 0.01, no neighbors within a factor of 5L
5. 2 galaxies with L >=0.01 within a factor of 5L, no others
6. More than 2 galaxies with L >= 0.01 within a factor of 5L
7. At least 1 galaxy with L >= 0.1 where that galaxy is a designated group member (Based on Tully 2015 2MASS group catalog)


Comes from:
picklePilot.py v 1.0 01/03/2018

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

from matplotlib import rc


###########################################################################

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
    
    
def add_to_isolated(Lya_v, e_Lya_v, Lya_W, e_Lya_W, Na, e_Na, b, e_b, W, e_W, target, z_target, RA_target, Dec_target):
    Lya_vs = isolated['Lya_vs']
    Lya_vs.append(float(Lya_v))

    e_Lya_vs = isolated['e_Lya_vs']
    e_Lya_vs.append(float(e_Lya_v))

    Lya_Ws = isolated['Lya_Ws']
    Lya_Ws.append(float(Lya_W))
    
    e_Lya_Ws = isolated['e_Lya_Ws']
    e_Lya_Ws.append(float(e_Lya_s))

    Nas = isolated['Nas']
    Nas.append(float(Na))
    
    e_Nas = isolated['e_Nas']
    e_Nas.append(float(e_Na))
    
    bs = isolated['bs']
    bs.append(float(b))

    e_bs = isolated['e_bs']
    e_bs.append(float(e_b))

    Ws = isolated['Ws']
    Ws.append(float(Ws))
    
    e_Ws = isolated['e_Ws']
    e_Ws.append(float(e_W))

    targets = isolated['targets']
    targets.append(target)

    z_targets = isolated['z_targets']
    z_targets.append(float(z_target))

    RA_targets = isolated['RA_targets']
    RA_targets.append(float(RA_target))
    
    Dec_targets = isolated['Dec_targets']
    Dec_targets.append(float(Dec_target))
    
    pass
    


def add_to_L_isolated(Lya_v, e_Lya_v, Lya_W, e_Lya_W, Na, e_Na, b, e_b, W, e_W, target, z_target, RA_target, Dec_target):
    Lya_vs = L_isolated['Lya_vs']
    Lya_vs.append(float(Lya_v))

    e_Lya_vs = L_isolated['e_Lya_vs']
    e_Lya_vs.append(float(e_Lya_v))

    Lya_Ws = L_isolated['Lya_Ws']
    Lya_Ws.append(float(Lya_W))
    
    e_Lya_Ws = L_isolated['e_Lya_Ws']
    e_Lya_Ws.append(float(e_Lya_s))

    Nas = L_isolated['Nas']
    Nas.append(float(Na))
    
    e_Nas = L_isolated['e_Nas']
    e_Nas.append(float(e_Na))
    
    bs = L_isolated['bs']
    bs.append(float(b))

    e_bs = L_isolated['e_bs']
    e_bs.append(float(e_b))

    Ws = L_isolated['Ws']
    Ws.append(float(Ws))
    
    e_Ws = L_isolated['e_Ws']
    e_Ws.append(float(e_W))

    targets = L_isolated['targets']
    targets.append(target)

    z_targets = L_isolated['z_targets']
    z_targets.append(float(z_target))

    RA_targets = L_isolated['RA_targets']
    RA_targets.append(float(RA_target))
    
    Dec_targets = L_isolated['Dec_targets']
    Dec_targets.append(float(Dec_target))
    
    pass
    

    
    
def add_to_L_associated_isolated(Lya_v, e_Lya_v, Lya_W, e_Lya_W, Na, e_Na, b, e_b, W, e_W, target, z_target,\
    RA_target, Dec_target, galaxyName, RA_galaxy, Dec_galaxy, impact, azimuth, pa, inclination,\
    adjustedInc, l, l_cus, Rvir, cus, diam, morph, vhel, vcorr, bestDist, e_bestDist,\
    group_num, group_mem, group_dist, Lstar, e_Lstar, B_mag):


    Lya_vs = L_associated_isolated['Lya_vs']
    Lya_vs.append(float(Lya_v))

    e_Lya_vs = L_associated_isolated['e_Lya_vs']
    e_Lya_vs.append(float(e_Lya_v))

    Lya_Ws = L_associated_isolated['Lya_Ws']
    Lya_Ws.append(float(Lya_W))
    
    e_Lya_Ws = L_associated_isolated['e_Lya_Ws']
    e_Lya_Ws.append(float(e_Lya_s))

    Nas = L_associated_isolated['Nas']
    Nas.append(float(Na))
    
    e_Nas = L_associated_isolated['e_Nas']
    e_Nas.append(float(e_Na))
    
    bs = L_associated_isolated['bs']
    bs.append(float(b))

    e_bs = L_associated_isolated['e_bs']
    e_bs.append(float(e_b))

    Ws = L_associated_isolated['Ws']
    Ws.append(float(Ws))
    
    e_Ws = L_associated_isolated['e_Ws']
    e_Ws.append(float(e_W))

    targets = L_associated_isolated['targets']
    targets.append(target)

    z_targets = L_associated_isolated['z_targets']
    z_targets.append(float(z_target))

    RA_targets = L_associated_isolated['RA_targets']
    RA_targets.append(float(RA_target))
    
    Dec_targets = L_associated_isolated['Dec_targets']
    Dec_targets.append(float(Dec_target))
    
    galaxyNames = L_associated_isolated['galaxyNames']
    galaxyNames.append(galaxyName)

    RA_galaxies = L_associated_isolated['RA_galaxies']
    RA_galaxies.append(float(RA_galaxy))
    
    Dec_galaxies = L_associated_isolated['Dec_galaxies']
    Dec_galaxies.append(float(Dec_galaxy))
    
    impacts = L_associated_isolated['impacts']
    impacts.append(float(impact))
    
    azimuths = L_associated_isolated['azimuths']
    azimuths.append(float(azimuth))
    
    pas = L_associated_isolated['pas']
    pas.append(float(pa))
    
    inclinations = L_associated_isolated['inclinations']
    inclinations.append(float(inclination))
    
    adjustedIncs = L_associated_isolated['adjustedIncs']
    adjustedIncs.append(float(adjustedInc))
    
    ls = L_associated_isolated['ls']
    ls.append(float(l))

    l_cuss = L_associated_isolated['l_cuss']
    l_cuss.append(float(l_cus))

    Rvirs = L_associated_isolated['Rvirs']
    Rvirs.append(float(Rvir))
    
    cuss = L_associated_isolated['cuss']
    cuss.append(float(cus))

    diams = L_associated_isolated['diams']
    diams.append(float(diam))
    
    morphs = L_associated_isolated['morphs']
    morphs.append(morph)
    
    vhels = L_associated_isolated['vhels']
    vhels.append(float(vhel))
    
    vcorrs = L_associated_isolated['vcorrs']
    vcorrs.append(float(vcorr))
    
    bestDists = L_associated_isolated['bestDists']
    bestDists.append(float(bestDist))
    
    e_bestDists = L_associated_isolated['e_bestDists']
    e_bestDists.append(float(e_bestDist))
    
    group_nums = L_associated_isolated['group_nums']
    group_nums.append(float(group_num))
    
    group_mems = L_associated_isolated['group_mems']
    group_mems.append(float(group_mem))

    group_dists = L_associated_isolated['group_dists']
    group_dists.append(float(group_dist))
    
    Lstars = L_associated_isolated['Lstars']
    Lstars.append(float(Lstar))
    
    e_Lstars = L_associated_isolated['e_Lstars']
    e_Lstars.append(float(e_Lstar))
    
    B_mags = L_associated_isolated['B_mags']
    B_mags.append(float(B_mag))
    
    pass
    



def add_to_L_associated(Lya_v, e_Lya_v, Lya_W, e_Lya_W, Na, e_Na, b, e_b, W, e_W, target, z_target,\
    RA_target, Dec_target, galaxyName, RA_galaxy, Dec_galaxy, impact, azimuth, pa, inclination,\
    adjustedInc, l, l_cus, Rvir, cus, diam, morph, vhel, vcorr, bestDist, e_bestDist,\
    group_num, group_mem, group_dist, Lstar, e_Lstar, B_mag):


    Lya_vs = L_associated['Lya_vs']
    Lya_vs.append(float(Lya_v))

    e_Lya_vs = L_associated['e_Lya_vs']
    e_Lya_vs.append(float(e_Lya_v))

    Lya_Ws = L_associated['Lya_Ws']
    Lya_Ws.append(float(Lya_W))
    
    e_Lya_Ws = L_associated['e_Lya_Ws']
    e_Lya_Ws.append(float(e_Lya_s))

    Nas = L_associated['Nas']
    Nas.append(float(Na))
    
    e_Nas = L_associated['e_Nas']
    e_Nas.append(float(e_Na))
    
    bs = L_associated['bs']
    bs.append(float(b))

    e_bs = L_associated['e_bs']
    e_bs.append(float(e_b))

    Ws = L_associated['Ws']
    Ws.append(float(Ws))
    
    e_Ws = L_associated['e_Ws']
    e_Ws.append(float(e_W))

    targets = L_associated['targets']
    targets.append(target)

    z_targets = L_associated['z_targets']
    z_targets.append(float(z_target))

    RA_targets = L_associated['RA_targets']
    RA_targets.append(float(RA_target))
    
    Dec_targets = L_associated['Dec_targets']
    Dec_targets.append(float(Dec_target))
    
    galaxyNames = L_associated['galaxyNames']
    galaxyNames.append(galaxyName)

    RA_galaxies = L_associated['RA_galaxies']
    RA_galaxies.append(float(RA_galaxy))
    
    Dec_galaxies = L_associated['Dec_galaxies']
    Dec_galaxies.append(float(Dec_galaxy))
    
    impacts = L_associated['impacts']
    impacts.append(float(impact))
    
    azimuths = L_associated['azimuths']
    azimuths.append(float(azimuth))
    
    pas = L_associated['pas']
    pas.append(float(pa))
    
    inclinations = L_associated['inclinations']
    inclinations.append(float(inclination))
    
    adjustedIncs = L_associated['adjustedIncs']
    adjustedIncs.append(float(adjustedInc))
    
    ls = L_associated['ls']
    ls.append(float(l))

    l_cuss = L_associated['l_cuss']
    l_cuss.append(float(l_cus))

    Rvirs = L_associated['Rvirs']
    Rvirs.append(float(Rvir))
    
    cuss = L_associated['cuss']
    cuss.append(float(cus))

    diams = L_associated['diams']
    diams.append(float(diam))
    
    morphs = L_associated['morphs']
    morphs.append(morph)
    
    vhels = L_associated['vhels']
    vhels.append(float(vhel))
    
    vcorrs = L_associated['vcorrs']
    vcorrs.append(float(vcorr))
    
    bestDists = L_associated['bestDists']
    bestDists.append(float(bestDist))
    
    e_bestDists = L_associated['e_bestDists']
    e_bestDists.append(float(e_bestDist))
    
    group_nums = L_associated['group_nums']
    group_nums.append(float(group_num))
    
    group_mems = L_associated['group_mems']
    group_mems.append(float(group_mem))

    group_dists = L_associated['group_dists']
    group_dists.append(float(group_dist))
    
    Lstars = L_associated['Lstars']
    Lstars.append(float(Lstar))
    
    e_Lstars = L_associated['e_Lstars']
    e_Lstars.append(float(e_Lstar))
    
    B_mags = L_associated['B_mags']
    B_mags.append(float(B_mag))
    
    pass




def add_to_L_nonassociated(Lya_v, e_Lya_v, Lya_W, e_Lya_W, Na, e_Na, b, e_b, W, e_W, target, z_target,\
    RA_target, Dec_target, galaxyName, RA_galaxy, Dec_galaxy, impact, azimuth, pa, inclination,\
    adjustedInc, l, l_cus, Rvir, cus, diam, morph, vhel, vcorr, bestDist, e_bestDist,\
    group_num, group_mem, group_dist, Lstar, e_Lstar, B_mag):


    Lya_vs = L_nonassociated['Lya_vs']
    Lya_vs.append(float(Lya_v))

    e_Lya_vs = L_nonassociated['e_Lya_vs']
    e_Lya_vs.append(float(e_Lya_v))

    Lya_Ws = L_nonassociated['Lya_Ws']
    Lya_Ws.append(float(Lya_W))
    
    e_Lya_Ws = L_nonassociated['e_Lya_Ws']
    e_Lya_Ws.append(float(e_Lya_s))

    Nas = L_nonassociated['Nas']
    Nas.append(float(Na))
    
    e_Nas = L_nonassociated['e_Nas']
    e_Nas.append(float(e_Na))
    
    bs = L_nonassociated['bs']
    bs.append(float(b))

    e_bs = L_nonassociated['e_bs']
    e_bs.append(float(e_b))

    Ws = L_nonassociated['Ws']
    Ws.append(float(Ws))
    
    e_Ws = L_nonassociated['e_Ws']
    e_Ws.append(float(e_W))

    targets = L_nonassociated['targets']
    targets.append(target)

    z_targets = L_nonassociated['z_targets']
    z_targets.append(float(z_target))

    RA_targets = L_nonassociated['RA_targets']
    RA_targets.append(float(RA_target))
    
    Dec_targets = L_nonassociated['Dec_targets']
    Dec_targets.append(float(Dec_target))
    
    galaxyNames = L_nonassociated['galaxyNames']
    galaxyNames.append(galaxyName)

    RA_galaxies = L_nonassociated['RA_galaxies']
    RA_galaxies.append(float(RA_galaxy))
    
    Dec_galaxies = L_nonassociated['Dec_galaxies']
    Dec_galaxies.append(float(Dec_galaxy))
    
    impacts = L_nonassociated['impacts']
    impacts.append(float(impact))
    
    azimuths = L_nonassociated['azimuths']
    azimuths.append(float(azimuth))
    
    pas = L_nonassociated['pas']
    pas.append(float(pa))
    
    inclinations = L_nonassociated['inclinations']
    inclinations.append(float(inclination))
    
    adjustedIncs = L_nonassociated['adjustedIncs']
    adjustedIncs.append(float(adjustedInc))
    
    ls = L_nonassociated['ls']
    ls.append(float(l))

    l_cuss = L_nonassociated['l_cuss']
    l_cuss.append(float(l_cus))

    Rvirs = L_nonassociated['Rvirs']
    Rvirs.append(float(Rvir))
    
    cuss = L_nonassociated['cuss']
    cuss.append(float(cus))

    diams = L_nonassociated['diams']
    diams.append(float(diam))
    
    morphs = L_nonassociated['morphs']
    morphs.append(morph)
    
    vhels = L_nonassociated['vhels']
    vhels.append(float(vhel))
    
    vcorrs = L_nonassociated['vcorrs']
    vcorrs.append(float(vcorr))
    
    bestDists = L_nonassociated['bestDists']
    bestDists.append(float(bestDist))
    
    e_bestDists = L_nonassociated['e_bestDists']
    e_bestDists.append(float(e_bestDist))
    
    group_nums = L_nonassociated['group_nums']
    group_nums.append(float(group_num))
    
    group_mems = L_nonassociated['group_mems']
    group_mems.append(float(group_mem))

    group_dists = L_nonassociated['group_dists']
    group_dists.append(float(group_dist))
    
    Lstars = L_nonassociated['Lstars']
    Lstars.append(float(Lstar))
    
    e_Lstars = L_nonassociated['e_Lstars']
    e_Lstars.append(float(e_Lstar))
    
    B_mags = L_nonassociated['B_mags']
    B_mags.append(float(B_mag))
    
    pass





def add_to_L_two(Lya_v, e_Lya_v, Lya_W, e_Lya_W, Na, e_Na, b, e_b, W, e_W, target, z_target,\
    RA_target, Dec_target, galaxyName, RA_galaxy, Dec_galaxy, impact, azimuth, pa, inclination,\
    adjustedInc, l, l_cus, Rvir, cus, diam, morph, vhel, vcorr, bestDist, e_bestDist,\
    group_num, group_mem, group_dist, Lstar, e_Lstar, B_mag):

 
    Lya_vs = L_two['Lya_vs']
    Lya_vs.append(float(Lya_v))

    e_Lya_vs = L_two['e_Lya_vs']
    e_Lya_vs.append(float(e_Lya_v))

    Lya_Ws = L_two['Lya_Ws']
    Lya_Ws.append(float(Lya_W))
    
    e_Lya_Ws = L_two['e_Lya_Ws']
    e_Lya_Ws.append(float(e_Lya_s))

    Nas = L_two['Nas']
    Nas.append(float(Na))
    
    e_Nas = L_two['e_Nas']
    e_Nas.append(float(e_Na))
    
    bs = L_two['bs']
    bs.append(float(b))

    e_bs = L_two['e_bs']
    e_bs.append(float(e_b))

    Ws = L_two['Ws']
    Ws.append(float(Ws))
    
    e_Ws = L_two['e_Ws']
    e_Ws.append(float(e_W))

    targets = L_two['targets']
    targets.append(target)

    z_targets = L_two['z_targets']
    z_targets.append(float(z_target))

    RA_targets = L_two['RA_targets']
    RA_targets.append(float(RA_target))
    
    Dec_targets = L_two['Dec_targets']
    Dec_targets.append(float(Dec_target))
    
    galaxyNames = L_two['galaxyNames']
    galaxyNames.append(galaxyName)

    RA_galaxies = L_two['RA_galaxies']
    RA_galaxies.append(float(RA_galaxy))
    
    Dec_galaxies = L_two['Dec_galaxies']
    Dec_galaxies.append(float(Dec_galaxy))
    
    impacts = L_two['impacts']
    impacts.append(float(impact))
    
    azimuths = L_two['azimuths']
    azimuths.append(float(azimuth))
    
    pas = L_two['pas']
    pas.append(float(pa))
    
    inclinations = L_two['inclinations']
    inclinations.append(float(inclination))
    
    adjustedIncs = L_two['adjustedIncs']
    adjustedIncs.append(float(adjustedInc))
    
    ls = L_two['ls']
    ls.append(float(l))

    l_cuss = L_two['l_cuss']
    l_cuss.append(float(l_cus))

    Rvirs = L_two['Rvirs']
    Rvirs.append(float(Rvir))
    
    cuss = L_two['cuss']
    cuss.append(float(cus))

    diams = L_two['diams']
    diams.append(float(diam))
    
    morphs = L_two['morphs']
    morphs.append(morph)
    
    vhels = L_two['vhels']
    vhels.append(float(vhel))
    
    vcorrs = L_two['vcorrs']
    vcorrs.append(float(vcorr))
    
    bestDists = L_two['bestDists']
    bestDists.append(float(bestDist))
    
    e_bestDists = L_two['e_bestDists']
    e_bestDists.append(float(e_bestDist))
    
    group_nums = L_two['group_nums']
    group_nums.append(float(group_num))
    
    group_mems = L_two['group_mems']
    group_mems.append(float(group_mem))

    group_dists = L_two['group_dists']
    group_dists.append(float(group_dist))
    
    Lstars = L_two['Lstars']
    Lstars.append(float(Lstar))
    
    e_Lstars = L_two['e_Lstars']
    e_Lstars.append(float(e_Lstar))
    
    B_mags = L_two['B_mags']
    B_mags.append(float(B_mag))
    
    pass

    
    

def add_to_L_two_plus(Lya_v, e_Lya_v, Lya_W, e_Lya_W, Na, e_Na, b, e_b, W, e_W, target, z_target,\
    RA_target, Dec_target, galaxyName, RA_galaxy, Dec_galaxy, impact, azimuth, pa, inclination,\
    adjustedInc, l, l_cus, Rvir, cus, diam, morph, vhel, vcorr, bestDist, e_bestDist,\
    group_num, group_mem, group_dist, Lstar, e_Lstar, B_mag):

 
    Lya_vs = L_two_plus['Lya_vs']
    Lya_vs.append(float(Lya_v))

    e_Lya_vs = L_two_plus['e_Lya_vs']
    e_Lya_vs.append(float(e_Lya_v))

    Lya_Ws = L_two_plus['Lya_Ws']
    Lya_Ws.append(float(Lya_W))
    
    e_Lya_Ws = L_two_plus['e_Lya_Ws']
    e_Lya_Ws.append(float(e_Lya_s))

    Nas = L_two_plus['Nas']
    Nas.append(float(Na))
    
    e_Nas = L_two_plus['e_Nas']
    e_Nas.append(float(e_Na))
    
    bs = L_two_plus['bs']
    bs.append(float(b))

    e_bs = L_two_plus['e_bs']
    e_bs.append(float(e_b))

    Ws = L_two_plus['Ws']
    Ws.append(float(Ws))
    
    e_Ws = L_two_plus['e_Ws']
    e_Ws.append(float(e_W))

    targets = L_two_plus['targets']
    targets.append(target)

    z_targets = L_two_plus['z_targets']
    z_targets.append(float(z_target))

    RA_targets = L_two_plus['RA_targets']
    RA_targets.append(float(RA_target))
    
    Dec_targets = L_two_plus['Dec_targets']
    Dec_targets.append(float(Dec_target))
    
    galaxyNames = L_two_plus['galaxyNames']
    galaxyNames.append(galaxyName)

    RA_galaxies = L_two_plus['RA_galaxies']
    RA_galaxies.append(float(RA_galaxy))
    
    Dec_galaxies = L_two_plus['Dec_galaxies']
    Dec_galaxies.append(float(Dec_galaxy))
    
    impacts = L_two_plus['impacts']
    impacts.append(float(impact))
    
    azimuths = L_two_plus['azimuths']
    azimuths.append(float(azimuth))
    
    pas = L_two_plus['pas']
    pas.append(float(pa))
    
    inclinations = L_two_plus['inclinations']
    inclinations.append(float(inclination))
    
    adjustedIncs = L_two_plus['adjustedIncs']
    adjustedIncs.append(float(adjustedInc))
    
    ls = L_two_plus['ls']
    ls.append(float(l))

    l_cuss = L_two_plus['l_cuss']
    l_cuss.append(float(l_cus))

    Rvirs = L_two_plus['Rvirs']
    Rvirs.append(float(Rvir))
    
    cuss = L_two_plus['cuss']
    cuss.append(float(cus))

    diams = L_two_plus['diams']
    diams.append(float(diam))
    
    morphs = L_two_plus['morphs']
    morphs.append(morph)
    
    vhels = L_two_plus['vhels']
    vhels.append(float(vhel))
    
    vcorrs = L_two_plus['vcorrs']
    vcorrs.append(float(vcorr))
    
    bestDists = L_two_plus['bestDists']
    bestDists.append(float(bestDist))
    
    e_bestDists = L_two_plus['e_bestDists']
    e_bestDists.append(float(e_bestDist))
    
    group_nums = L_two_plus['group_nums']
    group_nums.append(float(group_num))
    
    group_mems = L_two_plus['group_mems']
    group_mems.append(float(group_mem))

    group_dists = L_two_plus['group_dists']
    group_dists.append(float(group_dist))
    
    Lstars = L_two_plus['Lstars']
    Lstars.append(float(Lstar))
    
    e_Lstars = L_two_plus['e_Lstars']
    e_Lstars.append(float(e_Lstar))
    
    B_mags = L_two_plus['B_mags']
    B_mags.append(float(B_mag))
    
    pass

    


def add_to_L_group(Lya_v, e_Lya_v, Lya_W, e_Lya_W, Na, e_Na, b, e_b, W, e_W, target, z_target,\
    RA_target, Dec_target, galaxyName, RA_galaxy, Dec_galaxy, impact, azimuth, pa, inclination,\
    adjustedInc, l, l_cus, Rvir, cus, diam, morph, vhel, vcorr, bestDist, e_bestDist,\
    group_num, group_mem, group_dist, Lstar, e_Lstar, B_mag):

 
    Lya_vs = L_group['Lya_vs']
    Lya_vs.append(float(Lya_v))

    e_Lya_vs = L_group['e_Lya_vs']
    e_Lya_vs.append(float(e_Lya_v))

    Lya_Ws = L_group['Lya_Ws']
    Lya_Ws.append(float(Lya_W))
    
    e_Lya_Ws = L_group['e_Lya_Ws']
    e_Lya_Ws.append(float(e_Lya_s))

    Nas = L_group['Nas']
    Nas.append(float(Na))
    
    e_Nas = L_group['e_Nas']
    e_Nas.append(float(e_Na))
    
    bs = L_group['bs']
    bs.append(float(b))

    e_bs = L_group['e_bs']
    e_bs.append(float(e_b))

    Ws = L_group['Ws']
    Ws.append(float(Ws))
    
    e_Ws = L_group['e_Ws']
    e_Ws.append(float(e_W))

    targets = L_group['targets']
    targets.append(target)

    z_targets = L_group['z_targets']
    z_targets.append(float(z_target))

    RA_targets = L_group['RA_targets']
    RA_targets.append(float(RA_target))
    
    Dec_targets = L_group['Dec_targets']
    Dec_targets.append(float(Dec_target))
    
    galaxyNames = L_group['galaxyNames']
    galaxyNames.append(galaxyName)

    RA_galaxies = L_group['RA_galaxies']
    RA_galaxies.append(float(RA_galaxy))
    
    Dec_galaxies = L_group['Dec_galaxies']
    Dec_galaxies.append(float(Dec_galaxy))
    
    impacts = L_group['impacts']
    impacts.append(float(impact))
    
    azimuths = L_group['azimuths']
    azimuths.append(float(azimuth))
    
    pas = L_group['pas']
    pas.append(float(pa))
    
    inclinations = L_group['inclinations']
    inclinations.append(float(inclination))
    
    adjustedIncs = L_group['adjustedIncs']
    adjustedIncs.append(float(adjustedInc))
    
    ls = L_group['ls']
    ls.append(float(l))

    l_cuss = L_group['l_cuss']
    l_cuss.append(float(l_cus))

    Rvirs = L_group['Rvirs']
    Rvirs.append(float(Rvir))
    
    cuss = L_group['cuss']
    cuss.append(float(cus))

    diams = L_group['diams']
    diams.append(float(diam))
    
    morphs = L_group['morphs']
    morphs.append(morph)
    
    vhels = L_group['vhels']
    vhels.append(float(vhel))
    
    vcorrs = L_group['vcorrs']
    vcorrs.append(float(vcorr))
    
    bestDists = L_group['bestDists']
    bestDists.append(float(bestDist))
    
    e_bestDists = L_group['e_bestDists']
    e_bestDists.append(float(e_bestDist))
    
    group_nums = L_group['group_nums']
    group_nums.append(float(group_num))
    
    group_mems = L_group['group_mems']
    group_mems.append(float(group_mem))

    group_dists = L_group['group_dists']
    group_dists.append(float(group_dist))
    
    Lstars = L_group['Lstars']
    Lstars.append(float(Lstar))
    
    e_Lstars = L_group['e_Lstars']
    e_Lstars.append(float(e_Lstar))
    
    B_mags = L_group['B_mags']
    B_mags.append(float(B_mag))
    
    pass


    
    
def main():
    # correlation options
    maxSep = 500.
    agnSeparation = 4000.
    minVcorr = 450.
    minSize = 0.
    max_deltav = 400.
    min_likelihood = 0.01
    rigor = 5


    # assuming 'theFile' contains one name per line, read the file
    
    if getpass.getuser() == 'frenchd':
#         filename = '/Users/frenchd/Research/inclination/git_inclination/maps/LG_correlation_combined5_14_edit.csv'
#         pickleFilename = '/Users/frenchd/Research/inclination/git_inclination/picklePilot_plusSALT_14.p'

        targetlist_filename = '/Users/frenchd/Research/correlation/TARGETLIST_10_17_17_TOTAL.csv'

        filename = '/Users/frenchd/Research/inclination/git_inclination/correlatedTargetList_5_29_18_measurements.csv'
        
        # pickle files
        isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/isolated.p'
        L_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_isolated.p'
        L_associated_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated.p'
        L_associated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated.p'
        L_nonassociated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated.p'
        L_two_filename = '/Users/frenchd/Research/inclination/git_inclination/L_two.p'
        L_two_plus_filename = '/Users/frenchd/Research/inclination/git_inclination/L_two_plus.p'
        L_group_filename = '/Users/frenchd/Research/inclination/git_inclination/L_group.p'



    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    
    theFile = open(filename,'rU')
    pickleFile = open(pickleFilename,'wt')
    
    reader = csv.DictReader(theFile)
    
    # open the pickle files
    isolated_pickle_file = open(pickleFilename,'wt')
    L_isolated_pickle_file = open(pickleFilename,'wt')
    L_associated_isolated_file = open(pickleFilename,'wt')
    L_associated_file = open(pickleFilename,'wt')
    L_nonassociated_file = open(pickleFilename,'wt')
    L_two_file = open(pickleFilename,'wt')
    L_two_plus_file = open(pickleFilename,'wt')
    L_group_file = open(pickleFilename,'wt')

    
    # here are the dictionaries to be pickled. each one will the lists which follow
    isolated = {}
    L_isolated = {}
    L_associated_isolated = {}
    L_associated = {}
    L_nonassociated = {}
    L_two = {}
    L_two_plus = {}
    L_group = {}
    
    
    # isolated: add empty lists
    isolated['Lya_vs'] = []
    isolated['e_Lya_vs'] = []
    isolated['Lya_Ws'] = []
    isolated['e_Lya_Ws'] = []
    isolated['Nas'] = []
    isolated['e_Nas'] = []
    isolated['bs'] = []
    isolated['e_bs'] = []
    isolated['Ws'] = []
    isolated['e_Ws'] = []
    isolated['targets'] = []
    isolated['target_zs'] = []
    isolated['RAs'] = []
    isolated['Decs'] = []

    
    
    # L_isolated: add empty lists
    L_isolated['Lya_vs'] = []
    L_isolated['e_Lya_vs'] = []
    L_isolated['Lya_Ws'] = []
    L_isolated['e_Lya_Ws'] = []
    L_isolated['Nas'] = []
    L_isolated['e_Nas'] = []
    L_isolated['bs'] = []
    L_isolated['e_bs'] = []
    L_isolated['Ws'] = []
    L_isolated['e_Ws'] = []
    L_isolated['targets'] = []
    L_isolated['target_zs'] = []
    L_isolated['RAs'] = []
    L_isolated['Decs'] = []
    
    
    # L_associated_isolated: add empty lists
    L_associated_isolated['Lya_vs'] = []
    L_associated_isolated['e_Lya_vs'] = []
    L_associated_isolated['Lya_Ws'] = []
    L_associated_isolated['e_Lya_Ws'] = []
    L_associated_isolated['Nas'] = []
    L_associated_isolated['e_Nas'] = []
    L_associated_isolated['bs'] = []
    L_associated_isolated['e_bs'] = []
    L_associated_isolated['Ws'] = []
    L_associated_isolated['e_Ws'] = []
    L_associated_isolated['targets'] = []
    L_associated_isolated['z_targets'] = []
    L_associated_isolated['RA_targets'] = []
    L_associated_isolated['Dec_targets'] = []
    L_associated_isolated['galaxyNames'] = []
    L_associated_isolated['RA_galaxies'] = []
    L_associated_isolated['Dec_galaxies'] = []
    L_associated_isolated['impacts'] = []
    L_associated_isolated['azimuths'] = []
    L_associated_isolated['pas'] = []
    L_associated_isolated['inclinations'] = []
    L_associated_isolated['adjustedIncs'] = []
    L_associated_isolated['ls'] = []
    L_associated_isolated['l_cuss'] = []
    L_associated_isolated['Rvirs'] = []
    L_associated_isolated['cuss'] = []
    L_associated_isolated['diams'] = []
    L_associated_isolated['morphs'] = []
    L_associated_isolated['vhels'] = []
    L_associated_isolated['vcorrs'] = []
    L_associated_isolated['bestDists'] = []
    L_associated_isolated['e_bestDists'] = []
    L_associated_isolated['group_nums'] = []
    L_associated_isolated['group_mems'] = []
    L_associated_isolated['group_dists'] = []
    L_associated_isolated['Lstars'] = []
    L_associated_isolated['e_Lstars'] = []
    L_associated_isolated['B_mags'] = []

    
    # L_associated: add empty lists
    L_associated['Lya_vs'] = []
    L_associated['e_Lya_vs'] = []
    L_associated['Lya_Ws'] = []
    L_associated['e_Lya_Ws'] = []
    L_associated['Nas'] = []
    L_associated['e_Nas'] = []
    L_associated['bs'] = []
    L_associated['e_bs'] = []
    L_associated['Ws'] = []
    L_associated['e_Ws'] = []
    L_associated['targets'] = []
    L_associated['z_targets'] = []
    L_associated['RA_targets'] = []
    L_associated['Dec_targets'] = []
    L_associated['galaxyNames'] = []
    L_associated['RA_galaxies'] = []
    L_associated['Dec_galaxies'] = []
    L_associated['impacts'] = []
    L_associated['azimuths'] = []
    L_associated['pas'] = []
    L_associated['inclinations'] = []
    L_associated['adjustedIncs'] = []
    L_associated['ls'] = []
    L_associated['l_cuss'] = []
    L_associated['Rvirs'] = []
    L_associated['cuss'] = []
    L_associated['diams'] = []
    L_associated['morphs'] = []
    L_associated['vhels'] = []
    L_associated['vcorrs'] = []
    L_associated['bestDists'] = []
    L_associated['e_bestDists'] = []
    L_associated['group_nums'] = []
    L_associated['group_mems'] = []
    L_associated['group_dists'] = []
    L_associated['Lstars'] = []
    L_associated['e_Lstars'] = []
    L_associated['B_mags'] = []
    
    
    # L_associated: add empty lists
    L_two['Lya_vs'] = []
    L_two['e_Lya_vs'] = []
    L_two['Lya_Ws'] = []
    L_two['e_Lya_Ws'] = []
    L_two['Nas'] = []
    L_two['e_Nas'] = []
    L_two['bs'] = []
    L_two['e_bs'] = []
    L_two['Ws'] = []
    L_two['e_Ws'] = []
    L_two['targets'] = []
    L_two['z_targets'] = []
    L_two['RA_targets'] = []
    L_two['Dec_targets'] = []
    L_two['galaxyNames'] = []
    L_two['RA_galaxies'] = []
    L_two['Dec_galaxies'] = []
    L_two['impacts'] = []
    L_two['azimuths'] = []
    L_two['pas'] = []
    L_two['inclinations'] = []
    L_two['adjustedIncs'] = []
    L_two['ls'] = []
    L_two['l_cuss'] = []
    L_two['Rvirs'] = []
    L_two['cuss'] = []
    L_two['diams'] = []
    L_two['morphs'] = []
    L_two['vhels'] = []
    L_two['vcorrs'] = []
    L_two['bestDists'] = []
    L_two['e_bestDists'] = []
    L_two['group_nums'] = []
    L_two['group_mems'] = []
    L_two['group_dists'] = []
    L_two['Lstars'] = []
    L_two['e_Lstars'] = []
    L_two['B_mags'] = []
    
    
    # L_associated: add empty lists
    L_two_plus['Lya_vs'] = []
    L_two_plus['e_Lya_vs'] = []
    L_two_plus['Lya_Ws'] = []
    L_two_plus['e_Lya_Ws'] = []
    L_two_plus['Nas'] = []
    L_two_plus['e_Nas'] = []
    L_two_plus['bs'] = []
    L_two_plus['e_bs'] = []
    L_two_plus['Ws'] = []
    L_two_plus['e_Ws'] = []
    L_two_plus['targets'] = []
    L_two_plus['z_targets'] = []
    L_two_plus['RA_targets'] = []
    L_two_plus['Dec_targets'] = []
    L_two_plus['galaxyNames'] = []
    L_two_plus['RA_galaxies'] = []
    L_two_plus['Dec_galaxies'] = []
    L_two_plus['impacts'] = []
    L_two_plus['azimuths'] = []
    L_two_plus['pas'] = []
    L_two_plus['inclinations'] = []
    L_two_plus['adjustedIncs'] = []
    L_two_plus['ls'] = []
    L_two_plus['l_cuss'] = []
    L_two_plus['Rvirs'] = []
    L_two_plus['cuss'] = []
    L_two_plus['diams'] = []
    L_two_plus['morphs'] = []
    L_two_plus['vhels'] = []
    L_two_plus['vcorrs'] = []
    L_two_plus['bestDists'] = []
    L_two_plus['e_bestDists'] = []
    L_two_plus['group_nums'] = []
    L_two_plus['group_mems'] = []
    L_two_plus['group_dists'] = []
    L_two_plus['Lstars'] = []
    L_two_plus['e_Lstars'] = []
    L_two_plus['B_mags'] = []


    # L_associated: add empty lists
    L_group['Lya_vs'] = []
    L_group['e_Lya_vs'] = []
    L_group['Lya_Ws'] = []
    L_group['e_Lya_Ws'] = []
    L_group['Nas'] = []
    L_group['e_Nas'] = []
    L_group['bs'] = []
    L_group['e_bs'] = []
    L_group['Ws'] = []
    L_group['e_Ws'] = []
    L_group['targets'] = []
    L_group['z_targets'] = []
    L_group['RA_targets'] = []
    L_group['Dec_targets'] = []
    L_group['galaxyNames'] = []
    L_group['RA_galaxies'] = []
    L_group['Dec_galaxies'] = []
    L_group['impacts'] = []
    L_group['azimuths'] = []
    L_group['pas'] = []
    L_group['inclinations'] = []
    L_group['adjustedIncs'] = []
    L_group['ls'] = []
    L_group['l_cuss'] = []
    L_group['Rvirs'] = []
    L_group['cuss'] = []
    L_group['diams'] = []
    L_group['morphs'] = []
    L_group['vhels'] = []
    L_group['vcorrs'] = []
    L_group['bestDists'] = []
    L_group['e_bestDists'] = []
    L_group['group_nums'] = []
    L_group['group_mems'] = []
    L_group['group_dists'] = []
    L_group['Lstars'] = []
    L_group['e_Lstars'] = []
    L_group['B_mags'] = []



##########################################################################################
    # grab the target coordinates
    targetFile = open(targetlist_filename,'rU')
    targetReader = csv.DictReader(targetFile)
    
    target_coords = {}
    for t in targetReader:
        name = t['targetName']
        ra = t['degreesRA']
        dec = t['degreesDec']
        
        if not target_coords.has_key(target_coords):
            target_coords[name] = {'RA':ra,'Dec':dec}
    
    targetFile.close()
##########################################################################################
##########################################################################################
    # now the full data set
    
    for i in reader:
        target = i['target']
        
        target_ra = target_coords[target]['RA']
        target_dec = target_coords[target]['Dec']

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

        correlation = correlateTarget(target, maxSep, agnSeparation, minVcorr, minSize, slow=False, searchAll=False, returnDict=True):

        # if no galaxies are returned, add it to the isolated list
        if len(correlation) == 0:
            add_to_isolated(Lya_v, e_Lya_v, Lya_W, e_Lya_W, Na, e_Na, b, e_b, W, e_W, target, target_z, RA, Dec)
            
        else:
            candidates = []
            candidates_cus = []
            others = []
            for c in correlation:
                # unpack everything
                target = c['target']
                z_target = c['z_target']
                RA_target = c['RA_target']
                Dec_target = c['Dec_target']
                galaxyName c['Name']
                RA_galaxy = c['RAdeg']
                Dec_galaxy = c['DEdeg']
                impact = c['impact_parameter']
                azimuths = c['azimuth']
                pa = c['PA']
                inc = c['inc']
                adjustedInc c['adjustedIncs']
        
                R_vir = c['R_vir']
                MajDiam = c['MajDiam']
                MType = c['MType']
                Vhel = c['Vhel']
                vcorr = c['vcorr']
                bestDist = c['bestDist']
                e_bestDist = c['e_bestDist']
                group_num = c['group_num']
                group_mem = c['group_mem']
                group_dist = c['group_dist']
                Lstar = c['Lstars']
                e_Lstar = c['e_Lstars']
                B_mag = c['B_mags']
                
                # velocity difference = V_absorber - V_sys
                vel_dif = Lya_v - Vhel
    
                # check if the line is too far from a galaxy
                if abs(Lya_v - Vhel) > max_deltav:
                    # add to isolated if too far
                    add_to_isolated(Lya_v, e_Lya_v, Lya_W, e_Lya_W, Na, e_Na, b, e_b, W, e_W, target, z_target, RA_target, Dec_target)
                
                else:
                    # if there's a galaxy within max_deltav, calculate likelihood values
                
                    # try this "sphere of influence" value instead
                    cus = MajDiam**1.5

                    # first for the virial radius
                    likelihood = math.exp(-(impact/R_vir)**2) * math.exp(-(vel_dif/200.)**2)
                        
                    # then for the second 'virial like' m15 radius
                    likelihood_cus = math.exp(-(impact/m15)**2) * math.exp(-(vel_dif/200.)**2)                
                
#                     if rVir>= impact:
#                         likelihood = likelihood*2

                    # add likelihoods and stuff to the dictionary entries
                    correlation[c]['l'] = likelihood
                    correlation[c]['l_cus'] = likelihood_cus
                    correlation[c]['cus'] = cus
                    
                    # if they make the cut, add them to the candidate lists
                    if likelihood >= min_likelihood:
                        candidates.append([likelihood,correlation[c]])
                        
                    else:
                        others.append([likelihood,correlation[c]])
                    
                    else:
                        # otherwise this is an L_isolated absorber
                
            # now candidates have been added to that list and isolated have been taken care of.
            # handle the candidate list:
            #
            # if there are no candidates:
            if len(candidates) == 0:
                
                # then this is an L_isolated absorber
                # just grab the first 'other' galaxy and take the absorber info from it
                l, galaxy_info = others[0]
                add_to_L_isolated(
                galaxy_info['Lya_v'],\
                galaxy_info['e_Lya_v'],\
                galaxy_info['Lya_W'],\
                galaxy_info['e_Lya_W'],\
                galaxy_info['Na'],\
                galaxy_info['e_Na'],\
                galaxy_info['b'],\
                galaxy_info['e_b'],\
                galaxy_info['W'],\
                galaxy_info['e_W'],\
                galaxy_info['target'],\
                galaxy_info['z_target'],\
                galaxy_info['RA_target'],\
                galaxy_info['Dec_target'])
            
            
            # first sort by likelihood:
            candidates.sort(reverse=True)
            others.sort(reverse=True)
            
            # if there's only one entry in both, then it's an L_associated_isolated
            if len(candidates) == 1 and len(correlation) == 1:
                likelihood, galaxy_info = candidates[0]
                
                add_to_L_associated_isolated(galaxy_info['Lya_v'],\
                galaxy_info['e_Lya_v'],\
                galaxy_info['Lya_W'],\
                galaxy_info['e_Lya_W'],\
                galaxy_info['Na'],\
                galaxy_info['e_Na'],\
                galaxy_info['b'],\
                galaxy_info['e_b'],\
                galaxy_info['W'],\
                galaxy_info['e_W'],\
                galaxy_info['target'],\
                galaxy_info['z_target'],\
                galaxy_info['RA_target'],\
                galaxy_info['Dec_target'],\
                galaxy_info['galaxyName'],\
                galaxy_info['RA_galaxy'],\
                galaxy_info['Dec_galaxy'],\
                galaxy_info['impact'],\
                galaxy_info['azimuth'],\
                galaxy_info['pa'],\
                galaxy_info['inclination'],\
                galaxy_info['adjustedInc'],\
                galaxy_info['l'],\
                galaxy_info['l_cus'],\
                galaxy_info['Rvir'],\
                galaxy_info['cus'],\
                galaxy_info['diam'],\
                galaxy_info['morph'],\
                galaxy_info['vhel'],\
                galaxy_info['vcorr'],\
                galaxy_info['bestDist'],\
                galaxy_info['e_bestDist'],\
                galaxy_info['group_num'],\
                galaxy_info['group_mem'],\
                galaxy_info['group_dist'],\
                galaxy_info['Lstar'],\
                galaxy_info['e_Lstar'],\
                galaxy_info['B_mag'])
        
            # if there's only one with L > min_L, but others exist
            elif len(candidates) == 1 and len(correlation) > 1:
                likelihood, galaxy_info = candidates[0]
                
                # check to see if any of the L < min_L are within *rigor* of L
                others_within_rigor = []
                for o in others:
                    l, other_info = o
                    if l*rigor >= likelihood:
                        others_within_rigor.append(other_info)
                        
                # L_associated, none within L_other * rigor of this galaxy's L
                if len(others_within_rigor) == 0:
                    add_to_L_associated(galaxy_info['Lya_v'],\
                    galaxy_info['e_Lya_v'],\
                    galaxy_info['Lya_W'],\
                    galaxy_info['e_Lya_W'],\
                    galaxy_info['Na'],\
                    galaxy_info['e_Na'],\
                    galaxy_info['b'],\
                    galaxy_info['e_b'],\
                    galaxy_info['W'],\
                    galaxy_info['e_W'],\
                    galaxy_info['target'],\
                    galaxy_info['z_target'],\
                    galaxy_info['RA_target'],\
                    galaxy_info['Dec_target'],\
                    galaxy_info['galaxyName'],\
                    galaxy_info['RA_galaxy'],\
                    galaxy_info['Dec_galaxy'],\
                    galaxy_info['impact'],\
                    galaxy_info['azimuth'],\
                    galaxy_info['pa'],\
                    galaxy_info['inclination'],\
                    galaxy_info['adjustedInc'],\
                    galaxy_info['l'],\
                    galaxy_info['l_cus'],\
                    galaxy_info['Rvir'],\
                    galaxy_info['cus'],\
                    galaxy_info['diam'],\
                    galaxy_info['morph'],\
                    galaxy_info['vhel'],\
                    galaxy_info['vcorr'],\
                    galaxy_info['bestDist'],\
                    galaxy_info['e_bestDist'],\
                    galaxy_info['group_num'],\
                    galaxy_info['group_mem'],\
                    galaxy_info['group_dist'],\
                    galaxy_info['Lstar'],\
                    galaxy_info['e_Lstar'],\
                    galaxy_info['B_mag'])
                
                else:
                    # there exist galaxies within *rigor* of L, but with L_other < L_min
                    # first add the candidate info:
                    add_to_L_nonassociated(galaxy_info['Lya_v'],\
                    galaxy_info['e_Lya_v'],\
                    galaxy_info['Lya_W'],\
                    galaxy_info['e_Lya_W'],\
                    galaxy_info['Na'],\
                    galaxy_info['e_Na'],\
                    galaxy_info['b'],\
                    galaxy_info['e_b'],\
                    galaxy_info['W'],\
                    galaxy_info['e_W'],\
                    galaxy_info['target'],\
                    galaxy_info['z_target'],\
                    galaxy_info['RA_target'],\
                    galaxy_info['Dec_target'],\
                    galaxy_info['galaxyName'],\
                    galaxy_info['RA_galaxy'],\
                    galaxy_info['Dec_galaxy'],\
                    galaxy_info['impact'],\
                    galaxy_info['azimuth'],\
                    galaxy_info['pa'],\
                    galaxy_info['inclination'],\
                    galaxy_info['adjustedInc'],\
                    galaxy_info['l'],\
                    galaxy_info['l_cus'],\
                    galaxy_info['Rvir'],\
                    galaxy_info['cus'],\
                    galaxy_info['diam'],\
                    galaxy_info['morph'],\
                    galaxy_info['vhel'],\
                    galaxy_info['vcorr'],\
                    galaxy_info['bestDist'],\
                    galaxy_info['e_bestDist'],\
                    galaxy_info['group_num'],\
                    galaxy_info['group_mem'],\
                    galaxy_info['group_dist'],\
                    galaxy_info['Lstar'],\
                    galaxy_info['e_Lstar'],\
                    galaxy_info['B_mag'])
                    
                    # then add the others' info
                    for o in others_within_rigor:
                        add_to_L_nonassociated(o['Lya_v'],\
                        o['e_Lya_v'],\
                        o['Lya_W'],\
                        o['e_Lya_W'],\
                        o['Na'],\
                        o['e_Na'],\
                        o['b'],\
                        o['e_b'],\
                        o['W'],\
                        o['e_W'],\
                        o['target'],\
                        o['z_target'],\
                        o['RA_target'],\
                        o['Dec_target'],\
                        o['galaxyName'],\
                        o['RA_galaxy'],\
                        o['Dec_galaxy'],\
                        o['impact'],\
                        o['azimuth'],\
                        o['pa'],\
                        o['inclination'],\
                        o['adjustedInc'],\
                        o['l'],\
                        o['l_cus'],\
                        o['Rvir'],\
                        o['cus'],\
                        o['diam'],\
                        o['morph'],\
                        o['vhel'],\
                        o['vcorr'],\
                        o['bestDist'],\
                        o['e_bestDist'],\
                        o['group_num'],\
                        o['group_mem'],\
                        o['group_dist'],\
                        o['Lstar'],\
                        o['e_Lstar'],\
                        o['B_mag'])
                        
                        
            # if there's more than one with L > min_L
            else:
                likelihood, galaxy_info = candidates[0]
                
                candidates_within_rigor = []
                
                # first compare to other candidates
                for c in candidates:
                    l, candidate_info = c
                    
                    if l * rigor >= top_likelihood:
                        candidates_within_rigor.append(candidate_info)
                        
                for o in others:
                    l, other_info = o
                    
                    if l * rigor >= top_likelihood:
                        candidates_within_rigor.append(other_info)
                    
                # if none are close enough, add this one to L_associated
                # len == 1 means the top candidate is the only one in the list
                if len(candidates_within_rigor) == 1:
                    add_to_L_associated(galaxy_info['Lya_v'],\
                    galaxy_info['e_Lya_v'],\
                    galaxy_info['Lya_W'],\
                    galaxy_info['e_Lya_W'],\
                    galaxy_info['Na'],\
                    galaxy_info['e_Na'],\
                    galaxy_info['b'],\
                    galaxy_info['e_b'],\
                    galaxy_info['W'],\
                    galaxy_info['e_W'],\
                    galaxy_info['target'],\
                    galaxy_info['z_target'],\
                    galaxy_info['RA_target'],\
                    galaxy_info['Dec_target'],\
                    galaxy_info['galaxyName'],\
                    galaxy_info['RA_galaxy'],\
                    galaxy_info['Dec_galaxy'],\
                    galaxy_info['impact'],\
                    galaxy_info['azimuth'],\
                    galaxy_info['pa'],\
                    galaxy_info['inclination'],\
                    galaxy_info['adjustedInc'],\
                    galaxy_info['l'],\
                    galaxy_info['l_cus'],\
                    galaxy_info['Rvir'],\
                    galaxy_info['cus'],\
                    galaxy_info['diam'],\
                    galaxy_info['morph'],\
                    galaxy_info['vhel'],\
                    galaxy_info['vcorr'],\
                    galaxy_info['bestDist'],\
                    galaxy_info['e_bestDist'],\
                    galaxy_info['group_num'],\
                    galaxy_info['group_mem'],\
                    galaxy_info['group_dist'],\
                    galaxy_info['Lstar'],\
                    galaxy_info['e_Lstar'],\
                    galaxy_info['B_mag'])
                    
                # if one other is within *rigor*, add to L_two list
                elif len(candidates_within_rigor) == 2:
                    for galaxy_info in candidates_within_rigor:
                        add_to_L_two(galaxy_info['Lya_v'],\
                        galaxy_info['e_Lya_v'],\
                        galaxy_info['Lya_W'],\
                        galaxy_info['e_Lya_W'],\
                        galaxy_info['Na'],\
                        galaxy_info['e_Na'],\
                        galaxy_info['b'],\
                        galaxy_info['e_b'],\
                        galaxy_info['W'],\
                        galaxy_info['e_W'],\
                        galaxy_info['target'],\
                        galaxy_info['z_target'],\
                        galaxy_info['RA_target'],\
                        galaxy_info['Dec_target'],\
                        galaxy_info['galaxyName'],\
                        galaxy_info['RA_galaxy'],\
                        galaxy_info['Dec_galaxy'],\
                        galaxy_info['impact'],\
                        galaxy_info['azimuth'],\
                        galaxy_info['pa'],\
                        galaxy_info['inclination'],\
                        galaxy_info['adjustedInc'],\
                        galaxy_info['l'],\
                        galaxy_info['l_cus'],\
                        galaxy_info['Rvir'],\
                        galaxy_info['cus'],\
                        galaxy_info['diam'],\
                        galaxy_info['morph'],\
                        galaxy_info['vhel'],\
                        galaxy_info['vcorr'],\
                        galaxy_info['bestDist'],\
                        galaxy_info['e_bestDist'],\
                        galaxy_info['group_num'],\
                        galaxy_info['group_mem'],\
                        galaxy_info['group_dist'],\
                        galaxy_info['Lstar'],\
                        galaxy_info['e_Lstar'],\
                        galaxy_info['B_mag'])
                
                elif len(candidates_within_rigor) > 2:
                    for galaxy_info in candidates_within_rigor:
                        add_to_L_two_plus(galaxy_info['Lya_v'],\
                        galaxy_info['e_Lya_v'],\
                        galaxy_info['Lya_W'],\
                        galaxy_info['e_Lya_W'],\
                        galaxy_info['Na'],\
                        galaxy_info['e_Na'],\
                        galaxy_info['b'],\
                        galaxy_info['e_b'],\
                        galaxy_info['W'],\
                        galaxy_info['e_W'],\
                        galaxy_info['target'],\
                        galaxy_info['z_target'],\
                        galaxy_info['RA_target'],\
                        galaxy_info['Dec_target'],\
                        galaxy_info['galaxyName'],\
                        galaxy_info['RA_galaxy'],\
                        galaxy_info['Dec_galaxy'],\
                        galaxy_info['impact'],\
                        galaxy_info['azimuth'],\
                        galaxy_info['pa'],\
                        galaxy_info['inclination'],\
                        galaxy_info['adjustedInc'],\
                        galaxy_info['l'],\
                        galaxy_info['l_cus'],\
                        galaxy_info['Rvir'],\
                        galaxy_info['cus'],\
                        galaxy_info['diam'],\
                        galaxy_info['morph'],\
                        galaxy_info['vhel'],\
                        galaxy_info['vcorr'],\
                        galaxy_info['bestDist'],\
                        galaxy_info['e_bestDist'],\
                        galaxy_info['group_num'],\
                        galaxy_info['group_mem'],\
                        galaxy_info['group_dist'],\
                        galaxy_info['Lstar'],\
                        galaxy_info['e_Lstar'],\
                        galaxy_info['B_mag'])
            
            
            # now deal with group galaxies
            if len(candidates) >= 1:
                for c in candidates:
                    l, galaxy_info = c
                    group_num = galaxy_info['group_num']
                    
                    if group_num != -99:
                        add_to_L_group(galaxy_info['Lya_v'],\
                        galaxy_info['e_Lya_v'],\
                        galaxy_info['Lya_W'],\
                        galaxy_info['e_Lya_W'],\
                        galaxy_info['Na'],\
                        galaxy_info['e_Na'],\
                        galaxy_info['b'],\
                        galaxy_info['e_b'],\
                        galaxy_info['W'],\
                        galaxy_info['e_W'],\
                        galaxy_info['target'],\
                        galaxy_info['z_target'],\
                        galaxy_info['RA_target'],\
                        galaxy_info['Dec_target'],\
                        galaxy_info['galaxyName'],\
                        galaxy_info['RA_galaxy'],\
                        galaxy_info['Dec_galaxy'],\
                        galaxy_info['impact'],\
                        galaxy_info['azimuth'],\
                        galaxy_info['pa'],\
                        galaxy_info['inclination'],\
                        galaxy_info['adjustedInc'],\
                        galaxy_info['l'],\
                        galaxy_info['l_cus'],\
                        galaxy_info['Rvir'],\
                        galaxy_info['cus'],\
                        galaxy_info['diam'],\
                        galaxy_info['morph'],\
                        galaxy_info['vhel'],\
                        galaxy_info['vcorr'],\
                        galaxy_info['bestDist'],\
                        galaxy_info['e_bestDist'],\
                        galaxy_info['group_num'],\
                        galaxy_info['group_mem'],\
                        galaxy_info['group_dist'],\
                        galaxy_info['Lstar'],\
                        galaxy_info['e_Lstar'],\
                        galaxy_info['B_mag'])


##########################################################################################
##########################################################################################
##########################################################################################
    pickle.dump(isolated, isolated_pickle_file)
    pickle.dump(L_isolated, L_isolated_pickle_file)
    pickle.dump(L_associated_isolated, L_associated_isolated_file)
    pickle.dump(L_associated, L_associated_file)
    pickle.dump(L_nonassociated, L_nonassociated_file)
    pickle.dump(L_two, L_two_file)
    pickle.dump(L_two_plus, L_two_plus_file)
    pickle.dump(L_group, L_group_file)

    
    isolated_pickle_file.close()
    L_isolated_pickle_file.close()
    L_associated_isolated_file.close()
    L_associated_file.close()
    L_nonassociated_file.close()
    L_two_file.close()
    L_two_plus_file.close()
    L_group_file.close()

    theFile.close()
    
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    