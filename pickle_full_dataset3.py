#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)


$Id: pickle_full_dataset3.py v 2.3 06/21/2018

Create a pickle file with already correlated results

v2.1: Also pickle everything all in one, and one where all the likelihoods are summed for
each line (06/11/18)

v2.3: Add Lstar_sdss. Add more options for L value separation, such as the old method
of multiplying by 2 when impact <= rvir, etc.


Output bins:
1. Isolated (no galaxies within 500 kpc, 400 km/s)
2. No galaxies with L >= L_min
3. 1 galaxy with L >= L_min, isolated (no galaxies within 500 kpc, 400 km/s)
4. 1 galaxy with L >= L_min, no neighbors within a factor of rigor*L
5. 2 galaxies with L >=L_min within a factor of rigor*L, no others
6. More than 2 galaxies with L >= L_min within a factor of rigor*L
7. At least 1 galaxy with L >= L_min where that galaxy is a designated group member (Based on Tully 2015 2MASS group catalog)
8. All lines
9. L_summed, which includes all lines with L > L_min with the list of all nearby galaxy L values


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
    
    
    
    
def add_to_all(Lya_v, e_Lya_v, Lya_W, e_Lya_W, Na, e_Na, b, e_b, W, e_W, target, z_target, RA_target, Dec_target):
    Lya_vs = all['Lya_vs']
    Lya_vs.append(float(Lya_v))

    e_Lya_vs = all['e_Lya_vs']
    e_Lya_vs.append(float(e_Lya_v))

    Lya_Ws = all['Lya_Ws']
    Lya_Ws.append(float(Lya_W))
    
    e_Lya_Ws = all['e_Lya_Ws']
    e_Lya_Ws.append(float(e_Lya_W))

    Nas = all['Nas']
    Nas.append(float(Na))
    
    e_Nas = all['e_Nas']
    e_Nas.append(float(e_Na))
    
    bs = all['bs']
    bs.append(float(b))

    e_bs = all['e_bs']
    e_bs.append(float(e_b))

    Ws = all['Ws']
    Ws.append(float(W))
    
    e_Ws = all['e_Ws']
    e_Ws.append(float(e_W))

    targets = all['targets']
    targets.append(target)

    z_targets = all['z_targets']
    z_targets.append(float(z_target))

    RA_targets = all['RA_targets']
    RA_targets.append(float(RA_target))
    
    Dec_targets = all['Dec_targets']
    Dec_targets.append(float(Dec_target))
    
    pass
    
    
    
def add_to_isolated(Lya_v, e_Lya_v, Lya_W, e_Lya_W, Na, e_Na, b, e_b, W, e_W, target, z_target, RA_target, Dec_target):
    Lya_vs = isolated['Lya_vs']
    Lya_vs.append(float(Lya_v))

    e_Lya_vs = isolated['e_Lya_vs']
    e_Lya_vs.append(float(e_Lya_v))

    Lya_Ws = isolated['Lya_Ws']
    Lya_Ws.append(float(Lya_W))
    
    e_Lya_Ws = isolated['e_Lya_Ws']
    e_Lya_Ws.append(float(e_Lya_W))

    Nas = isolated['Nas']
    Nas.append(float(Na))
    
    e_Nas = isolated['e_Nas']
    e_Nas.append(float(e_Na))
    
    bs = isolated['bs']
    bs.append(float(b))

    e_bs = isolated['e_bs']
    e_bs.append(float(e_b))

    Ws = isolated['Ws']
    Ws.append(float(W))
    
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
    e_Lya_Ws.append(float(e_Lya_W))

    Nas = L_isolated['Nas']
    Nas.append(float(Na))
    
    e_Nas = L_isolated['e_Nas']
    e_Nas.append(float(e_Na))
    
    bs = L_isolated['bs']
    bs.append(float(b))

    e_bs = L_isolated['e_bs']
    e_bs.append(float(e_b))

    Ws = L_isolated['Ws']
    Ws.append(float(W))
    
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
    RA_target, Dec_target, Name, RAdeg, DEdeg, impact, azimuth, PA, inc,\
    adjustedInc, l, l_cus, R_vir, cus, MajDiam, MType, Vhel, vcorr, bestDist, e_bestDist,\
    group_num, group_mem, group_dist, Lstar_med, e_Lstar_med, Lstar_sdss, e_Lstar_sdss,\
    Bmag, Bmag_sdss):


    Lya_vs = L_associated_isolated['Lya_vs']
    Lya_vs.append(float(Lya_v))

    e_Lya_vs = L_associated_isolated['e_Lya_vs']
    e_Lya_vs.append(float(e_Lya_v))

    Lya_Ws = L_associated_isolated['Lya_Ws']
    Lya_Ws.append(float(Lya_W))
    
    e_Lya_Ws = L_associated_isolated['e_Lya_Ws']
    e_Lya_Ws.append(float(e_Lya_W))

    Nas = L_associated_isolated['Nas']
    Nas.append(float(Na))
    
    e_Nas = L_associated_isolated['e_Nas']
    e_Nas.append(float(e_Na))
    
    bs = L_associated_isolated['bs']
    bs.append(float(b))

    e_bs = L_associated_isolated['e_bs']
    e_bs.append(float(e_b))

    Ws = L_associated_isolated['Ws']
    Ws.append(float(W))
    
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
    
    Names = L_associated_isolated['Names']
    Names.append(Name)

    RA_galaxies = L_associated_isolated['RA_galaxies']
    RA_galaxies.append(float(RAdeg))
    
    Dec_galaxies = L_associated_isolated['Dec_galaxies']
    Dec_galaxies.append(float(DEdeg))
    
    impacts = L_associated_isolated['impacts']
    impacts.append(float(impact))
    
    azimuths = L_associated_isolated['azimuths']
    azimuths.append(float(azimuth))
    
    PAs = L_associated_isolated['PAs']
    PAs.append(float(PA))
    
    incs = L_associated_isolated['incs']
    incs.append(float(inc))
    
    adjustedIncs = L_associated_isolated['adjustedIncs']
    adjustedIncs.append(float(adjustedInc))
    
    ls = L_associated_isolated['ls']
    ls.append(float(l))

    l_cuss = L_associated_isolated['l_cuss']
    l_cuss.append(float(l_cus))

    R_virs = L_associated_isolated['R_virs']
    R_virs.append(float(R_vir))
    
    cuss = L_associated_isolated['cuss']
    cuss.append(float(cus))

    MajDiams = L_associated_isolated['MajDiams']
    MajDiams.append(float(MajDiam))
    
    MTypes = L_associated_isolated['MTypes']
    MTypes.append(MType)
    
    Vhels = L_associated_isolated['Vhels']
    Vhels.append(float(Vhel))
    
    vcorrs = L_associated_isolated['vcorrs']
    vcorrs.append(float(vcorr))
    
    bestDists = L_associated_isolated['bestDists']
    bestDists.append(float(bestDist))
    
    e_bestDists = L_associated_isolated['e_bestDists']
    e_bestDists.append(float(e_bestDist))
    
    group_nums = L_associated_isolated['group_nums']
    group_nums.append(group_num)
    
    group_mems = L_associated_isolated['group_mems']
    group_mems.append(group_mem)

    group_dists = L_associated_isolated['group_dists']
    group_dists.append(float(group_dist))
    
    Lstar_meds = L_associated_isolated['Lstar_meds']
    Lstar_meds.append(float(Lstar_med))
    
    e_Lstar_meds = L_associated_isolated['e_Lstar_meds']
    e_Lstar_meds.append(float(e_Lstar_med))

    Lstar_sdsss = L_associated_isolated['Lstar_sdsss']
    Lstar_sdsss.append(float(Lstar_sdss))
    
    e_Lstar_sdsss = L_associated_isolated['e_Lstar_sdsss']
    e_Lstar_sdsss.append(float(e_Lstar_sdss))

    Bmags = L_associated_isolated['Bmags']
    Bmags.append(float(Bmag))
    
    Bmag_sdsss = L_associated_isolated['Bmag_sdsss']
    Bmag_sdsss.append(float(Bmag_sdss))
    
    pass
    



def add_to_L_associated(Lya_v, e_Lya_v, Lya_W, e_Lya_W, Na, e_Na, b, e_b, W, e_W, target, z_target,\
    RA_target, Dec_target, Name, RAdeg, DEdeg, impact, azimuth, PA, inc,\
    adjustedInc, l, l_cus, R_vir, cus, MajDiam, MType, Vhel, vcorr, bestDist, e_bestDist,\
    group_num, group_mem, group_dist, Lstar_med, e_Lstar_med, Lstar_sdss, e_Lstar_sdss,\
    Bmag, Bmag_sdss):


    Lya_vs = L_associated['Lya_vs']
    Lya_vs.append(float(Lya_v))

    e_Lya_vs = L_associated['e_Lya_vs']
    e_Lya_vs.append(float(e_Lya_v))

    Lya_Ws = L_associated['Lya_Ws']
    Lya_Ws.append(float(Lya_W))
    
    e_Lya_Ws = L_associated['e_Lya_Ws']
    e_Lya_Ws.append(float(e_Lya_W))

    Nas = L_associated['Nas']
    Nas.append(float(Na))
    
    e_Nas = L_associated['e_Nas']
    e_Nas.append(float(e_Na))
    
    bs = L_associated['bs']
    bs.append(float(b))

    e_bs = L_associated['e_bs']
    e_bs.append(float(e_b))

    Ws = L_associated['Ws']
    Ws.append(float(W))
    
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
    
    Names = L_associated['Names']
    Names.append(Name)

    RA_galaxies = L_associated['RA_galaxies']
    RA_galaxies.append(float(RAdeg))
    
    Dec_galaxies = L_associated['Dec_galaxies']
    Dec_galaxies.append(float(DEdeg))
    
    impacts = L_associated['impacts']
    impacts.append(float(impact))
    
    azimuths = L_associated['azimuths']
    azimuths.append(float(azimuth))
    
    PAs = L_associated['PAs']
    PAs.append(float(PA))
    
    incs = L_associated['incs']
    incs.append(float(inc))
    
    adjustedIncs = L_associated['adjustedIncs']
    adjustedIncs.append(float(adjustedInc))
    
    ls = L_associated['ls']
    ls.append(float(l))

    l_cuss = L_associated['l_cuss']
    l_cuss.append(float(l_cus))

    R_virs = L_associated['R_virs']
    R_virs.append(float(R_vir))
    
    cuss = L_associated['cuss']
    cuss.append(float(cus))

    MajDiams = L_associated['MajDiams']
    MajDiams.append(float(MajDiam))
    
    MTypes = L_associated['MTypes']
    MTypes.append(MType)
    
    Vhels = L_associated['Vhels']
    Vhels.append(float(Vhel))
    
    vcorrs = L_associated['vcorrs']
    vcorrs.append(float(vcorr))
    
    bestDists = L_associated['bestDists']
    bestDists.append(float(bestDist))
    
    e_bestDists = L_associated['e_bestDists']
    e_bestDists.append(float(e_bestDist))
    
    group_nums = L_associated['group_nums']
    group_nums.append(group_num)
    
    group_mems = L_associated['group_mems']
    group_mems.append(group_mem)

    group_dists = L_associated['group_dists']
    group_dists.append(float(group_dist))
    
    Lstar_meds = L_associated['Lstar_meds']
    Lstar_meds.append(float(Lstar_med))
    
    e_Lstar_meds = L_associated['e_Lstar_meds']
    e_Lstar_meds.append(float(e_Lstar_med))
    
    Lstar_sdsss = L_associated['Lstar_sdsss']
    Lstar_sdsss.append(float(Lstar_sdss))
    
    e_Lstar_sdsss = L_associated['e_Lstar_sdsss']
    e_Lstar_sdsss.append(float(e_Lstar_sdss))
    
    Bmags = L_associated['Bmags']
    Bmags.append(float(Bmag))

    Bmag_sdsss = L_associated['Bmag_sdsss']
    Bmag_sdsss.append(float(Bmag_sdss))
    
    pass




def add_to_L_nonassociated(Lya_v, e_Lya_v, Lya_W, e_Lya_W, Na, e_Na, b, e_b, W, e_W, target, z_target,\
    RA_target, Dec_target, Name, RAdeg, DEdeg, impact, azimuth, PA, inc,\
    adjustedInc, l, l_cus, R_vir, cus, MajDiam, MType, Vhel, vcorr, bestDist, e_bestDist,\
    group_num, group_mem, group_dist, Lstar_med, e_Lstar_med, Lstar_sdss, e_Lstar_sdss,
    Bmag, Bmag_sdss, other):


    Lya_vs = L_nonassociated['Lya_vs']
    Lya_vs.append(float(Lya_v))

    e_Lya_vs = L_nonassociated['e_Lya_vs']
    e_Lya_vs.append(float(e_Lya_v))

    Lya_Ws = L_nonassociated['Lya_Ws']
    Lya_Ws.append(float(Lya_W))
    
    e_Lya_Ws = L_nonassociated['e_Lya_Ws']
    e_Lya_Ws.append(float(e_Lya_W))

    Nas = L_nonassociated['Nas']
    Nas.append(float(Na))
    
    e_Nas = L_nonassociated['e_Nas']
    e_Nas.append(float(e_Na))
    
    bs = L_nonassociated['bs']
    bs.append(float(b))

    e_bs = L_nonassociated['e_bs']
    e_bs.append(float(e_b))

    Ws = L_nonassociated['Ws']
    Ws.append(float(W))
    
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
    
    Names = L_nonassociated['Names']
    Names.append(Name)

    RA_galaxies = L_nonassociated['RA_galaxies']
    RA_galaxies.append(float(RAdeg))
    
    Dec_galaxies = L_nonassociated['Dec_galaxies']
    Dec_galaxies.append(float(DEdeg))
    
    impacts = L_nonassociated['impacts']
    impacts.append(float(impact))
    
    azimuths = L_nonassociated['azimuths']
    azimuths.append(float(azimuth))
    
    PAs = L_nonassociated['PAs']
    PAs.append(float(PA))
    
    incs = L_nonassociated['incs']
    incs.append(float(inc))
    
    adjustedIncs = L_nonassociated['adjustedIncs']
    adjustedIncs.append(float(adjustedInc))
    
    ls = L_nonassociated['ls']
    ls.append(float(l))

    l_cuss = L_nonassociated['l_cuss']
    l_cuss.append(float(l_cus))

    R_virs = L_nonassociated['R_virs']
    R_virs.append(float(R_vir))
    
    cuss = L_nonassociated['cuss']
    cuss.append(float(cus))

    MajDiams = L_nonassociated['MajDiams']
    MajDiams.append(float(MajDiam))
    
    MTypes = L_nonassociated['MTypes']
    MTypes.append(MType)
    
    Vhels = L_nonassociated['Vhels']
    Vhels.append(float(Vhel))
    
    vcorrs = L_nonassociated['vcorrs']
    vcorrs.append(float(vcorr))
    
    bestDists = L_nonassociated['bestDists']
    bestDists.append(float(bestDist))
    
    e_bestDists = L_nonassociated['e_bestDists']
    e_bestDists.append(float(e_bestDist))
    
    group_nums = L_nonassociated['group_nums']
    group_nums.append(group_num)
    
    group_mems = L_nonassociated['group_mems']
    group_mems.append(group_mem)

    group_dists = L_nonassociated['group_dists']
    group_dists.append(float(group_dist))
    
    Lstar_meds = L_nonassociated['Lstar_meds']
    Lstar_meds.append(float(Lstar_med))
    
    e_Lstar_meds = L_nonassociated['e_Lstar_meds']
    e_Lstar_meds.append(float(e_Lstar_med))
    
    Lstar_sdsss = L_nonassociated['Lstar_sdsss']
    Lstar_sdsss.append(float(Lstar_sdss))
    
    e_Lstar_sdsss = L_nonassociated['e_Lstar_sdsss']
    e_Lstar_sdsss.append(float(e_Lstar_sdss))
    
    Bmags = L_nonassociated['Bmags']
    Bmags.append(float(Bmag))
    
    Bmag_sdsss = L_nonassociated['Bmag_sdsss']
    Bmag_sdsss.append(float(Bmag_sdss))
    
    others = L_nonassociated['others']
    others.append(other)

    pass


def add_to_L_two(Lya_v, e_Lya_v, Lya_W, e_Lya_W, Na, e_Na, b, e_b, W, e_W, target, z_target,\
    RA_target, Dec_target, Name, RAdeg, DEdeg, impact, azimuth, PA, inc,\
    adjustedInc, l, l_cus, R_vir, cus, MajDiam, MType, Vhel, vcorr, bestDist, e_bestDist,\
    group_num, group_mem, group_dist, Lstar_med, e_Lstar_med, Lstar_sdss, e_Lstar_sdss,\
    Bmag, Bmag_sdss, other):

 
    Lya_vs = L_two['Lya_vs']
    Lya_vs.append(float(Lya_v))

    e_Lya_vs = L_two['e_Lya_vs']
    e_Lya_vs.append(float(e_Lya_v))

    Lya_Ws = L_two['Lya_Ws']
    Lya_Ws.append(float(Lya_W))
    
    e_Lya_Ws = L_two['e_Lya_Ws']
    e_Lya_Ws.append(float(e_Lya_W))

    Nas = L_two['Nas']
    Nas.append(float(Na))
    
    e_Nas = L_two['e_Nas']
    e_Nas.append(float(e_Na))
    
    bs = L_two['bs']
    bs.append(float(b))

    e_bs = L_two['e_bs']
    e_bs.append(float(e_b))

    Ws = L_two['Ws']
    Ws.append(float(W))
    
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
    
    Names = L_two['Names']
    Names.append(Name)

    RA_galaxies = L_two['RA_galaxies']
    RA_galaxies.append(float(RAdeg))
    
    Dec_galaxies = L_two['Dec_galaxies']
    Dec_galaxies.append(float(DEdeg))
    
    impacts = L_two['impacts']
    impacts.append(float(impact))
    
    azimuths = L_two['azimuths']
    azimuths.append(float(azimuth))
    
    PAs = L_two['PAs']
    PAs.append(float(PA))
    
    incs = L_two['incs']
    incs.append(float(inc))
    
    adjustedIncs = L_two['adjustedIncs']
    adjustedIncs.append(float(adjustedInc))
    
    ls = L_two['ls']
    ls.append(float(l))

    l_cuss = L_two['l_cuss']
    l_cuss.append(float(l_cus))

    R_virs = L_two['R_virs']
    R_virs.append(float(R_vir))
    
    cuss = L_two['cuss']
    cuss.append(float(cus))

    MajDiams = L_two['MajDiams']
    MajDiams.append(float(MajDiam))
    
    MTypes = L_two['MTypes']
    MTypes.append(MType)
    
    Vhels = L_two['Vhels']
    Vhels.append(float(Vhel))
    
    vcorrs = L_two['vcorrs']
    vcorrs.append(float(vcorr))
    
    bestDists = L_two['bestDists']
    bestDists.append(float(bestDist))
    
    e_bestDists = L_two['e_bestDists']
    e_bestDists.append(float(e_bestDist))
    
    group_nums = L_two['group_nums']
    group_nums.append(group_num)
    
    group_mems = L_two['group_mems']
    group_mems.append(group_mem)

    group_dists = L_two['group_dists']
    group_dists.append(float(group_dist))
    
    Lstar_meds = L_two['Lstar_meds']
    Lstar_meds.append(float(Lstar_med))
    
    e_Lstar_meds = L_two['e_Lstar_meds']
    e_Lstar_meds.append(float(e_Lstar_med))
    
    Lstar_sdsss = L_two['Lstar_sdsss']
    Lstar_sdsss.append(float(Lstar_sdss))
    
    e_Lstar_sdsss = L_two['e_Lstar_sdsss']
    e_Lstar_sdsss.append(float(e_Lstar_sdss))
    
    Bmags = L_two['Bmags']
    Bmags.append(float(Bmag))
    
    Bmag_sdsss = L_two['Bmag_sdsss']
    Bmag_sdsss.append(float(Bmag_sdss))
    
    others = L_two['others']
    others.append(other)
    
    pass

    
    

def add_to_L_three_plus(Lya_v, e_Lya_v, Lya_W, e_Lya_W, Na, e_Na, b, e_b, W, e_W, target, z_target,\
    RA_target, Dec_target, Name, RAdeg, DEdeg, impact, azimuth, PA, inc,\
    adjustedInc, l, l_cus, R_vir, cus, MajDiam, MType, Vhel, vcorr, bestDist, e_bestDist,\
    group_num, group_mem, group_dist, Lstar_med, e_Lstar_med, Lstar_sdss, e_Lstar_sdss, \
    Bmag, Bmag_sdss, other):

 
    Lya_vs = L_three_plus['Lya_vs']
    Lya_vs.append(float(Lya_v))

    e_Lya_vs = L_three_plus['e_Lya_vs']
    e_Lya_vs.append(float(e_Lya_v))

    Lya_Ws = L_three_plus['Lya_Ws']
    Lya_Ws.append(float(Lya_W))
    
    e_Lya_Ws = L_three_plus['e_Lya_Ws']
    e_Lya_Ws.append(float(e_Lya_W))

    Nas = L_three_plus['Nas']
    Nas.append(float(Na))
    
    e_Nas = L_three_plus['e_Nas']
    e_Nas.append(float(e_Na))
    
    bs = L_three_plus['bs']
    bs.append(float(b))

    e_bs = L_three_plus['e_bs']
    e_bs.append(float(e_b))

    Ws = L_three_plus['Ws']
    Ws.append(float(W))
    
    e_Ws = L_three_plus['e_Ws']
    e_Ws.append(float(e_W))

    targets = L_three_plus['targets']
    targets.append(target)

    z_targets = L_three_plus['z_targets']
    z_targets.append(float(z_target))

    RA_targets = L_three_plus['RA_targets']
    RA_targets.append(float(RA_target))
    
    Dec_targets = L_three_plus['Dec_targets']
    Dec_targets.append(float(Dec_target))
    
    Names = L_three_plus['Names']
    Names.append(Name)

    RA_galaxies = L_three_plus['RA_galaxies']
    RA_galaxies.append(float(RAdeg))
    
    Dec_galaxies = L_three_plus['Dec_galaxies']
    Dec_galaxies.append(float(DEdeg))
    
    impacts = L_three_plus['impacts']
    impacts.append(float(impact))
    
    azimuths = L_three_plus['azimuths']
    azimuths.append(float(azimuth))
    
    PAs = L_three_plus['PAs']
    PAs.append(float(PA))
    
    incs = L_three_plus['incs']
    incs.append(float(inc))
    
    adjustedIncs = L_three_plus['adjustedIncs']
    adjustedIncs.append(float(adjustedInc))
    
    ls = L_three_plus['ls']
    ls.append(float(l))

    l_cuss = L_three_plus['l_cuss']
    l_cuss.append(float(l_cus))

    R_virs = L_three_plus['R_virs']
    R_virs.append(float(R_vir))
    
    cuss = L_three_plus['cuss']
    cuss.append(float(cus))

    MajDiams = L_three_plus['MajDiams']
    MajDiams.append(float(MajDiam))
    
    MTypes = L_three_plus['MTypes']
    MTypes.append(MType)
    
    Vhels = L_three_plus['Vhels']
    Vhels.append(float(Vhel))
    
    vcorrs = L_three_plus['vcorrs']
    vcorrs.append(float(vcorr))
    
    bestDists = L_three_plus['bestDists']
    bestDists.append(float(bestDist))
    
    e_bestDists = L_three_plus['e_bestDists']
    e_bestDists.append(float(e_bestDist))
    
    group_nums = L_three_plus['group_nums']
    group_nums.append(group_num)
    
    group_mems = L_three_plus['group_mems']
    group_mems.append(group_mem)

    group_dists = L_three_plus['group_dists']
    group_dists.append(float(group_dist))
    
    Lstar_meds = L_three_plus['Lstar_meds']
    Lstar_meds.append(float(Lstar_med))
    
    e_Lstar_meds = L_three_plus['e_Lstar_meds']
    e_Lstar_meds.append(float(e_Lstar_med))
    
    Lstar_sdsss = L_three_plus['Lstar_sdsss']
    Lstar_sdsss.append(float(Lstar_sdss))
    
    e_Lstar_sdsss = L_three_plus['e_Lstar_sdsss']
    e_Lstar_sdsss.append(float(e_Lstar_sdss))
    
    Bmags = L_three_plus['Bmags']
    Bmags.append(float(Bmag))
    
    Bmag_sdsss = L_three_plus['Bmag_sdsss']
    Bmag_sdsss.append(float(Bmag_sdss))
    
    others = L_three_plus['others']
    others.append(other)
    
    pass




def add_to_L_group(Lya_v, e_Lya_v, Lya_W, e_Lya_W, Na, e_Na, b, e_b, W, e_W, target, z_target,\
    RA_target, Dec_target, Name, RAdeg, DEdeg, impact, azimuth, PA, inc,\
    adjustedInc, l, l_cus, R_vir, cus, MajDiam, MType, Vhel, vcorr, bestDist, e_bestDist,\
    group_num, group_mem, group_dist, Lstar_med, e_Lstar_med, Lstar_sdss, e_Lstar_sdss, \
    Bmag, Bmag_sdss):

 
    Lya_vs = L_group['Lya_vs']
    Lya_vs.append(float(Lya_v))

    e_Lya_vs = L_group['e_Lya_vs']
    e_Lya_vs.append(float(e_Lya_v))

    Lya_Ws = L_group['Lya_Ws']
    Lya_Ws.append(float(Lya_W))
    
    e_Lya_Ws = L_group['e_Lya_Ws']
    e_Lya_Ws.append(float(e_Lya_W))

    Nas = L_group['Nas']
    Nas.append(float(Na))
    
    e_Nas = L_group['e_Nas']
    e_Nas.append(float(e_Na))
    
    bs = L_group['bs']
    bs.append(float(b))

    e_bs = L_group['e_bs']
    e_bs.append(float(e_b))

    Ws = L_group['Ws']
    Ws.append(float(W))
    
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
    
    Names = L_group['Names']
    Names.append(Name)

    RA_galaxies = L_group['RA_galaxies']
    RA_galaxies.append(float(RAdeg))
    
    Dec_galaxies = L_group['Dec_galaxies']
    Dec_galaxies.append(float(DEdeg))
    
    impacts = L_group['impacts']
    impacts.append(float(impact))
    
    azimuths = L_group['azimuths']
    azimuths.append(float(azimuth))
    
    PAs = L_group['PAs']
    PAs.append(float(PA))
    
    incs = L_group['incs']
    incs.append(float(inc))
    
    adjustedIncs = L_group['adjustedIncs']
    adjustedIncs.append(float(adjustedInc))
    
    ls = L_group['ls']
    ls.append(float(l))

    l_cuss = L_group['l_cuss']
    l_cuss.append(float(l_cus))

    R_virs = L_group['R_virs']
    R_virs.append(float(R_vir))
    
    cuss = L_group['cuss']
    cuss.append(float(cus))

    MajDiams = L_group['MajDiams']
    MajDiams.append(float(MajDiam))
    
    MTypes = L_group['MTypes']
    MTypes.append(MType)
    
    Vhels = L_group['Vhels']
    Vhels.append(float(Vhel))
    
    vcorrs = L_group['vcorrs']
    vcorrs.append(float(vcorr))
    
    bestDists = L_group['bestDists']
    bestDists.append(float(bestDist))
    
    e_bestDists = L_group['e_bestDists']
    e_bestDists.append(float(e_bestDist))
    
    group_nums = L_group['group_nums']
    group_nums.append(group_num)
    
    group_mems = L_group['group_mems']
    group_mems.append(group_mem)

    group_dists = L_group['group_dists']
    group_dists.append(float(group_dist))
    
    Lstar_meds = L_group['Lstar_meds']
    Lstar_meds.append(float(Lstar_med))
    
    e_Lstar_meds = L_group['e_Lstar_meds']
    e_Lstar_meds.append(float(e_Lstar_med))
    
    Lstar_sdsss = L_group['Lstar_sdsss']
    Lstar_sdsss.append(float(Lstar_sdss))
    
    e_Lstar_sdsss = L_group['e_Lstar_sdsss']
    e_Lstar_sdsss.append(float(e_Lstar_sdss))
    
    Bmags = L_group['Bmags']
    Bmags.append(float(Bmag))
    
    Bmag_sdsss = L_group['Bmag_sdsss']
    Bmag_sdsss.append(float(Bmag_sdss))
    
    pass




def add_to_L_summed(Lya_v, e_Lya_v, Lya_W, e_Lya_W, Na, e_Na, b, e_b, W, e_W, target, z_target,\
    RA_target, Dec_target, Name, RAdeg, DEdeg, impact, azimuth, PA, inc,\
    adjustedInc, l, l_cus, R_vir, cus, MajDiam, MType, Vhel, vcorr, bestDist, e_bestDist,\
    group_num, group_mem, group_dist, Lstar_med, e_Lstar_med, Lstar_sdss, e_Lstar_sdss,\
    Bmag, Bmag_sdss, summed_l, summed_l_cus):

 
    Lya_vs = L_summed['Lya_vs']
    Lya_vs.append(float(Lya_v))

    e_Lya_vs = L_summed['e_Lya_vs']
    e_Lya_vs.append(float(e_Lya_v))

    Lya_Ws = L_summed['Lya_Ws']
    Lya_Ws.append(float(Lya_W))
    
    e_Lya_Ws = L_summed['e_Lya_Ws']
    e_Lya_Ws.append(float(e_Lya_W))

    Nas = L_summed['Nas']
    Nas.append(float(Na))
    
    e_Nas = L_summed['e_Nas']
    e_Nas.append(float(e_Na))
    
    bs = L_summed['bs']
    bs.append(float(b))

    e_bs = L_summed['e_bs']
    e_bs.append(float(e_b))

    Ws = L_summed['Ws']
    Ws.append(float(W))
    
    e_Ws = L_summed['e_Ws']
    e_Ws.append(float(e_W))

    targets = L_summed['targets']
    targets.append(target)

    z_targets = L_summed['z_targets']
    z_targets.append(float(z_target))

    RA_targets = L_summed['RA_targets']
    RA_targets.append(float(RA_target))
    
    Dec_targets = L_summed['Dec_targets']
    Dec_targets.append(float(Dec_target))
    
    Names = L_summed['Names']
    Names.append(Name)

    RA_galaxies = L_summed['RA_galaxies']
    RA_galaxies.append(float(RAdeg))
    
    Dec_galaxies = L_summed['Dec_galaxies']
    Dec_galaxies.append(float(DEdeg))
    
    impacts = L_summed['impacts']
    impacts.append(float(impact))
    
    azimuths = L_summed['azimuths']
    azimuths.append(float(azimuth))
    
    PAs = L_summed['PAs']
    PAs.append(float(PA))
    
    incs = L_summed['incs']
    incs.append(float(inc))
    
    adjustedIncs = L_summed['adjustedIncs']
    adjustedIncs.append(float(adjustedInc))
    
    ls = L_summed['ls']
    ls.append(float(l))

    l_cuss = L_summed['l_cuss']
    l_cuss.append(float(l_cus))

    R_virs = L_summed['R_virs']
    R_virs.append(float(R_vir))
    
    cuss = L_summed['cuss']
    cuss.append(float(cus))

    MajDiams = L_summed['MajDiams']
    MajDiams.append(float(MajDiam))
    
    MTypes = L_summed['MTypes']
    MTypes.append(MType)
    
    Vhels = L_summed['Vhels']
    Vhels.append(float(Vhel))
    
    vcorrs = L_summed['vcorrs']
    vcorrs.append(float(vcorr))
    
    bestDists = L_summed['bestDists']
    bestDists.append(float(bestDist))
    
    e_bestDists = L_summed['e_bestDists']
    e_bestDists.append(float(e_bestDist))
    
    group_nums = L_summed['group_nums']
    group_nums.append(group_num)
    
    group_mems = L_summed['group_mems']
    group_mems.append(group_mem)

    group_dists = L_summed['group_dists']
    group_dists.append(float(group_dist))
    
    Lstar_meds = L_summed['Lstar_meds']
    Lstar_meds.append(float(Lstar_med))
    
    e_Lstar_meds = L_summed['e_Lstar_meds']
    e_Lstar_meds.append(float(e_Lstar_med))
    
    Lstar_sdsss = L_summed['Lstar_sdsss']
    Lstar_sdsss.append(float(Lstar_sdss))
    
    e_Lstar_sdsss = L_summed['e_Lstar_sdsss']
    e_Lstar_sdsss.append(float(e_Lstar_sdss))
    
    Bmags = L_summed['Bmags']
    Bmags.append(float(Bmag))
    
    Bmag_sdsss = L_summed['Bmag_sdsss']
    Bmag_sdsss.append(float(Bmag_sdss))
    
    summed_ls = L_summed['summed_ls']
    summed_ls.append(summed_l)
    
    summed_l_cuss = L_summed['summed_l_cuss']
    summed_l_cuss.append(summed_l_cus)
    
    pass
    
    
    
    
def main():
    # correlation options
    maxSep = 500.
    agnSeparation = 4000.
    minVcorr = 450.
    minSize = 0.
    max_deltav = 400.
    min_likelihood = 0.005
    rigor = 5
    v_norm = 150.
    l_norm = 1.
    
    # sort based on likelihood cus instead of the regular one?
    use_likelihood_cus = False
    
    # double l if impact <= 1 R_vir?
    double_l_within_rvir = False
    

    # assuming 'theFile' contains one name per line, read the file
    
    if getpass.getuser() == 'frenchd':
#         filename = '/Users/frenchd/Research/inclination/git_inclination/maps/LG_correlation_combined5_14_edit.csv'
#         pickleFilename = '/Users/frenchd/Research/inclination/git_inclination/picklePilot_plusSALT_14.p'

#         targetlist_filename = '/Users/frenchd/Research/correlation/TARGETLIST_10_17_17_TOTAL.csv'
        targetlist_filename = '/Users/frenchd/Research/correlation/TARGETLIST_06_06_18_TOTAL.csv'

#         filename = '/Users/frenchd/Research/inclination/git_inclination/correlatedTargetList_5_29_18_measurements.csv'
#         filename = '/Users/frenchd/Research/inclination/git_inclination/targets/correlatedTargetList_5_29_18_measurements_copy_short.csv'
#         filename = '/Users/frenchd/Research/inclination/git_inclination/targets/correlatedTargetList_5_29_18_measurements_copy.csv'
        filename = '/Users/frenchd/Research/inclination/git_inclination/targets/correlatedTargetList_5_29_18_measurements_copy.csv'

        # pickle files
        all_filename = '/Users/frenchd/Research/inclination/git_inclination/all7_min005_v150.p'
        isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/isolated7_min005_v150.p'
        L_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_isolated7_min005.p'
        L_associated_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated7_min005_v150.p'
        L_associated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated7_min005_v150.p'
        L_nonassociated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated7_min005_v150.p'
        L_two_filename = '/Users/frenchd/Research/inclination/git_inclination/L_two7_min005_v150.p'
        L_three_plus_filename = '/Users/frenchd/Research/inclination/git_inclination/L_three_plus7_min005_v150.p'
        L_group_filename = '/Users/frenchd/Research/inclination/git_inclination/L_group7_min005_v150.p'
        L_summed_filename = '/Users/frenchd/Research/inclination/git_inclination/L_summed7_min005_v150.p'


    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    
    theFile = open(filename,'rU')    
    reader = csv.DictReader(theFile)
    
    # open the pickle files
    all_file = open(all_filename,'wt')
    isolated_file = open(isolated_filename,'wt')
    L_isolated_file = open(L_isolated_filename,'wt')
    L_associated_isolated_file = open(L_associated_isolated_filename,'wt')
    L_associated_file = open(L_associated_filename,'wt')
    L_nonassociated_file = open(L_nonassociated_filename,'wt')
    L_two_file = open(L_two_filename,'wt')
    L_three_plus_file = open(L_three_plus_filename,'wt')
    L_group_file = open(L_group_filename,'wt')
    L_summed_file = open(L_summed_filename,'wt')


    # all: add empty lists
    all['Lya_vs'] = []
    all['e_Lya_vs'] = []
    all['Lya_Ws'] = []
    all['e_Lya_Ws'] = []
    all['Nas'] = []
    all['e_Nas'] = []
    all['bs'] = []
    all['e_bs'] = []
    all['Ws'] = []
    all['e_Ws'] = []
    all['targets'] = []
    all['z_targets'] = []
    all['RA_targets'] = []
    all['Dec_targets'] = []
  
  
    
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
    isolated['z_targets'] = []
    isolated['RA_targets'] = []
    isolated['Dec_targets'] = []

    
    
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
    L_isolated['z_targets'] = []
    L_isolated['RA_targets'] = []
    L_isolated['Dec_targets'] = []
    
    
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
    L_associated_isolated['Names'] = []
    L_associated_isolated['RA_galaxies'] = []
    L_associated_isolated['Dec_galaxies'] = []
    L_associated_isolated['impacts'] = []
    L_associated_isolated['azimuths'] = []
    L_associated_isolated['PAs'] = []
    L_associated_isolated['incs'] = []
    L_associated_isolated['adjustedIncs'] = []
    L_associated_isolated['ls'] = []
    L_associated_isolated['l_cuss'] = []
    L_associated_isolated['R_virs'] = []
    L_associated_isolated['cuss'] = []
    L_associated_isolated['MajDiams'] = []
    L_associated_isolated['MTypes'] = []
    L_associated_isolated['Vhels'] = []
    L_associated_isolated['vcorrs'] = []
    L_associated_isolated['bestDists'] = []
    L_associated_isolated['e_bestDists'] = []
    L_associated_isolated['group_nums'] = []
    L_associated_isolated['group_mems'] = []
    L_associated_isolated['group_dists'] = []
    L_associated_isolated['Lstar_meds'] = []
    L_associated_isolated['e_Lstar_meds'] = []
    L_associated_isolated['Lstar_sdsss'] = []
    L_associated_isolated['e_Lstar_sdsss'] = []
    L_associated_isolated['Bmags'] = []
    L_associated_isolated['Bmag_sdsss'] = []

    
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
    L_associated['Names'] = []
    L_associated['RA_galaxies'] = []
    L_associated['Dec_galaxies'] = []
    L_associated['impacts'] = []
    L_associated['azimuths'] = []
    L_associated['PAs'] = []
    L_associated['incs'] = []
    L_associated['adjustedIncs'] = []
    L_associated['ls'] = []
    L_associated['l_cuss'] = []
    L_associated['R_virs'] = []
    L_associated['cuss'] = []
    L_associated['MajDiams'] = []
    L_associated['MTypes'] = []
    L_associated['Vhels'] = []
    L_associated['vcorrs'] = []
    L_associated['bestDists'] = []
    L_associated['e_bestDists'] = []
    L_associated['group_nums'] = []
    L_associated['group_mems'] = []
    L_associated['group_dists'] = []
    L_associated['Lstar_meds'] = []
    L_associated['e_Lstar_meds'] = []
    L_associated['Lstar_sdsss'] = []
    L_associated['e_Lstar_sdsss'] = []
    L_associated['Bmags'] = []
    L_associated['Bmag_sdsss'] = []

    
    
    # L_nonassociated: add empty lists
    L_nonassociated['Lya_vs'] = []
    L_nonassociated['e_Lya_vs'] = []
    L_nonassociated['Lya_Ws'] = []
    L_nonassociated['e_Lya_Ws'] = []
    L_nonassociated['Nas'] = []
    L_nonassociated['e_Nas'] = []
    L_nonassociated['bs'] = []
    L_nonassociated['e_bs'] = []
    L_nonassociated['Ws'] = []
    L_nonassociated['e_Ws'] = []
    L_nonassociated['targets'] = []
    L_nonassociated['z_targets'] = []
    L_nonassociated['RA_targets'] = []
    L_nonassociated['Dec_targets'] = []
    L_nonassociated['Names'] = []
    L_nonassociated['RA_galaxies'] = []
    L_nonassociated['Dec_galaxies'] = []
    L_nonassociated['impacts'] = []
    L_nonassociated['azimuths'] = []
    L_nonassociated['PAs'] = []
    L_nonassociated['incs'] = []
    L_nonassociated['adjustedIncs'] = []
    L_nonassociated['ls'] = []
    L_nonassociated['l_cuss'] = []
    L_nonassociated['R_virs'] = []
    L_nonassociated['cuss'] = []
    L_nonassociated['MajDiams'] = []
    L_nonassociated['MTypes'] = []
    L_nonassociated['Vhels'] = []
    L_nonassociated['vcorrs'] = []
    L_nonassociated['bestDists'] = []
    L_nonassociated['e_bestDists'] = []
    L_nonassociated['group_nums'] = []
    L_nonassociated['group_mems'] = []
    L_nonassociated['group_dists'] = []
    L_nonassociated['Lstar_meds'] = []
    L_nonassociated['e_Lstar_meds'] = []
    L_nonassociated['Lstar_sdsss'] = []
    L_nonassociated['e_Lstar_sdsss'] = []
    L_nonassociated['Bmags'] = []
    L_nonassociated['Bmag_sdsss'] = []
    L_nonassociated['others'] = []

    
    
    # L_two: add empty lists
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
    L_two['Names'] = []
    L_two['RA_galaxies'] = []
    L_two['Dec_galaxies'] = []
    L_two['impacts'] = []
    L_two['azimuths'] = []
    L_two['PAs'] = []
    L_two['incs'] = []
    L_two['adjustedIncs'] = []
    L_two['ls'] = []
    L_two['l_cuss'] = []
    L_two['R_virs'] = []
    L_two['cuss'] = []
    L_two['MajDiams'] = []
    L_two['MTypes'] = []
    L_two['Vhels'] = []
    L_two['vcorrs'] = []
    L_two['bestDists'] = []
    L_two['e_bestDists'] = []
    L_two['group_nums'] = []
    L_two['group_mems'] = []
    L_two['group_dists'] = []
    L_two['Lstar_meds'] = []
    L_two['e_Lstar_meds'] = []
    L_two['Lstar_sdsss'] = []
    L_two['e_Lstar_sdsss'] = []
    L_two['Bmags'] = []
    L_two['Bmag_sdsss'] = []
    L_two['others'] = []

    
    # L_three_plus: add empty lists
    L_three_plus['Lya_vs'] = []
    L_three_plus['e_Lya_vs'] = []
    L_three_plus['Lya_Ws'] = []
    L_three_plus['e_Lya_Ws'] = []
    L_three_plus['Nas'] = []
    L_three_plus['e_Nas'] = []
    L_three_plus['bs'] = []
    L_three_plus['e_bs'] = []
    L_three_plus['Ws'] = []
    L_three_plus['e_Ws'] = []
    L_three_plus['targets'] = []
    L_three_plus['z_targets'] = []
    L_three_plus['RA_targets'] = []
    L_three_plus['Dec_targets'] = []
    L_three_plus['Names'] = []
    L_three_plus['RA_galaxies'] = []
    L_three_plus['Dec_galaxies'] = []
    L_three_plus['impacts'] = []
    L_three_plus['azimuths'] = []
    L_three_plus['PAs'] = []
    L_three_plus['incs'] = []
    L_three_plus['adjustedIncs'] = []
    L_three_plus['ls'] = []
    L_three_plus['l_cuss'] = []
    L_three_plus['R_virs'] = []
    L_three_plus['cuss'] = []
    L_three_plus['MajDiams'] = []
    L_three_plus['MTypes'] = []
    L_three_plus['Vhels'] = []
    L_three_plus['vcorrs'] = []
    L_three_plus['bestDists'] = []
    L_three_plus['e_bestDists'] = []
    L_three_plus['group_nums'] = []
    L_three_plus['group_mems'] = []
    L_three_plus['group_dists'] = []
    L_three_plus['Lstar_meds'] = []
    L_three_plus['e_Lstar_meds'] = []
    L_three_plus['Lstar_sdsss'] = []
    L_three_plus['e_Lstar_sdsss'] = []
    L_three_plus['Bmags'] = []
    L_three_plus['Bmag_sdsss'] = []
    L_three_plus['others'] = []


    # L_group: add empty lists
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
    L_group['Names'] = []
    L_group['RA_galaxies'] = []
    L_group['Dec_galaxies'] = []
    L_group['impacts'] = []
    L_group['azimuths'] = []
    L_group['PAs'] = []
    L_group['incs'] = []
    L_group['adjustedIncs'] = []
    L_group['ls'] = []
    L_group['l_cuss'] = []
    L_group['R_virs'] = []
    L_group['cuss'] = []
    L_group['MajDiams'] = []
    L_group['MTypes'] = []
    L_group['Vhels'] = []
    L_group['vcorrs'] = []
    L_group['bestDists'] = []
    L_group['e_bestDists'] = []
    L_group['group_nums'] = []
    L_group['group_mems'] = []
    L_group['group_dists'] = []
    L_group['Lstar_meds'] = []
    L_group['e_Lstar_meds'] = []
    L_group['Lstar_sdsss'] = []
    L_group['e_Lstar_sdsss'] = []
    L_group['Bmags'] = []
    L_group['Bmag_sdsss'] = []


    # L_summed: add empty lists
    L_summed['Lya_vs'] = []
    L_summed['e_Lya_vs'] = []
    L_summed['Lya_Ws'] = []
    L_summed['e_Lya_Ws'] = []
    L_summed['Nas'] = []
    L_summed['e_Nas'] = []
    L_summed['bs'] = []
    L_summed['e_bs'] = []
    L_summed['Ws'] = []
    L_summed['e_Ws'] = []
    L_summed['targets'] = []
    L_summed['z_targets'] = []
    L_summed['RA_targets'] = []
    L_summed['Dec_targets'] = []
    L_summed['Names'] = []
    L_summed['RA_galaxies'] = []
    L_summed['Dec_galaxies'] = []
    L_summed['impacts'] = []
    L_summed['azimuths'] = []
    L_summed['PAs'] = []
    L_summed['incs'] = []
    L_summed['adjustedIncs'] = []
    L_summed['ls'] = []
    L_summed['l_cuss'] = []
    L_summed['R_virs'] = []
    L_summed['cuss'] = []
    L_summed['MajDiams'] = []
    L_summed['MTypes'] = []
    L_summed['Vhels'] = []
    L_summed['vcorrs'] = []
    L_summed['bestDists'] = []
    L_summed['e_bestDists'] = []
    L_summed['group_nums'] = []
    L_summed['group_mems'] = []
    L_summed['group_dists'] = []
    L_summed['Lstar_meds'] = []
    L_summed['e_Lstar_meds'] = []
    L_summed['Lstar_sdsss'] = []
    L_summed['e_Lstar_sdsss'] = []
    L_summed['Bmags'] = []
    L_summed['Bmag_sdsss'] = []
    L_summed['summed_ls'] = []
    L_summed['summed_l_cuss'] = []


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
    # now the full data set
    total = 1506
    counter = 0
    stopCount = 5000
    print
    print 'Starting loop!'
    print
    for i in reader:
        counter +=1
        sys.stdout.write("\r Percent Complete: {0} / {1}".format(counter, total))
        sys.stdout.flush()
        
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
        
        proceed = False
        if isNumber(Lya_v):
            if bfind(identified.lower(), 'lya') or bfind(identified.lower(), 'maybe'):
                proceed = True
                
                if not isNumber(e_Lya_v):
                    e_Lya_v = 1
                
        try:
            RA_target = target_coords[target]['RA']
            Dec_target = target_coords[target]['Dec']
        except Exception, e:
            proceed = False
            print
            print "Target not found: ",target
            
        if counter >= stopCount:
            proceed = False

        if proceed:
            # add this to the all set
            add_to_all(Lya_v,\
                        e_Lya_v,\
                        Lya_W,\
                        e_Lya_W,\
                        Na,\
                        e_Na,\
                        b,\
                        e_b,\
                        W,\
                        e_W,\
                        target,\
                        z_target,\
                        RA_target,\
                        Dec_target)
            
            # do the correlation
            correlation = correlateSingle.correlateTarget(target, maxSep, agnSeparation, minVcorr, minSize, slow=False, searchAll=False, returnDict=True)

            # if no galaxies are returned, add it to the isolated list
            if len(correlation) == 0:
                add_to_isolated(Lya_v, e_Lya_v, Lya_W, e_Lya_W, Na, e_Na, b, e_b, W, e_W, target, z_target, RA_target, Dec_target)
            
            else:
                candidates = []
                candidates_cus = []
                others = []
            
                summed_likelihoods = []
                summed_likelihood_cuss = []
                
                
                for c in correlation:
                    # add line stuff to the galaxy dictionaries
                    correlation[c]['identified'] = identified
                    correlation[c]['Lya_v'] = Lya_v
                    correlation[c]['e_Lya_v'] = e_Lya_v
                    correlation[c]['Lya_W'] = Lya_W
                    correlation[c]['e_Lya_W'] = e_Lya_W
                    correlation[c]['Na'] = Na
                    correlation[c]['e_Na'] = e_Na
                    correlation[c]['b'] = b
                    correlation[c]['e_b'] = e_b
                    correlation[c]['W'] = W
                    correlation[c]['e_W'] = e_W

                    
                    # unpack everything
                    target = correlation[c]['target']
                    z_target = correlation[c]['z_target']
                    RA_target = correlation[c]['RA_target']
                    Dec_target = correlation[c]['Dec_target']
                    Name = correlation[c]['Name']
                    RAdeg = correlation[c]['RAdeg']
                    DEdeg = correlation[c]['DEdeg']
                    impact = correlation[c]['impact']
                    azimuths = correlation[c]['azimuth']
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

                    Vhel = float(Vhel)
                    Lya_v = float(Lya_v)
                    
                
                    # velocity difference = V_absorber - V_sys
                    vel_dif = Lya_v - Vhel
    
                    # check if the line is too far from a galaxy
                    if abs(Lya_v - Vhel) <= max_deltav:                
                        # if there's a galaxy within max_deltav, calculate likelihood values
                
                        if isNumber(MajDiam) and MajDiam != -99.99 and MajDiam != '-99.99':
                            # "sphere of influence" value for likelihood_custom
                            MajDiam = float(MajDiam)
                            impact = float(impact)
                            R_vir = float(R_vir)
                            
                            cus = MajDiam**1.5

                            # first for the virial radius
                            likelihood = l_norm * math.exp(-(impact/R_vir)**2) * math.exp(-(vel_dif/v_norm)**2)
                        
                            # then for the second 'virial like' m15 radius
                            likelihood_cus = l_norm * math.exp(-(impact/cus)**2) * math.exp(-(vel_dif/v_norm)**2)
                
                            # sort based on which likelihood metric?
                            if use_likelihood_cus:
                                l_used = likelihood_cus
                            else:
                                l_used = likelihood
                            
                            # multiply the likelihood by two if within 1 R_vir?
                            if double_l_within_rvir:
                                if impact <= R_vir:
                                    l_used = l_used * 2

                            # add likelihoods and stuff to the dictionary entries
                            correlation[c]['l'] = likelihood
                            correlation[c]['l_cus'] = likelihood_cus
                            correlation[c]['cus'] = cus
                    
                            # if they make the cut, add them to the candidate lists
                            if l_used >= min_likelihood:
                                candidates.append([l_used,correlation[c]])
                        
                            else:
                                others.append([l_used,correlation[c]])
                                
                            
                            # add to the summed list:
                            summed_likelihoods.append(likelihood)
                            summed_likelihood_cuss.append(likelihood_cus)
                                
                                
                        else:
                            likelihood = -99.99
                            others.append([likelihood, correlation[c]])
                
                
                
                # now candidates have been added to that list and isolated have been taken care of.
                # first add to the summed list:
                add_to_L_summed(correlation[c]['Lya_v'],\
                                correlation[c]['e_Lya_v'],\
                                correlation[c]['Lya_W'],\
                                correlation[c]['e_Lya_W'],\
                                correlation[c]['Na'],\
                                correlation[c]['e_Na'],\
                                correlation[c]['b'],\
                                correlation[c]['e_b'],\
                                correlation[c]['W'],\
                                correlation[c]['e_W'],\
                                correlation[c]['target'],\
                                correlation[c]['z_target'],\
                                correlation[c]['RA_target'],\
                                correlation[c]['Dec_target'],\
                                correlation[c]['Name'],\
                                correlation[c]['RAdeg'],\
                                correlation[c]['DEdeg'],\
                                correlation[c]['impact'],\
                                correlation[c]['azimuth'],\
                                correlation[c]['PA'],\
                                correlation[c]['inc'],\
                                correlation[c]['adjustedInc'],\
                                -99,\
                                -99,\
                                correlation[c]['R_vir'],\
                                -99,\
                                correlation[c]['MajDiam'],\
                                correlation[c]['MType'],\
                                correlation[c]['Vhel'],\
                                correlation[c]['vcorr'],\
                                correlation[c]['bestDist'],\
                                correlation[c]['e_bestDist'],\
                                correlation[c]['group_num'],\
                                correlation[c]['group_mem'],\
                                correlation[c]['group_dist'],\
                                correlation[c]['Lstar_med'],\
                                correlation[c]['e_Lstar_med'],\
                                correlation[c]['Lstar_sdss'],\
                                correlation[c]['e_Lstar_sdss'],\
                                correlation[c]['Bmag'],\
                                correlation[c]['Bmag_sdss'],\
                                summed_likelihoods,\
                                summed_likelihood_cuss)
                
                
                # next handle the candidate list:
                
                # first sort by likelihood:
                candidates.sort(reverse=True)
                others.sort(reverse=True)
                
                # if there are no galaxies, add to isolated
                if len(candidates) == 0 and len(others) == 0:
                    add_to_isolated(Lya_v,\
                                    e_Lya_v,\
                                    Lya_W,\
                                    e_Lya_W,\
                                    Na,\
                                    e_Na,\
                                    b,\
                                    e_b,\
                                    W,\
                                    e_W,\
                                    target,\
                                    z_target,\
                                    RA_target,\
                                    Dec_target)
                
                
                # if there are no candidates, regardless of how many others are there:
                elif len(candidates) == 0:
                    # then this is an L_isolated absorber
                    # just grab the first 'other' galaxy and take the absorber info from it
                    l, galaxy_info = others[0]
                    add_to_L_isolated(galaxy_info['Lya_v'],\
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
            
            
                # if there's only one entry in both, then it's an L_associated_isolated
                elif len(candidates) == 1 and len(others) == 0:
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
                                                galaxy_info['Name'],\
                                                galaxy_info['RAdeg'],\
                                                galaxy_info['DEdeg'],\
                                                galaxy_info['impact'],\
                                                galaxy_info['azimuth'],\
                                                galaxy_info['PA'],\
                                                galaxy_info['inc'],\
                                                galaxy_info['adjustedInc'],\
                                                galaxy_info['l'],\
                                                galaxy_info['l_cus'],\
                                                galaxy_info['R_vir'],\
                                                galaxy_info['cus'],\
                                                galaxy_info['MajDiam'],\
                                                galaxy_info['MType'],\
                                                galaxy_info['Vhel'],\
                                                galaxy_info['vcorr'],\
                                                galaxy_info['bestDist'],\
                                                galaxy_info['e_bestDist'],\
                                                galaxy_info['group_num'],\
                                                galaxy_info['group_mem'],\
                                                galaxy_info['group_dist'],\
                                                galaxy_info['Lstar_med'],\
                                                galaxy_info['e_Lstar_med'],\
                                                galaxy_info['Lstar_sdss'],\
                                                galaxy_info['e_Lstar_sdss'],\
                                                galaxy_info['Bmag'],\
                                                galaxy_info['Bmag_sdss'])
        
                # if there's only one with L > min_L, but others exist
                elif len(candidates) == 1 and len(others) >= 1:
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
                                            galaxy_info['Name'],\
                                            galaxy_info['RAdeg'],\
                                            galaxy_info['DEdeg'],\
                                            galaxy_info['impact'],\
                                            galaxy_info['azimuth'],\
                                            galaxy_info['PA'],\
                                            galaxy_info['inc'],\
                                            galaxy_info['adjustedInc'],\
                                            galaxy_info['l'],\
                                            galaxy_info['l_cus'],\
                                            galaxy_info['R_vir'],\
                                            galaxy_info['cus'],\
                                            galaxy_info['MajDiam'],\
                                            galaxy_info['MType'],\
                                            galaxy_info['Vhel'],\
                                            galaxy_info['vcorr'],\
                                            galaxy_info['bestDist'],\
                                            galaxy_info['e_bestDist'],\
                                            galaxy_info['group_num'],\
                                            galaxy_info['group_mem'],\
                                            galaxy_info['group_dist'],\
                                            galaxy_info['Lstar_med'],\
                                            galaxy_info['e_Lstar_med'],\
                                            galaxy_info['Lstar_sdss'],\
                                            galaxy_info['e_Lstar_sdss'],\
                                            galaxy_info['Bmag'],\
                                            galaxy_info['Bmag_sdss'])
                
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
                                            galaxy_info['Name'],\
                                            galaxy_info['RAdeg'],\
                                            galaxy_info['DEdeg'],\
                                            galaxy_info['impact'],\
                                            galaxy_info['azimuth'],\
                                            galaxy_info['PA'],\
                                            galaxy_info['inc'],\
                                            galaxy_info['adjustedInc'],\
                                            galaxy_info['l'],\
                                            galaxy_info['l_cus'],\
                                            galaxy_info['R_vir'],\
                                            galaxy_info['cus'],\
                                            galaxy_info['MajDiam'],\
                                            galaxy_info['MType'],\
                                            galaxy_info['Vhel'],\
                                            galaxy_info['vcorr'],\
                                            galaxy_info['bestDist'],\
                                            galaxy_info['e_bestDist'],\
                                            galaxy_info['group_num'],\
                                            galaxy_info['group_mem'],\
                                            galaxy_info['group_dist'],\
                                            galaxy_info['Lstar_med'],\
                                            galaxy_info['e_Lstar_med'],\
                                            galaxy_info['Lstar_sdss'],\
                                            galaxy_info['e_Lstar_sdss'],\
                                            galaxy_info['Bmag'],\
                                            galaxy_info['Bmag_sdss'],\
                                            others_within_rigor)
                        

                # if there's more than one with L > min_L
                else:
                    likelihood, galaxy_info = candidates.pop(0)
                
                    candidates_within_rigor = []
                
                    # first compare to other candidates
                    for c in candidates:
                        l, candidate_info = c
                    
                        if l * rigor >= likelihood:
                            candidates_within_rigor.append(candidate_info)
                        
                    for o in others:
                        l, other_info = o
                    
                        if l * rigor >= likelihood:
                            candidates_within_rigor.append(other_info)
                    
                    # if none are close enough, add this one to L_associated
                    if len(candidates_within_rigor) == 0:
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
                                            galaxy_info['Name'],\
                                            galaxy_info['RAdeg'],\
                                            galaxy_info['DEdeg'],\
                                            galaxy_info['impact'],\
                                            galaxy_info['azimuth'],\
                                            galaxy_info['PA'],\
                                            galaxy_info['inc'],\
                                            galaxy_info['adjustedInc'],\
                                            galaxy_info['l'],\
                                            galaxy_info['l_cus'],\
                                            galaxy_info['R_vir'],\
                                            galaxy_info['cus'],\
                                            galaxy_info['MajDiam'],\
                                            galaxy_info['MType'],\
                                            galaxy_info['Vhel'],\
                                            galaxy_info['vcorr'],\
                                            galaxy_info['bestDist'],\
                                            galaxy_info['e_bestDist'],\
                                            galaxy_info['group_num'],\
                                            galaxy_info['group_mem'],\
                                            galaxy_info['group_dist'],\
                                            galaxy_info['Lstar_med'],\
                                            galaxy_info['e_Lstar_med'],\
                                            galaxy_info['Lstar_sdss'],\
                                            galaxy_info['e_Lstar_sdss'],\
                                            galaxy_info['Bmag'],\
                                            galaxy_info['Bmag_sdss'])
                    
                    # if one other is within *rigor*, add to L_two list
                    elif len(candidates_within_rigor) == 1:
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
                                    galaxy_info['Name'],\
                                    galaxy_info['RAdeg'],\
                                    galaxy_info['DEdeg'],\
                                    galaxy_info['impact'],\
                                    galaxy_info['azimuth'],\
                                    galaxy_info['PA'],\
                                    galaxy_info['inc'],\
                                    galaxy_info['adjustedInc'],\
                                    galaxy_info['l'],\
                                    galaxy_info['l_cus'],\
                                    galaxy_info['R_vir'],\
                                    galaxy_info['cus'],\
                                    galaxy_info['MajDiam'],\
                                    galaxy_info['MType'],\
                                    galaxy_info['Vhel'],\
                                    galaxy_info['vcorr'],\
                                    galaxy_info['bestDist'],\
                                    galaxy_info['e_bestDist'],\
                                    galaxy_info['group_num'],\
                                    galaxy_info['group_mem'],\
                                    galaxy_info['group_dist'],\
                                    galaxy_info['Lstar_med'],\
                                    galaxy_info['e_Lstar_med'],\
                                    galaxy_info['Lstar_sdss'],\
                                    galaxy_info['e_Lstar_sdss'],\
                                    galaxy_info['Bmag'],\
                                    galaxy_info['Bmag_sdss'],\
                                    candidates_within_rigor)
                
                    elif len(candidates_within_rigor) > 2:
                        add_to_L_three_plus(galaxy_info['Lya_v'],\
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
                                        galaxy_info['Name'],\
                                        galaxy_info['RAdeg'],\
                                        galaxy_info['DEdeg'],\
                                        galaxy_info['impact'],\
                                        galaxy_info['azimuth'],\
                                        galaxy_info['PA'],\
                                        galaxy_info['inc'],\
                                        galaxy_info['adjustedInc'],\
                                        galaxy_info['l'],\
                                        galaxy_info['l_cus'],\
                                        galaxy_info['R_vir'],\
                                        galaxy_info['cus'],\
                                        galaxy_info['MajDiam'],\
                                        galaxy_info['MType'],\
                                        galaxy_info['Vhel'],\
                                        galaxy_info['vcorr'],\
                                        galaxy_info['bestDist'],\
                                        galaxy_info['e_bestDist'],\
                                        galaxy_info['group_num'],\
                                        galaxy_info['group_mem'],\
                                        galaxy_info['group_dist'],\
                                        galaxy_info['Lstar_med'],\
                                        galaxy_info['e_Lstar_med'],\
                                        galaxy_info['Lstar_sdss'],\
                                        galaxy_info['e_Lstar_sdss'],\
                                        galaxy_info['Bmag'],\
                                        galaxy_info['Bmag_sdss'],\
                                        candidates_within_rigor)
            
            
                # now deal with group galaxies
                if len(candidates) >= 1:
                    for c in candidates:
                        l, galaxy_info = c
                        group_num = str(galaxy_info['group_num'])
                    
                        if group_num != '-99':
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
                            galaxy_info['Name'],\
                            galaxy_info['RAdeg'],\
                            galaxy_info['DEdeg'],\
                            galaxy_info['impact'],\
                            galaxy_info['azimuth'],\
                            galaxy_info['PA'],\
                            galaxy_info['inc'],\
                            galaxy_info['adjustedInc'],\
                            galaxy_info['l'],\
                            galaxy_info['l_cus'],\
                            galaxy_info['R_vir'],\
                            galaxy_info['cus'],\
                            galaxy_info['MajDiam'],\
                            galaxy_info['MType'],\
                            galaxy_info['Vhel'],\
                            galaxy_info['vcorr'],\
                            galaxy_info['bestDist'],\
                            galaxy_info['e_bestDist'],\
                            galaxy_info['group_num'],\
                            galaxy_info['group_mem'],\
                            galaxy_info['group_dist'],\
                            galaxy_info['Lstar_med'],\
                            galaxy_info['e_Lstar_med'],\
                            galaxy_info['Lstar_sdss'],\
                            galaxy_info['e_Lstar_sdss'],\
                            galaxy_info['Bmag'],\
                            galaxy_info['Bmag_sdss'])
                            
                            break


##########################################################################################
##########################################################################################
##########################################################################################

    pickle.dump(isolated, isolated_file)
    pickle.dump(L_isolated, L_isolated_file)
    pickle.dump(L_associated_isolated, L_associated_isolated_file)
    pickle.dump(L_associated, L_associated_file)
    pickle.dump(L_nonassociated, L_nonassociated_file)
    pickle.dump(L_two, L_two_file)
    pickle.dump(L_three_plus, L_three_plus_file)
    pickle.dump(L_group, L_group_file)
    pickle.dump(L_summed, L_summed_file)
    pickle.dump(all, all_file)

    
    isolated_file.close()
    L_isolated_file.close()
    L_associated_isolated_file.close()
    L_associated_file.close()
    L_nonassociated_file.close()
    L_two_file.close()
    L_three_plus_file.close()
    L_group_file.close()
    L_summed_file.close()
    all_file.close()


    theFile.close()
    print
    print 'Done!'
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    