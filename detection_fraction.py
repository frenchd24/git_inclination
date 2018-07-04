#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)


$Id: detection_fraction.py v 1.0 07/01/2018

Determine the detection fraction. Correlate the sightlines with galaxies, and figure out 
for what level of likelihood an absorber is always found, etc.

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
    maxSep = 1000.
    agnSeparation = 4000.
    minVcorr = 450.
    minSize = 0.
    max_deltav = 400.
    rigor = 5
    
    Lstar_min = 1.0
    
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
        gtPickleFilename = '/Users/frenchd/Research/GT_update2/pickleGT_filteredAll.p'

#         filename = '/Users/frenchd/Research/inclination/git_inclination/correlatedTargetList_5_29_18_measurements.csv'
#         filename = '/Users/frenchd/Research/inclination/git_inclination/targets/correlatedTargetList_5_29_18_measurements_copy_short.csv'
#         filename = '/Users/frenchd/Research/inclination/git_inclination/targets/correlatedTargetList_5_29_18_measurements_copy.csv'
        filename = '/Users/frenchd/Research/inclination/git_inclination/targets/correlatedTargetList_5_29_18_measurements_copy.csv'

        # pickle files
        detection_fraction_filename = '/Users/frenchd/Research/inclination/git_inclination/detection_fraction_lstarmin_1.p'

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

##########################################################################################
    # now for likelihood thresholds
    # lists of \Delta v for galaxies within 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 0.75 L of a target
    dv_l001 = []
    dv_l005 = []
    dv_l01 = []
    dv_l05 = []
    dv_l1 = []
    dv_l5 = []
    dv_l75 = []
    
    # count of detections within each likelihood window
    dv400_l001_det = 0
    dv400_l005_det = 0
    dv400_l01_det = 0
    dv400_l05_det = 0
    dv400_l1_det = 0
    dv400_l5_det = 0
    dv400_l75_det = 0
    
    
    # count of non-detections within each likelihood window
    dv400_l001_non = 0
    dv400_l005_non = 0
    dv400_l01_non = 0
    dv400_l05_non = 0
    dv400_l1_non = 0
    dv400_l5_non = 0
    dv400_l75_non = 0
    
    
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
        
        proceed = False
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


        # do the correlation
        correlation = correlateSingle.correlateTarget(target, maxSep, agnSeparation, minVcorr, minSize, slow=False, searchAll=False, returnDict=True)
            
        for c in correlation:
            proceed = False
            
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
            
            if str(Lstar_med) == '-99.99':
                Lstar = Lstar_sdss
            else:
                Lstar = Lstar_med
                
            if float(Lstar) >= Lstar_min:
                proceed = True


            if proceed:
                # if there's a galaxy within 500 kpc, search for a line 
                if float(impact) <= 1000:
                    detection = False
                    for Lya_v in Lya_vs:
                        if isNumber(Lya_v):
                            dv = float(Lya_v) - Vhel
                            dv_1000.append(dv)
                        
                            if abs(dv) <= 400:
                                detection = True
                            
                    if detection:
                        dv400_imp1000_det +=1
                    else:
                        dv400_imp1000_non +=1
            
                # if there's a galaxy within 500 kpc, search for a line 
                if float(impact) <= 750:
                    detection = False
                    for Lya_v in Lya_vs:
                        if isNumber(Lya_v):
                            dv = float(Lya_v) - Vhel
                            dv_750.append(dv)
                        
                            if abs(dv) <= 400:
                                detection = True
                            
                    if detection:
                        dv400_imp750_det +=1
                    else:
                        dv400_imp750_non +=1
            
                # if there's a galaxy within 500 kpc, search for a line 
                if float(impact) <= 500:
                    detection = False
                    for Lya_v in Lya_vs:
                        if isNumber(Lya_v):
                            dv = float(Lya_v) - Vhel
                            dv_500.append(dv)
                        
                            if abs(dv) <= 400:
                                detection = True
                            
                    if detection:
                        dv400_imp500_det +=1
                    else:
                        dv400_imp500_non +=1
            
                # if there's a galaxy within 400 kpc, search for a line 
                if float(impact) <= 400:
                    detection = False
                    for Lya_v in Lya_vs:
                        if isNumber(Lya_v):
                            dv = float(Lya_v) - Vhel
                            dv_400.append(dv)
            
                            if abs(dv) <= 400:
                                detection = True
                            
                    if detection:
                        dv400_imp400_det +=1
                    else:
                        dv400_imp400_non +=1
            
                # if there's a galaxy within 300 kpc, search for a line 
                if float(impact) <= 300:
                    detection = False
                    for Lya_v in Lya_vs:
                        if isNumber(Lya_v):
                            dv = float(Lya_v) - Vhel
                            dv_300.append(dv)
                        
                            if abs(dv) <= 400:
                                detection = True
                            
                    if detection:
                        dv400_imp300_det +=1
                    else:
                        dv400_imp300_non +=1

                # if there's a galaxy within 200 kpc, search for a line 
                if float(impact) <= 200:
                    detection = False
                    for Lya_v in Lya_vs:
                        if isNumber(Lya_v):
                            dv = float(Lya_v) - Vhel
                            dv_200.append(dv)

                            if abs(dv) <= 400:
                                detection = True
                            
                    if detection:
                        dv400_imp200_det +=1
                    else:
                        dv400_imp200_non +=1

                # if there's a galaxy within 100 kpc, search for a line 
                if float(impact) <= 100:
                    detection = False
                    for Lya_v in Lya_vs:
                        if isNumber(Lya_v):
                            dv = float(Lya_v) - Vhel
                            dv_100.append(dv)
                        
                            if abs(dv) <= 400:
                                detection = True
                            
                    if detection:
                        dv400_imp100_det +=1
                    else:
                        dv400_imp100_non +=1

                # if there's a galaxy within 50 kpc, search for a line 
                if float(impact) <= 50:
                    detection = False
                    for Lya_v in Lya_vs:
                        if isNumber(Lya_v):
                            dv = float(Lya_v) - Vhel
                            dv_50.append(dv)
            
                            if abs(dv) <= 400:
                                detection = True
                            
                    if detection:
                        dv400_imp50_det +=1
                    else:
                        dv400_imp50_non +=1
                    
                # if there's a galaxy within 25 kpc, search for a line 
                if float(impact) <= 25:
                    detection = False
                    for Lya_v in Lya_vs:
                        if isNumber(Lya_v):
                            dv = float(Lya_v) - Vhel
                            dv_25.append(dv)
            
                            if abs(dv) <= 400:
                                detection = True
                            
                    if detection:
                        dv400_imp25_det +=1
                    else:
                        dv400_imp25_non +=1
                        
                        print 'non-detection within 25: {0} - {1}'.format(target, Name)
                    
            
    ##########################################################################################            
                # now do it for likelihood

                if isNumber(MajDiam) and MajDiam != -99.99 and MajDiam != '-99.99':
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

            
                    # if the likelihood is greater than 0.001, see if there's a corresponding
                    if l_used >= 0.001:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_l001.append(dv)
                        
                                if abs(dv) <= 400:
                                    detection = True
                            
                        if detection:
                            dv400_l001_det +=1
                        else:
                            dv400_l001_non +=1
            
                    # if the likelihood is greater than 0.005, see if there's a corresponding
                    if l_used >= 0.005:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_l005.append(dv)
                        
                                if abs(dv) <= 400:
                                    detection = True
                            
                        if detection:
                            dv400_l005_det +=1
                        else:
                            dv400_l005_non +=1


                    # if the likelihood is greater than 0.01, see if there's a corresponding
                    if l_used >= 0.01:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_l01.append(dv)
                        
                                if abs(dv) <= 400:
                                    detection = True
                            
                        if detection:
                            dv400_l01_det +=1
                        else:
                            dv400_l01_non +=1
                        
                        
                        
                    # if the likelihood is greater than 0.05, see if there's a corresponding
                    if l_used >= 0.05:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_l05.append(dv)
                        
                                if abs(dv) <= 400:
                                    detection = True
                            
                        if detection:
                            dv400_l05_det +=1
                        else:
                            dv400_l05_non +=1
                        
                        
                    # if the likelihood is greater than 0.1, see if there's a corresponding
                    if l_used >= 0.1:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_l1.append(dv)
                        
                                if abs(dv) <= 400:
                                    detection = True
                            
                        if detection:
                            dv400_l1_det +=1
                        else:
                            dv400_l1_non +=1
                        
                        
                    # if the likelihood is greater than 0.5, see if there's a corresponding
                    if l_used >= 0.5:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_l5.append(dv)
                        
                                if abs(dv) <= 400:
                                    detection = True
                            
                        if detection:
                            dv400_l5_det +=1
                        else:
                            dv400_l5_non +=1
                        
                        
                    # if the likelihood is greater than 0.75, see if there's a corresponding
                    if l_used >= 0.75:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_l75.append(dv)
                        
                                if abs(dv) <= 400:
                                    detection = True
                            
                        if detection:
                            dv400_l75_det +=1
                        else:
                            dv400_l75_non +=1
                        
                        
                        
                        
    ##########################################################################################            
                # now do it for imp/R_vir
                
                    l_used = impact/R_vir
            
                    # if the imp/Rvir is greater than 0.25, see if there's a corresponding line
                    if l_used >= 0.25:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_rvir025.append(dv)
                        
                                if abs(dv) <= 400:
                                    detection = True
                            
                        if detection:
                            dv400_rvir025_det +=1
                        else:
                            dv400_rvir025_non +=1
            
                    # if the imp/Rvir is greater than 0.5, see if there's a corresponding line
                    if l_used >= 0.5:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_rvir05.append(dv)
                        
                                if abs(dv) <= 400:
                                    detection = True
                            
                        if detection:
                            dv400_rvir05_det +=1
                        else:
                            dv400_rvir05_non +=1


                    # if the imp/Rvir is greater than 0.75, see if there's a corresponding line
                    if l_used >= 0.75:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_rvir075.append(dv)
                        
                                if abs(dv) <= 400:
                                    detection = True
                            
                        if detection:
                            dv400_rvir075_det +=1
                        else:
                            dv400_rvir075_non +=1
                        
                        
                        
                    # if the imp/Rvir is greater than 1.0, see if there's a corresponding line
                    if l_used >= 1.0:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_rvir1.append(dv)
                        
                                if abs(dv) <= 400:
                                    detection = True
                            
                        if detection:
                            dv400_rvir1_det +=1
                        else:
                            dv400_rvir1_non +=1
                        
                        
                    # if the imp/Rvir is greater than 1.5, see if there's a corresponding line
                    if l_used >= 1.5:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_rvir15.append(dv)
                        
                                if abs(dv) <= 400:
                                    detection = True
                            
                        if detection:
                            dv400_rvir15_det +=1
                        else:
                            dv400_rvir15_non +=1
                        
                        
                    # if the imp/Rvir is greater than 2.0, see if there's a corresponding line
                    if l_used >= 2.0:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_rvir2.append(dv)
                        
                                if abs(dv) <= 400:
                                    detection = True
                            
                        if detection:
                            dv400_rvir2_det +=1
                        else:
                            dv400_rvir2_non +=1
                        
                        
                    # if the imp/Rvir is greater than 1.5, see if there's a corresponding line
                    if l_used >= 2.5:
                        detection = False
                        for Lya_v in Lya_vs:
                            if isNumber(Lya_v):
                                dv = float(Lya_v) - Vhel
                                dv_rvir25.append(dv)
                        
                                if abs(dv) <= 400:
                                    detection = True
                            
                        if detection:
                            dv400_rvir25_det +=1
                        else:
                            dv400_rvir25_non +=1


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
        print 'Detection fraction for 0.001 L: ', float(dv400_l001_det) / (dv400_l001_det   + dv400_l001_non)
        print 'Detection fraction for 0.005 L: ', float(dv400_l005_det) / (dv400_l005_det   + dv400_l005_non)
        print 'Detection fraction for 0.01 L: ', float(dv400_l01_det)   / (dv400_l01_det    + dv400_l01_non)
        print 'Detection fraction for 0.05 L: ', float(dv400_l05_det)   / (dv400_l05_det    + dv400_l05_non)
        print 'Detection fraction for 0.1 L: ', float(dv400_l1_det)     / (dv400_l1_det     + dv400_l1_non)
        print 'Detection fraction for 0.5 L: ', float(dv400_l5_det)     / (dv400_l5_det     + dv400_l5_non)
        print 'Detection fraction for 0.75 L: ', float(dv400_l75_det)   / (dv400_l75_det    + dv400_l75_non)
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
    #     print 'dv_l001: ',dv_l001
        print
        print
    
    except Exception, e:
        print 'Error: ',e
        print
        
    
    full_dict = {}
    
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



    # now for likelihood thresholds
    full_dict['dv_l001'] = dv_l001
    full_dict['dv_l005'] = dv_l005
    full_dict['dv_l01'] = dv_l01
    full_dict['dv_l05'] = dv_l05
    full_dict['dv_l1'] = dv_l1
    full_dict['dv_l5'] = dv_l5
    full_dict['dv_l75'] = dv_l75
    
    # now for likelihood detections
    full_dict['dv400_l001_det'] = dv400_l001_det
    full_dict['dv400_l005_det'] = dv400_l005_det
    full_dict['dv400_l01_det'] = dv400_l01_det
    full_dict['dv400_l05_det'] = dv400_l05_det
    full_dict['dv400_l1_det'] = dv400_l1_det
    full_dict['dv400_l5_det'] = dv400_l5_det
    full_dict['dv400_l75_det'] = dv400_l75_det
    
    # now for likelihood non-detections
    full_dict['dv400_l001_non'] = dv400_l001_non
    full_dict['dv400_l005_non'] = dv400_l005_non
    full_dict['dv400_l01_non'] = dv400_l01_non
    full_dict['dv400_l05_non'] = dv400_l05_non
    full_dict['dv400_l1_non'] = dv400_l1_non
    full_dict['dv400_l5_non'] = dv400_l5_non
    full_dict['dv400_l75_non'] = dv400_l75_non
    
    
    
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
    
    
    
##########################################################################################
##########################################################################################
##########################################################################################
    
    pickle.dump(full_dict, detection_fraction)
    detection_fraction.close()

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
    