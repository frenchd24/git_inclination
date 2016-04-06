#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plotImpact_vir.py, v 1.0 03/16/2016

Plots histograms of impact parameter, and various versions of normalized impact parameters
also plots impact parameter vs R_vir

v1: separated from plotImpactHist_Diam2.py, include a new function to plot Wakker & Savage
    2009 data as well. (03/16/16)

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
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# ## for Palatino and other serif fonts use:
# #rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)

    

###########################################################################

    
def main():
    # assuming 'theFile' contains one name per line, read the file
    
    
    if getpass.getuser() == 'David':
        pickleFilename = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/pilotData2.p'
        resultsFilename = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_8_edit2.csv'
        saveDirectory = '/Users/David/Research_Documents/inclination/git_inclination/pilot_paper_code/plots2/'
        WS09data = '/Users/David/Research_Documents/inclination/git_inclination/WS2009_lya_data.tsv'

    elif getpass.getuser() == 'frenchd':
        pickleFilename = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/pilotData2.p'
        resultsFilename = '/usr/users/frenchd/inclination/git_inclination/LG_correlation_combined5_8_edit2.csv'
        saveDirectory = '/usr/users/frenchd/inclination/git_inclination/pilot_paper_code/plots2/'
        WS09data = '/usr/users/frenchd/inclination/git_inclination/WS2009_lya_data.tsv'

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # use the old pickle file to get the full galaxy dataset info
    pickleFile = open(pickleFilename,'rU')
    fullDict = pickle.load(pickleFile)
    
    pickleFile.close()
    
    
    # save each plot?
    save = False
    
    results = open(resultsFilename,'rU')
    reader = csv.DictReader(results)
    
    WS = open(WS09data,'rU')
    WSreader = csv.DictReader(WS,delimiter=';')
    
    virInclude = False
    cusInclude = False
    finalInclude = True
    
    # if match, then the includes in the file have to MATCH the includes above. e.g., if 
    # virInclude = False, cusInclude = True, finalInclude = False, then only systems
    # matching those three would be included. Otherwise, all cusInclude = True would be included
    # regardless of the others
    match = False
    
    # all the lists to be used for associated lines
    lyaVList = []
    lyaWList = []
    lyaErrList = []
    naList = []
    bList = []
    impactList = []
    azList = []
    incList = []
    fancyIncList = []
    cosIncList = []
    cosFancyIncList = []
    paList = []
    vcorrList = []
    majList = []
    difList = []
    envList = []
    morphList = []
    m15List = []
    virList = []
    likeList = []
    likem15List = []
    
    
    # WS lists
    WSvcorr = []
    WSdiam = []
    WSimpact =[]
    WSew = []
    WSvel = []
    WSlya = []
    WSvel_dif = []
    WSvir = []
    WSlike = []
    
    l_min = 0.001

    
    for w in WSreader:
        vcorr = w['HV']
        diam = w['Diam']
        rho = w['rho']
        ew = w['EWLya']
        vel = w['LyaVel']
        lya = w['Lya']
        
        if lya == 'Lya  ' and isNumber(diam) and isNumber(ew) and isNumber(rho):
            if float(rho) <=500.0:
                # this is a single galaxy association
                vir = calculateVirialRadius(float(diam))
                
                vel_dif = float(vcorr) - float(vel)
    
                # try this "sphere of influence" value instead
                m15 = float(diam)**1.5

                # first for the virial radius
                likelihood = math.exp(-(float(rho)/vir)**2) * math.exp(-(vel_dif/200.)**2)
                
                if vir>= float(rho):
                    likelihood = likelihood*2
                    
                # then for the second 'virial like' m15 radius
                likelihoodm15 = math.exp(-(float(rho)/m15)**2) * math.exp(-(vel_dif/200.)**2)
                
                if m15>= float(rho):
                    likelihoodm15 = likelihoodm15*2
                    
                if likelihood <= likelihoodm15:
                    likelihood = likelihoodm15
                    
                WSlike.append(likelihood)
                
#                 l_min=0
                
                if likelihood >= l_min:
                
                    WSvcorr.append(float(vcorr))
                    WSdiam.append(float(diam))
                    WSvir.append(vir)
                    WSimpact.append(float(rho))
                    WSew.append(float(ew))
                    WSvel.append(float(vel))
                    WSlya.append(lya)
                    WSvel_dif.append(vel_dif)
    
    
    for l in reader:
        include_vir = eval(l['include_vir'])
        include_cus = eval(l['include_custom'])
        include = eval(l['include'])
        
        go = False
        if match:
            if virInclude == include_vir and cusInclude == include_cus:
                go = True
            else:
                go = False
                
        else:
            if virInclude and include_vir:
                go = True
                
            elif cusInclude and include_cus:
                go = True
                
            elif finalInclude and include:
                go = True
            
            else:
                go = False
        
        if go:
            AGNra_dec = eval(l['degreesJ2000RA_DecAGN'])
            galaxyRA_Dec = eval(l['degreesJ2000RA_DecGalaxy'])
            lyaV = l['Lya_v']
            lyaW = l['Lya_W'].partition('pm')[0]
            lyaW_err = l['Lya_W'].partition('pm')[2]
            env = l['environment']
            galaxyName = l['galaxyName']
            impact = l['impactParameter (kpc)']
            galaxyDist = l['distGalaxy (Mpc)']
            pa = l['positionAngle (deg)']
            RC3pa = l['RC3pa (deg)']
            morph = l['morphology']
            vcorr = l['vcorrGalaxy (km/s)']
            maj = l['majorAxis (kpc)']
            minor = l['minorAxis (kpc)']
            inc = l['inclination (deg)']
            az = l['azimuth (deg)']
            b = l['b'].partition('pm')[0]
            b_err = l['b'].partition('pm')[2]
            na = eval(l['Na'].partition(' pm ')[0])
            na_err = eval(l['Na'].partition(' pm ')[2])
            likelihood = l['likelihood']
            likelihoodm15 = l['likelihood_1.5']
            virialRadius = l['virialRadius']
            m15 = l['d^1.5']
            vel_diff = l['vel_diff']
            
            if isNumber(RC3pa) and not isNumber(pa):
                pa = RC3pa
            
            if isNumber(inc):
                cosInc = cos(float(inc) * pi/180.)
                
                if isNumber(maj) and isNumber(minor):
                    q0 = 0.2
                    fancyInc = calculateFancyInclination(maj,minor,q0)
                    cosFancyInc = cos(fancyInc * pi/180)
                else:
                    fancyInc = -99
                    cosFancyInc = -99
            else:
                cosInc = -99
                inc = -99
                fancyInc = -99
                cosFancyInc = -99
            
            # all the lists to be used for associated lines
            lyaVList.append(float(lyaV))
            lyaWList.append(float(lyaW))
            lyaErrList.append(float(lyaW_err))
            naList.append(na)
            bList.append(float(b))
            impactList.append(float(impact))
            azList.append(az)
            incList.append(float(inc))
            fancyIncList.append(fancyInc)
            cosIncList.append(cosInc)
            cosFancyIncList.append(cosFancyInc)
            paList.append(pa)
            vcorrList.append(vcorr)
            majList.append(maj)
            difList.append(float(vel_diff))
            envList.append(float(env))
            morphList.append(morph)
            m15List.append(m15)
            virList.append(virialRadius)
            likeList.append(likelihood)
            likem15List.append(likelihoodm15)

    results.close()
    WS.close()
            
    # lists for the full galaxy dataset
    allPA = fullDict['allPA']
    allInclinations = fullDict['allInclinations']
    allCosInclinations = fullDict['allCosInclinations']
    allFancyInclinations = fullDict['allFancyInclinations']
    allCosFancyInclinations = fullDict['allCosFancyInclinations']
    
    total = 0
    totalNo = 0
    totalYes = 0
    totalIsolated = 0
    totalGroup = 0
    


##########################################################################################
##########################################################################################
    # plot a histogram of the Wakker & Savage 2009 data likelihoods
    #
    
    plot_like_WS2009 = False
    save = False
    
    if plot_like_WS2009:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        labelg = 'WS2009 Absorber'
        alpha = 0.85

#         bins = arange(0,0.01,0.0005)
        bins = 15
        
        plot2 = hist(WSlike,bins = bins,alpha=alpha)
            
#         legend(scatterpoints=1,prop={'size':12},loc=2)
#         ax.grid(b=None,which='major',axis='both')

        if save:
            savefig('{0}/hist(WS2009_like).pdf'.format(saveDirectory),format='pdf')
        else:
            show()

##########################################################################################
##########################################################################################
    # plot impact parameters of absorbers as a function of virial radius of the associated
    # galaxy, add in the Wakker & Savage 2009 data 
    #
    
    plotImpact_vs_virial_WS2009 = True
    save = True
    
    if plotImpact_vs_virial_WS2009:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Redshifted Absorber'
        labelb = "Blueshifted Absorber"
        labelg = 'WS09 Absorber'
        alpha = 0.85

        for d,i,v in zip(difList,impactList,virList):
            # check if all the values are okay
            print 'd: ',d
            if isNumber(d) and isNumber(i) and isNumber(v):
                if d!=-99 and i!=-99 and v!=-99:
#                     print 'd: ',d
                    if d>0:
                        # galaxy is 'behind' absorber, so GAS = blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(v,i,c='Blue',s=50,label= labelb,alpha=alpha)
                    if d<0:
                        # gas is RED shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(v,i,c='Red',s=50,label= labelr,alpha=alpha)
                
                    plot1 = scatter(v,i,c=color,s=50,alpha=alpha)
        
        plot2 = scatter(WSvir,WSimpact,s=50,alpha=alpha,c='Green',label=labelg)
            
        xlabel(r'$\rm R_{vir}$ (kpc)')
        ylabel(r'$\rm \rho$ (kpc)')
        legend(scatterpoints=1,prop={'size':12},loc=1)
        ax.grid(b=None,which='major',axis='both')
        xlim(-1,350)
        ylim(-10,500)

        if save:
            savefig('{0}/impact(virial)_WS09_lmin_cut.pdf'.format(saveDirectory),format='pdf')
        else:
            show()

##########################################################################################
##########################################################################################
    # plot impact parameters of absorbers as a function of virial radius of the associated
    # galaxy
    #
    
    plotImpact_vs_virial = False
    save = False
    
    if plotImpact_vs_virial:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        alpha = 0.85

        for d,i,v in zip(difList,impactList,virList):
            # check if all the values are okay
            print 'd: ',d
            if isNumber(d) and isNumber(i) and isNumber(v):
                if d!=-99 and i!=-99 and v!=-99:
#                     print 'd: ',d
                    if d>0:
                        # galaxy is 'behind' absorber, so GAS = blue shifted
                        color = 'Blue'
                        if countb == 0:
                            countb +=1
                            plotb = ax.scatter(v,i,c='Blue',s=50,label= labelb,alpha=alpha)
                    if d<0:
                        # gas is RED shifted compared to galaxy
                        color = 'Red'
                        if countr == 0:
                            countr +=1
                            plotr = ax.scatter(v,i,c='Red',s=50,label= labelr,alpha=alpha)
                
                    plot1 = scatter(v,i,c=color,s=50,alpha=alpha)
            
        xlabel(r'$\rm R_{vir}$ (kpc)')
        ylabel(r'$\rm \rho$ (kpc)')
        legend(scatterpoints=1,prop={'size':12},loc=2)
        ax.grid(b=None,which='major',axis='both')

        if save:
            savefig('{0}/impact(virial).pdf'.format(saveDirectory),format='pdf')
        else:
            show()
            
            
##########################################################################################
##########################################################################################  

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


if __name__=="__main__":
    # do the work
    main()
    