#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)


$Id:  plot_impact_Rvir.py, v 1.8 06/11/18

Plot R_vir vs impact parameter


Based on:
$Id:  plotImpact_vir_final.py, v 1.7 01/05/18)

Plots impact parameter vs R_vir in a number of ways


v1: separated from plotImpactHist_Diam2.py, include a new function to plot Wakker & Savage
    2009 data as well. (03/16/16)
    
v1.1: remake plots with v_hel instead of vcorr (4/21/16)

v1.2: remake plots with new large galaxy sample (7/13/16) -> /plots4/

v1.3: add ability to limit results based on environment and/or likelihood (7/14/16)

v1.4: minor updates to formatting for new pilot paper sample (large galaxies only)
        - (8/05/16)
        
v1.5: more minor updates for LG_correlation_combined5_11_25cut_edit4.csv -> /plots5/
    (9/26/16)
    
v1.6: include line indicating upper envelope on impact/R_vir = 2.14597 for L=0.01
    (10/07/16)
    
v1.7: update for the first referee report (12/16/16)
    - larger marker points (60), make open markers for imp/R_vir <=1 systems
    
v1.8 update for AAS_2017 (01/05/18)
'''

import sys
import os
import csv
from scipy import stats

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

fontScale = 16
rc('text', usetex=True)
rc('font', size=16, family='serif', weight='normal')
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


def errors(a):
    # return the standard error in the mean for the input array
    return stats.sem(a)
    


def main():
    # Plot impact parameters of absorbers as a function of virial radius of the associated
    # galaxy, include median histograms    
    plot_impact_virial_median = True
    plot_impact_virial_median_save = True
    

    # some colors
    color_blue = '#436bad'      # french blue
    color_red = '#ec2d01'     # tomato red

    if getpass.getuser() == 'frenchd':

#         pickleFilename = '/Users/frenchd/Research/inclination/git_inclination/picklePilot_plusSALT_14.p'
#         gtPickleFilename = '/Users/frenchd/Research/inclination/git_inclination/pickleGT.p'
#         saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/plotting_code/figs'

#         gtPickleFilename = '/Users/frenchd/Research/inclination/git_inclination/pickleGT_filteredAll.p'
        gtPickleFilename = '/Users/frenchd/Research/GT_update2/pickleGT_filteredAll.p'

        saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/plotting_code/figs/'
        
        isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/isolated2.p'
        L_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_isolated2.p'
        L_associated_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated2.p'
        L_associated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated2.p'
        L_nonassociated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated2.p'
        L_two_filename = '/Users/frenchd/Research/inclination/git_inclination/L_two2.p'
        L_two_plus_filename = '/Users/frenchd/Research/inclination/git_inclination/L_two_plus2.p'
        L_group_filename = '/Users/frenchd/Research/inclination/git_inclination/L_group2.p'


    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # pickle file for the whole galaxy table:
    gtPickleFile = open(gtPickleFilename,'rU')
    gtDict = pickle.load(gtPickleFile)
    gtPickleFile.close()
    
    
    # open all the pickle files
    isolated_file = open(isolated_filename,'r')
    L_isolated_file = open(L_isolated_filename,'r')
    L_associated_isolated_file = open(L_associated_isolated_filename,'r')
    L_associated_file = open(L_associated_filename,'r')
    L_nonassociated_file = open(L_nonassociated_filename,'r')
    L_two_file = open(L_two_filename,'r')
    L_two_plus_file = open(L_two_plus_filename,'r')
    L_group_file = open(L_group_filename,'r')

    # unload the data from them
    isolated = pickle.load(isolated_file)
    L_isolated = pickle.load(L_isolated_file)
    L_associated_isolated = pickle.load(L_associated_isolated_file)
    L_associated = pickle.load(L_associated_file)
    L_nonassociated = pickle.load(L_nonassociated_file)
    L_two = pickle.load(L_two_file)
    L_two_plus = pickle.load(L_two_plus_file)
    L_group = pickle.load(L_group_file)
    
    # close the files
    isolated_file.close()
    L_isolated_file.close()
    L_associated_isolated_file.close()
    L_associated_file.close()
    L_nonassociated_file.close()
    L_two_file.close()
    L_two_plus_file.close()
    L_group_file.close()
    
    
    # save each plot?
    save = False
    
#     results = open(resultsFilename,'rU')
#     reader = csv.DictReader(results)
    
#     WS = open(WS09data,'rU')
#     WSreader = csv.DictReader(WS,delimiter=';')
    
    maxEnv = 3000
    minL = 0.001
    
    # if match, then the includes in the file have to MATCH the includes above. e.g., if 
    # virInclude = False, cusInclude = True, finalInclude = False, then only systems
    # matching those three would be included. Otherwise, all cusInclude = True would be included
    # regardless of the others
    match = False
    
    dataSet = L_associated_isolated
    
    Lya_vs = dataSet['Lya_vs']
    e_Lya_vs = dataSet['e_Lya_vs']
    Lya_Ws = dataSet['Lya_Ws']
    e_Lya_Ws = dataSet['e_Lya_Ws']
    Nas = dataSet['Nas']
    e_Nas = dataSet['e_Nas']
    bs = dataSet['bs']
    e_bs = dataSet['e_bs']
    Ws = dataSet['Ws']
    e_Ws = dataSet['e_Ws']
    targets = dataSet['targets']
    z_targets = dataSet['z_targets']
    RA_targets = dataSet['RA_targets']
    Dec_targets = dataSet['Dec_targets']
    Names = dataSet['Names']
    RA_galaxies = dataSet['RA_galaxies']
    Dec_galaxies = dataSet['Dec_galaxies']
    impacts = dataSet['impacts']
    azimuths = dataSet['azimuths']
    PAs = dataSet['PAs']
    incs = dataSet['incs']
    adjustedIncs = dataSet['adjustedIncs']
    ls = dataSet['ls']
    l_cuss = dataSet['l_cuss']
    R_virs = dataSet['R_virs']
    cuss = dataSet['cuss']
    MajDiams = dataSet['MajDiams']
    MTypes = dataSet['MTypes']
    Vhels = dataSet['Vhels']
    vcorrs = dataSet['vcorrs']
    bestDists = dataSet['bestDists']
    e_bestDists = dataSet['e_bestDists']
    group_nums = dataSet['group_nums']
    group_mems = dataSet['group_mems']
    group_dists = dataSet['group_dists']
    Lstar_meds = dataSet['Lstar_meds']
    e_Lstar_meds = dataSet['e_Lstar_meds']
    Bmags = dataSet['Bmags']



    majorAxisL = gtDict['majorAxis']
    incL = gtDict['inc']
    adjustedIncL = gtDict['adjustedInc']
    paL = gtDict['PA']
    BmagL = gtDict['Bmag']
#     Bmag_sdssL = gtDict['Bmag_sdss']
    RID_medianL = gtDict['RID_median']
    RID_meanL = gtDict['RID_mean']
    RID_stdL = gtDict['RID_std']
    VhelL = gtDict['Vhel']
    RAdegL = gtDict['RAdeg']
    DEdegL = gtDict['DEdeg']
    NameL= gtDict['Name']
    
    allPA = paL
    allInclinations = []
    allAdjustedIncs = []
    allCosInclinations = []

#     print 'type: ',type(incL)
    for i in incL:
        if i != -99:
            i = float(i)
            allInclinations.append(i)
            
            i2 = pi/180. * i
            cosi2 = cos(i)
            allCosInclinations.append(cosi2)
            
    allFancyInclinations = []
    allCosFancyCosInclinations = []
    for i in adjustedIncL:
        if i != -99:
            i = float(i)

            allAdjustedIncs.append(i)
            
            i2 = pi/180. * i
            cosi2 = cos(i)
            allCosFancyCosInclinations.append(cosi2)
            
    allDiameter = majorAxisL

    print 'finished with this shit'

    total = 0
    totalNo = 0
    totalYes = 0
    totalIsolated = 0
    totalGroup = 0


            
##########################################################################################
##########################################################################################
    
    if plot_impact_virial_median:
        fig = figure(figsize=(7.8,5.9))
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        alpha = 0.7
        alphaInside = 0.7
        markerSize = 60
        errorAlpha = 0.15
        
        plotErrors = False
        
        labelr = r'$\rm Redshifted ~Absorber$'
        labelb = r'$\rm Blueshifted ~Absorber$'
        bSymbol = 'D'
        rSymbol = 'o'
        
        rImpact = []
        rVir = []
        bImpact = []
        bVir = []
        allImpact = []
        allVir = []
        
        for Lya_v, Vhel, i, r, w, e in zip(Lya_vs, Vhels, impacts, R_virs, Lya_Ws, e_Lya_Ws):
            vel_dif = float(Lya_v) - float(Vhel)

            yVal = float(i)
            xVal = float(r)
            allImpact.append(yVal)
            allVir.append(xVal)

            if vel_dif >=0:
                # gas is RED shifted compared to galaxy
                color = color_red
                symbol = rSymbol
                
                if float(i) > float(r):
                    # impact parameter > virial radius
                    a = alpha
                    fc = color
                    ec = 'black'
            
                if float(i) <= float(r):
                    # impact parameter <= virial radius
                    a = alphaInside
                    fc = 'none'
                    ec = color                        
                
                rImpact.append(yVal)
                rVir.append(xVal)
                
                if countr == 0:
                    countr +=1
                    plotr = ax.scatter(xVal,
                                        yVal,
                                        marker=symbol,
                                        c=color_red,
                                        s=markerSize,
                                        facecolor=fc,
                                        edgecolor=ec,
                                        alpha=a)

            else:
                # galaxy is 'behind' absorber, so GAS = blue shifted
                color = color_blue
                symbol = bSymbol

                if float(i) > float(r):
                    # impact parameter > virial radius
                    a = alpha
                    fc = color
                    ec = 'black'
            
                if float(i) <= float(r):
                    # impact parameter <= virial radius
                    a = alphaInside
                    fc = 'none'
                    ec = color                        
                
                bImpact.append(yVal)
                bVir.append(xVal)
                
                if countb == 0:
                    countb +=1
                    plotb = ax.scatter(xVal,
                                        yVal,
                                        marker=symbol,
                                        c=color_blue,
                                        facecolor=fc,
                                        edgecolor=ec,
                                        s=markerSize,
                                        alpha=a)
                                                    

            plot1 = scatter(xVal,
                            yVal,
                            marker=symbol,
                            c=color,
                            s=markerSize,
                            facecolor=fc,
                            edgecolor=ec,
                            alpha=a)
        
        
        # x-axis
        majorLocator   = MultipleLocator(50)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(25)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        
        # y axis
        majorLocator   = MultipleLocator(100)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(50)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        
        binSize = 50
        bins = arange(0,400,binSize)
           

        # redshifted ones
        bin_means,edges,binNumber = stats.binned_statistic(array(rVir), array(rImpact),\
        statistic='mean', bins=bins)
        
        bin_errors,edges_e,binNumber_e = stats.binned_statistic(array(rVir), array(rImpact), \
        statistic=lambda y: errors(y), bins=bins)
        
        bin_std,edges_std,binNumber_std = stats.binned_statistic(array(rVir), array(rImpact), \
        statistic=lambda y: std(y), bins=bins)
        
        print 'bin_means,edges,binNumber : ',bin_means,edges,binNumber
        print
        print 'bin_errors, edges_e,binNumber_e : ',bin_errors,edges_e,binNumber_e
        print
        print 'bin_std,edges_std,binNumber_std : ',bin_std,edges_std,binNumber_std
        
        # the mean
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        
        # the errors
        left_e,right_e = edges_e[:-1],edges_e[1:]        
        X_e = array([left_e,right_e]).T.flatten()
        Y_e = array([nan_to_num(bin_errors),nan_to_num(bin_errors)]).T.flatten()
        
        yErrorsTop = Y + Y_e
        yErrorsBot = Y - Y_e

        if plotErrors:
            plot(X_e,yErrorsBot, ls='solid',color='red',lw=1,alpha=errorAlpha)
            plot(X_e,yErrorsTop, ls='solid',color='red',lw=1,alpha=errorAlpha)
            fill_between(X_e, yErrorsBot, yErrorsTop, facecolor='red', interpolate=True,alpha=errorAlpha)        
        
        plot(X,Y, c='red',ls='dotted',lw=2.1,alpha=alpha+0.2,label=r'$\rm Mean ~Redshifted ~EW$')
            
            
        # blueshifted ones
        bin_means,edges,binNumber = stats.binned_statistic(array(bVir), array(bImpact), \
        statistic='mean', bins=bins)
        
        bin_errors,edges_e,binNumber_e = stats.binned_statistic(array(bVir), array(bImpact), \
        statistic=lambda y: errors(y), bins=bins)
        
        bin_std,edges_std,binNumber_std = stats.binned_statistic(array(bVir), array(bImpact), \
        statistic=lambda y: std(y), bins=bins)
        
        print 'bin_means,edges,binNumber : ',bin_means,edges,binNumber
        print
        print 'bin_errors, edges_e,binNumber_e : ',bin_errors,edges_e,binNumber_e
        print
        print 'bin_std,edges_std,binNumber_std : ',bin_std,edges_std,binNumber_std
        
        # the mean
        left,right = edges[:-1],edges[1:]        
        X = array([left,right]).T.flatten()
        Y = array([nan_to_num(bin_means),nan_to_num(bin_means)]).T.flatten()
        
        # the errors
        left_e,right_e = edges_e[:-1],edges_e[1:]        
        X_e = array([left_e,right_e]).T.flatten()
        Y_e = array([nan_to_num(bin_errors),nan_to_num(bin_errors)]).T.flatten()
        
        yErrorsTop = Y + Y_e
        yErrorsBot = Y - Y_e

        if plotErrors:
            plot(X_e,yErrorsBot, ls='solid',color='blue',lw=1,alpha=errorAlpha)
            plot(X_e,yErrorsTop, ls='solid',color='blue',lw=1,alpha=errorAlpha)
            fill_between(X_e, yErrorsBot, yErrorsTop, facecolor='blue', interpolate=True,alpha=errorAlpha)
        
        plot(X,Y, c='blue',ls='dashed',lw=1.5,alpha=alpha+0.1,label=r'$\rm Mean ~Blueshifted ~EW$')
            

            
        # add an upper envelope line at impact/R_vir = 2.14
        envX = [150,350]
        envY = [321.895,751.088]
        plot(envX,envY, c='black',ls='dashed',lw=2.5,alpha=0.9,label=r'$\rm Max ~ \rho/R_{vir}~ cutoff$')
            
        xlabel(r'$\rm R_{vir} ~[kpc]$')
        ylabel(r'$\rm Impact ~Parameter ~[kpc]$')
        legend(scatterpoints=1,prop={'size':14},loc=1,fancybox=True)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,650)
#         xlim(150,350)
        xlim(0,350)

        tight_layout()
        
        print
        print
        print "allImpact: ",allImpact
        print "allVir: ",allVir
        print
        print 'Done'

        if plot_impact_virial_median_save:
            savefig('{0}/impact(virial)_{1}_median.pdf'.format(saveDirectory, binSize),\
            format='pdf',bbox_inches='tight')
        else:
            show()
            

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


if __name__=="__main__":
    # do the work
    main()
    