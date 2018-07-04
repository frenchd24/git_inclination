#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_azimuth_map.py, v 1.0 07/01/18


Based on:
$Id:  plotAzimuthMap2.py, v 2.1 1/20/17

Plot absorption properties on a map showing their distribution around a central galaxy

v2: Separate blue and red into separate plots, so now it's a plane of 6. Size points by 
    EW (10/06/2016)
    
    
v2.1: Use open symbols for rho <= R_vir to be consistent with the other plots. -> /plots6/
    (1/20/17)

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
from matplotlib.patches import Ellipse


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

fontScale = 18
rc('text', usetex=True)
rc('font', size=18, family='serif', weight='normal')
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

    
def main():
    # plot dopplar b parameter as a function of impact w/o splitting, add other sets if you want
    plot_az_map_plus = True
    plot_az_map_plus_save = True

    
    # plot_number = 1 for just the isolated sample, =2 adds the associated, =3 adds two+
    # =4 adds groups with 2 or more members
    plot_number = 1

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
        
        isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/isolated8.p'
        L_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_isolated8.p'
        L_associated_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated8.p'
        L_associated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated8.p'
        L_nonassociated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated8.p'
        L_two_filename = '/Users/frenchd/Research/inclination/git_inclination/L_two8.p'
        L_three_plus_filename = '/Users/frenchd/Research/inclination/git_inclination/L_three_plus8.p'
        L_group_filename = '/Users/frenchd/Research/inclination/git_inclination/L_group8.p'

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
    L_three_plus_file = open(L_three_plus_filename,'r')
    L_group_file = open(L_group_filename,'r')

    # unload the data from them
    isolated = pickle.load(isolated_file)
    L_isolated = pickle.load(L_isolated_file)
    L_associated_isolated = pickle.load(L_associated_isolated_file)
    L_associated = pickle.load(L_associated_file)
    L_nonassociated = pickle.load(L_nonassociated_file)
    L_two = pickle.load(L_two_file)
    L_three_plus = pickle.load(L_three_plus_file)
    L_group = pickle.load(L_group_file)
    
    # close the files
    isolated_file.close()
    L_isolated_file.close()
    L_associated_isolated_file.close()
    L_associated_file.close()
    L_nonassociated_file.close()
    L_two_file.close()
    L_three_plus_file.close()
    L_group_file.close()
    
    # which dataset to use for plotting?
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


########################################################################################
    
    if plot_az_map_plus:
        
        # azimuths for each inclination interval
        azBlue_1 = []
        azRed_1 = []
        azBlue_2 = []
        azRed_2 = []
        azBlue_3 = []
        azRed_3 = []
        
        # scaled sizes for EW-scaled markers
        size_1b = []
        size_1r = []
        size_2b = []
        size_2r = []
        size_3b = []
        size_3r = []
        
        # new coordinates of absorption for plotting (w/ galaxy at (0,0))
        yList_1r = []
        xList_1r = []
        yList_2r = []
        xList_2r = []
        yList_3r = []
        xList_3r = []
        
        yList_1b = []
        xList_1b = []
        yList_2b = []
        xList_2b = []
        yList_3b = []
        xList_3b = []
        
        # save rho/R_vir also
        rR_1r = []
        rR_2r = []
        rR_3r = []
        
        rR_1b = []
        rR_2b = []
        rR_3b = []
        
        
        L_associated_Lya_vs = L_associated['Lya_vs']
        L_associated_Lya_Ws = L_associated['Lya_Ws']
        L_associated_impacts = L_associated['impacts']
        L_associated_azimuths = L_associated['azimuths']
        L_associated_adjustedIncs = L_associated['adjustedIncs']
        L_associated_R_virs = L_associated['R_virs']
        L_associated_Vhels = L_associated['Vhels']

#         impacts2 = list(impacts) + list(L_associated_impacts)
#         azimuths2 = list(azimuths) + list(L_associated_azimuths)
#         adjustedIncs2 = list(adjustedIncs) + list(L_associated_adjustedIncs)
#         R_virs2 = list(R_virs) + list(L_associated_R_virs)
#         Lya_vs2 = list(Lya_vs) + list(L_associated_Lya_vs)
#         Lya_Ws2 = list(Lya_Ws) + list(L_associated_Lya_Ws)
#         Vhels2 = list(Vhels) + list(L_associated_Vhels)

        impacts2 = list(impacts)
        azimuths2 = list(azimuths)
        adjustedIncs2 = list(adjustedIncs)
        R_virs2 = list(R_virs)
        Lya_vs2 = list(Lya_vs)
        Lya_Ws2 = list(Lya_Ws)
        Vhels2 = list(Vhels)
        
        
        # set a minimum EW limit if you want
        EW_min = 100.

        impacts3 = []
        azimuths3 = []
        adjustedIncs3 = []
        R_virs3 = []
        Lya_vs3 = []
        Vhels3 = []
        Lya_Ws3 = []
        for i, a, fInc, vir, Lya_v, Vhel, ew in zip(impacts2, azimuths2, adjustedIncs2, R_virs2, Lya_vs2, Vhels2, Lya_Ws2):
            if float(ew) >= EW_min:
                impacts3.append(i)
                azimuths3.append(a)
                adjustedIncs3.append(fInc)
                R_virs3.append(vir)
                Lya_vs3.append(Lya_v)
                Vhels3.append(Vhel)
                Lya_Ws3.append(ew)
        
                
        largestEW = max(Lya_Ws3)
        smallestEW = min(Lya_Ws3)
        maxSize = 600
        minSize = 20
                
                
        # calculate the position on the sky for each absorption feature wrt to the galaxy
        for i, a, fInc, vir, Lya_v, Vhel, ew in zip(impacts3, azimuths3, adjustedIncs3, R_virs3, Lya_vs3, Vhels3, Lya_Ws3):
            dif = float(Lya_v) - float(Vhel)
        
            if float(a) >= 0. and float(fInc)>=0. and float(fInc)<=90.:
            
                # y coordinate
                y = (float(i)/float(vir)) * sin((a*pi)/180.)
            
                # x coordinate
                x = (float(i)/float(vir)) * cos(a*pi/180.)
                
                # rho / Rvir
                rhoRvir = float(i)/float(vir)
                
                # new size for the marker point
                newSize = ((float(ew) - smallestEW)/(largestEW - smallestEW)) * (maxSize - minSize) + minSize
                    
                if fInc <=40:
                    if dif >0:
                        #blue absorber
                        azBlue_1.append(a)
                        xList_1b.append(x)
                        yList_1b.append(y)
                        size_1b.append(newSize)
                        rR_1b.append(rhoRvir)

                    else:
                        azRed_1.append(a)
                        xList_1r.append(x)
                        yList_1r.append(y)
                        size_1r.append(newSize)
                        rR_1r.append(rhoRvir)

                if fInc > 40 and fInc <=65:
                    if dif >0:
                        #blue absorber
                        azBlue_2.append(a)
                        yList_2b.append(y)
                        xList_2b.append(x)
                        size_2b.append(newSize)
                        rR_2b.append(rhoRvir)

                    else:
                        azRed_2.append(a)
                        yList_2r.append(y)
                        xList_2r.append(x)
                        size_2r.append(newSize)
                        rR_2r.append(rhoRvir)

                if fInc > 65:
                    if dif >0:
                        #blue absorber
                        azBlue_3.append(a)
                        yList_3b.append(y)
                        xList_3b.append(x)
                        size_3b.append(newSize)
                        rR_3b.append(rhoRvir)

                    else:
                        azRed_3.append(a)
                        yList_3r.append(y)
                        xList_3r.append(x)
                        size_3r.append(newSize)
                        rR_3r.append(rhoRvir)

            else:
                print 'float(a) <0: ',r,d,i,a,fInc

        # calculate the average red vs blue azimuth line
        blueAvg1 = mean(azBlue_1)
        print 'blueAvg1: ',blueAvg1
        redAvg1 = mean(azRed_1)
        print 'redAvg1: ',redAvg1
        print
        blueAvg2 = mean(azBlue_2)
        print 'blueAvg2: ',blueAvg2
        redAvg2 = mean(azRed_2)
        print 'redAvg2: ',redAvg2
        print
        blueAvg3 = mean(azBlue_3)
        print 'blueAvg3: ',blueAvg3
        redAvg3 = mean(azRed_3)
        print 'redAvg3: ',redAvg3
        print
        
        xyBlueAvg1 = (500.* cos(blueAvg1 * pi/180.), 500.* sin(blueAvg1 * pi/180.))
        xyRedAvg1 = (500.* cos(redAvg1 * pi/180.), 500.* sin(redAvg1 * pi/180.))
        print 'xyBlueAvg1: ',xyBlueAvg1
        print 'xyRedAvg1: ',xyRedAvg1
        xyBlueAvg2 = (500.* cos(blueAvg2 * pi/180.), 500.* sin(blueAvg2 * pi/180.))
        xyRedAvg2 = (500.* cos(redAvg2 * pi/180.), 500.* sin(redAvg2 * pi/180.))
        print 'xyBlueAvg2: ',xyBlueAvg2
        print 'xyRedAvg2: ',xyRedAvg2
        xyBlueAvg3 = (500.* cos(blueAvg3 * pi/180.), 500.* sin(blueAvg3 * pi/180.))
        xyRedAvg3 = (500.* cos(redAvg3 * pi/180.), 500.* sin(redAvg3 * pi/180.))
        print 'xyBlueAvg3: ',xyBlueAvg3
        print 'xyRedAvg3: ',xyRedAvg3

##########################################################################################

        # plot the distributions
        fig = figure(figsize=(14.5,11))
        alpha = 0.6
        bSymbol = 'D'
        rSymbol = 'o'
    
        ax1 = fig.add_subplot(2,3,1)
        subplots_adjust(left=None,bottom=None,right=None,top=None,wspace=0.001,hspace=0.001)

        # Low inclinations = circular mostly
        e = Ellipse(xy=(0,0), width=1.0, height=0.9, angle=0)
#         legend(scatterpoints=1,prop={'size':15},loc=2,fancybox=True)
                    
        # no transparency
        e.set_alpha(0.3)
        ax1.add_artist(e)
        e.set_facecolor('black')
        e.set_edgecolor('black')
        label1 = r'$\rm Inc \leq 40$'
        
#         y coordinate
#         y = (float(i)/float(vir)) * sin((a*pi)/180.)
#     
#         x coordinate
#         x = (float(i)/float(vir)) * cos(a*pi/180.)
        
        for x, y, a, s,r in zip(xList_1b,yList_1b,azBlue_1,size_1b,rR_1b):
            # rho/R_vir = y / sin((a*pi)/180.) as above
            rhoRvir = y / sin((a*pi)/180.)
            print 'rhoRvir: {0} vs rB_: {1} '.format(rhoRvir,r)
            
            if r <=1:
                fc = 'none'
                ec = color_blue
                ax1.scatter(x,y,c=color_blue,marker=bSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)
            else:
                fc = color_blue
                ec = 'black'
                ax1.scatter(x,y,c=color_blue,marker=bSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)

#         ax1.scatter(xList_1r,yList_1r,c=color_red,marker=rSymbol,alpha=alpha,s=50)
    
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', alpha=1, facecolor='none')

        # place a text box in upper right in axes coords
        ax1.text(0.715, 0.93, label1, transform=ax1.transAxes, fontsize=15,verticalalignment='top', bbox=props)

        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax1.yaxis.set_major_locator(majorLocator)
        ax1.yaxis.set_major_formatter(majorFormatter)
        ax1.yaxis.set_minor_locator(minorLocator)
        plt.setp(ax1.get_xticklabels(), visible=False)
    
        # y axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax1.yaxis.set_major_locator(majorLocator)
        ax1.yaxis.set_major_formatter(majorFormatter)
        ax1.yaxis.set_minor_locator(minorLocator)

        xlim(0,2.0)
        ylim(0,2.5)
        ylabel(r'$\rm \rho / R_{vir}$')

#         ax1.get_yaxis().set_tick_params(which='both', direction='out')

##########################################################################################
        ax2 = fig.add_subplot(2,3,2,sharey=ax1)
        
        # medium inclination = ellipse
        e = Ellipse(xy=(0,0), width=1.0, height=0.6, angle=0)
                    
        # no transparency
        e.set_alpha(0.3)
        ax2.add_artist(e)
        e.set_facecolor('black')
        e.set_edgecolor('black')
        label2 = r'$\rm 40 < Inc \leq 65$'
        
        for x, y, a, s, r in zip(xList_2b,yList_2b,azBlue_2,size_2b,rR_2b):
            # rho/R_vir = y / sin((a*pi)/180.) as above
            rhoRvir = y / sin((a*pi)/180.)
            print 'rhoRvir: {0} vs rB_: {1} '.format(rhoRvir,r)
            
            if r <=1:
                fc = 'none'
                ec = color_blue
                ax2.scatter(x,y,c=color_blue,marker=bSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)
            else:
                fc = color_blue
                ec = 'black'
                ax2.scatter(x,y,c=color_blue,marker=bSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)
        
#         ax2.scatter(xList_2b,yList_2b,c=color_blue,marker=bSymbol,alpha=alpha,s=size_2b)
    
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', alpha=1, facecolor='none')

        # place a text box in upper right in axes coords
        ax2.text(0.58, 0.93, label2, transform=ax2.transAxes, fontsize=15,verticalalignment='top', bbox=props)
    
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax2.yaxis.set_major_locator(majorLocator)
        ax2.yaxis.set_major_formatter(majorFormatter)
        ax2.yaxis.set_minor_locator(minorLocator)
        plt.setp(ax2.get_xticklabels(), visible=False)

        # y axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax2.yaxis.set_major_locator(majorLocator)
        ax2.yaxis.set_major_formatter(majorFormatter)
        ax2.yaxis.set_minor_locator(minorLocator)
        plt.setp(ax2.get_xticklabels(), visible=False)

        xlim(0,2.0)
        ylim(0,2.5)
        plt.setp(ax2.get_yticklabels(), visible=False)
        xlabel(r'$\rm \rho / R_{vir}$')

#         ax2.yaxis.tick_right()
#         ax2.get_yaxis().set_tick_params(which='both', direction='out')

##########################################################################################
        ax3 = fig.add_subplot(2,3,3,sharey=ax1)
        
        # plot the flat galaxy line
        e = Ellipse(xy=(0,0), width=1.0, height=0.3, angle=0)
                    
        # no transparency
        e.set_alpha(0.3)
        ax3.add_artist(e)
        e.set_facecolor('black')
        e.set_edgecolor('black')
        label3 = r'$\rm Inc > 65$'
        
        for x, y, a, s, r in zip(xList_3b,yList_3b,azBlue_3,size_3b,rR_3b):
            # rho/R_vir = y / sin((a*pi)/180.) as above
            rhoRvir = y / sin((a*pi)/180.)
            print 'rhoRvir: {0} vs rB_: {1} '.format(rhoRvir,r)
            
            if r <=1:
                fc = 'none'
                ec = color_blue
                ax3.scatter(x,y,c=color_blue,marker=bSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)
            else:
                fc = color_blue
                ec = 'black'
                ax3.scatter(x,y,c=color_blue,marker=bSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)
        
#         ax3.scatter(xList_3b,yList_3b,c=color_blue,marker=bSymbol,alpha=alpha,s=size_3b)
    
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', alpha=1, facecolor='none')

        # place a text box in upper right in axes coords
        ax3.text(0.715, 0.93, label3, transform=ax3.transAxes, fontsize=15,verticalalignment='top', bbox=props)
    
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax3.yaxis.set_major_locator(majorLocator)
        ax3.yaxis.set_major_formatter(majorFormatter)
        ax3.yaxis.set_minor_locator(minorLocator)
        plt.setp(ax3.get_xticklabels(), visible=False)

        # y axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax3.yaxis.set_major_locator(majorLocator)
        ax3.yaxis.set_major_formatter(majorFormatter)
        ax3.yaxis.set_minor_locator(minorLocator)
        plt.setp(ax3.get_yticklabels(), visible=False)

        xlim(0,2.0)
        ylim(0,2.5)
#             xlabel(r'$\rm \rho / R_{vir}$')


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

        ax4 = fig.add_subplot(2,3,4,sharex=ax1)

        # Low inclinations = circular mostly
        e = Ellipse(xy=(0,0), width=1.0, height=0.9, angle=0)
#         legend(scatterpoints=1,prop={'size':15},loc=2,fancybox=True)
                    
        # no transparency
        e.set_alpha(0.3)
        ax4.add_artist(e)
        e.set_facecolor('black')
        e.set_edgecolor('black')
        label1 = r'$\rm Inc \leq 40$'

        for x, y, a, s, r in zip(xList_1r,yList_1r,azRed_1,size_1r,rR_1r):
            # rho/R_vir = y / sin((a*pi)/180.) as above
            rhoRvir = y / sin((a*pi)/180.)
            print 'rhoRvir: {0} vs rB_: {1} '.format(rhoRvir,r)
            
            if r <=1:
                fc = 'none'
                ec = color_red
                ax4.scatter(x,y,c=color_red,marker=rSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)
            else:
                fc = color_red
                ec = 'black'
                ax4.scatter(x,y,c=color_red,marker=rSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)

#         ax4.scatter(xList_1r,yList_1r,c=color_red,marker=rSymbol,alpha=alpha,s=size_1r)
    
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', alpha=1, facecolor='none')

        # place a text box in upper right in axes coords
        ax4.text(0.715, 0.93, label1, transform=ax4.transAxes, fontsize=15,verticalalignment='top', bbox=props)

        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax4.yaxis.set_major_locator(majorLocator)
        ax4.yaxis.set_major_formatter(majorFormatter)
        ax4.yaxis.set_minor_locator(minorLocator)
        ax4.set_xticks([0.0,0.5,1.0,1.5])
    
        # y axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax4.yaxis.set_major_locator(majorLocator)
        ax4.yaxis.set_major_formatter(majorFormatter)
        ax4.yaxis.set_minor_locator(minorLocator)
        ax4.set_yticks([0.0,0.5,1.0,1.5,2.0])

        xlim(0,2.0)
        ylim(0,2.5)
        ylabel(r'$\rm \rho / R_{vir}$')

#         ax1.get_yaxis().set_tick_params(which='both', direction='out')

##########################################################################################
        ax5 = fig.add_subplot(2,3,5,sharey=ax4,sharex=ax2)
        
        # medium inclination = ellipse
        e = Ellipse(xy=(0,0), width=1.0, height=0.6, angle=0)
                    
        # no transparency
        e.set_alpha(0.3)
        ax5.add_artist(e)
        e.set_facecolor('black')
        e.set_edgecolor('black')
        label5 = r'$\rm 40 < Inc \leq 65$'
#         ax5.scatter(xList_2b,yList_2b,c=color_blue,marker=bSymbol,alpha=alpha,s=50)

        for x, y, a, s, r in zip(xList_2r,yList_2r,azRed_2,size_2r,rR_2r):
            # rho/R_vir = y / sin((a*pi)/180.) as above
            rhoRvir = y / sin((a*pi)/180.)
            print 'rhoRvir: {0} vs rB_: {1} '.format(rhoRvir,r)
            
            if rhoRvir <=1:
                fc = 'none'
                ec = color_red
                ax5.scatter(x,y,c=color_red,marker=rSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)
            else:
                fc = color_red
                ec = 'black'
                ax5.scatter(x,y,c=color_red,marker=rSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)

#         ax5.scatter(xList_2r,yList_2r,c=color_red,marker=rSymbol,alpha=alpha,s=size_2r)
    
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', alpha=1, facecolor='none')

        # place a text box in upper right in axes coords
        ax5.text(0.58, 0.93, label5, transform=ax5.transAxes, fontsize=15,verticalalignment='top', bbox=props)
    
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax5.yaxis.set_major_locator(majorLocator)
        ax5.yaxis.set_major_formatter(majorFormatter)
        ax5.yaxis.set_minor_locator(minorLocator)
        ax5.set_xticks([0.0,0.5,1.0,1.5])

        # y axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax5.yaxis.set_major_locator(majorLocator)
        ax5.yaxis.set_major_formatter(majorFormatter)
        ax5.yaxis.set_minor_locator(minorLocator)
        plt.setp(ax5.get_yticklabels(), visible=False)
        ax5.set_yticks([0.0,0.5,1.0,1.5,2.0])

        xlim(0,2.0)
        ylim(0,2.5)
        xlabel(r'$\rm \rho / R_{vir}$')

#         ax2.yaxis.tick_right()
#         ax2.get_yaxis().set_tick_params(which='both', direction='out')

##########################################################################################
        ax6 = fig.add_subplot(2,3,6,sharex=ax3,sharey=ax5)
        # plot the flat galaxy line
        e = Ellipse(xy=(0,0), width=1.0, height=0.3, angle=0)
                    
        # no transparency
        e.set_alpha(0.3)
        ax6.add_artist(e)
        e.set_facecolor('black')
        e.set_edgecolor('black')
        label6 = r'$\rm Inc > 65$'
#         ax6.scatter(xList_3b,yList_3b,c=color_blue,marker=bSymbol,alpha=alpha,s=50)

        for x, y, a, s, r in zip(xList_3r,yList_3r,azRed_3,size_3r,rR_3r):
            # rho/R_vir = y / sin((a*pi)/180.) as above
            rhoRvir = y / sin((a*pi)/180.)
            print 'rhoRvir: {0} vs rB_: {1} '.format(rhoRvir,r)

            if r <=1:
                fc = 'none'
                ec = color_red
                ax6.scatter(x,y,c=color_red,marker=rSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)
            else:
                fc = color_red
                ec = 'black'
                ax6.scatter(x,y,c=color_red,marker=rSymbol,alpha=alpha,s=s,\
                facecolor=fc,edgecolor=ec)

#         ax6.scatter(xList_3r,yList_3r,c=color_red,marker=rSymbol,alpha=alpha,s=size_3r)
    
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', alpha=1, facecolor='none')

        # place a text box in upper right in axes coords
        ax6.text(0.715, 0.93, label6, transform=ax6.transAxes, fontsize=15,verticalalignment='top', bbox=props)
    
        # x-axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax6.yaxis.set_major_locator(majorLocator)
        ax6.yaxis.set_major_formatter(majorFormatter)
        ax6.yaxis.set_minor_locator(minorLocator)
        ax6.set_xticks([0.0,0.5,1.0,1.5,2.0])
    
        # y axis
        majorLocator   = MultipleLocator(0.5)
        majorFormatter = FormatStrFormatter(r'$\rm %0.1f$')
        minorLocator   = MultipleLocator(0.25)
        ax6.yaxis.set_major_locator(majorLocator)
        ax6.yaxis.set_major_formatter(majorFormatter)
        ax6.yaxis.set_minor_locator(minorLocator)
        ax6.set_yticks([0.0,0.5,1.0,1.5,2.0])
        plt.setp(ax6.get_yticklabels(), visible=False)

        xlim(0,2.0)
        ylim(0,2.5)
#             xlabel(r'$\rm \rho / R_{vir}$')

        if plot_az_map_plus_save:
            savefig('{0}/azimuthMap_separate2_isolated_EWmin_{1}.pdf'.format(saveDirectory, EW_min),\
            format='pdf',bbox_inches='tight')
        else:
            show()


###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    