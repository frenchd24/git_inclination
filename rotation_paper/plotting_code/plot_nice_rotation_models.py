#!/Users/frenchd/anaconda2/bin/python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_nice_rotation_models.py, v1.0 07/12/18

Plot NFW and regular rotation models all in one in a publication quality plot

v1.1:
Re-plot with inward-facing tick marks

'''


import sys
import os
import csv
import time
from copy import deepcopy


from pylab import *
# import atpy
# from math import *
from utilities import *
from scipy import stats
import getpass
import math
import pickle
import json
import io
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

from scipy.optimize import curve_fit
import matplotlib.ticker as ticker

from scipy.interpolate import interp1d
from scipy import interpolate
from scipy import interpolate, optimize
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline


from matplotlib import rc
fontScale = 18
rc('text', usetex=True)
rc('font', size=fontScale,family='serif',weight='medium')
rc('xtick.major',size=8,width=0.6)
rc('xtick.minor',size=5,width=0.6)
rc('ytick.major',size=8,width=0.6)
rc('ytick.minor',size=5,width=0.6)
rc('xtick',labelsize = fontScale)
rc('ytick',labelsize = fontScale)
rc('axes',labelsize = fontScale)
rc('xtick',labelsize = fontScale)
rc('ytick',labelsize = fontScale)
# rc('font', weight = 450)
# rc('axes',labelweight = 'bold')
rc('axes',linewidth = 1,labelweight='normal')
rc('axes',titlesize='small')
rc('xtick',direction='in')
rc('ytick',direction='in')


'''
========================================================
'''



def NFW(r,v200,c,r200):

    x = r/r200
    top = (np.log(1 + c*x) - c*x / (1 + c*x))

    bottom = x * (np.log(1 + c) - c / (1 + c))

    vr = v200 * np.sqrt(top/bottom)

    return vr



def plot_cylinder(p0,p1,R):
    #   I totally jacked this code from here:
    #
    #   '''
    #     Created on Sun Oct  2 18:33:10 2016
    # 
    #     Modified from https://stackoverflow.com/questions/38076682/how-to-add-colors-to-each-individual-face-of-a-cylinder-using-matplotlib
    #     to add "end caps" and to undo fancy coloring.
    # 
    #     @author: astrokeat
    #   '''

    #axis and radius
#     p0 = np.array([1, 3, 2]) #point at one end
#     p1 = np.array([8, 5, 9]) #point at other end
#     R = 5

    #vector in direction of axis
    v = p1 - p0

    #find magnitude of vector
    mag = norm(v)

    #unit vector in direction of axis
    v = v / mag

    #make some vector not in the same direction as v
    not_v = np.array([1, 0, 0])
    if (v == not_v).all():
        not_v = np.array([0, 1, 0])

    #make vector perpendicular to v
    n1 = np.cross(v, not_v)
    #normalize n1
    n1 /= norm(n1)

    #make unit vector perpendicular to v and n1
    n2 = np.cross(v, n1)

    #surface ranges over t from 0 to length of axis and 0 to 2*pi
    t = np.linspace(0, mag, 2)
    theta = np.linspace(0, 2 * np.pi, 100)
    rsample = np.linspace(0, R, 2)

    #use meshgrid to make 2d arrays
    t, theta2 = np.meshgrid(t, theta)

    rsample,theta = np.meshgrid(rsample, theta)

    #generate coordinates for surface
    # "Tube"
    X, Y, Z = [p0[i] + v[i] * t + R * np.sin(theta2) * n1[i] + R * np.cos(theta2) *       n2[i] for i in [0, 1, 2]]
    # "Bottom"
    X2, Y2, Z2 = [p0[i] + rsample[i] * np.sin(theta) * n1[i] + rsample[i] * np.cos(theta) * n2[i] for i in [0, 1, 2]]
    # "Top"
    X3, Y3, Z3 = [p0[i] + v[i]*mag + rsample[i] * np.sin(theta) * n1[i] + rsample[i] * np.cos(theta) * n2[i] for i in [0, 1, 2]]

    return (X, Y, Z), (X2, Y2, Z2), (X3, Y3, Z3)




def main():
    # open the data files
    directory = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/SALT/NGC3633/'

    cyl_model_filename = '{0}/NGC3633-RX_J1121.2+0326_model_NFWFalse.p'.format(directory)
    NFW_model_filename = '{0}/NGC3633-RX_J1121.2+0326_model_NFWTrue.p'.format(directory)

    cyl_model_file = open(cyl_model_filename,'r')
    cyl_model = pickle.load(cyl_model_file)

    NFW_model_file = open(NFW_model_filename,'r')
    NFW_model = pickle.load(NFW_model_file)

    cyl_model_file.close()
    NFW_model_file.close()


    # which thing to plot?
    plot_velocity = False
    plot_3Dmodel = False
    plot_3Dmodel_movie = False
    
    plot_NFW_fit = True
    plot_rotation_curve = True
    
    movie_directory = '/Users/frenchd/Research/inclination/git_inclination/thesis/DMF_thesis/NGC3633_movie/'

    
##########################################################################################
    # get the data from the JSON files
##########################################################################################
    # fields in JSON file are:
    #
    # 'name': galaxy name
    # 'vsys_published': Vhel as found on NED
    # 'dvsys_published': error in Vhel as found in NED
    # 'inclination': inclination
    # 'di': inclination error
    # 'centerTrace': center of spectrum
    # 'distance': best distance to galaxy
    # 'vsys_measured': measured Vhel
    # 'vsys_measured_err': measured error in Vhel
    # 'left_vrot_incCorrected_avg': average left wing velocity, corrected for inc
    # 'left_vrot_incCorrected_avg_err': error in average left wing velocity corrected for inc
    # 'right_vrot_incCorrected_avg': average right wing velocity, corrected for inc
    # 'right_vrot_incCorrected_avg_err':  error in average right wing velocity corrected for inc
    # 'left_vrot_avg': average left wing velocity (redshift subtracted)
    # 'left_vrot_avg_err': error in average left wing velocity (redshift subtracted)
    # 'right_vrot_avg': average right wing velocity (redshift subtracted)
    # 'right_vrot_avg_err': error in average right wing velocity (redshift subtracted)
    # 'vrot_vals': observed velocities (redshift but not inc corrected)
    # 'vrot_errs': errors in observed velocities (redshift but not inc corrected)
    # 'vrot_incCorrected_vals': inclination corrected velocities (redshift subtracted)
    # 'vrot_incCorrected_errs': errors in inclination corrected velocities
    # 'vrot_observed': observed velocities (Vhel + rotation)
    # 'agn': include any information about AGN here
    # 'xVals': physical (kpc) x axis along the slit

    directory = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/rot_curves/'
    galaxyName = 'NGC3633'
    
#     filename = 'CGCG039-137-summary4.json'
#     filename = 'ESO343-G014-summary4.json'
#     filename = 'RFGC3781-summary4.json'
#     filename = 'IC5325-summary4.json'
#     filename = 'MCG-03-58-009-summary4.json'
#     filename = 'NGC1566-summary4.json'
#     filename = 'NGC3513-summary4.json'
#     filename = 'NGC3633-summary4.json'
#     filename = 'NGC4536-summary4.json'
#     filename = 'NGC4939-summary4.json'
#     filename = 'NGC5364-summary4.json'

    filename = '{0}-summary4.json'.format(galaxyName)

    
    with open(directory+filename) as data_file:
        data = json.load(data_file)    
        
        vrot_vals = data['vrot_vals']
        vrot_incCorrected_vals = data['vrot_incCorrected_vals']
        vrot_incCorrected_errs = data['vrot_incCorrected_errs']
        
        right_vrot_incCorrected_avg = data['right_vrot_incCorrected_avg']
        right_vrot_incCorrected_avg_err = data['right_vrot_incCorrected_avg_err']

        left_vrot_incCorrected_avg = data['left_vrot_incCorrected_avg']
        left_vrot_incCorrected_avg_err = data['left_vrot_incCorrected_avg_err']


        xVals = data['xVals']
        inc = data['inclination']
        vsys_measured = data['vsys_measured']
        vsys_measured_err = data['vsys_measured_err']

        
        RA_galaxy = data['RAdeg']
        Dec_galaxy = data['DEdeg']
        dist = data['dist']
        majDiam = data['majDiam']
        inc = data['inclination']
        PA = data['PA']
        agn = data['agn']
        
        R_vir = calculateVirialRadius(majDiam)





##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
    if plot_velocity:

        # get the data
        cyl_model_intersect_v_list = cyl_model['intersect_v_list']
        cyl_model_v_proj_list = cyl_model['v_proj_list']

        NFW_model_intersect_v_list = NFW_model['intersect_v_list']
        NFW_model_v_proj_list = NFW_model['v_proj_list']
        

        # do the plotting
#         fig = plt.figure(figsize=(6.7,7.7))
        fig = plt.figure(figsize=(8.0, 6.7))

        ax = fig.add_subplot(1,1,1)


        # first cylindrical
        ax.plot(cyl_model_intersect_v_list,
                    cyl_model_v_proj_list, 
                    color='black',
                    lw=2,
                    linestyle='dashed',
                    markersize=0,
                    label=r'$\rm Cylindrical$')
           
        # next NFW     
        ax.plot(NFW_model_intersect_v_list, 
                    NFW_model_v_proj_list, 
                    color='green',
                    lw=2,
                    markersize=0,
                    label = r'$\rm NFW$')
                
                
        ylabel(r'$\rm Projected~Rotation~Vel. ~[km~s^{-1}]$')
        xlabel(r'$\rm Vel.~Along~Sightline ~[km~s^{-1}]$')
            
        # x-axis
        majorLocator   = MultipleLocator(10)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
    
        # y-axis
        majorLocator   = MultipleLocator(20)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(10)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
    
    
        leg = ax.legend(scatterpoints=1,prop={'size':14},loc='lower right',fancybox=True)
    #         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(-160, 0)
        xlim(-30, 30)

        savefig('{0}/NGC3633-RX_J1121.2+0326_model_plot3.pdf'.format(directory),format='pdf',bbox_inches='tight')


##########################################################################################
##########################################################################################

    if plot_NFW_fit:
        # do the plotting
        fig = plt.figure(figsize=(8.0, 6.7))
        fig = plt.figure(figsize=(7.7,5.7))
        ax = fig.add_subplot(1,1,1)
        
        
        '''
        This function makes a nice looking plot showing the NFW fit and data
    
        xData - the x values for the observed rotation curve 
        yData - the y values for the observed rotation curve
        popt - the fit values
        x_lim - the maximum x value to plot to
    
        returns the figure object, to be shown or saved
        '''

        xData = []
        yData = []
        yErrs = []
        yData_abs = []
        xData_abs = []
        
        
        for v, e, x in zip(vrot_incCorrected_vals, vrot_incCorrected_errs, xVals):
            xData.append(x)
            yData.append(v)
            yErrs.append(e)
            
            
            # fold it over for NFW plotting
            xData_abs.append(abs(x))
            yData_abs.append(abs(v))
            
        # turn them into numpy arrays
        xData_abs = np.array(xData_abs)
        yData_abs = np.array(yData_abs)
        yErrs = np.array(yErrs)

    
        v200 = 111.91
        c = 21.4
        r200 = 60.03
    
        x_fit = linspace(0, round(R_vir,0)+10, num=1000)
    
        scatter(xData, yData, color='black', s=40, lw=0, label = r'$\rm Data$')
    #     plot(x_fit, NFW(x_fit, *popt), 'r-',label='fit: a={0}, rho={1}'.format(*popt))
        plot(x_fit, NFW(x_fit, *popt), 'r-', color='green', alpha=0.7, lw=2, \
        label='NFW Fit: V200={0}, c={1}, R200={2}'.format(round(v200,2),round(c,2),round(r200,2)))
    
        legend(scatterpoints=1,prop={'size':14},loc='lower right',fancybox=True)
    
        xlim(0, x_lim)
        ylim(0, round(np.nanmax(yData),-1) + 15)

        # x-axis
        majorLocator   = MultipleLocator(25)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(5)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        # y axis
        majorLocator   = MultipleLocator(50)
        majorFormatter = FormatStrFormatter(r'$\rm %d$')
        minorLocator   = MultipleLocator(10)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)

        xlabel(r'$\rm R ~[kpc]$')
        ylabel(r'$\rm \emph{v}_{{rot}} ~[km s^{{-1}}]$')
    
    
        leg = ax.legend(scatterpoints=1,prop={'size':14},loc='lower right',fancybox=True)
    #         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(-160, 0)
        xlim(-30, 30)

        savefig('{0}/NGC3633-NFW_fit_Rvir10.pdf'.format(directory),format='pdf',bbox_inches='tight')


##########################################################################################
##########################################################################################
##########################################################################################


    # plot the model figure cylinder thing
    if plot_3Dmodel:
        # first get the data
        z = NFW_model['z_sightline']
        y = NFW_model['y_sightline']
        x = NFW_model['x_sightline']
    
        # now the cylinder
        intersect = NFW_model['intersect']
        R = NFW_model['R']
        p0 = NFW_model['p0']
        p1 = NFW_model['p1']

        fig = plt.figure(figsize=(7.0,8.7))
        ax = fig.add_subplot(1,1,1,projection='3d')
        
    
        print 'x: ',x
        print 'y: ',y
        print 'z: ',z
        print
        
        plotExtent = 600
        
        # plot the sightline
        ax.plot(x, y, z, color='black', alpha = 0.6, lw=2)

        # intersect point
        print 'intersect: ',intersect[0],intersect[1],intersect[2]

    #     ax.plot([0,v[0]], [0,v[1]], [0,v[2]], color='green',lw=plotExtent/100)
    #     ax.plot([0,v_90[0]], [0,v_90[1]], [0,v_90[2]], color='purple',lw=plotExtent/100)
    #     ax.plot([intersect[0],v_90[0]], [intersect[1],v_90[1]], [intersect[2],v_90[2]], color='purple',lw=plotExtent/100)


        # put a star on the intersect
    #         planePoint_end = [-1.18639357e-01, 6.80095455e+02, -2.46470324e+02]
    #         planePoint_end2 = [1.18630006e-01, -6.79357841e+02, 2.48330210e+02]
    #         ax.plot([planePoint_end[0]],[planePoint_end[1]],[planePoint_end[2]],color='red',marker='*',lw=0)
    #         ax.plot([planePoint_end2[0]],[planePoint_end2[1]],[planePoint_end2[2]],color='green',marker='*',lw=0)

        ax.plot([intersect[0]],[intersect[1]],[intersect[2]],color='red',marker='*',lw=0)
        
    ##########################################################################################
    ##########################################################################################
        # plot the cylinder

        tube,bottom,top = plot_cylinder(p0,p1,R)

        X, Y, Z = tube
        X2, Y2, Z2 = bottom
        X3, Y3, Z3 = top

        alphaTube = 0.15
        alphaBottom = 0.3
        alphaTop = 0.5
    
        colorTube = 'blue'
        colorBottom = 'red'
        colorTop = 'blue'

    #     ax=plt.subplot(111, projection='3d')
        ax.plot_surface(X, Y, Z, color=colorTube, alpha = alphaTube)
        ax.plot_surface(X2, Y2, Z2, color=colorBottom, alpha = alphaBottom)
        ax.plot_surface(X3, Y3, Z3, color=colorTop, alpha = alphaTop)
    
        ax.set_xlim(-plotExtent, plotExtent)
        ax.set_ylim(-plotExtent, plotExtent)
        ax.set_zlim(-plotExtent, plotExtent)

        # reverse the RA axis so negative is on the right
    #     ax = plt.gca()
        ax.invert_xaxis()

        # rotate the plot
    #         ax.view_init(elev=10., azim=5)
    #         ax.view_init(elev=15., azim=20)
        ax.view_init(elev=5., azim=5)

        yticks((-600, -400, -200, 0, 200, 400, 600), (600, 400, 200, 0, -200, -400, -600))
#         xticks((-600, -400, -200, 0, 200, 400, 600), (-600, -400, -200, 0, 200, 400, 600))
        xticks((-600, -300, 0, 300, 600), (-600, -300, 0, 300, 600))

    
        [t.set_va('center') for t in ax.get_yticklabels()]
        [t.set_ha('center') for t in ax.get_yticklabels()]
    
        [t.set_va('bottom') for t in ax.get_xticklabels()]
        [t.set_ha('left') for t in ax.get_xticklabels()]
    
        [t.set_va('center') for t in ax.get_zticklabels()]
        [t.set_ha('right') for t in ax.get_zticklabels()]
    
        ax.grid(True)
        ax.xaxis.pane.set_edgecolor('black')
        ax.yaxis.pane.set_edgecolor('black')
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False

        ax.xaxis._axinfo['tick']['inward_factor'] = 0
        ax.xaxis._axinfo['tick']['outward_factor'] = 0.6
        ax.yaxis._axinfo['tick']['inward_factor'] = 0
        ax.yaxis._axinfo['tick']['outward_factor'] = 0.6
        ax.zaxis._axinfo['tick']['inward_factor'] = 0
        ax.zaxis._axinfo['tick']['outward_factor'] = 0.6
        ax.zaxis._axinfo['tick']['outward_factor'] = 0.6
    
        ax.set_xlabel(r'$\rm z ~ [kpc]$')
        ax.set_ylabel(r'$\rm R.A.~ [kpc]$')
        ax.set_zlabel(r'$\rm Dec.~ [kpc]$')
    
        x_label = r'$\rm z ~ [kpc]$'
        y_label = r'$\rm R.A.~ [kpc]$'
        z_label = r'$\rm Dec.~ [kpc]$'
    
        z_label = 'Dec. [kpc]'
    #         z_label = 'abcde'
    
        ax.set_zlabel(z_label, rotation=0, fontsize=14, labelpad=40)
    #         ax.set_xlabel(x_label, rotation=0, fontsize=14, labelpad=40)
    #         ax.set_ylabel(y_label, rotation=0, fontsize=14, labelpad=40)

    #         ax.xaxis.set_label_coords(10.0, -200.02)
    
        ax.xaxis.labelpad=18.0
        ax.yaxis.labelpad=1.0
        ax.zaxis.labelpad=28.0

        tight_layout()

        savefig('{0}/NGC3633-RX_J1121.2+0326_3Dmodel_plot4.pdf'.format(directory),format='pdf',bbox_inches='tight')



    # plot the model figure cylinder thing as a rotating movie
    if plot_3Dmodel_movie:
        steps = 100
    
        for i in arange(steps):
            i +=1
            
            color_blue = '#436bad'      # french blue
            color_red = '#ec2d01'     # tomato red

            color_purple = '#7570b3'
            color_purple2 = '#984ea3'
            color_purple3 = '#7570b3'
            color_purple4 = '#810f7c'

            color_green = '#1b9e77'
            color_orange = '#d95f02'
            color_pink = '#e7298a'
            color_lime = '#66a61e'
            color_yellow = '#e6ab02'
            color_brown = '#a6761d'
            color_coal = '#666666'
        
        
            # first get the data
            z = NFW_model['z_sightline']
            y = NFW_model['y_sightline']
            x = NFW_model['x_sightline']
        
            len_step = len(z)/steps
    
            # now the cylinder
            intersect = NFW_model['intersect']
            R = NFW_model['R']
            p0 = NFW_model['p0']
            p1 = NFW_model['p1']

            fig = plt.figure(figsize=(7.0,8.7))
            ax = fig.add_subplot(1,1,1,projection='3d')
        
    
            print 'x: ',x
            print 'y: ',y
            print 'z: ',z
            print
        
            plotExtent = 600
        
            # plot the sightline
#             ax.plot(x, y, z, color='black', alpha = 0.6, lw=2)
            ax.plot(x[:i*len_step], y[:i*len_step], z[:i*len_step], color='black', alpha = 0.7, lw=4)

            # intersect point
            print 'intersect: ',intersect[0],intersect[1],intersect[2]

        #     ax.plot([0,v[0]], [0,v[1]], [0,v[2]], color='green',lw=plotExtent/100)
        #     ax.plot([0,v_90[0]], [0,v_90[1]], [0,v_90[2]], color='purple',lw=plotExtent/100)
        #     ax.plot([intersect[0],v_90[0]], [intersect[1],v_90[1]], [intersect[2],v_90[2]], color='purple',lw=plotExtent/100)


            # put a star on the intersect
        #         planePoint_end = [-1.18639357e-01, 6.80095455e+02, -2.46470324e+02]
        #         planePoint_end2 = [1.18630006e-01, -6.79357841e+02, 2.48330210e+02]
        #         ax.plot([planePoint_end[0]],[planePoint_end[1]],[planePoint_end[2]],color='red',marker='*',lw=0)
        #         ax.plot([planePoint_end2[0]],[planePoint_end2[1]],[planePoint_end2[2]],color='green',marker='*',lw=0)

            
            ax.plot([intersect[0]],[intersect[1]],[intersect[2]], color=color_yellow, marker='*', lw=0)
        
        ##########################################################################################
        ##########################################################################################
            # plot the cylinder

            tube,bottom,top = plot_cylinder(p0,p1,R)

            X, Y, Z = tube
            X2, Y2, Z2 = bottom
            X3, Y3, Z3 = top

            alphaTube = 0.15
            alphaBottom = 0.8
            alphaTop = 0.9   

#             colorTube = 'blue'
#             colorBottom = 'red'
#             colorTop = 'blue'
    
            colorTube = color_purple3
            
            colorBottom = color_red
            colorTop = color_blue


        #     ax=plt.subplot(111, projection='3d')
            ax.plot_surface(X, Y, Z, color=colorTube, alpha = alphaTube)
            ax.plot_surface(X2, Y2, Z2, color=colorBottom, alpha = alphaBottom)
            ax.plot_surface(X3, Y3, Z3, color=colorTop, alpha = alphaTop)
    
            ax.set_xlim(-plotExtent, plotExtent)
            ax.set_ylim(-plotExtent, plotExtent)
            ax.set_zlim(-plotExtent, plotExtent)

            # reverse the RA axis so negative is on the right
        #     ax = plt.gca()
            ax.invert_xaxis()

            # rotate the plot
        #         ax.view_init(elev=10., azim=5)
        #         ax.view_init(elev=15., azim=20)
#             ax.view_init(elev=5., azim=5)

            
            ax.view_init(elev=10., azim=3.6*i)


            yticks((-600, -400, -200, 0, 200, 400, 600), (600, 400, 200, 0, -200, -400, -600))
    #         xticks((-600, -400, -200, 0, 200, 400, 600), (-600, -400, -200, 0, 200, 400, 600))
            xticks((-600, -300, 0, 300, 600), (-600, -300, 0, 300, 600))

    
            [t.set_va('center') for t in ax.get_yticklabels()]
            [t.set_ha('center') for t in ax.get_yticklabels()]
    
            [t.set_va('bottom') for t in ax.get_xticklabels()]
            [t.set_ha('left') for t in ax.get_xticklabels()]
    
            [t.set_va('center') for t in ax.get_zticklabels()]
            [t.set_ha('right') for t in ax.get_zticklabels()]
    
            ax.grid(True)
            ax.xaxis.pane.set_edgecolor('black')
            ax.yaxis.pane.set_edgecolor('black')
            ax.xaxis.pane.fill = False
            ax.yaxis.pane.fill = False
            ax.zaxis.pane.fill = False

            ax.xaxis._axinfo['tick']['inward_factor'] = 0
            ax.xaxis._axinfo['tick']['outward_factor'] = 0.6
            ax.yaxis._axinfo['tick']['inward_factor'] = 0
            ax.yaxis._axinfo['tick']['outward_factor'] = 0.6
            ax.zaxis._axinfo['tick']['inward_factor'] = 0
            ax.zaxis._axinfo['tick']['outward_factor'] = 0.6
            ax.zaxis._axinfo['tick']['outward_factor'] = 0.6
    
            ax.set_xlabel(r'$\rm z ~ [kpc]$')
            ax.set_ylabel(r'$\rm R.A.~ [kpc]$')
            ax.set_zlabel(r'$\rm Dec.~ [kpc]$')
    
            x_label = r'$\rm z ~ [kpc]$'
            y_label = r'$\rm R.A.~ [kpc]$'
            z_label = r'$\rm Dec.~ [kpc]$'
    
            z_label = 'Dec. [kpc]'
        #         z_label = 'abcde'
    
            ax.set_zlabel(z_label, rotation=0, fontsize=14, labelpad=40)
        #         ax.set_xlabel(x_label, rotation=0, fontsize=14, labelpad=40)
        #         ax.set_ylabel(y_label, rotation=0, fontsize=14, labelpad=40)

        #         ax.xaxis.set_label_coords(10.0, -200.02)
    
            ax.xaxis.labelpad=18.0
            ax.yaxis.labelpad=1.0
            ax.zaxis.labelpad=28.0

            tight_layout()

            savefig('{0}/{1}.jpg'.format(movie_directory, i), dpi=200, format='jpg', bbox_inches='tight')



    
if __name__ == '__main__':
    main()
    