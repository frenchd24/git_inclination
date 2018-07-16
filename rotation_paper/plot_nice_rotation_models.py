#!/Users/frenchd/anaconda2/bin/python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  plot_nice_rotation_models.py, v1.0 07/12/18

Plot NFW and regular rotation models all in one in a publication quality plot


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
fontScale = 14
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


'''
========================================================
'''

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
    plot_3Dmodel = True


    if plot_velocity:

        # get the data
        cyl_model_intersect_v_list = cyl_model['intersect_v_list']
        cyl_model_v_proj_list = cyl_model['v_proj_list']

        NFW_model_intersect_v_list = NFW_model['intersect_v_list']
        NFW_model_v_proj_list = NFW_model['v_proj_list']
        

        # do the plotting
        fig = plt.figure(figsize=(6.7,7.7))
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
                
                
        ylabel(r'$\rm V_{proj} ~[km~s^{-1}]$')
        xlabel(r'$\rm Intersect ~[km~s^{-1}]$')
            
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
    
    
        leg = ax.legend(scatterpoints=1,prop={'size':12},loc='lower right',fancybox=True)
    #         leg.get_frame().set_alpha(0.5)

        ax.grid(b=None,which='major',axis='both')
        ylim(-160, 0)
        xlim(-30, 30)

        savefig('{0}/NGC3633-RX_J1121.2+0326_model_plot.pdf'.format(directory),format='pdf',bbox_inches='tight')


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

        savefig('{0}/NGC3633-RX_J1121.2+0326_3Dmodel_plot3.pdf'.format(directory),format='pdf',bbox_inches='tight')


    
if __name__ == '__main__':
    main()
    