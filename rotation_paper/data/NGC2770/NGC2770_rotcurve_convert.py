#!/Users/frenchd/anaconda2/bin/python

'''
By David M. French (frenchd@astro.wisc.edu)

Test

'''

import sys
# import os
import csv
import numpy as np
import getpass
from utilities import *
# from pylab import *

import math
import matplotlib.pyplot as plt

from scipy import stats


def cot(x):
    # cotangent
    return 1/math.tan(x)
    
    
def csc(x):
    # cosecant
    return 1/math.sin(x)
    

def inclination_error(v,dv,i,di):
    # calculates the quadrature error in the final velocity value
    # w = observed wavelength of Halpha line center
    # dw = error in wavelength of line center
    # v = systemic velocity
    # dv = error in systemic velocity

    i = i*math.pi/180.
    di = di*math.pi/180
    
    # wavelength term
    incTerm = (v * cot(i) * csc(i) * di)**2

    # v_sys term
    vsysTerm = (csc(i) * dv)**2

    totalError = math.sqrt(vsysTerm + incTerm)

    return totalError


def main():
    # filename = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/NGC2770/NGC2770_rotation_curve.csv'
    # outfile = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/NGC2770/NGC2770_rotation_curve3.csv'

    filename = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/NGC2770/NGC2770_rotation_curve_redo.csv'
    # outfile = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/NGC2770/NGC2770_rotation_curve_redo2.csv'
    outfile = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/NGC2770/NGC2770_rotation_curve_redo4.csv'

    inFile = open(filename,'rU')
    reader = csv.DictReader(inFile)

    fieldnames = ('R (kpc)','vel', 'err')

    outFile = open(outfile,'wt')
    writer = csv.DictWriter(outFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)

    dist = 29.2
    dist_err = 4.18
    vsys = 1948.
    vsys_err = 4
    inclination = 80.
    inc_err = 2.

    # list of axis ratios b/a
    ratios = [0.28, 0.26, 0.26, 0.29, 0.3, 0.2, 0.23, 0.2, 0.23]
    incs = []
    for r in ratios:
        q0 = 0.2
        inc = calculateFancyInclination(1,r,q0)
        incs.append(inc)
    
    print 'mean inclination: ',np.mean(incs)
    print 'median inclination: ',np.median(incs)
    print 'std inclination: ',np.std(incs)
    print 'sem inclination: ',stats.sem(incs)
    print
    
    
    
#     inclination = np.mean(incs)
#     inc_err = stats.sem(incs)
    inc_err = round(np.std(incs),0)
    print 'new inc_err: ',inc_err

    for i in reader:
    #     r = i['R (arcsec)']
        r = i['R (arcmin)']
        v = i['vel']
    
        v2 = float(v) - vsys
    
        v3 = v2/np.sin(inclination*np.pi/180.)
    
    #     r2 = float(r)
        r2 = float(r) * 60.

        r2_lin, r2_lin = calculateLinearDiameters(r2, r2, dist)
        
#         err = inclination_error(vsys,vsys_err, inclination, inc_err)
        err = inclination_error(v2, vsys_err, inclination, inc_err)
        print 'err_form: ', err
        
#         err = (v2/np.sin(inclination*np.pi/180.) - v2/np.sin((inclination-inc_err) * np.pi/180.))
#         print 'err: ',err

                
        # write info to file
        
        objectInfoList = [round(r2_lin,4),v3, round(err,2)]
        row = dict((f,o) for f,o in zip(fieldnames,objectInfoList))
        writer.writerow(row)
    
    print 'Done.'
    
if __name__ == '__main__':
    main()
    
    