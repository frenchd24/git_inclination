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
#     filename = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/NGC3067/NGC3067_rotation_curve.csv'
#     outfile = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/NGC3067/NGC3067_rotation_curve2.csv'

    filename = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/NGC3067/NGC3067_rubin_rotation_curve_redo2.csv'
    outfile = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/NGC3067/NGC3067_rubin_rotation_curve_redo4.csv'

    inFile = open(filename,'rU')
    reader = csv.DictReader(inFile)

    fieldnames = ('R (kpc)','vel', 'err')

    outFile = open(outfile,'wt')
    writer = csv.DictWriter(outFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)

    dist = 20.4
    dist_err = 2.68
    vsys = 1465.
    vsys_err = 5.
    inclination = 71.
    inc_err = 2.
    
    vsys_err +=2

    # list of axis ratios b/a
    ratios = [0.35, 0.35, 0.36, 0.36, 0.38, 0.38, 0.38, 0.5, 0.57]
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


    for i in reader:
        r = float(i['R (kpc)'])
        v = float(i['vel'])
        
        err = inclination_error(v, vsys_err, inclination, inc_err)
        
        # write info to file
        objectInfoList = [round(r, 4), v, round(err,2)]
        row = dict((f,o) for f,o in zip(fieldnames,objectInfoList))
        writer.writerow(row)
    
    
    print 'Done.'
    
if __name__ == '__main__':
    main()