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

    filename = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/UGC04238/UGC04238_rotation_curve_redo.csv'
    outfile = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/UGC04238/UGC04238_rotation_curve_redo3.csv'

    inFile = open(filename,'rU')
    reader = csv.DictReader(inFile)

    fieldnames = ('R (kpc)','vel','err')

    outFile = open(outfile,'wt')
    writer = csv.DictWriter(outFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)
    
    dist = 22.4
    dist_err = 4.21
    vsys = 1544.0
    vsys_err = 7.0
#     inclination = 62.0 -> UGC04238_rotation_curve_redo2.csv
    inclination = 75.0 # -> UGC04238_rotation_curve_redo3.csv

    inc_err = 2.0
    
    vsys_err +=2

    # list of axis ratios b/a
    ratios = [0.570, 0.320, 0.320, 0.310, 0.590, 0.570]
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
    inc_err = np.std(incs)

    for i in reader:
    #     r = i['R (arcsec)']
        r = float(i['R (kpc)'])
        v = float(i['vel'])
        e = float(i['err'])
        
        r = float(r)
        v = float(v)
        
        err_i = inclination_error(v, vsys_err, inclination, inc_err)
        print 'v: ',v
        print 'err_i: ',err_i
        
        err = np.sqrt(e**2 + err_i**2)
        print 'final err: ',err
        
        # write info to file
        objectInfoList = [round(r, 4), v, round(err,2)]
        row = dict((f,o) for f,o in zip(fieldnames,objectInfoList))
        writer.writerow(row)
    
    print 'Done.'
    
if __name__ == '__main__':
    main()
    