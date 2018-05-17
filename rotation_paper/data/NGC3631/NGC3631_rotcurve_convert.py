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


# filename = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/NGC2770/NGC2770_rotation_curve.csv'
# outfile = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/NGC2770/NGC2770_rotation_curve3.csv'

filename = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/NGC3631/NGC3631_rotation_curve.csv'
outfile = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/NGC3631/NGC3631_rotation_curve2.csv'

inFile = open(filename,'rU')
reader = csv.DictReader(inFile)

fieldnames = ('R (kpc)','vel')

outFile = open(outfile,'wt')
writer = csv.DictWriter(outFile, fieldnames=fieldnames)
headers = dict((n,n) for n in fieldnames)
writer.writerow(headers)

dist = 8.56
vsys = 1156.
inclination = 17.

for i in reader:
#     r = i['R (arcsec)']
    r = i['R (arcmin)']
    v = i['vel']
    
#     v2 = float(v) - vsys
    v3 = v
    
#     v3 = v2/np.sin(inclination*np.pi/180.)
    
#     r2 = float(r)
    r2 = float(r) * 60.

    r2_lin, r2_lin = calculateLinearDiameters(r2, r2, dist)
        
                
    # write info to file
    objectInfoList = [round(r2_lin,4),v3]
    row = dict((f,o) for f,o in zip(fieldnames,objectInfoList))
    writer.writerow(row)
    
print 'Done.'
    