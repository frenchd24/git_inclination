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


filename = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/NGC3067/NGC3067_rotation_curve.csv'
outfile = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/NGC3067/NGC3067_rotation_curve2.csv'

inFile = open(filename,'rU')
reader = csv.DictReader(inFile)

fieldnames = ('R (kpc)','vel')

outFile = open(outfile,'wt')
writer = csv.DictWriter(outFile, fieldnames=fieldnames)
headers = dict((n,n) for n in fieldnames)
writer.writerow(headers)

dist = 20.4

for i in reader:
    r = i['R (arcsec)']
#     r = i['R (arcmin)']
    v = i['vel']
    
#     r2 = float(r) * 60.
    r2 = float(r)

    r2_lin, r2_lin = calculateLinearDiameters(r2, r2, dist)
        
                
    # write info to file
    objectInfoList = [round(r2_lin,4),v]
    row = dict((f,o) for f,o in zip(fieldnames,objectInfoList))
    writer.writerow(row)
    
print 'Done.'
    