#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: print_morph_distribution.py, v 1.0 09/23/2016

A program to print out all types of morphologies in the galaxy table


'''

import sys
import os
import csv
# import string
import warnings
import numpy
# import atpy
import getpass
from utilities import *
import math


# from astropy.io.votable import parse,tree

# from vo.table import parse
# import vo.tree

###########################################################################

    
def main():
    # check which computer we're on, and grab the appropriate file
    
    user = getpass.getuser()
    if user == 'David':
        filename = '/Users/David/Research_Documents/gt/NewGalaxyTable5.csv'
    
    elif user == 'frenchd':
        filename = '/usr/users/frenchd/gt/NewGalaxyTable5.csv'
    
    else:
        print 'User not recognised: ',user
        print 'Exiting.'
        sys.exit()
        
    # open the file and read out the table
    theFile = open(filename,'rU')
    tableReader = csv.DictReader(theFile)

    morphDict = {}
    for line in tableReader:
        morphology = line['morphology']
        
        if morphDict.has_key(morphology):
            i = morphDict[morphology]
            i +=1
            morphDict[morphology] = i
        else:
            morphDict[morphology] = 1
    
    
    keys = morphDict.keys()
    values = morphDict.values()
    
    both = zip(values,keys)
    both.sort(reverse=True)
    
#     sortedValues, sortedKeys = both
    
    for i in both:
        print i    
    
    
    theFile.close()
    print "Done."
    print
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
