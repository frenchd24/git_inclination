#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: createTargetTable_tex.py, v 1.0 07/20/2016

Create table 1 for the pilot paper -> list of targets with ra and dec, z, program ID, 
grating, obs ID, obs date, texp and s/n.

"""

# from __future__ import division
import optparse
import sys
import os
# import tempfile
import csv
import string
from math import *
import numpy
import getpass
import correlateSingle7 as correlateSingle
from utilities import *

    
################################################################

def main(opts):
    # This function reformats Bart's file
    
    
    # newest target file
    filename = '/usr/users/frenchd/correlation/TARGETLIST_6_23_16_reduced.csv'
    
    
    results = open(resultsFilename,'rU')
    reader = csv.DictReader(results)

    fieldnames = ('galaxyName',\
    'AGNname',\
    'degreesJ2000RA_DecGalaxy',\
    'degreesJ2000RA_DecAGN',\
    'impactParameter (kpc)',\
    'vcorrGalaxy (km/s)',\
    'radialVelocity (km/s)',\
    'distGalaxy (Mpc)',\
    'linDiameters (kpc)',\
    'inclination (deg)',\
    'positionAngle (deg)',\
    'azimuth (deg)',\
    'RC3flag',\
    'RC3type',\
    'RC3inc (deg)',\
    'RC3pa (deg)',\
    'morphology',\
    'AGNz',\
    'galaxyRedshift',\
    'groups_dist_std (Mpc)')
        
    
        
    # sort by velocity
    result.sort(reverse=True)
    
    for object in result:
        size, rest = object
        AGNz = rest['AGNredshift']
        AGNra_dec = rest['degreesJ2000RA_DecAGN']
        
        nameList.append(rest['galaxyName'])
        ra_decList.append(rest['degreesJ2000RA_DecGalaxy'])
        separationList.append(rest['impactParameter (kpc)'])
        vcorrList.append(rest['vcorrGalaxy (km/s)'])
        v_helList.append(rest['radialVelocity (km/s)'])
        diameterList.append(rest['linDiameters (kpc)'])
        morphList.append(rest['morphology'])
        
        if opts.verbose:
            for f,i in zip(fieldnames,rest):
                print '{0}: {1}'.format(f,i)
            print
            print
    

    if isNumber(AGNz):
        AGNvel = float(AGNz)*3.0e5
    else:
        AGNvel = 'x'
    
    print '####################################################################'
    print 'Summary: '
    print
    print 'AGN: ',opts.AGNname,', AGN position: ',AGNra_dec,', AGN vel: ',AGNvel
    print

    summaryList = []
    summaryList.append(('galaxyname','ra & dec','separation','v_hel','diameter (kpc)','morphology'))
    for object in result:
        size,rest = object
        summaryList.append((rest['galaxyName'],\
        rest['degreesJ2000RA_DecGalaxy'],\
        rest['impactParameter (kpc)'],\
        rest['radialVelocity (km/s)'],\
        rest['linDiameters (kpc)'],\
        rest['morphology']))

    padding = 4
    widths =[\
    max(len(str(d)) for d in nameList) + padding,\
    max(len(str(d)) for d in ra_decList) + padding,\
    max(len(str(d)) for d in separationList) + padding,\
    max(len(str(d)) for d in vcorrList) + padding,\
    max(len(str(d)) for d in diameterList) + padding,\
    max(len(str(d)) for d in morphList) + padding]

    for row in summaryList:
        print "".join(str(i).ljust(width) for i,width in zip(row,widths))
    
    print
    print 'Finished. {0} total correlated galaxies found.'.format(len(result))
    
else:
    print 'Could not find the AGN, {0}, in the search table.'.format(opts.AGNname)
        


if __name__=="__main__":
    main()
