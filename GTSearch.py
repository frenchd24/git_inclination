#!/usr/bin/env python

'''
By David M. French (frenchd@astro.wisc.edu)

$Id: GTSearch.py, v 1.1 03/09/2016

A program to grab galaxy data from the master galaxy table (v1.0 09/12/2013)

v1.1: updated for NewGalaxyTable5.csv , and moved to /usr/users/frenchd/gt/


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

def absoluteMag_noExtinc(m,d):
    # m is apparent magnitude, d is distance in Mpc
    M = float(m) - 5*math.log10((float(d)*10**6)/10)
    return M


def absoluteMag(m,d,e):
    # m is apparent magnitude, d is distance in Mpc, e is extinction E(B-V)
    M = float(m) - 5*math.log10((float(d)*10**6)/10) - 3.1*float(e)
    return M


def printOutInfo(line,t):
    # grab all the information in the table for an object and return it all if t == f, 
    # or just the basics if t == b
    
    preferredName = line['preferredName']
    oldName = line['oldName']
    redshift = line['redshift']
    degreesJ2000 = eval(line['degreesJ2000RA_Dec'])
    J2000 = eval(line['J2000RA_Dec'])
    gLongLat = eval(line['galacticLong_Lat'])
    radVel = line['radialVelocity (km/s)']
    vcorr = line['vcorr (km/s)']
    distvcorr = line['distvcorr (Mpc)']
    rid = eval(line['rIndependentDistMean_sd_min_max (Mpc)'])
    bestDist = line['Best Distance (Mpc)']
    angDiameter = line['angDiameters (arcsec)']
    linDiameter = line['linDiameters (kpc)']
    inclination = line['inclination (deg)']
    positionAngle = line['positionAngle (deg)']
    diameterKey = line['diameterKey']
    RC3flag = line['RC3flag']
    RC3type = line['RC3type']
    RC3inc = line['RC3inc (deg)']
    RC3pa = line['RC3pa (deg)']
    Groups = eval(line['Groups_Dist_std (Mpc)'])
    groupsInfo = line['groupsInfo']
    morphology = line['morphology']
    distanceIndicator = line['distanceIndicator']
    luminosityClass = line['luminosityClass']
    EBminusV = line['EBminusV']
    EBminusV_new = line['E(B-V)_new']
    B_median = line['B_median']
    B_sdss_median = line['B_sdss_median']
    B_median_Lstar = line['B_median_Lstar']
    B_max_Lstar = line['B_max_Lstar']
    B_sdss_median_Lstar = line['B_sdss_median_Lstar']
    B_sdss_max_Lstar = line['B_sdss_median_Lstar']
    Lstar = line['Lstar']
    photometry = eval(line['photometry'])
    altNames = eval(line['alternativeNames'])
    
    if isNumber(inclination):
        bestInc = inclination
    else:
        bestInc = RC3inc
    
    if isNumber(positionAngle):
        bestPA = positionAngle
    else:
        bestPA = RC3pa
    
    if isNumber(B_median):
        bestB = B_median
    else:
        bestB = B_sdss_median
    
    if isNumber(B_median_Lstar):
        bestLstar = B_median_Lstar
    else:
        bestLstar = B_sdss_median_Lstar
    
    # if t == 'b':
    print
    print 'Preferred Name: ',preferredName
    print 'Equatorial Coordinates: ',degreesJ2000
    print 'redshift: ', redshift
    print 'vcorr: ',vcorr
    print 'Best Distance: ',bestDist
    print 'Linear diameter: ',linDiameter,' kpc'
    print 'Angular diameter: ',angDiameter, ' (arcsec)'
    print 'Inclination: ',bestInc
    print 'Position angle: ',bestPA
    print 'Morphology: ',morphology
    print
    print 'Best L* estimates: '
    print 'B_median: ',bestB
    print 'Lstar: ',bestLstar
    print
    print 'All estimates: '
    print 'B_median: ',B_median
    print 'B_sdss_median: ',B_sdss_median
    print 'B_median_Lstar: ', B_median_Lstar
    print 'B_max_Lstar: ', B_max_Lstar
    print 'B_sdss_median_Lstar: ', B_sdss_median_Lstar
    print 'B_sdss_median_Lstar: ', B_sdss_max_Lstar
    print
    print 'All: '
    for i in eval(Lstar):
        print 'Lstar: ', i
        
    if t == 'f':
        # print some additional info
        print
        print '------ Additional Info -------'
        print 'Old name: ',oldName
        print 'Equatorial Coordiantes: ',J2000
        print 'Galactic Coordinates: ',gLongLat
        print 'Radial Velocity: ',radVel
        print 'Redshift independent dist (Mean_sd_min_max): ',rid
        print 'RC3flag: ',RC3flag
        print 'RC3type: ',RC3type
        print 'RC3 Inclination: ',RC3inc
        print 'RC3 Position Angle: ',RC3pa
        print 'Groups_distances_std: ',Groups
        print 'Distance Indicator: ',distanceIndicator
        print 'Luminosity Class: ',luminosityClass
        print 'EBminusV_new: ',EBminusV_new
        print 'Photometry measurments: '
        for m in photometry:
            print '\t ',m
        
        print
        print 'Alternative names: '
        for a in altNames:
            print '\t ',a
            
        print
        
    return
    
    
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
    
    # ask for the galaxy name
    quit = False
    while not quit:
        theFile.seek(0)
        tableReader = csv.DictReader(theFile)
        
        response = raw_input("Galaxy name: ")
        response = response.replace(' ','')
        response = response.upper()
    
        # search the table for this name
        for line in tableReader:
            preferredName = line['preferredName'].upper()
            oldName = line['oldName'].upper()
            alternatives = eval(line['alternativeNames'])
            found = False
            if response == preferredName or response == oldName:
#                 typeResponse = 'n'
                typeResponse = 'f'

                while typeResponse != 'b' and typeResponse != 'f':
                    typeResponse = raw_input("Return basic data (b) or full (f)?: ")
                    
                dummy = printOutInfo(line,typeResponse)
                found = True
                break
        
            else:
                for name in alternatives:
                    name = name.upper()
                    if name == response:
                        typeResponse = 'n'
                        while typeResponse != 'b' and typeResponse != 'f':
                            typeResponse = raw_input("Return basic data (b) or full (f)?: ")
                            
                        printOutInfo(line,typeResponse)
                        found = True
                        break

            if found:
                break
                            
        if not found:
            quitAnswer = 'x'
            while quitAnswer != 'n' and quitAnswer != 'y':
                quitAnswer = raw_input('Sorry, this galaxy was not found, search again? (y/n) ')
            if quitAnswer == 'n':
                quit = True
            else:
                quit = False
        
        if found:
            againAnswer = 'x'
            while againAnswer != 'y' and againAnswer != 'n':
                againAnswer = raw_input("Search again? (y/n) ")
            if againAnswer == 'y':
                quit = False
            if againAnswer == 'n':
                quit = True
                break
            
 
    theFile.close()
    print "Done."
    print
    
###############################################################################

if __name__=="__main__":
    # do the work
    main()
