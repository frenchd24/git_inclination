#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: trim_correlation_list.py, v 1.0 03/04/2016

trim NGT5-TG6_500Correlation_500cutoff in RA and DEC to match telescope limits

NOT FINISHED. DOES NOTHING RIGHT NOW!!!


"""

# from __future__ import division
import sys
import os
# import tempfile
import csv
import string
from math import *
import numpy
import correlateSingle6 as correlateSingle
from utilities import *

    
################################################################

# def isNumber(s):
#     try:
#         float(s)
#         return True
#     except Exception,e:
#         return False
        
def isNumber_andGreater(s,n):
    try:
        s=float(s)
        if s >= n:
            return True
        else:
            return False
    except Exception,e:
        return False



def main():
    # This function reformats Bart's file
    
    # hubble constant used throughout
    hubbleC = 71.0
    
    # open the files
    
    filename = '/usr/users/frenchd/gt/NGT5-TG6_500Correlation_500cutoff.csv'
    file = open(filename,'rU')        
    outFilename = '/usr/users/frenchd/gt/NGT5-TG6_500Correlation_500cutoff_SALT.csv'
    
    # read in the csv files as dictionary csv files
    reader = csv.DictReader(file)
    
    # note: decimalJ2000RA_Dec = old 'J2000RA_Dec', whereas now J2000RA_Dec refers to coordinates in HH:MM:SS and DD:MM:SS format
    fieldnames = ('AGNname',\
    'galaxyName',\
    'degreesJ2000RA_DecAGN',\
    'degreesJ2000RA_DecGalaxy',\
    'impactParameter (kpc)',\
    'vcorrGalaxy (km/s)',\
    'distGalaxy (Mpc)',\
    'majorAxis (kpc)',\
    'minorAxis (kpc)',\
    'Lstar',\
    'inclination (deg)',\
    'positionAngle (deg)',\
    'RC3pa (deg)',\
    'azimuth (deg)',\
    'RC3type',\
    'morphology',\
    'galaxyRedshift',\
    'AGNredshift',\
    'spectrumStatus',\
    'groups_dist_std (Mpc)')
    
    writerOutFile = open(outFilename,'wt')
#     writer = csv.DictWriter(writerOutFile, fieldnames=fieldnames)
#     headers = dict((n,n) for n in fieldnames)
#     writer.writerow(headers)

    endGroup = True
    counter = 0
    groupCounter = 0

    # RA limits
    maxRA = 0
    minRA = 0
    
    # Dec limits
    maxDec = 0
    minDec = 0

    
    for galaxyRow in reader:
        count +=1
        percentComplete = round((float(count)/2484.)*100,3)
        sys.stdout.write('Percent complete: {0}\r'.format(percentComplete))
        
        galaxyName = galaxyRow['galaxyName']
#         AGNname = galaxyRow['AGNname']
        galaxyPosition = eval(str(galaxyRow['degreesJ2000RA_DecGalaxy']))
#         AGNposition = eval(str(galaxyRow['degreesJ2000RA_DecAGN']))
#         separation = galaxyRow['impactParameter (kpc)']
#         galaxyVcorr = galaxyRow['vcorrGalaxy (km/s)']
#         galaxyDist = galaxyRow['distGalaxy (Mpc)']
#         major,minor = eval(str(galaxyRow['linDiameters (kpc)']))
#         lstar = galaxyRow['Lstar']
#         inclination = galaxyRow['inclination (deg)']
#         positionAngle = galaxyRow['positionAngle (deg)']
#         azimuth = galaxyRow['azimuth (deg)']
#         RC3flag = galaxyRow['RC3flag']
#         RC3type = galaxyRow['RC3type']
#         RC3inc = galaxyRow['RC3inc (deg)']
#         RC3pa = galaxyRow['RC3pa (deg)']
#         morphology = galaxyRow['morphology']
#         AGNz = galaxyRow['AGNredshift']
#         galaxyRedshift = galaxyRow['galaxyRedshift']
#         group = eval(str(galaxyRow['groups_dist_std (Mpc)']))

        # decide how good each galaxy is based on info available

        ra_g,dec_g = galaxyPosition
        ra_agn,dec_agn = AGNposition
        if isNumber(ra_g):
            ra_g = str(round(ra_g,5))
            dec_g = str(round(dec_g,5))
    
        if isNumber(ra_agn):
            ra_agn = str(round(float(ra_agn),5))
            dec_agn = str(round(float(dec_agn),5))
    
        objectInfoList = [AGNname,\
        galaxyName,\
        (ra_agn,dec_agn),\
        (ra_g,dec_g),\
        separation,\
        galaxyVcorr,\
        galaxyDist,\
        major,minor,\
        lstar,\
        inclination,\
        positionAngle,\
        RC3pa,\
        azimuth,\
        RC3type,\
        morphology,\
        galaxyRedshift,\
        AGNz,\
        spectrumStatus,\
        group]


        # make list ordered by major axis size
        objectInfoList2 = [major,objectInfoList]
        fullTargetList.append(objectInfoList2)
                
        # sort the list according to galaxy size
        length = len(fullTargetList)
        fullTargetList.sort(reverse=True)
    
    # finally, write all this to file
    for object in fullTargetList:
        size, rest = object
        row = dict((f,o) for f,o in zip(fieldnames,rest))
        writer.writerow(row)
                                        
                                        
    file.close()
    writerOutFile.close()
    
    print
    print 'Complete!'
    print 'Output: {0}'.format(outFilename)
    print


if __name__=="__main__":
    main()
