#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: corrolateAGN15.py, v 15.1 03/04/2016

Corrolate AGN targets with galaxies, producing a sorted correlation list

Updated version of: correlateAGN7_new.py.

This version correlates with the current full targetlist. Precedence is given to 
galaxies with known pa, inclination, and size. Sorting is by size (largest first) 

v8: Include AGN S/N and galaxy photometry info

v10: Sort entire list by galaxy diameter, include AGN redshift, and only include AGN
with nonzero redshift and S/N >10

v10.1: round off long floats

v11.2: add 'spectrumStatus' column to correlation table

v12: use correlateSingle3.py module throughout. No longer run from commandline with optparse.
    this version made LG_correlation.csv
    
v13: correlate with updated table: NewGalaxyTable3.csv. try to make azimuth work again.
    makes: '/usr/users/frenchd/gt/NGT3-TG6_2000Correlation_full_500cutoff2.csv'
    
    do a b<150kpc correlation for Bart for HST proposal - 04/07/15
    
    
v14: try to drastically speed this up -> making NGT3-TG6_2000Correlation_full_500cutoff4.csv,
    NGT3-TG6_2000Correlation_full_500cutoff5.csv

v15: same, this time making NGT3-TG6_2000Correlation_full_500cutoff6.csv
    - (actually said: corrolateAGN14.py, v 14.0 04/07/2015)
    
v15.1 same as above, now making NGT5-TG6_500Correlation_500cutoff.csv (03/04/16)
    - started using correlateSingle6 here

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

# def convertRAandDec(ra,dec):
#     # this function converts ra and dec in degrees to HH:MM:SS and DD:MM:SS
#     
#     raHours,raMinutes = str(((float(ra)*24)/360)).split('.')
#     raMinutes,raSeconds = str(float('0.'+raMinutes) *60).split('.')
#     raSeconds = float('0.'+raSeconds) *60
# 
#     decDegrees,decMinutes = str(dec).split('.')
#     decMinutes,decSeconds = str(float('0.'+decMinutes)*60).split('.')
#     decSeconds = float('0.'+decSeconds)*60
#     
#     return (raHours,raMinutes,round(raSeconds,2)),(decDegrees,decMinutes,round(decSeconds,2))


# def calculatevcorr(ra,dec,velocity):
#     rav = 186.7833
#     decv = 12.9333
#     vcorr = velocity + 300*(sin(radians(dec)) * sin(radians(decv)) + cos(radians(dec))*cos(radians(decv)) * cos(radians(ra-rav)))
#     return vcorr
    
# def returnLinDiameters(major,minor,distance):
#     newMajor = tan(radians(float(major)/60))*distance
#     newMinor = tan(radians(float(minor)/60))*distance
#     return (newMajor,newMinor)
#     
# def calculateInclination(major,minor):
#     print 'minor/major: ',minor,'/',major
#     print 'float(minor)/float(major): ',float(minor),'/',float(major)
#     inclination = acos(float(minor)/float(major))*180/pi
#     return inclination

# def calculateAngularSeparation(ra1,dec1,ra2,dec2,dist):
#     angSep = acos(cos(dec1)*cos(dec2) + sin(dec1)*sin(dec2)*cos(ra1-ra2))
#     distSep = angSep*dist*1000
#     return distSep
    
# def calculateAngularSeparation(ra1,dec1,ra2,dec2,dist):
#     angSep = acos(sin(dec1)*sin(dec2) + cos(dec1)*cos(dec2)*cos(ra1-ra2))
#     distSep = angSep*dist*1000
#     return distSep



def main():
    # This function reformats Bart's file
    
    # hubble constant used throughout
    hubbleC = 71.0
    
    # open the files
    
    AGNfilename = '/usr/users/frenchd/correlation/TARGETLISTupdate_6.csv'
    galaxyFilename = '/usr/users/frenchd/gt/NewGalaxyTable5.csv'
    galaxyFile = open(galaxyFilename,'rU')        
    AGNfile = open(AGNfilename,'rU')
    outFilename = '/usr/users/frenchd/gt/NGT5-TG6_500Correlation_500cutoff.csv'
    
    # read in the csv files as dictionary csv files
    galaxyReader = csv.DictReader(galaxyFile)
    AGNReader = csv.DictReader(AGNfile)
    
    
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
    writer = csv.DictWriter(writerOutFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)

    endGroup = True
    counter = 0
    groupCounter = 0

    # maximum separation between AGN target and galaxy
    maxSep = 500.0
    
    # how far behind the galaxy must the AGN be in km/s?
    agnSeparation = 0
    
    # minimum S/N for AGN target - not implemented for correlateTarget_fast
    minimumSN = 0

    # minimum vcorr for each galaxy
    minimumVcorr = 500.0
    
    # minimum size of galaxy to include. 'half' = top 50% largest correlated galaxies,
    # 'quarter' = top 25% largest, otherwise a number represents minimum radius
    # 
    minSize = 10.0
    
    # True if you only want to include the largest correlated galaxy to each AGN target
    largestOnly = False
    
    # maximum distance for galaxy (Mpc)
    maxDist = False

    count = 0
    AGNlist = []
    for AGNrow in AGNReader:
        # instantiate list to be populated by info for targets within maxSep
        # this will later be used to make minSize cuts to the correlated galaxy list

        # these are the newest field names for the new targetlist
        AGNname = AGNrow['AGNname']
        AGNra_dec = eval(AGNrow['degreesJ2000RA_DecAGN'])
        AGNz = AGNrow['z']
        spectrumStatus = AGNrow['spectrumStatus']
        AGNlist.append([AGNname,AGNra_dec,AGNz,spectrumStatus])

    
    galaxyList = []
    for galaxyRow in galaxyReader:
        galaxyName = galaxyRow['preferredName']
        galaxyPosition = eval(str(galaxyRow['degreesJ2000RA_Dec']))
        galaxyDist = galaxyRow['Best Distance (Mpc)']
        group = eval(str(galaxyRow['Groups_Dist_std (Mpc)']))
        major,minor = eval(str(galaxyRow['linDiameters (kpc)']))
        galaxyVcorr = galaxyRow['vcorr (km/s)']
        morphology = galaxyRow['morphology']
        RC3flag = galaxyRow['RC3flag']
        RC3inc = galaxyRow['RC3inc (deg)']
        RC3type = galaxyRow['RC3type']
        RC3pa = galaxyRow['RC3pa (deg)']
        positionAngle = galaxyRow['positionAngle (deg)']
        inclination = galaxyRow['inclination (deg)']
        galaxyRedshift = galaxyRow['redshift']
        Lstar_bmed = galaxyRow['B_median_Lstar']
        Lstar_sdss = galaxyRow['B_sdss_median_Lstar']
        
        if isNumber(Lstar_bmed):
            lstar = float(Lstar_bmed)
        elif isNumber(Lstar_sdss):
            lstar = float(Lstar_sdss)
        else:
            lstar = 'x'
        
        if isNumber(inclination):
            inclination = round(eval(inclination),1)

        # decide how good each galaxy is based on info available
    
        ra_g,dec_g = galaxyPosition
        if isNumber(ra_g):
            ra_g = str(round(ra_g,5))
            dec_g = str(round(dec_g,5))
        
        if isNumber(galaxyVcorr):
            galaxyVcorr = round(float(galaxyVcorr),2)
        
        if isNumber(major):
            major = round(float(major),2)
    
        if isNumber(minor):
            minor = round(float(minor),2)
        
        if isNumber(inclination):
            inclination = int(float(inclination))
        
        if isNumber(galaxyDist) and isNumber(galaxyVcorr):
        
            goSize = True
            if minSize:
                if isNumber(major):
                    if float(major) >= float(minSize):
                        goSize = True
                    else:
                        goSize = False
            
            goDist = True
            if maxDist:
                if float(galaxyDist) <= float(maxDist):
                    goDist = True
                else:
                    goDist = False
                    
            if float(galaxyVcorr) >= minimumVcorr and goSize and goDist:
                galaxyList.append([galaxyName,\
                galaxyPosition,\
                galaxyDist,\
                group,\
                major,\
                minor,\
                galaxyVcorr,\
                morphology,\
                RC3flag,\
                RC3type,\
                RC3inc,\
                RC3pa,\
                positionAngle,\
                inclination,\
                galaxyRedshift,\
                lstar])
                    
    
    print len(galaxyList)
    print len(AGNlist)
    # now loop through the AGN list
    
    for agn in AGNlist:
        AGNname,AGNra_dec,AGNz,spectrumStatus = agn
        correlation = correlateSingle.correlateTarget_fast(agn, galaxyList, maxSep, agnSeparation, minimumVcorr, minSize, slow=False)
        galaxyInfo = correlation[AGNname]
        
        fullTargetList = []
        count +=1
        percentComplete = round((float(count)/504.)*100,3)
        sys.stdout.write('Percent complete: {0}\r'.format(percentComplete))
        
        if galaxyInfo:
            for row in galaxyInfo:
                vcorr, galaxyRow = row
                galaxyName = galaxyRow['galaxyName']
                AGNname = galaxyRow['AGNname']
                galaxyPosition = eval(str(galaxyRow['degreesJ2000RA_DecGalaxy']))
                AGNposition = eval(str(galaxyRow['degreesJ2000RA_DecAGN']))
                separation = galaxyRow['impactParameter (kpc)']
                galaxyVcorr = galaxyRow['vcorrGalaxy (km/s)']
                galaxyDist = galaxyRow['distGalaxy (Mpc)']
                major,minor = eval(str(galaxyRow['linDiameters (kpc)']))
                lstar = galaxyRow['Lstar']
                inclination = galaxyRow['inclination (deg)']
                positionAngle = galaxyRow['positionAngle (deg)']
                azimuth = galaxyRow['azimuth (deg)']
                RC3flag = galaxyRow['RC3flag']
                RC3type = galaxyRow['RC3type']
                RC3inc = galaxyRow['RC3inc (deg)']
                RC3pa = galaxyRow['RC3pa (deg)']
                morphology = galaxyRow['morphology']
                AGNz = galaxyRow['AGNredshift']
                galaxyRedshift = galaxyRow['galaxyRedshift']
                group = eval(str(galaxyRow['groups_dist_std (Mpc)']))

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
                                        
                                        
    galaxyFile.close()
    AGNfile.close()
    writerOutFile.close()
    
    print
    print 'Complete!'
    print 'Output: {0}'.format(outFilename)
    print


if __name__=="__main__":
    main()
