#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: output_SALT_list.py, v 1.1 03/08/2016

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
    
    filename = '/usr/users/frenchd/gt/NGT5-TG6_500Correlation_500cutoff_sorted.csv'
    file = open(filename,'rU')        
    outFilename = '/usr/users/frenchd/gt/NGT5-TG6_500Correlation_500cutoff_SALT_sorted_20cut2_full.csv'
    outFilename_SALT = '/usr/users/frenchd/gt/NGT5-TG6_500Correlation_500cutoff_SALT_sorted_20cut2.csv'

    
    # read in the csv files as dictionary csv files
    reader = csv.DictReader(file)
    
    writerOutFile = open(outFilename,'wt')
    writerOutFile_SALT = open(outFilename_SALT,'wt')


    # RA region to avoid
    maxRA = 11
    minRA = 10
    
    # Dec limits
    maxDec = 11.
    minDec = -76.
    
    minSize = 20.

    count =-1
    entryDict = {}
    fullTargetList = []
    for galaxyRow in reader:
        count +=1
        percentComplete = round((float(count)/2484.)*100,3)
        sys.stdout.write('Percent complete: {0}\r'.format(percentComplete))
        
        galaxyName = galaxyRow['galaxyName']
        AGNname = galaxyRow['AGNname']
        galaxyPosition = eval(str(galaxyRow['degreesJ2000RA_DecGalaxy']))
        AGNposition = eval(str(galaxyRow['degreesJ2000RA_DecAGN']))
        separation = galaxyRow['impactParameter (kpc)']
        galaxyVcorr = galaxyRow['vcorrGalaxy (km/s)']
        galaxyDist = galaxyRow['distGalaxy (Mpc)']
        major = eval(galaxyRow['majorAxis (kpc)'])
        minor = eval(galaxyRow['minorAxis (kpc)'])
        lstar = galaxyRow['Lstar']
        inclination = galaxyRow['inclination (deg)']
        positionAngle = galaxyRow['positionAngle (deg)']
        azimuth = galaxyRow['azimuth (deg)']
#         RC3flag = galaxyRow['RC3flag']
        RC3type = galaxyRow['RC3type']
#         RC3inc = galaxyRow['RC3inc (deg)']
        RC3pa = galaxyRow['RC3pa (deg)']
        morphology = galaxyRow['morphology']
        AGNz = galaxyRow['AGNredshift']
        galaxyRedshift = galaxyRow['galaxyRedshift']
        group = eval(str(galaxyRow['groups_dist_std (Mpc)']))
        

        # convert ra and dec to sexagesimal
        ra_g,dec_g = galaxyPosition
        ra_g2, dec_g2 = convertRAandDec(ra_g,dec_g,'sexagesimal')

        go = False
        if float(dec_g) >=minDec:
            if float(dec_g) <= maxDec:
                if float(ra_g) >= maxRA or float(ra_g) <= minRA:
                    if float(major) >= minSize:
                        go = True
        
        if go:
            print
            print galaxyName
            print ra_g2, dec_g2
                   
            line = '"{0}" {1} {2} {3} {4} {5} {6} {7} {8}\n'.format(galaxyName,'Galaxy',ra_g2[0],ra_g2[1], ra_g2[2],dec_g2[0],dec_g2[1],dec_g2[2],2000.0)

            if entryDict.has_key(galaxyName):
                # i is [count,line], so add to the count
                i = entryDict[galaxyName]
                count = int(i[0]) +1
                entryDict[galaxyName] = [count,line]
                
                # update the list of AGN near this target
                agns, gList = targetList[0],targetList[1]
                agns.append(AGNname)
                targetList = [agns,gList]

            # otherwise make a new dictionary entry for this galaxy
            else:
                # name, 'n' in dictionary is now associated with list containing 'p' photometry value
                entryDict[galaxyName] = [1,line]
                
                targetList = [[AGNname],\
                [galaxyName,\
                galaxyPosition,\
                separation,\
                galaxyVcorr,\
                galaxyDist,\
                major,\
                minor,\
                lstar,\
                inclination,\
                positionAngle,\
                azimuth,\
                RC3type,\
                RC3pa,\
                morphology,\
                galaxyRedshift,\
                group]]
            
            fullTargetList.append(targetList)

    names = entryDict.keys()
    lines = entryDict.values()
        
    lines.sort(reverse=True)
    
    print '----- stats ------'
    for s in lines:
        print s
        
        count,line = s[0],s[1]
        writerOutFile_SALT.write(line)

    
##########################################################################################
    # write all the stuff to file

    fieldnames = ('AGNname',\
    'galaxyName',\
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
    'groups_dist_std (Mpc)')
    
    writerOutFile = open(outFilename,'wt')
    writer = csv.DictWriter(writerOutFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)


    for object in fullTargetList:
        row = dict((f,o) for f,o in zip(fieldnames,object))
        writer.writerow(row)

                                        
    file.close()
    writerOutFile.close()
    
    print
    print 'Complete!'
    print 'Output: {0}'.format(outFilename)
    print


if __name__=="__main__":
    main()
