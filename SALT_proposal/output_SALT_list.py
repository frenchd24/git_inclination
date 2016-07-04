#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: output_SALT_list.py, v 1.0 03/08/2016

trim NGT5-TG6_500Correlation_500cutoff in RA and DEC to match telescope limits

Made NGT5-TG6_500Correlation_500cutoff_SALT_sorted_15cut.csv for updated SALT proposal (05/03/16) 

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
import getpass

    
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


def calculateLikelihood(impact, diameter):
    # return likelihood
    
    rVir = calculateVirialRadius(float(diameter))
    l = math.exp(-(float(impact)/rVir)**2)
    
    return l
    

def main():
    # This function reformats Bart's file
    
    # hubble constant used throughout
    hubbleC = 71.0
    
    # open the files
    
    if getpass.getuser() == 'David':
        filename = '/Users/David/Research_Documents/gt/NGT5-TG6_500Correlation_500cutoff_sorted.csv'
        file = open(filename,'rU')
        outFilename = '/Users/David/Research_Documents/inclination/git_inclination/SALT_proposal/NGT5-TG6_500Correlation_500cutoff_SALT_sorted_15cut_full.csv'
        outFilename_SALT = '/Users/David/Research_Documents/inclination/git_inclination/SALT_proposal/NGT5-TG6_500Correlation_500cutoff_SALT_sorted_15cut.csv'
        
    elif getpass.getuser() == 'frenchd':
        filename = '/usr/users/frenchd/inclination/git_inclination/SALT_proposal/NGT5-TG6_500Correlation_500cutoff_sorted.csv'
        file = open(filename,'rU')
        outFilename = '/usr/users/frenchd/inclination/git_inclination/SALT_proposal/NGT5-TG6_500Correlation_500cutoff_SALT_sorted_15cut_full.csv'
        outFilename_SALT = '/usr/users/frenchd/inclination/git_inclination/SALT_proposal/NGT5-TG6_500Correlation_500cutoff_SALT_sorted_15cut.csv'

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # read in the csv files as dictionary csv files
    reader = csv.DictReader(file)
    
    writerOutFile = open(outFilename,'wt')
    writerOutFile_SALT = open(outFilename_SALT,'wt')


    # RA region to avoid
    maxRA = 11.
    minRA = 10.
    
    # Dec limits
    maxDec = 11.
    minDec = -76.
    
    # minimum size for a galaxy to be considered (in kpc)
    minSize = 15.
    
    # currently not used
    l_min = 0.001
    
    # minimum separation between the galaxy and the AGN (in km/s)
    minSep = 4000
    
    # list of AGN targets that aren't any good
#     avoidList = ['UM228',\
#     'SDSSJ143004.07+022213.3',\
#     'SDSSJ130524.30+035731.0',\
#     'SDSSJ124423.37+021540.4',\
#     'SDSSJ123647.80+060048.0',\
#     'SDSSJ123426.80+072411.0',\
#     'SDSSJ122018.43+064119.6',\
#     'SDSSJ121640.60+071224.0',\
#     'SBS0335-052',\
#     'RX_J1232.5+0603',\
#     'RX_J1223.2+0922',\
#     'RX_J1218.8+1015',\
#     'RBS1089',\
#     'QSO0246-3050',\
#     'PHL1444',\
#     'PG1115+080',\
#     'NVSSJ152511-171436',\
#     'NGC7552',\
#     'NGC4696',\
#     'NGC4593',\
#     'NGC3783',\
#     'MS0244.6-3020',\
#     'LBQS1230-0015',\
#     'LBQS1222+0901',\
#     'LBQS0302-0019',\
#     'IRAS11598-0112',\
#     'HE0241-3043',\
#     'ESO031-G08',\
#     '4C10.34',\
#     '2dFGRS_S394Z150']

# included in avoidList because of no data available:
#     'SDSSJ123647.80+060048.0',\
#    'SDSSJ123426.80+072411.0',\
#    'SDSSJ122018.43+064119.6'
#    'RX_J1232.5+0603',\
#    'RX_J1223.2+0922',\
#    'RX_J1218.8+1015',\
#    'RBS1089',\
#    'QSO0246-3050',\
#    'PHL1444',\
#    'PG1115+080',\
#    'NGC4696',\
#    'MS0244.6-3020',\
#    'LBQS1230-0015',\
#    'LBQS1222+0901',\
#    'LBQS0302-0019',\
#    'IRAS11598-0112',\
#    'HE0241-3043',\
#    '4C10.34',\

# don't know why this is included:     'RX_J1429.6+0321',\




    avoidList = ['SDSSJ143004.07+022213.3',\
    'SDSSJ124423.37+021540.4',\
    'SDSSJ123647.80+060048.0',\
    'SDSSJ123426.80+072411.0',\
    'SDSSJ122018.43+064119.6',\
    'SBS0335-052',\
    'RX_J1232.5+0603',\
    'RX_J1223.2+0922',\
    'RX_J1218.8+1015',\
    'RBS1089',\
    'QSO0246-3050',\
    'PHL1444',\
    'PG1115+080',\
    'NVSSJ152511-171436',\
    'NGC7552',\
    'NGC4696',\
    'NGC4593',\
    'NGC3783',\
    'MS0244.6-3020',\
    'LBQS1230-0015',\
    'LBQS1222+0901',\
    'LBQS0302-0019',\
    'IRAS11598-0112',\
    'HE0241-3043',\
    'ESO031-G08',\
    '4C10.34',\
    '2dFGRS_S394Z150',\
    '2dFGRSS394Z150',\
    'RX_J1429.6+0321',\
    'LBQS1206+1052',\
    'QSO1502-4154',\
    'SDSSJ125013.50+073441.5']

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
        
        # calculate the separation between the galaxy and the AGN
        AGNvel = float(AGNz)* 3*10**5
        sep = AGNvel - float(galaxyVcorr)

        go = False
        goAGN = True
        
        # does this AGN suck?
        for a in avoidList:
            if AGNname == a:
                goAGN = False
                break
                
        if goAGN and sep >= minSep:
            if float(dec_g) >=minDec:
                if float(dec_g) <= maxDec:
                    # avoid the 10-11 HR range in RA for SALT
                    if float(ra_g2[0]) >= maxRA or float(ra_g2[0]) < minRA:
                        print '{0} >= {1} or {2} < {3}'.format(float(ra_g2[0]),maxRA,float(ra_g2[0]),minRA)
                        if float(major) >= minSize:
                            go = True
        
        if go:
            print
            print galaxyName
            print ra_g2, dec_g2
                   
            line = '"{0}" {1} {2} {3} {4} {5} {6} {7} {8}\n'.format(galaxyName,'Galaxy',ra_g2[0],ra_g2[1], ra_g2[2],dec_g2[0],dec_g2[1],dec_g2[2],2000.0)

            if isNumber(major):
                likelihood = calculateLikelihood(separation,major)
            else:
                likelihood = -99

            if entryDict.has_key(galaxyName):
                # i is [likelihood,line], so see if this one is any better than the last
                i = entryDict[galaxyName]
                old_likelihood = i[0]
                
                # if the new association has a higher likelihood, use it instead
                if likelihood > old_likelihood:
                    entryDict[galaxyName] = [likelihood,line]
                else:
                    pass

            # otherwise make a new dictionary entry for this galaxy
            else:
                # name, 'n' in dictionary is now associated with list containing 'p' photometry value
                entryDict[galaxyName] = [likelihood,line]
            
            targetList = [galaxyName,\
            AGNname,\
            AGNposition,\
            galaxyPosition,\
            separation,\
            galaxyVcorr,\
            galaxyDist,\
            major,\
            minor,\
            lstar,\
            likelihood,\
            inclination,\
            positionAngle,\
            azimuth,\
            RC3type,\
            RC3pa,\
            morphology,\
            AGNz,\
            galaxyRedshift,\
            group]
            
            fullTargetList.append([likelihood,targetList])

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

    fieldnames = ('galaxyName',\
    'AGNname',\
    'degreesJ2000RA_DecAGN',\
    'degreesJ2000RA_DecGalaxy',\
    'impactParameter (kpc)',\
    'vcorrGalaxy (km/s)',\
    'distGalaxy (Mpc)',\
    'majorAxis (kpc)',\
    'minorAxis (kpc)',\
    'Lstar',\
    'likelihood',\
    'inclination (deg)',\
    'positionAngle (deg)',\
    'azimuth (deg)',\
    'RC3type',\
    'RC3pa (deg)',\
    'morphology',\
    'AGNredshift',\
    'galaxyRedshift',\
    'groups_dist_std (Mpc)')
    
    writerOutFile = open(outFilename,'wt')
    writer = csv.DictWriter(writerOutFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)

    fullTargetList.sort(reverse=True)

    for object in fullTargetList:
        likelihood, rest = object
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
