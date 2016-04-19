#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: make_SALT_blocks.py, v 1.0 04/17/2016

Auto-fill SALT observing blocks

"""

import sys
import csv
import string
from math import *
import numpy
from utilities import *
# from astropy.io.votable import parse,tree
import xml.etree.ElementTree as ET
import shutil


##########################################################################################

def main():

    template = '/Users/David/Research_Documents/inclination/git_inclination/SALT_proposal/SALT_targets/target_blocks/target_template.xml'
    datafile = '/Users/David/Research_Documents/inclination/git_inclination/SALT_proposal/SALT_targets_all_phase2.csv'
    
    data = open(datafile,'rU')
    reader = csv.DictReader(data)
    
    limit = 3
    count = 0
    for r in reader:
        count +=1
        if limit >=count:
            galaxyName = r['galaxyName']
            type = r['target type']
            ra = r['right ascension']
            dec = r['declination']
            equinox = r['equinox']
            bandpass = r['bandpass']
            mag = r['minimum magnitude']
            ranking = r['ranking']
            obsT = r['observing time']
            lunar = r['maximum lunar phase']
            expT = r['exposure time']
            sn = r['s/n']
            morph = r['morphology']
            redshift = r['galaxyRedshift']
            slit = r['slit width']
            grating = r['grating']
            camera = r['camera station']
            angle = r['grating angle']
            halpha = r['h-alpha']
            impact = r['impactParameter (kpc)']
            vcorr = r['vcorrGalaxy (km/s)']
            major = r['majorAxis (kpc)']
            diam = r['angDiameter (arcsec)']
            like = r['likelihood']
            lya  = r['possible_lya']
            inc = r['inclination (deg)']
            pa = r['PA']
            az = r['azimuth (deg)']
        
            # new template file name
            newTemp = '/Users/David/Research_Documents/inclination/git_inclination/SALT_proposal/SALT_targets/target_blocks/{0}_template.xml'.format(galaxyName)
            shutil.copy(template,newTemp)
        
            # now make a new block file
        
            tree = ET.parse(newTemp)
            root = tree.getroot()
        
            root[0].text = galaxyName 
            # root[1] = blockcode
            # root[2] = priority
            root[3].text = ranking
            # root[4] = transparency
            root[5][1].text = lunar
            
            tree.write(newTemp)
        
        else:
            break
        
#     tree.close()
    
    
if __name__ == '__main__':
    main()