#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: createTargetTable_tex.py, v 1.0 02/22/2018


Based on: createTargetTable_tex.py, v 1.0 07/20/2016

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

from astropy.io import ascii
from astropy import table
import pickle

################################################################

def main():
    # This function reformats Bart's file
    
    if getpass.getuser() == 'frenchd':
        target_filename = '/Users/frenchd/Research/correlation/TARGETLIST_06_06_18_TOTAL.csv'
#         resultsName = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/salt_galaxy_sightlines_cut.csv'
#         resultsName = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/salt_galaxy_sightlines_cut_plus_ancillary_fits.csv'
#         outName = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/table2_2.txt'
#         outName2 = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/table2_2_alt.txt'

#         results_filename = '/Users/frenchd/Research/inclination/git_inclination/targets/correlatedTargetList_5_29_18_measurements_copy.csv'

        all_filename = '/Users/frenchd/Research/inclination/git_inclination/all8_double.p'

        outName = '/Users/frenchd/Research/inclination/git_inclination/thesis_final_chapter/table1_4.txt'
        
        target_SN_filename = '/Users/frenchd/Research/inclination/git_inclination/thesis_final_chapter/target_SN.csv'

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # Bart's TARGETLIST
    target_file = open(target_filename,'rU')
    target_reader = csv.DictReader(target_file)
    
    # Bart's SN measurements for these targets
    target_SN_file = open(target_SN_filename, 'rU')
    SN_reader = csv.DictReader(target_SN_file)
    
    output = open(outName,'wt')
    
    all_file = open(all_filename,'r')
    all = pickle.load(all_file)
    all_file.close()
    
    # here is the data from the pickle files
    targets = all['targets']
    z_targets = all['z_targets']
    RA_targets = all['RA_targets']
    Dec_targets = all['Dec_targets']
    Lya_vs = all['Lya_vs']
    Lya_Ws = all['Lya_Ws']
    Lya_bs = all['bs']
    
    # put all the measurements in a dictionary
    target_measurement_dict = {}
    for target, lya_v, lya_w, lya_b in zip(targets, Lya_vs, Lya_Ws, Lya_bs):
        if target_measurement_dict.has_key(target):
            others = target_measurement_dict[target]
            others.append({'Lya_v':lya_v, 'Lya_W':lya_w, 'Lya_b':lya_b})
            
        else:
            target_measurement_dict[target] = [{'Lya_v':lya_v,
                                                'Lya_W':lya_w,
                                                'Lya_b':lya_b}]
    
    targets.sort()
        
    summaryList = []
    summaryList.append(('Target',\
    'R.A.',\
    'Dec.',\
    'z',\
    'Program',\
    'T_exp',\
    'SN',\
    'Lya_v',\
    'Lya_W',\
    'Lya_b'))
    
    nameList = ['Target']
    raList = ['R.A.']
    decList = ['Dec.']
    zList = ['z']
    programList = ['Program']
#     gratingList = ['Grating']
#     obsIDList = ['Obs ID']
#     obsDateList = ['Obs Date']
    texpList = ['T_exp [ks]']
#     snList = ['S/N [1238]']

    SNList = ['S/N [1238']
    lya_vList = ['Lya_v']
    lya_WList = ['Lya_W']
    lya_bList = ['Lya_b']
    
    used_dict = {}
    target_dict = {}
    target_SN_dict = {}
    
    # get the SN measurements
    for sn in SN_reader:
        name = sn['targetname']
        sn = sn['SN']
        
        target_SN_dict[name] = sn
    
    
    
    for l in target_reader:
        target = l['targetName']
        
        if not target_dict.has_key(target) and target_SN_dict.has_key(target):
            ra = l['ra']
            dec = l['dec']
            z = l['z']
    #                 sn = l['SNratio']
            programID = l['programID']
    #                 programPI = l['programPI']
    #                 instrument = l['instrument']
    #                 grating = l['grating']
            texp = l['exposureTime']
    #                 sn = l['SNexpected']
    
            SN = target_SN_dict[target]
            
            target_dict[target] = {'ra':ra,
                                    'dec':dec,
                                    'z':z,
                                    'programID':programID,
                                    'texp':texp,
                                    'SN':SN}
                
    used_targets = []
    for target in targets:
        if target not in used_targets:
            measurements = target_measurement_dict[target]
        
            for measurement in measurements:
                ra = target_dict[target]['ra']
                dec = target_dict[target]['dec']
                z = target_dict[target]['z']
                programID = target_dict[target]['programID']
                texp = target_dict[target]['texp']
                SN = target_dict[target]['SN']
    
                Lya_v = measurement['Lya_v']
                Lya_W = measurement['Lya_W']
                Lya_b = measurement['Lya_b']

        #             ra = eval(ra.replace('(','').replace(')','').replace(',',' '))
                ra = eval(ra)
                dec = eval(dec)

                ra1 = int(ra[0])
                ra2 = int(ra[1])
                ra3 = ra[2]
    
                dec1 = int(dec[0])
                dec2 = int(dec[1])
                dec3 = dec[2]
                
                ra1 = str(ra1)
                ra2 = str(ra2)
                ra3 = str(ra3)
                dec1 = str(dec1)
                dec2 = str(dec2)
                dec3 = str(dec3)
        
                if len(ra1) == 1:
                    ra1 = '0'+ra1
            
                if len(ra2) == 1:
                    ra2 = '0'+ra2
            
                if len(ra3) == 3:
                    num1, decimal, num2 = ra3
                    ra3 = '0' + num1 + decimal + num2
        
        
                if bfind(dec1,'-'):
                    if len(dec1) == 2:
                        minus, number = dec1
                        dec1 = '$' + minus + '$' + '0' + number
                    else:
                        minus, number1, number2 = dec1
                        dec1 = '$' + minus + '$' + number1 + number2
                else:
                    if len(dec1) == 1:
                        dec1 = '$+$' + '0' + dec1
                    else:
                        dec1 = '$+$' + dec1
        
        
                if len(dec2) == 1:
                    dec2 = '0' + dec2
        
        
                if len(dec3) == 3:
                    num1, decimal, num2 = dec3
                    dec3 = '0' + num1 + decimal + num2
            
    
                ra_str = ra1 + ' ' + ra2 + ' ' + ra3
                dec_str = dec1 + ' ' + dec2 + ' ' + dec3
        
                target_table = target.replace('_','\_')

                summaryList.append((target_table,\
                ra_str,\
                dec_str,\
                str(z),\
                programID,\
                str(texp),
                str(SN),
                str(Lya_v),
                str(Lya_W),
                str(Lya_b)))
        
                nameList.append(target_table)
                raList.append(ra)
                decList.append(dec)
                zList.append(z)
                programList.append(programID)
                texpList.append(texp)
                lya_vList.append(Lya_v)
                lya_WList.append(Lya_W)
                lya_bList.append(Lya_b)
                SNList.append(SN)
                
                used_targets.append(target)
    

    padding = 4
    widths =[\
    max(len(str(d)) for d in nameList) + padding,\
    max(len(str(d)) for d in raList) + padding,\
    max(len(str(d)) for d in decList) + padding,\
    max(len(str(d)) for d in zList) + padding,\
    max(len(str(d)) for d in programList) + padding,\
    max(len(str(d)) for d in texpList) + padding,\
    max(len(str(d)) for d in SNList) + padding,\
    max(len(str(d)) for d in lya_vList) + padding,\
    max(len(str(d)) for d in lya_WList) + padding,\
    max(len(str(d)) for d in lya_bList) + padding]



    for row in summaryList:
        output.write("".join(str(i + '  &').ljust(width) for i,width in zip(row,widths))+'\\\\\n')


    target_file.close()
    target_SN_file.close()
    output.close()

if __name__=="__main__":
    main()
