#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: convert_correlatedList.py, v 1.0 01/22/2018

Convert Bart's correlation table (correlatedList_12_17_17.rtf) to a csv format with
only the AGN info.



Comes from: makeAGNcsv4.py, v 4.0 03/17/2014

A modification from makeHalosCSV.py to instead create a dictionary csv from Bart's AGN
targetlist. From makeHalosCSV.py version 4:

    Make Halos.txt into a dictionary key csv file (and convert coordinates to degrees)

    vs2: Make Halos2.txt into a dictionary key csv file (and convert coordinates to degrees)

    vs3: Make halos3.txt into a dictionary key csv file (and convert coordinates to degrees)

    vs4: Make NewFilamentAGN.txt into a dictionary key csv file (etc)
    
    
This program successfully made: agn_targetlistCSV.csv on 9/16/2013

This program successfully made: newTARGETLIST2.csv on 02/27/14
This program successfully made: newTARGETLIST2_alt.csv on 02/28/14
This program successfully made: finalTARGETLIST.csv on 03/17/14

"""

from __future__ import division
import optparse
import sys
import os
import tempfile
import csv
import string
from math import *
from utilities import *
import re

    
################################################################



def main():
    # This function reformats Bart's file

    hubbleC = 71
    

#     AGNFile = open('/Users/frenchd/Research/inclination/git_inclination/targets/correlatedList_12_17_17 copy.rtf','rU')
#     outFilename = '/Users/frenchd/Research/inclination/git_inclination/targets/correlatedTargetList_12_17_17.csv'

    AGNFile = open('/Users/frenchd/Research/inclination/git_inclination/targets/correlatedList_5_29_18_data_edit.txt','rU')
    outFilename = '/Users/frenchd/Research/inclination/git_inclination/targets/correlatedTargetList_5_29_18_data.csv'

    AGNFilelines = AGNFile.readlines()
    
    fieldnames = ('target', 'identified', 'v_limits', 'e_v_center', 'Lya_v', 'Lya_W',\
    'e_Lya_W', 'Na', 'e_Na', 'b', 'e_b', 'W', 'e_W', 'comment', 'z_target', 'RAdeg_target',\
    'DEdeg_target','type', 'glon', 'glat', 'z', 'program PI', 'program #', 'obsid',\
    'SN(N=1238,C=1548,O=1031)', 'calver', 'nightonly', 'IGMPARS', 'FPN')
    
    writerOutFile = open(outFilename,'wt')
    writer = csv.DictWriter(writerOutFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)

    start  = False
    separator = '-------------'
    d = {}
    lastThing = False
    print 'starting...'
    counter = 0
    for row in AGNFilelines:
        # test for start of the good stuff
        
#         print 'row: ',row
        if not start:
            if bfind(row,separator):
                start = True
                print 'Start!'
                
        if start:
            if bfind(row,separator):
                # this is the start of a new target
                lastThing = separator
                print 'found first separator'
                
            else:
                if lastThing == separator:
                    row = row.strip()
                    
                    # replace all white space with a single ','
#                     row = re.sub("\s+", ",", row.strip())

                    # replace multiple whitespace with a a single one
                    row = re.sub("\s+", " ", row.strip())

#                     # split the row
#                     sr = row.split(',')
#                     print 'splitrow: ',sr

                    # split the row
                    sr = row.split(' ')
                    print 'splitrow: ',sr
        
                    target = sr[0]
                    
                    if target != "Classification" and len(sr) >10:
                        type = sr[1]
                        glon = sr[2]
                        glat = sr[3]
                        z = sr[4]
                        PI = sr[5].replace('PI:','')
                        progNum = str(sr[6])
                        obsid = sr[7]
                        SN = sr[8]
                        calver = sr[9]
                        nightonly = sr[10].replace('\n','')
                        IGMPARS = sr[11].replace('\n','')
                        FPN = sr[12].replace('\n','')
                        
                        
                        print 'progNum: ',progNum
        
                
                        objectInfoList = [target,\
                        '',\
                        '',\
                        '',\
                        '',\
                        '',\
                        '',\
                        '',\
                        '',\
                        '',\
                        '',\
                        '',\
                        '',\
                        '',\
                        '',\
                        '',\
                        '',\
                        type,
                        glon,\
                        glat,\
                        z,\
                        PI,\
                        progNum,\
                        obsid,\
                        SN,\
                        calver,\
                        nightonly,
                        IGMPARS,
                        FPN]
                    
                        print objectInfoList
        
                        row = dict((f,o) for f,o in zip(fieldnames,objectInfoList))
                        writer.writerow(row)
                        
                        counter +=1
                        lastThing = row
                    
    print 'Count: ',counter
    
    AGNFile.close()
    writerOutFile.close()

if __name__=="__main__":
    main()
