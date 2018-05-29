#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: merge_correlatedLists.py, v 1.0 05/29/2018

Combine the measured correlatedList with the new one 
(i.e., put measurements in the new list)


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

    measuredList_filename = '/Users/frenchd/Research/inclination/git_inclination/targets/correlatedTargetList_12_17_17_finished.csv'
    newList_filename = '/Users/frenchd/Research/inclination/git_inclination/targets/correlatedTargetList_5_29_18_data.csv'
    out_filename = '/Users/frenchd/Research/inclination/git_inclination/targets/correlatedTargetList_5_29_18_measurements.csv'
    
    fieldnames = ('target', 'identified', 'v_limits', 'e_v_center', 'Lya_v', 'Lya_W',\
    'e_Lya_W', 'Na', 'e_Na', 'b', 'e_b', 'W', 'e_W', 'comment', 'z_target', 'RAdeg_target',\
    'DEdeg_target','type', 'glon', 'glat', 'z', 'program PI', 'program #', 'obsid',\
    'SN(N=1238,C=1548,O=1031)', 'calver', 'nightonly', 'IGMPARS', 'FPN')
    
    measuredList_file = open(measuredList_filename,'rU')
    measured_reader = csv.DictReader(measuredList_file)
    
    newList_file = open(newList_filename,'rU')
    new_reader = csv.DictReader(newList_file)
    
    writerOutFile = open(out_filename,'wt')
    writer = csv.DictWriter(writerOutFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)


    measurement_dict = {}
    for l in measured_reader:
        target = l['target']
        identified = l['identified']
        v_limits = l['v_limits']
        e_v_center = l['e_v_center']
        Lya_v = l['Lya_v']
        Lya_W = l['Lya_W']
        e_Lya_W = l['e_Lya_W']
        Na = l['Na']
        e_Na = l['e_Na']
        b = l['b']
        e_b = l['e_b']
        W = l['W']
        e_W = l['e_W']
        comment = l['comment']
        z_target = l['z_target']
        
        if measurement_dict.has_key(target):
            measurements = measurement_dict[target]
            measurements.append([identified, v_limits, e_v_center, Lya_v, Lya_W, e_Lya_W,
                                Na, e_Na, b, e_b, W, e_W, comment, z_target])

        else:
            measurement_dict[target] = [[identified,\
            v_limits,\
            e_v_center,\
            Lya_v,\
            Lya_W,\
            e_Lya_W,\
            Na,\
            e_Na,\
            b,\
            e_b,\
            W,\
            e_W,\
            comment,\
            z_target]]
            

    print 'measurement_dict: ', measurement_dict
    print
    
    for i in new_reader:
        target = i['target']
        RAdeg_target = i['RAdeg_target']
        DEdeg_target = i['DEdeg_target']
        type = i['type']
        glon = i['glon']
        glat = i['glat']
        z = i['z']
        programPI = i['program PI']
        programNum = i['program #']
        obsid = i['obsid']
        SN = i['SN(N=1238,C=1548,O=1031)']
        calver = i['calver']
        nightonly = i['nightonly']
        IGMPARS = i['IGMPARS']
        FPN= i['FPN']
        

        if measurement_dict.has_key(target):
            measurements = measurement_dict[target]
            for m in measurements:
#                 print 'm: ',m
                identified = m[0]
                v_limits = m[1]
                e_v_center = m[2]
                Lya_v = m[3]
                Lya_W = m[4]
                e_Lya_W = m[5]
                Na = m[6]
                e_Na = m[7]
                b = m[8]
                e_b  = m[9]
                W = m[10]
                e_W = m[11]
                comment = m[12]
                z_target = m[13]
                
                outputList = [target,\
                identified,\
                v_limits,\
                e_v_center,\
                Lya_v,\
                Lya_W,\
                e_Lya_W,\
                Na,\
                e_Na,\
                b,\
                e_b ,\
                W,\
                e_W,\
                comment,\
                z_target,\
                RAdeg_target,\
                DEdeg_target,\
                type,\
                glon,\
                glat,\
                z,\
                programPI,\
                programNum,\
                obsid,\
                SN,\
                calver,\
                nightonly,
                IGMPARS,
                FPN]
                
                lineRow = dict((f,o) for f,o in zip(fieldnames,outputList))
                writer.writerow(lineRow)
                               
        else:
            identified = i['identified']
            v_limits = i['v_limits']
            e_v_center = i['e_v_center']
            Lya_v = i['Lya_v']
            Lya_W = i['Lya_W']
            e_Lya_W = i['e_Lya_W']
            Na = i['Na']
            e_Na = i['e_Na']
            b = i['b']
            e_b  = i['e_b']
            W = i['W']
            e_W = i['e_W']
            comment = i['comment']
            z_target = i['z_target']
            
            
        outputList = [target,\
        identified,\
        v_limits,\
        e_v_center,\
        Lya_v,\
        Lya_W,\
        e_Lya_W,\
        Na,\
        e_Na,\
        b,\
        e_b ,\
        W,\
        e_W,\
        comment,\
        z_target,\
        RAdeg_target,\
        DEdeg_target,\
        type,\
        glon,\
        glat,\
        z,\
        programPI,\
        programNum,\
        obsid,\
        SN,\
        calver,\
        nightonly,
        IGMPARS,
        FPN]
        
        lineRow = dict((f,o) for f,o in zip(fieldnames,outputList))
        writer.writerow(lineRow)
    
    newList_file.close()
    measuredList_file.close()
    writerOutFile.close()

if __name__=="__main__":
    main()
