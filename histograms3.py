#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  histograms3.py, v 3.0 05/06/2015

Plot some stuff for the 100largest initial results

Make plots for AAS winter 2014 poster. Uses LG_correlation_combined2.csv file

Updated for the pilot paper (05/06/15)

'''

import sys
import os
import csv

from pylab import *
# import atpy
from math import *
from utilities import *
import getpass

# from astropy.io.votable import parse,tree

# from vo.table import parse
# import vo.tree


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# ## for Palatino and other serif fonts use:
# #rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)



    
def parseGalaxyNames(nameList):
    # format galaxy names for the url
    
    newNameList = []
    for name in nameList:
        nname = name.strip()
        nname = urllib.quote_plus(nname)
        nname = nname.replace('\n','')
        newNameList.append(nname)
    return newNameList
    

def createCSVTable(outFile,fieldnames):
    # creates and returns a DictReader object populated with header names
    
    writer = csv.DictWriter(outFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)
    return writer


def decidePreferredName(listOfAlternatives,oldName):
    # decides which name to use as the first name, returns that name
    
    print 
    bestName = None
    found = False
    
    # priority:
    first = ['NGC ','IC ','UGC ']
    second = ['MRK ','MCG ','ISO ','SBS ']
    third = ['CGCG ','IRAS ','RXJ ','FGC ','KUG ','PGC ','SDSS ','VCC ']
    fourth = ['2MASS ','2DF ','6DF ','HIPASS ','2MASX ']
    
    for f in first:
        if oldName.find(f) != -1:
            bestName = oldName
            found = True
            break
    
    if not found:
        for name in listOfAlternatives:
            for f in first:
                if name.find(f) != -1 and name.find('[') == -1:
                    bestName = name
                    found = True
                    break
                    
            if not found:
                for s in second:
                    if name.find(s) != -1 and name.find('[') == -1:
                        bestName = name
                        found = True            
                        break
                        
            if not found:
                for t in third:
                    if name.find(t) != -1 and name.find('[') == -1:
                        bestName = name
                        found = True
                        break
    
            if not found:
                for q in fourth:
                    if name.find(q) != -1 and name.find('[') == -1:
                        bestName = name
                        found = True
                        break
            if not found:
                bestName = oldName
    
    found2 = False
    if not found:
        for name in listOfAlternatives:
            if str(name) == oldName:
                bestName = name
                found2 = True
                found = True
                break
    
    if not found2:
        bestName = listOfAlternatives[0]
        
    print 'bestName: ',bestName
#   if bestName != oldName:
#       listOfAlternatives.remove(bestName)
#       listOfAlternatives.append(oldName)
    return bestName


    
def returnLinDiameters(major,minor,distance):
    # input major and minor in arcsec, distance in Mpc
    # outputs major and minor in kpc
    newMajor = math.tan(math.radians(float(major)))*(distance*1000)
    newMinor = math.tan(math.radians(float(minor)))*(distance*1000)
    return (newMajor,newMinor)
    
    

def returnAngDiameters(major,minor,distance):
    # input distances in mpc, major and minor is in kpc
    # outputs angular diameters in arcsec
    newMajor = math.atan((float(major)/1000)/float(distance))*(1/3600)
    newMinor = math.atan((float(minor)/1000)/float(distance))*(1/3600)
    return (newMajor,newMinor)
    
    

###########################################################################

    
def main():
    # assuming 'theFile' contains one name per line, read the file
    
    
    if getpass.getuser() == 'David':
        galaxyFilename = '/Users/David/Research_Documents/NewGalaxyTable3.csv'
        filename = '/Users/David/Research_Documents/inclination/LG_correlation_combined3.csv'
        saveDirectory = '/Users/David/Research_Documents/inclination/pilot_paper/figures/'

    elif getpass.getuser() == 'frenchd':
        galaxyFilename = '/usr/users/frenchd/gt/NewGalaxyTable3.csv'
        filename = '/usr/users/frenchd/inclination/LG_correlation_combined2.csv'
        saveDirectory = '/usr/users/inclination/pilot_paper/figures'

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    
    theFile = open(filename,'rU')
    galaxyFile = open(galaxyFilename,'rU')
    
    reader = csv.DictReader(theFile)
    galaxyReader = csv.DictReader(galaxyFile)
    
    
    
    # save each plot?
    save = False
    
    # create and write output table header and fieldnames
#     fieldnames = ('preferredName','oldName','J2000RA_Dec','alternativeNames')
#     writerOutFile = open(outname,'wt')
#     writer = createCSVTable(writerOutFile,fieldnames)
    
    
    # overall structure: fullDict is a dictionary with all the lines and their data in it
    # separated into 'associated' and 'ambiguous' as the two keys. Associated contains
    # all the lists of data for lines associated with a galaxy. Ambiguous contains all
    # the lists of data for lines not unambiguously associated (could be many galaxies
    # or none)
    
    fullDict = {}
    
    # all the lists to be used for associated lines
    lyaVList = []
    lyaWList = []
    lyaErrorList = []
    naList = []
    bList = []
    impactList = []
    azList = []
    newAzList = []
    incList = []
    fancyIncList = []
    cosIncList = []
    fancyCosIncList = []
    paList = []
    vcorrList = []
    majList = []
    difList = []
    combList = []
    envList = []
    morphList = []
    galaxyNameList = []
    
    # all the lists to be used for ambiguous lines
    lyaVListAmb = []
    lyaWListAmb = []
    lyaErrorListAmb = []
    naListAmb = []
    bListAmb = []
    impactListAmb = []
    azListAmb = []
    newAzListAmb = []
    incListAmb = []
    fancyIncListAmb = []
    cosIncListAmb = []
    fancyCosIncListAmb = []
    paListAmb = []
    vcorrListAmb = []
    majListAmb = []
    difListAmb = []
    combListAmb = []
    envListAmb = []
    morphListAmb = []
    galaxyNameListAmb = []
    
    
    total = 0
    totalNo = 0
    totalYes = 0
    totalIsolated = 0
    totalGroup = 0
    
       
#     associated = array()
#     ambiguous = array()
    
    
    for line in reader:
        #grab all the values
        print 'line: ',line
        AGNra_dec = eval(line['degreesJ2000RA_DecAGN'])
        galaxyRA_Dec = eval(line['degreesJ2000RA_DecGalaxy'])
        lyaV = line['Lya_v']
        lyaW = line['Lya_W']
        env = line['environment']
        galaxyName = line['galaxyName']
        impact = line['impactParameter (kpc)']
        galaxyDist = line['distGalaxy (Mpc)']
        pa = line['positionAngle (deg)']
        morph = line['morphology']
        vcorr = line['vcorrGalaxy (km/s)']
        maj = line['majorAxis (kpc)']
        min = line['minorAxis (kpc)']
        inc = line['inclination (deg)']
        az = line['corrected_az (deg)']
        b = line['b'].partition('pm')[0]
        include = line['include']


        # do some counting
        if include == 'yes' or include == 'no':
            total +=1
            if isNumber(env):
                if int(env) == 0:
                    totalIsolated +=1
                if int(env) >=2:
                    totalGroup +=1
            
        if include == 'yes':
            totalYes +=1
        if include == 'no':
            totalNo +=1
        
        # for all lines
        if isNumber(lyaV) and isNumber(b):
        
            # split up column density estimates and errors
            na = str(line['Na']).split()[0]
            na = float(eval(na.replace('E','e')))
            
            # split up equivalent width measurements and errors
            i = lyaW.find('pm')
            lyaW2 = float(lyaW[:i])
            lyaErr = float(str(lyaW.split('pm')[1]))
            
            b = float(b)
            
            if isNumber(vcorr):
                vcorr = float(vcorr)
                dif = float(lyaV) - vcorr
            else:
                vcorr = 'x'
                dif = 'x'
        
            # azimuth
            if isNumber(az):
                az = abs(float(az))
            else:
                az = -99
            
            # now compute a new one with the working azimuth code
            if isNumber(pa) and isNumber(galaxyDist):
                newAz = calculateAzimuth(galaxyRA_Dec[0], galaxyRA_Dec[1], AGNra_dec[0], AGNra_dec[1], galaxyDist, pa)
            else:
                newAz = -99
            
            
            # position angle
            if isNumber(pa):
                pa = float(pa)
            else:
                pa = -99
        
            # inclination
            if isNumber(inc):
                inc = float(inc)
                cosInc = cos(pi/180 *inc)
            else:
                inc = -99
                cosInc = -99
            
            # major axis diameter and fancy inclination calculation
            if isNumber(maj):
                maj = float(maj)
                if isNumber(min):
                    # minimum axis ratio: q0 = 0.13
                    q0 = 0.13
                    fancyInc = calculateFancyInclination(maj,min,q0)
                    
                    fancyCosInc = cos(pi/180 * float(fancyInc))
                else:
                    fancyInc = -99
                    fancyCosInc = -99
                    
            else:
                maj = -99
                fancyInc = -99
                fancyCosInc = -99
                
            
                
            # try to determine a morphology
            localType = 'x'
            morphology = morph.lower()
            m = morphology[:3]
            if morphology != 'x':
                if bfind(m,'s'):
                    if not bfind(m,'s0'):
                        # straight spiral type
                        localType = 's'
                    else:
                        # lenticular or S0 type
                        localType = 'e'
                elif bfind(m,'e') or bfind(m,'dwarf') or bfind(m,'pec'):
                    localType = 'e'
                elif bfind(m,'len'):
                    localType = 'e'
        
                elif bfind(m,'Ir') or bfind(m,'Im') or bfind(m,'I ') or bfind(m,'IA'):
                    localType = 's'
        
                else:
                    localType = 'x'
            
            # for associated lines
#             if galaxyName !='x' and isNumber(env) and isNumber(lyaV) and isNumber(na) and inc !=-99:
            if include == 'yes':
                lyaV = float(lyaV)
                lyaVList.append(lyaV)
                lyaWList.append(lyaW2)
                lyaErrorList.append(lyaErr)
                naList.append(na)
                bList.append(b)
                impactList.append(impact)
                azList.append(az)
                newAzList.append(newAz)
                incList.append(inc)
                fancyIncList.append(fancyInc)
                cosIncList.append(cosInc)
                fancyCosIncList.append(fancyCosInc)
                paList.append(pa)
                vcorrList.append(vcorr)
                majList.append(maj)
                difList.append(dif)
                envList.append(float(env))
                morphList.append(localType)
                galaxyNameList.append(galaxyName)

            # for ambiguous lines that have a measured galaxy velocity
            if include == 'no' and isNumber(dif):
                lyaVAmb = float(lyaV)
                lyaVListAmb.append(lyaV)
                lyaWListAmb.append(lyaW2)
                lyaErrorListAmb.append(lyaErr)
                naListAmb.append(na)
                bListAmb.append(b)
                impactListAmb.append(impact)
                azListAmb.append(az)
                newAzListAmb.append(newAz)
                incListAmb.append(inc)
                fancyIncListAmb.append(fancyInc)
                cosIncListAmb.append(cosInc)
                fancyCosIncListAmb.append(fancyCosInc)
                paListAmb.append(pa)
                vcorrListAmb.append(vcorr)
                majListAmb.append(maj)
                difListAmb.append(dif)
                envListAmb.append(float(env))
                morphListAmb.append(localType)
                galaxyNameListAmb.append(galaxyName)



    # do some stats:
    incList2 = []
    azList2 = []
    for i,a in zip(incList,azList):
        if i>=50:
            incList2.append(i)
        if a>50:
            azList2.append(a)
            
    print '{0} of {1} are inclined >=50%'.format(len(incList2),len(incList))
    print '{0} of {1} have azimuth >=50%'.format(len(azList2),len(azList))
    print 'a: ',a
    print 'len: ',len(incList)
    print
    print 'total yes: ',totalYes
    print 'total no: ',totalNo
    print 'total: ',total
    print
    print 'totalIsolated: ',totalIsolated
    print 'totalGroup: ',totalGroup


########################################################################################
########################################################################################
########################################################################################
########################################################################################


    # grab all inclinations and position angles in the galaxy dataset
    allInclinations = []
    allCosInclinations = []
    allFancyInclinations = []
    allPA = []
    for line in galaxyReader:
        diameters = eval(line['linDiameters (kpc)'])
        inc = line['inclination (deg)']
        pa = line['positionAngle (deg)']
        
        if isNumber(diameters[0]) and isNumber(inc):
            allCosInclinations.append(cos(pi/180 * round(float(inc),1)))
            allInclinations.append(round(float(inc),1))
            
            # computer fancy inclination with minimum: q0 = 0.13
            q0 = 0.13
            if isNumber(diameters[1]):
                fancyInc = calculateFancyInclination(diameters[0], diameters[1],q0)
            else:
                fancyInc = -99
                
            allFancyInclinations.append(fancyInc)
            
        if isNumber(pa):
            allPA.append(float(pa))
            
    galaxyFile.close()

########################################################################################

    # make a histogram of all the azimuth angles (old, hand measured ones)
    plotAzHist = False
    
    if plotAzHist:
        fig = figure()
        ax = fig.add_subplot(111)
        bins = [0,10,20,30,40,50,60,70,80,90]

        plot1 = hist(azList,bins=bins,histtype='bar')
        title('Distribution of old azimuths')
        xlabel('Azimuth (deg)')
        ylabel('Number')
        xlim(0,90)
        
        if save:
            savefig('{0}/hist(azimuth).pdf'.format(saveDirectory),format='pdf')
        else:
            show()
      
########################################################################################
 
    # make a histogram of all the azimuth angles (new, automatic ones)
    plotNewAzHist = True
    
    if plotNewAzHist:
        fig = figure()
        ax = fig.add_subplot(111)
        bins = [0,10,20,30,40,50,60,70,80,90]

        plot1 = hist(newAzList,bins=bins,histtype='bar')
        title("Distribution of new azimuths")
        xlabel('Azimuth (deg)')
        ylabel('Number')
        xlim(0,90)
        
        if save:
            savefig('{0}/hist(new_azimuth).pdf'.format(saveDirectory),format='pdf')
        else:
            show()
      
########################################################################################

    # make a histogram of the position angle distribution for both the associated galaxy
    # sample and the full galaxy table
    plotPAHist = False
    
    if plotPAHist:
        fig = figure()
        
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(211)
        bins = [0,10,20,30,40,50,60,70,80,90]
        
#         fig.subplots_adjust(left=0.01, bottom=0.01, right=0.01, top=0.01, wspace=0.01, hspace=0.01)
        fig.subplots_adjust(hspace=0.4)

        plot1 = hist(paList,bins=bins,histtype='bar')
        title('Absorber-associated galaxies')
        xlabel('Position Angle (deg)')
        ylabel('Number')
        
        fig.add_subplot(212)
        bins = [0,10,20,30,40,50,60,70,80,90]

        plot1 = hist(allPA,bins=bins,histtype='bar')
        title('Full Galaxy Sample')
        xlabel('Position Angle (deg)')
        ylabel('Number')
#         ylim(0,8000)
#         tight_layout()

        if save:
            savefig('{0}/hist(PA).pdf'.format(saveDirectory),format='pdf')
        else:
            show()

#########################################################################################
    
    # make a histogram of the distribution of impact parameters for associated galaxies
    plotImpactHist = False
    
    if plotImpactHist:
        fig = figure(figsize=(10,2))
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(111)
        
#         bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        plot1 = hist(array(impactList)/array(majList),bins=10,histtype='bar')
        
        title('Distribution of impact parameters')
        xlabel('Impact Parameter (kpc)')
        ylabel('Number')
        ax.tick_params(axis='x', labelsize=0)
        ax.tick_params(axis='y',labelsize=8)
#         ylim(0,5)

        if save:
            savefig('{0}/hist(Impact).pdf'.format(saveDirectory),format='pdf')
        else:
            show()

########################################################################################

    # make a histogram of the distribution of Lyalpha equivalent widths for both the 
    # associated and ambiguous samples
    plotLyaWHist_both = False
    
    if plotLyaWHist_both:
        fig = figure(figsize=(2,8))
        ax = fig.add_subplot(211)
        bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
    
        plot1 = hist(lyaWList,bins=10,histtype='bar',orientation = 'horizontal')
        
        ax = fig.add_subplot(212)
        plot1 = hist(lyaWListAmb,bins=10,histtype='bar',orientation = 'horizontal')
        
        title('Distribution of Lya W')
        xlabel(r'Equivalent Width ($\rm m\AA$)')
        xlabel('Number')
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y',labelsize=0)
#         xlim(0,11)
#         tight_layout()
        
        if save:
            savefig('{0}/hist(lyaW).pdf'.format(saveDirectory),format='pdf')
        else:
            show()

########################################################################################

    # I have no idea what this is or was supposed to do.
    plotInc_vs_dif = False
    
    if plotInc_vs_dif:
        fig = figure()
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(211)
#         bins = [0,10,20,30,40,50,60,70,80,90]
#         bins = [5,15,25,35,45,55,65,75,85,95]

    #     bins = [5,15,25,35,45,55,65,75,85]
#         bins = [0,15,30,45,60,75,90]
#         bins = [0,30,60,90]
        bins = [0,30,60,90]
        incBlue = []
        incRed = []
        for i,d,l in zip(incList,difList,lyaWList):
            if d>0:
                incBlue.append(i)
            else:
                incRed.append(i)
                
        n, bins, patches = hist(incBlue, bins)
        setp(patches, 'facecolor', 'blue', 'alpha', 0.5)               
#         plot1 = hist(incBlue,bins=bins,histtype='bar',c='blue',alpha=0.5)
        title('Blue Shifted Absorbers')
        xlabel('Inclination (deg) / W (AA)')
        ylabel('Number')
#         xlim(0,90)
        
        ax = fig.add_subplot(212)
#         bins = [5,15,25,35,45,55,65,75,85,95]
    #     bins = [5,15,25,35,45,55,65,75,85]
#         bins = [0,15,30,45,60,75,90]
    
        n, bins, patches = hist(incRed, bins)
        setp(patches, 'facecolor', 'red', 'alpha', 0.5)
#         plot2 = hist(incRed,bins=bins,histtype='bar',c='red',alpha=0.5)
        title('Red Shifted Absorbers')
        xlabel('Inclination (deg) / W (AA)')
        ylabel('Number')
#         ylim(0,1)
#         xlim(0,90)
        tight_layout()


        # give me the stats:
        incList2 = []
        azList2 = []
        for i,a in zip(incList,azList):
            if i>=50:
                incList2.append(i)
            if a>50:
                azList2.append(a)
        print '{0} of {1} are inclined >=50%'.format(len(incList2),len(incList))
        print '{0} of {1} have azimuth >=50%'.format(len(azList2),len(azList))
        print 'a: ',a
        print 'len: ',len(incList)
        
        if save:
            savefig('{0}/inc_vs_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()

########################################################################################

    # plot histograms of the associated galaxies' inclinations along with that of the full
    # galaxy set
    plotIncHist_full = False
    
    if plotIncHist_full:
        fig = figure()
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(211)
        bins = [0,10,20,30,40,50,60,70,80,90]
        subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)

    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        plot1 = hist(incList,bins=bins,histtype='bar')
        title('Absorber-associated galaxies inclinations')
        xlabel('Inclination (deg)')
        ylabel('Number')
#         ylim(0,5)
        
        ax = fig.add_subplot(212)
        bins = [0,10,20,30,40,50,60,70,80,90]
    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        plot1 = hist(allInclinations,bins=bins,histtype='bar')
        title('Full Galaxy Sample inclinations')
        xlabel('Inclination (deg)')
        ylabel('Number')
#         ylim(0,8000)
#         tight_layout()

        if save:
            savefig('{0}/hist(inclination).pdf'.format(saveDirectory),format='pdf')
        else:
            show()

########################################################################################

    # plot histograms of the cos(inclinations) for both associated galaxies and the 
    # full galaxy data set
    plotCosIncHist_full = False
    
    if plotCosIncHist_full:
        fig = figure()
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(211)
        bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)

    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        plot1 = hist(cosIncList,bins=bins,histtype='bar')
        title('Absorber-associated galaxies cos(inclination)')
        xlabel('Inclination (deg)')
        ylabel('Number')
#         ylim(0,5)

        ax = fig.add_subplot(212)
        bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        plot1 = hist(allCosInclinations,bins=bins,histtype='bar')
        title('Full galaxy sample cos(inclination)')
        xlabel('Inclination (deg)')
        ylabel('Number')
#         ylim(0,8000)
#         tight_layout()

        if save:
            savefig('{0}/hist(cos(inclination)).pdf'.format(saveDirectory),format='pdf')
        else:
            show()

########################################################################################

    # plot histograms of the fancy inclination for both associated galaxies and the 
    # full galaxy data set
    plotFancyIncHist_full = False
    
    if plotFancyIncHist_full:
        fig = figure()
#         subplots_adjust(hspace=0.200)
        ax = fig.add_subplot(211)
        bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)

    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        plot1 = hist(fancyIncList,bins=bins,histtype='bar')
        title('Absorber-associated galaxies fancy inclination')
        xlabel('Inclination (deg)')
        ylabel('Number')
#         ylim(0,5)

        ax = fig.add_subplot(212)
        bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
    #     bins = [5,15,25,35,45,55,65,75,85]
    #     bins = [0,15,30,45,60,75,90]
        plot1 = hist(allFancyInclinations,bins=bins,histtype='bar')
        title('Full galaxy sample fancy inclination')
        xlabel('Inclination (deg)')
        ylabel('Number')
#         ylim(0,8000)
#         tight_layout()

        if save:
            savefig('{0}/hist(fancy_inclination).pdf'.format(saveDirectory),format='pdf')
        else:
            show()

########################################################################################

    # cos(inclination) histograms for redshifted vs blueshifted distributions of absorbers
    plotCosIncDifHist_full = False
    
    if plotCosIncDifHist_full:
    
#         fig = figure(figsize=(12,10))
#         fig = figure()
#         subplots_adjust(hspace=0.200)

#         subplots_adjust(hspace=0.7)
#         ax = fig.add_subplot(311)
#         subplot(311)
        bins = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
        blue = []
        red = []
        
        blueLya = []
        blueLyaErr = []

        redLya = []
        redLyaErr = []

        
        for d,i,l,e in zip(difList,cosIncList,lyaWList,lyaErrorList):
            if d >0:
                # blue shifted absorber, but galaxy is REDSHIFTED
                print 'd: ',d
                blue.append(i)
                blueLya.append(l)
                blueLyaErr.append(e)
            else:
                # red shifted absorber, but galaxy is BLUESHIFTED
                red.append(i)
                redLya.append(l)
                redLyaErr.append(e)
        
#         plot1 = hist([red,blue],bins=bins,histtype='bar',color=['Red','blue'],alpha=0.5)
 #        title('Red shifted aborption: Galaxies')
#         xlabel('Inclination (deg)')
#         ylim(0,6)
#         ylabel('Number')
        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        subplots_adjust(top = 0.85,bottom = 0.25)
        
        hist(red,bins=bins,histtype='bar',color='red',alpha = 0.9)
        print 'average red: ',average(redLya)
        print 'median red: ', median(redLya)
        print 'average(error): ',average(redLyaErr)
        print 'median(error): ',median(redLyaErr)
        print
        title('Red shifted absorption: Galaxies')
        xlabel('Cos(inclination) = b/a')
#         ylim(0,6)
        ylabel('Number')
        
        if save:
            savefig('{0}/hist(cos(inclination))_dif_red.pdf'.format(saveDirectory),format='pdf')
        else:
            show()


#         fig.add_subplot(212)
#         subplot(212)
        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        subplots_adjust(top = 0.85,bottom=0.25)

        hist(blue,bins=bins,histtype='bar',color='Blue',alpha = 0.9)
        print 'average blue: ',average(blueLya)
        print 'median blue: ',median(blueLya)
        print 'avereage(error): ',average(blueLyaErr)
        print 'median(error): ',median(blueLyaErr)
        
        title('Blue shifted absorption: Galaxies')
        xlabel('Cos(inclination) = b/a')
#         ylim(0,5)
        ylabel('Number')
        
        if save:
            savefig('{0}/hist(cos(inclination))_dif_blue.pdf'.format(saveDirectory),format='pdf')
        else:
            show()

#         ax = fig.add_subplot(312)
#         subplot(312)    
        fig = figure(figsize=(10,2))
        ax = fig.add_subplot(111)
        subplots_adjust(top = 0.85,bottom=0.25)

        hist(allCosInclinations,bins=bins,histtype='bar',color = 'green',alpha=0.9)
        title('Full galaxy sample inclinations')
        xlabel('Cos(inclination) = b/a')
        ylabel('Number')
#         tight_layout()

        if save:
            savefig('{0}/hist(cos(inclination))_fulldataset.pdf'.format(saveDirectory),format='pdf')
        else:
            show()

########################################################################################
########################################################################################

    # plot equivalent width as a function of galaxy diameter
    plotW_Diameter = False
    
    if plotW_Diameter:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,i,w in zip(difList,majList,lyaWList):
            count +=1
            if d>0:
                # galaxy is behind absorber, so gas is blue shifted
                color = 'Blue'
                if countb == 0:
                    countb +=1
                    plotb = ax.scatter(i,w,c='Blue',s=50,label= labelb)
            if d<0:
                # gas is red shifted compared to galaxy
                color = 'Red'
                if countr == 0:
                    countr +=1
                    plotr = ax.scatter(i,w,c='Red',s=50,label= labelr)
                    
            plot1 = scatter(i,w,c=color,s = 50)
        
        title('Equivalent width vs. galaxy diameter')
        xlabel('Major Axis (kpc)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1200)
        ax.legend(scatterpoints=1)
        
        if save:
            savefig('{0}/LyaW(diameter).pdf'.format(saveDirectory),format='pdf')
        else:
            show()

        
########################################################################################
    
    # plot equivalent width as a function of impact parameter, splitting up red and 
    # blue shifted absorption
    plotW_b = False
    
    if plotW_b:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,i,w in zip(difList,impactList,lyaWList):
            count +=1
            if d>0:
                # galaxy is behind absorber, so gas is blue shifted
                color = 'Blue'
                if countb == 0:
                    countb +=1
                    plotb = ax.scatter(i,w,c='Blue',s=50,label= labelb)
            if d<0:
                # gas is red shifted compared to galaxy
                color = 'Red'
                if countr == 0:
                    countr +=1
                    plotr = ax.scatter(i,w,c='Red',s=50,label= labelr)
                    
            plot1 = scatter(i,w,c=color,s = 50)
        
        title('W(impact parameter) for red and blue shifted absorption')
        xlabel('Impact Parameter (kpc)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1200)
        ax.legend(scatterpoints=1)
        
        if save:
            savefig('{0}/W(impact)_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
           
########################################################################################

    # plot equivalent width as a function of impact parameter/diameter, split between
    # red and blue shifted absorption
    plotW_b_diam= False
    
    if plotW_b_diam:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,i,w,m in zip(difList,impactList,lyaWList,majList):
            count +=1
            if d>0:
                # galaxy is behind absorber, so gas is blue shifted
                color = 'Blue'
                if countb == 0:
                    countb +=1
                    plotb = ax.scatter(i/m,w,c='Blue',s=50,label= labelb)
            if d<0:
                # gas is red shifted compared to galaxy
                color = 'Red'
                if countr == 0:
                    countr +=1
                    plotr = ax.scatter(i/m,w,c='Red',s=50,label= labelr)
                    
            plot1 = scatter(i/m,w,c=color,s = 50)
            
        # make the legend work properly
#         labelr = 'Red Shifted Absorber'
#         labelb = "Blue Shifted Absorber"
#         plotb = scatter(i[countb]/m[countb],w[countb],c='Blue',s=50,label= labelb)
#         plotr = scatter(i[countr]/m[countr],w[countr],c='Red',s=50,label= labelr)
        
        title('W(impact/diameter) for red and blue shifted absorption')
        xlabel('Impact Parameter / Diameter')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1200)
        xlim(-1,70)

        ax.legend(scatterpoints=1)
        
        if save:
            savefig('{0}/W(impact-diam)_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
        
########################################################################################

    # plot apparent column density as a function of impact parameter, split between red and
    # blue shifted absorption
    plotNaV_b = False
    
    if plotNaV_b:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,i,n in zip(difList,impactList,naList):
            if d>0:
                # galaxy is behind absorber, so gas is blue shifted
                color = 'Blue'
                if countb == 0:
                    countb +=1
                    plotb = ax.scatter(i,n,c='Blue',s=50,label= labelb)
            if d<0:
                # gas is red shifted compared to galaxy
                color = 'Red'
                if countr == 0:
                    countr +=1
                    plotr = ax.scatter(i,n,c='Red',s=50,label= labelr)
                
            plot1 = scatter(i,n,c=color,s=50)
            
        title('Apparent N(HI) vs impact parameter for red vs blue absorption')
        xlabel('Impact Parameter (kpc)')
        ylabel(r'Na(v) (cm$^{\rm -2}$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,5e14)
        xlim(0,500)
        
        if save:
            savefig('{0}/NaV(impact)_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
        
########################################################################################        

    # plot apparent column density as a function of impact parameter/diameter for red
    # and blue shifted absorption
    plotNaV_b_diam = False
    
    if plotNaV_b_diam:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,i,n,m in zip(difList,impactList,naList,majList):
            if d>0:
                # galaxy is behind absorber, so gas is blue shifted
                color = 'Blue'
                if countb == 0:
                    countb +=1
                    plotb = ax.scatter(i/m,n,c='Blue',s=50,label= labelb)
            if d<0:
                # gas is red shifted compared to galaxy
                color = 'Red'
                if countr == 0:
                    countr +=1
                    plotr = ax.scatter(i/m,n,c='Red',s=50,label= labelr)
                
            plot1 = scatter(i/m,n,c=color,s=50)
            
        title('Apparent N(HI) vs impact/diameter for red vs blue absorption')
        xlabel(r'Impact Parameter / Diameter')
        ylabel(r'Na(v) (cm$^{\rm -2}$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,5e14)
        xlim(-1,70)
        
        if save:
            savefig('{0}/NaV(impact/diameter)_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
             
#########################################################################################
    
    # plot equivalent width as a function of azimuth angle (old one) for red vs blue
    # shifted absorption
    plotW_Az = False
    
    if plotW_Az:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,a,w,m in zip(difList,azList,lyaWList,majList):
            if d>0:
                # galaxy is behind absorber, so gas is blue shifted
                color = 'Blue'
                if countb == 0:
                    countb +=1
                    plotb = ax.scatter(a,w,c='Blue',s=50,label= labelb)
            if d<0:
                # gas is red shifted compared to galaxy
                color = 'Red'
                if countr == 0:
                    countr +=1
                    plotr = ax.scatter(a,w,c='Red',s=50,label= labelr)
                
            plot1 = scatter(a,w,c=color,s=50)
            
        title('W(azimuth_old) for red vs blue shifted absorption')
        xlabel(r'Azimuth (deg)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(0,1200)
        xlim(0,90)

        if save:
            savefig('{0}/W(azimuth_old)_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()

########################################################################################

    # plot equivalent width as a function of azimuth normalized by galaxy size, separated
    # into red and blue shifted absorption samples
    plotW_Az_major = False
    
    if plotW_Az_major:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        
        # give some stats:
        lessThan45 = 0
        for a in azList:
            if a <=45:
                lessThan45 +=1
        print '{0}/{1} have az <= 45 degrees'.format(lessThan45,len(azList))
        print 'average, median azimuth: ',average(azList),', ',median(azList)
        
        for d,a,w,m in zip(difList,azList,lyaWList,majList):
            if d>0:
                # galaxy is behind absorber, so gas is blue shifted
                color = 'Blue'
                if countb == 0:
                    countb +=1
                    plotb = ax.scatter(a/m,w,c='Blue',s=50,label= labelb)
            if d<0:
                # gas is red shifted compared to galaxy
                color = 'Red'
                if countr == 0:
                    countr +=1
                    plotr = ax.scatter(a/m,w,c='Red',s=50,label= labelr)
                
            plot1 = scatter(a/m,w,c=color,s=50)
            
        title('W(azimuth_old/diameter) for red vs blue absorption')
        xlabel(r'Azimuth / Major Axis')
        ylabel(r'Equivalent Width ($\rm \AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
#         ylim(0,600)
#         xlim(0,90)

        if save:
            savefig('{0}/W(azimuth_old/diameter)_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()

########################################################################################

    # plot equivalent width as a function of inclination for red and blue shifted
    # absorption
    plotW_Inc = False
    
    if plotW_Inc:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,a,w,m in zip(difList,incList,lyaWList,majList):
            if d>0:
                # galaxy is behind absorber, so gas is blue shifted
                color = 'Blue'
                if countb == 0:
                    countb +=1
                    plotb = ax.scatter(a,w,c='Blue',s=50,label= labelb)
            if d<0:
                # gas is red shifted compared to galaxy
                color = 'Red'
                if countr == 0:
                    countr +=1
                    plotr = ax.scatter(a,w,c='Red',s=50,label= labelr)
                
            plot1 = scatter(a,w,c=color,s=50)
            
        title('W(inclination) for red vs blue shifted absorption')
        xlabel(r'Inclination (deg)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(0,90)
        
        if save:
            savefig('{0}/W(inclination)_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()

########################################################################################

    # plot equivalent width as a function of cos(inclination) for red and blue shifted
    # absorption
    plotW_CosInc = False
    
    if plotW_CosInc:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,a,w,m in zip(difList,cosIncList,lyaWList,majList):
            if d>0:
                # galaxy is behind absorber, so gas is blue shifted
                color = 'Blue'
                if countb == 0:
                    countb +=1
                    plotb = ax.scatter(a,w,c='Blue',s=50,label= labelb)
            if d<0:
                # gas is red shifted compared to galaxy
                color = 'Red'
                if countr == 0:
                    countr +=1
                    plotr = ax.scatter(a,w,c='Red',s=50,label= labelr)
                
            plot1 = scatter(a,w,c=color,s=50)
            
        title('W(cos(inclination)) for red vs blue shifted absorption')
        xlabel(r'Cos(inclination) = b/a')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(0,1)
        
        if save:
            savefig('{0}/W(cos(inclination))_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()

########################################################################################

    # plot equivalent width as a function of cos(inclination) with red and blue shifted
    # absorption represented by a color bar
    plotW_CosInc_colorbar= False
    
    if plotW_CosInc_colorbar:
        # colormap the velocity difference of the absorber
        averaged =[]

        blueMap = cm.Blues
        redMap = cm.Reds
        
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
#         labelr = 'Red Shifted Absorber'
#         labelb = "Blue Shifted Absorber"
        rdif = []
        rcosInc = []
        rlyaW = []
        rMaj = []
        bdif = []
        bcosInc = []
        blyaW = []
        bMaj = []
        for d,i,w,m in zip(difList,cosIncList,lyaWList,majList):
            if d>0:
                # galaxy is behind absorber, so gas is blue shifted
                color = 'Blue'
                if countb == 0:
                    bdif.append(d)
                    bcosInc.append(i)
                    blyaW.append(w)
                    bMaj.append(m)
#                     plotb = ax.scatter(a, w, cmap=blueMap, c=d, s=50, vmin=0, vmax=400)

            if d<0:
                # gas is red shifted compared to galaxy
                color = 'Red'
                if countr == 0:
#                     plotr = ax.scatter(a, w, cmap=redMap, c=d, s=50, vmin=0, vmax=400)
                    rdif.append(d)
                    rcosInc.append(i)
                    rlyaW.append(w)
                    rMaj.append(m)


        print
        print 'average, median redshifts: ',average(rdif),', ',median(rdif)
        print 'average, median blueshifts: ',average(bdif),', ',median(bdif)
        
        plotr = ax.scatter(rcosInc, rlyaW, cmap=redMap, c=rdif, s=50, vmin=0, vmax=400)
        norm = matplotlib.colors.Normalize(vmin = 0, vmax = 400)
        cbarRed = plt.colorbar(plotr,cmap=redMap,orientation='horizontal')
        cbarRed.set_label('vcorr - absorber velocity (km/s)')
        
        plotb = ax.scatter(bcosInc, blyaW, cmap=blueMap, c=bdif, s=50, vmin=0, vmax=400)
        cbarBlue = plt.colorbar(plotb,cmap=blueMap,orientation='vertical')
        cbarBlue.set_label('vcorr - absorber velocity (km/s)')
#       cbar.ax.set_yticklabels(ticks)
                    
        xlabel(r'Cos(inclination) = b/a')
        ylabel(r'Equivalent Width ($\rm \AA$)')
#         legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
#         ylim(0,600)

        if save:
            savefig('{0}/W(cos(inclination))_dif_colorbar.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
        
        
        
#         fig = figure()
#         ax = fig.add_subplot(211)
#         hist(rdif,color='red')
#        
#         ax = fig.add_subplot(212)
#         hist(bdif,color='blue')
#         show()
        

########################################################################################

    # equivalent width as a function of position angle for both red and blue shifted
    # absorption features
    plotPA = False
    
    if plotPA:
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
        labelr = 'Red Shifted Absorber'
        labelb = "Blue Shifted Absorber"
        for d,a,w,m in zip(difList,paList,lyaWList,majList):
            if d>0:
                # galaxy is behind absorber, so gas is blue shifted
                color = 'Blue'
                if countb == 0:
                    countb +=1
                    plotb = ax.scatter(a,w,c='Blue',s=50,label= labelb)
            if d<0:
                # gas is red shifted compared to galaxy
                color = 'Red'
                if countr == 0:
                    countr +=1
                    plotr = ax.scatter(a,w,c='Red',s=50,label= labelr)
                
            plot1 = scatter(a,w,c=color,s=50)
            
        title('W vs position angle for red vs blue shifted absorption')
        xlabel(r'Position Angle (deg)')
        ylabel(r'Equivalent Width ($\rm m\AA$)')
        legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
        ylim(-1,1200)
        xlim(0,180)
        
        if save:
            savefig('{0}/W(positionAngle)_dif.pdf'.format(saveDirectory),format='pdf')
        else:
            show()

########################################################################################

    # plot equivalent width as a function of cos(fancy-inclination) with red and blue shifted
    # absorption represented by a color bar
    plotW_FancyCosInc_colorbar= False
    
    if plotW_FancyCosInc_colorbar:
        # colormap the velocity difference of the absorber
        averaged =[]

        blueMap = cm.Blues
        redMap = cm.Reds
        
        fig = figure()
        ax = fig.add_subplot(111)
        countb = 0
        countr = 0
        count = -1
#         labelr = 'Red Shifted Absorber'
#         labelb = "Blue Shifted Absorber"
        rdif = []
        rcosInc = []
        rlyaW = []
        rMaj = []
        bdif = []
        bcosInc = []
        blyaW = []
        bMaj = []
        for d,i,w,m in zip(difList,fancyCosIncList,lyaWList,majList):
            if d>0:
                # galaxy is behind absorber, so gas is blue shifted
                color = 'Blue'
                if countb == 0:
                    bdif.append(d)
                    bcosInc.append(i)
                    blyaW.append(w)
                    bMaj.append(m)
#                     plotb = ax.scatter(a, w, cmap=blueMap, c=d, s=50, vmin=0, vmax=400)

            if d<0:
                # gas is red shifted compared to galaxy
                color = 'Red'
                if countr == 0:
#                     plotr = ax.scatter(a, w, cmap=redMap, c=d, s=50, vmin=0, vmax=400)
                    rdif.append(d)
                    rcosInc.append(i)
                    rlyaW.append(w)
                    rMaj.append(m)


        print
        print 'average, median redshifts: ',average(rdif),', ',median(rdif)
        print 'average, median blueshifts: ',average(bdif),', ',median(bdif)
        
        plotr = ax.scatter(rcosInc, rlyaW, cmap=redMap, c=rdif, s=50, vmin=0, vmax=400)
        norm = matplotlib.colors.Normalize(vmin = 0, vmax = 400)
        cbarRed = plt.colorbar(plotr,cmap=redMap,orientation='horizontal')
        cbarRed.set_label('vcorr - absorber velocity (km/s)')
        
        plotb = ax.scatter(bcosInc, blyaW, cmap=blueMap, c=bdif, s=50, vmin=0, vmax=400)
        cbarBlue = plt.colorbar(plotb,cmap=blueMap,orientation='vertical')
        cbarBlue.set_label('vcorr - absorber velocity (km/s)')
#       cbar.ax.set_yticklabels(ticks)
                    
        xlabel(r'Cos(inclination) = b/a')
        ylabel(r'Equivalent Width ($\rm \AA$)')
#         legend(scatterpoints=1)
        ax.grid(b=None,which='major',axis='both')
#         ylim(0,600)

        if save:
            savefig('{0}/W(cos(fancy-inclination))_dif_colorbar.pdf'.format(saveDirectory),format='pdf')
        else:
            show()
        
        
        
#         fig = figure()
#         ax = fig.add_subplot(211)
#         hist(rdif,color='red')
# 
#         ax = fig.add_subplot(212)
#         hist(bdif,color='blue')
#         show()
        

########################################################################################

# how to make side histograms:

    # the random data
#     x = array(azList)
#     y = array(lyaWList)
    
#     fig = figure()
#     ax = fig.add_subplot(211)    
#     plot1 = hist(x,bins=25,histtype='bar')
#     show()
    
#     print 'x:,',x
#     print 'y: ',y
#     print 'type: ',type(x)
# 
#     nullfmt   = NullFormatter()         # no labels
# 
#     # definitions for the axes 
#     left, width = 0.1, 0.65
#     bottom, height = 0.1, 0.65
#     bottom_h = left_h = left+width+0.02
# 
#     rect_scatter = [left, bottom, width, height]
#     rect_histx = [left, bottom_h, width, 0.2]
#     rect_histy = [left_h, bottom, 0.2, height]
# 
#     # start with a rectangular Figure
#     plt.figure(1, figsize=(8,8))
# 
#     axScatter = plt.axes(rect_scatter)
#     axHistx = plt.axes(rect_histx)
#     axHisty = plt.axes(rect_histy)
# 
#     # no labels
#     axHistx.xaxis.set_major_formatter(nullfmt)
#     axHisty.yaxis.set_major_formatter(nullfmt)
# 
#     # the scatter plot:
#     axScatter.scatter(x, y)
# 
#     # now determine nice limits by hand:
#     binwidth = 0.25
#     xmax = np.max(np.fabs(x))
#     ymax = np.max(np.fabs(y))
#     print 'ymax: ',ymax
#     print 'xmax: ',xmax
#     xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
# 
# #     print 'xymax',xymax
#     lim = ( int(xymax/binwidth) + 1) * binwidth
# 
#     axScatter.set_xlim( (0, xlim) )
#     axScatter.set_ylim( (0, ylim) )
# 
#     bins = arange(0, lim + binwidth, binwidth)
#     axHistx.hist(x, bins=bins)
#     axHisty.hist(y, bins=bins, orientation='horizontal')
# 
#     axHistx.set_xlim( axScatter.get_xlim() )
#     axHisty.set_ylim( axScatter.get_ylim() )
# 
#     plt.show()
               
               
               
    theFile.close()
    
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    