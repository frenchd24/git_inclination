#!/usr/bin/env python
"""
By David French (frenchd@astro.wisc.edu)

$Id: plotCorrelationMap13.py, v 13 12/06/17

This program takes in a list of AGN targets and generates an environment map (i.e. nearby
galaxies) for each. 


Adapted from:
'
    plotAGNCorrelationSquare9.py, v 9.0 02/12/2014

    Make a square plot showing the positions of correlated galaxies around a target AGN

    Starts with the older 'plotAGNCorrelationSquare6.py' and makes adjustments for the current
    data.

    v8: include inclinations in galaxy annotation tags

    v9: plots all targets in the input file, and does everything better faster
'

v1: Removes group vs field galaxy bits

v2: Draw ellipses to denote inclination and position angle at the same time

v3: Make nicer plots for poster (09/10/14)

v4: change the labeling for individual plots (09/11/14)

v5: adopt to make plots for general correlation tables

v6: include option to make skyplot as well (plot just the ra and dec of all objects)
    This was final AAS poster version, made targetmaps3
    
v7: More general version - 06/19/15

v7.1: Now saves into 'ambiguous' or 'associated' folders depending on the 'include' answer
    - 06/25/15

7.2: Reverse PA direction. Now PA is measured clockwise from N (up) to E (on right). This
    changes angle = float(PA[i]) to angle = float(-PA[i]) in plotting
    
    Also changed the name annotation slightly, to xy=(-60, size[i]/6) from 
    xy=(-65, size[i]/4), and unknown type marker size increased to size[i]*4
    
    E type galaxies now have no edge, spiral have black edges, unknown are diamonds
    
    - 07/07/15
    
    - unfortunately, this version was still named "plotCorrelationMap7_1.py".
    

7.3: Sort the map table by the quantity: abs(dif_velocity)*impact parameter/virial radius
    - low means the galaxy in question close in both velocity and distance. These get
     listed first.

7.4: Same as above, but tweaking the "likelihood" parameter to the following:
    abs(dif_velocity) * (impact parameter - virial radius)

7.5: Now the likelihood parameter is:
      abs(dif_velocity) * (impact parameter - virial radius) / (virial radius)**2
      
      if this value is >2 times the next largest, call it associated and put it in that folder
      
      11/6/15: (approximate date) Now use a gaussian form likelihood parameter, and take
      the velocity difference into account as well (normalized to 200 km/s)
        -looks like: exp((impact/virial)^2) * exp((vel_dif/200)^2) - Targetmaps12
      
      11/11/15: continued with more tweaks to try to get satellite galaxies to be rated
            lower. Targetmaps13 uses diam^1.4 instead of a virial radius
            
            Targetmaps14 uses diam^1.5
            
    11/15/15: Targetmaps15 uses diam^1.5, and multiplies likelihood by 2 if the virial
            radius is larger than the impact parameter
            
            Targetmaps16 requires 3* greater likelihood for associated tag
            
            Targetmaps17 back to standard virial radius, but with the above tweaks for 
            multiplying likelihood *2 and requiring 3* difference
            
    11/17/15: Targetmaps18 back to diam^1.5, and now with 5* difference required
    
    11/19/15: Targetmaps19 uses standard r_vir, but requires a 5* difference
    
    
v8.0: Computes both standard R_vir and diam^1.5 results. Adds the results automatically to
LG_correlation_combined5.csv, 

v9.0: velocity colormap in units of delta-v instead of absolute velocity units

v9.1: 'include' now requires a hard L limit (L>=0.001 etc)

v10: correlate with galaxy v_helio, not vcorr - using correlateSingle7 now
    - (03/24/16 - 4/13/16)
    
v11: moved to .../inclination/git_inclination/ so it can be updated with git, and 
    now including rc tex formating for the plots (05/12/2016)
    
v11.1: updates for the newest round of sightlines in LG_correlation_combined5_10.csv
       - Made LG_correlation_combined5_11.csv and targetmaps34 (07/06/16)
       
v11.2: minor formatting updates. (8/08/16)

v11.3: more minor updates - make the tick labels NOT bold (10/03/16)
    - make LG_correlation_combined5_11_25cut_edit5.csv and targetmaps37/

v11.3: more minor updates - fix the colorbar bolded ticks problem (10/10/16)
    - remake LG_correlation_combined5_11_25cut_edit5.csv and targetmaps37/
    
v11.4: removes the 'include_folder' for saving - saves everything in the same general
    directory, so you don't have to fuss around when doing single correlations
    (1/16/17)

v11.5: include the speed of light for converting from redshift, upgrade to 
    correlateSingle9.py - uses the newest galaxy table (7/21/17)
    
v12: this comes from plotCorrelationMap_single.py -> specifically made to plot single
    galaxy maps instead of from a big list, for the metals project (Taesun?)
    
    Adapting to make this more flexible. Now you can plot single or multiple from a file
    or list. You can also include both galaxies and AGN as both center targets AND/OR
    included in the environment field


v13: version12 incorrectly plots the orientation of galaxies. They need to be reflected
    over the RA axis, and the axis does also. This will match east to the left standard
    for astronomical images. (12/06/17)
    
v13.1: diamonds were plotted backwards (-PA), fixed this. (12/18/17)

"""

import sys
import os
import csv
import string
import math
import ast
from pylab import *
import correlateSingle11 as correlateSingle
from matplotlib.patches import Ellipse
from utilities import *
import getpass

from matplotlib import rc
fontScale = 18
rc('text', usetex=True)
rc('font', size=18,family='serif',weight='medium')
rc('xtick.major',size=8,width=0.6)
rc('xtick.minor',size=5,width=0.6)
rc('ytick.major',size=8,width=0.6)
rc('ytick.minor',size=5,width=0.6)
rc('xtick',labelsize = fontScale)
rc('ytick',labelsize = fontScale)
rc('axes',labelsize = fontScale)
rc('xtick',labelsize = fontScale)
rc('ytick',labelsize = fontScale)
# rc('font', weight = 450)
# rc('axes',labelweight = 'bold')
rc('axes',linewidth = 1,labelweight='normal')
rc('axes',titlesize='small')

    
################################################################


    
def buildTargetList(file,AGNheader,velocityHeader):
    # builds a list consisting of tuples: (AGNname, velocity, include)
    #
    # requires a column called 'galaxyName' in the file. Only includes targets for which
    # 'galaxyName' != 'x'
    #
    # - AGN name is the name of the target
    # - velocity is the velocity of the absorption feature centroid in km/s
    # - include is a boolean indicating whether the results should be saved into the
    #   'ambiguous' folder (False), or the 'associated' folder (True)
    
    f = open(file,'rU')
    reader = csv.DictReader(f)
    targetList = []
    for l in reader:
        AGNname = l[AGNheader]
        include = l['include']
        if include == 'yes':
            include = True
        else:
            include = False
        
        # velocityHeader is the name of the columns containing the center velocity, or 
        # it can also be just a number designating the center velocity
        if isNumber(velocityHeader):
            vel = velocityHeader
        else:
            vel = l[velocityHeader]
            
        galaxyName = l['galaxyName']
        if isNumber(vel) and galaxyName !='x':
            pair = (AGNname,int(vel), include)
            targetList.append(pair)
    
    f.close()
    return targetList


def buildFullTargetList(file,AGNheader,velocityHeader):
    # makes a targetmap for every absorption line, regardless of whether there is an
    # associated galaxy name or not, so long as the line velocity entry is a number
    #
    # builds a list consisting of tuples: (AGNname, velocity, include)
    #
    # - AGN name is the name of the target
    # - velocity is the velocity of the absorption feature centroid in km/s
    # - include is a boolean indicating whether the results should be saved into the
    #   'ambiguous' folder (False), or the 'associated' folder (True)
    
    f = open(file,'rU')
    reader = csv.DictReader(f)
    targetList = []
#     startT = 'QSO1500-4140'
#     startV = '9757'
#     goOn = False
    for l in reader:
        AGNname = l[AGNheader]
        
        
#         include = l['include']
#         if include == 'yes':
#             include = True
#         else:
#             include = False
        
        # velocityHeader is the name of the columns containing the center velocity, or 
        # it can also be just a number designating the center velocity
    
    
    
        # added
#         vel = l[velocityHeader]
#         if AGNname == startT:
#             if vel == startV:
#                 goOn = True
#         
#         if goOn:
        # added
        
        if isNumber(velocityHeader):
            vel = velocityHeader
        else:
            vel = l[velocityHeader]
   
        if isNumber(vel):
            include = True
            pair = (AGNname,int(vel),include)
            targetList.append(pair)
    
    f.close()
    return targetList
        

def main():
    # main function to create targetmaps around selected sightlines
    nullFloat = -99.99
    nullInt = -99
    nullStr = 'x'
        
    counter = 0
    AGNList = []
    rowList = []
    masterVirList = []
    masterCustomList = []
    
    # max impact parameter to use
    maxSep = 500
    
    # +/- galaxy velocity to search within around absorption velocity
    velocityWindow = 400
    
    # minimum galaxy velocity to include (False to ignore). Usually = 500
    minVcorr = 500
    
    # minimum galaxy size to include (False to ignore)
    minSize = False
    
    # minimum separation in km/s between the redshift of the AGN and the galaxy (False to ignore) 
#     agnSeparation = 4000.
    agnSeparation = False

    
    # include name tags on galaxies? They don't scale very well...
    includeNameTags = True
    
    # x and y name tag offset
#     yTagOffset = 6
#     xTagOffset = -60
    yTagOffset = 2
    xTagOffset = -10
    
    # name tag font size
    nameTagFont = 4.5
    
    # include a title on the plots?
    includeTitle = True
    
    # also make a plot of just real positions of galaxies on the sky in RA and Dec coords?
    includeSkyPlot = False
    
    # Save the map plots?
    saveMaps = True
    
    # Save the individual map plot tables?
    saveMapTables = True
    
    # Save the full results with "include" tags? This is the whole big correlation table
    # which looks like LG_correlation_combined5_11_25cut_edit4.csv
    saveResults = False
    
    # 2nd place galaxy likelihood * rigor <= 1st place galaxy for 'include'
    rigor = 5
    
    # hard limit for likelihood
    l_min = 0.01
    
    # bypass l_min for lone galaxies? (i.e. include lone galaxies no matter what likelihood is)
    loner = False
    
    # the speed of light
    c = 2.9979*10**5
    
    # maximum number of galaxies to plot on a single window. Just set it high to ignore
    maxPlotObjects = 5000
    
    # ticks
    xAxisMajorTicks = 200
    xAxisMinorTicks = 100
    yAxisMajorTicks = 200
    yAxisMinorTicks = 100
    velocityStepSize = 100
    
    # sort results into /associated/, ~/ambiguous/, and ~/nonassociated/ folders?
    # if True, these folders must already exist
    # if False, puts all the files into saveDirectory as set below
    sortIntoFolders = False
    
    # include AGN background targets as well?
    includeAGN = True
    
    # include name tags for AGN?
    includeAGNnameTags = True
    
    # size of stars for AGN
    AGNsize = 50

    # color map: 
    colmap = cm.RdBu_r
#     colmap = cm.inferno

    # which way to plot RA axis? RAeastLeft = True puts east to the left, matching
    # the standard. I.e., RA increases toward the left
    RAeastLeft = True

    
    # Location and name of targetfile, and where to save figures and tables. Set this 
    # even if you'll manually enter targets below
    user = getpass.getuser()
    if user == "David":
#         targetFile = '/Users/David/Research_Documents/inclination/git_inclination/LG_correlation_combined5_11_25cut_edit4.csv'
#         saveDirectory = '/Users/David/Research_Documents/iraf/NGC3633/'
#         outputFile = '/Users/David/Research_Documents/iraf/NGC3633_correlation.csv'

#         targetFile = '/Users/David/Research_Documents/metal_absorbers/met.dat'
#         saveDirectory = '/Users/David/Research_Documents/metal_absorbers/'
#         outputFile = '/Users/David/Research_Documents/metal_absorbers/metal_absorbers.csv'
        pass

    elif user == "frenchd":
#         targetFile = '/Users/frenchd/Research/fullListMaps/LG_correlation_combined5_11_25cut_edit4.csv'
#         saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/'
#         outputFile = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/data/test.csv'
#         saveDirectory = '/Users/frenchd/Research/fullListMaps/'
#         outputFile = '/Users/frenchd/Research/fullListMaps/fullListMaps.csv'

#         targetFile = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/include_maps/salt_sightlines_all.csv'
#         saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/include_maps/'
#         outputFile = '/Users/frenchd/Research/inclination/git_inclination/rotation_paper/include_maps/salt_sightlines_all_results.csv'
#         saveDirectory = '/Users/frenchd/Research/test/'
#         outputFile = '/Users/frenchd/Research/test/test.csv'

        targetFile = '/Users/frenchd/Research/inclination/git_inclination/LG_correlation_combined5_12_edit_plusSALTcut.csv'
#         saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/maps/'
#         outputFile = '/Users/frenchd/Research/inclination/git_inclination/maps/LG_correlation_combined5_14.csv'
        
        saveDirectory = '/Users/frenchd/Research/test/'
        outputFile = '/Users/frenchd/Research/test/test.csv'
    else:
        print "Unknown user: ",user
        sys.exit()
    
    # what are the column names in this file for the AGN name and absorption velocity?
    targetHeader = 'Target'
    velocityHeader = 'Lya_v'
#     z_targetHeader = 'z_target'
#     v_limitsHeader = 'v_limits'
    z_targetHeader = 'AGNredshift'
    v_limitsHeader = 'vlimits'
    Lya_WHeader = 'Lya_W'
    NaHeader = 'Na'
    bHeader = 'b'
    identifiedHeader = 'identified'
    AGN_coordsHeader = ('degreesJ2000RA_DecAGN')

    # targets from a file, use this:
    targets = buildFullTargetList(targetFile,targetHeader,velocityHeader)

    
    # or build up a custom list of AGN names and absorption velocities here:
#     targets = [('1H0419-577',0.003678*c,True),\
#     ('3C273.0',0.005277*c,True)]

#     targets = [('MCG-03-58-009',9015,True)]
#     targets = [('NGC0891',528,True)]
#     targets = [('NGC3633',2587,True)]

#     targets = [('NGC4939',3093,True)]
#     targets = [('RX_J1217.2+2749',1326,True)]
#     targets = [('NGC6140',910,True)]

    


#     targets = [('CGCG039-137',6918,True),\
#     ('ESO343-G014',9139,True),\
#     ('IC5325',1512,True),\
#     ('MCG-03-58-009',9015,True),\
#     ('NGC1566',1502,True),\
#     ('NGC3513',1204,True),\
#     ('NGC3633',2587,True),\
#     ('NGC3640',1298,True),\
#     ('NGC4536',1867,True),\
#     ('NGC4939',3093,True),\
#     ('NGC5364',1238,True),\
#     ('NGC5786',2975,True),\
#     ('UGC09760',2094,True)]
    
    
#     targets = [('MRK279',9294,True),\
#     ('PG0838+770',721,True),\
#     ('PG0838+770',2911,True),\
#     ('SDSSJ104335.90+115129.0',717,True),\
#     ('SDSSJ104335.90+115129.0',882,True),\
#     ('SDSSJ104335.90+115129.0',1030,True),\
#     ('MRK504',9706,True),\
#     ('2dFGRS_S393Z082',1241,True),\
#     ('US2816',2848,True),\
#     ('2E1530+1511',1953,True),\
#     ('2E1530+1511',1795,True),\
#     ('FBQSJ1134+2555',9552,True),\
#     ('FBQSJ1134+2555',6343,True),\
#     ('FBQSJ1134+2555',3070,True),\
#     ('RX_J2139.7+0246',9219,True),\
#     ('RX_J2139.7+0246',4181,True),\
#     ('RX_J2139.7+0246',4083,True),\
#     ('CSO1124',1653,True),\
#     ('HE0241-3043',1219,True),\
#     ('HE0241-3043',1310,True),\
#     ('CSO327',1812,True),\
#     ('RX_J1303.7+2633',8955,True),\
#     ('RX_J1303.7+2633',7853,True)]

#     targets = [('RX_J2139.7+0246',9219,True),\
#     ('RX_J2139.7+0246',4181,True),\
#     ('RX_J2139.7+0246',4083,True),\
#     ('CSO1124',1653,True),\
#     ('HE0241-3043',1219,True),\
#     ('HE0241-3043',1310,True),\
#     ('CSO327',1812,True),\
#     ('RX_J1303.7+2633',8955,True),\
#     ('RX_J1303.7+2633',7853,True)]

##########################################################################################
    # THINGS galaxies
    
#     targets = [('NGC2403',133,True),\
#     ('MCG-02-07-026',2102,True),\
#     ('NGC4789A',374,True),\
#     ('NGC0628',657,True),\
#     ('NGC0925',553,True),\
#     ('NGC2841',638,True),\
#     ('NGC2903',550,True),\
#     ('NGC2976',3,True),\
#     ('NGC3077',14,True),\
#     ('NGC3184',592,True),\
#     ('NGC3198',660,True),\
#     ('NGC3351',778,True),\
#     ('NGC3521',801,True),\
#     ('NGC3627',727,True),\
#     ('NGC4214',291,True),\
#     ('NGC4449',207,True),\
#     ('NGC4736',308,True),\
#     ('NGC4826',408,True),\
#     ('NGC5055',484,True),\
#     ('NGC5194',463,True),\
#     ('NGC5236',513,True),\
#     ('NGC5457',241,True),\
#     ('NGC6946',40,True),\
#     ('NGC7331',816,True),\
#     ('NGC7793',230,True)]

#     targets = [('NGC4579',1517,True)]
#     targets = [('NGC3198',660,True)]
#     targets = [('NGC4414',716,True)]
#     targets = [('NGC5907',667,True)]
#     targets = [('NGC0973',4855,True)]
#     targets = [('UGC04277',5459,True)]
#     targets = [('NGC5529',2875,True)]
#     targets = [('NGC4157',774,True)]
#     targets = [('NGC4565',1230,True)]
#     targets = [('NGC3982',1109,True)]
#     targets = [('NGC4527',1736,True)]
#     targets = [('NGC3079',1116,True)]
#     targets = [('UGC05272',513,True)]
#     targets = [('UGC06399',791,True)]
#     targets = [('UGC06446',645,True)]
#     targets = [('NGC3972',852,True)]
#     targets = [('NGC3985',948,True)]
#     targets = [('UGC07089',770,True)]
#     targets = [('NGC4218',975,True)]
#     targets = [('NGC4455',637,True)]
#     targets = [('NGC5033',875,True)]
#     targets = [('UGC09211',686,True)]
#     targets = [('UGC09211',686,True)]
#     targets = [('NGC3877',895,True)]
#     targets = [('NGC3893',967,True)]
#     targets = [('NGC3718',993,True)]
#     targets = [('NGC0891',528,True)]
#     targets = [('NGC4529',2536,True)]
#     targets = [('UGC04238',1544,True)]


#     targets = [('NGC2770',1947,True),\
#     ('NGC3432',616,True),\
#     ('NGC3666',1060,True),\
#     ('NGC3769',737,True),\
#     ('NGC3949',800,True),\
#     ('NGC4157',774,True),\
#     ('NGC4414',716,True),\
#     ('NGC4534',802,True),\
#     ('NGC5951',1780,True),\
#     ('NGC7741',750,True),\
#     ('NGC7817',2309,True),\
#     ('UGC05459',1112,True),\
#     ('UGC08146',670,True)]


#     targets = [('NGC4238',2762,True)]
#     targets = [('NGC3351',778,True)]
#     targets = [('NGC4254',2407,True)]
#     targets = [('NGC4559',807,True)]

#     targets = [('NGC3726',866,True)]
#     targets = [('NGC3067',1476,True)]
#     targets = [('PG1302-102',3447,True)]
#     targets = [('RX_J1142.5+2503',550,True)]
#     targets = [('2E1530+1511',1953,True)]
#     targets = [('MRK335',1954,True)]
#     targets = [('RX_J1236.0+2641',794,True),
#                 ('RX_J1236.0+2641',1009,True),
#                 ('RX_J1236.0+2641',1166,True),
#                 ('RX_J1236.0+2641',1254,True)]

#     targets = [('SDSSJ112448.30+531818.0',645,True),
#                 ('SDSSJ112448.30+531818.0',1156,True)]

#     targets = [('SBS1116+523',731,True)]

#     targets = [('CSO1208',874,True)]
#     targets = [('MRK876',939,True)]
#     targets = [('PG0804+761',1143,True)]

#     targets = [('SDSSJ112439.50+113117.0',1047,True)]
#     targets = [('SDSSJ112632.90+120437.0',1060,True)]
#     targets = [('MRK335',1954,True)]
#     targets = [('PG1259+593',670,True)]
#     targets = [('SDSSJ101622.60+470643.0',661,True)]
#     targets = [('RBS1503',667,True)]
#     targets = [('SBS1503+570',667,True)]
#     targets = [('SDSSJ112448.30+531818.0',1019,True), ('SDSSJ112448.30+531818.0',1141,True)]

#     targets = [('SDSSJ112632.90+120437.0',1060,True)]
#     targets = [('2E1530+1511',1953,True)]

#     targets = [('NGC2770',1947,True)]
#     targets = [('MRK876',939,True),\
#     ('MRK876',3478,True),\
#     ('MRK876',4508,True),\
#     ('MRK876',5036,True),\
#     ('MRK876',6037,True),\
#     ('MRK876',7005,True),\
#     ('MRK876',9895,True)]
    
    
    targets = [('1H0419-577',1075,True),\
            ('1H0419-577',1123,True),\
            ('1H0419-577',1188,True),\
            ('1H0419-577',1264,True),\
            ('1H0419-577',2020,True)]







    
    c = 0
    for i in targets:
        # find AGN environment using the imported version of correlateSingle
        targetName,center,include = i
        correlation = correlateSingle.correlateTarget(targetName, maxSep, agnSeparation, minVcorr, minSize, slow=False,searchAll=True)
        galaxyInfo = correlation[targetName]
        
#         print 'galaxyInfo: ',galaxyInfo
                
        if includeAGN:
            correlationAGN = correlateSingle.correlateGalaxy(targetName, maxSep, agnSeparation, minVcorr, minSize)
            
            if correlationAGN:
                AGNinfo = correlationAGN[targetName]
            else:
                AGNinfo = []
        
#         print '{0} = {1}'.format(targetName,len(galaxyInfo))
#         print

#         galaxyInfo.sort()

        # instantiate some lists for later
        galaxyNames = []
        separations = []
        positions = []
        plotPositionsRA = []
        plotPositionsDec = []
        plotRA = []
        plotDec = []
        plotAGNposition = []
        plotSizes = []
        plotVelocity = []
        PA = []
        inc = []
        typeList = []
        infoDict = {}
        
        AGNRAs = []
        AGNDecs = []
        plotPositionAGNRA = []
        plotPositionAGNDec = []
        AGNnames = []
        
        if includeAGN:
            for r in AGNinfo:
                vhel, AGNrow = r
                            
                AGNRA = AGNrow['AGNRA']
                AGNRAs.append(AGNRA)
                AGNDec = AGNrow['AGNDec']
                AGNDecs.append(AGNDec)
            
                galaxyDist = AGNrow['bestDist']
            
                targetRA = float(AGNrow['galaxyRA'])
                targetDec = float(AGNrow['galaxyDec'])
                
                AGNname = AGNrow['AGNname']
                AGNnames.append(AGNname)
            
    
                #find plot placement
                # calculate angular separations in ra and dec to determine positions on chart w.r.t. target
                AGNRA,AGNDec = float(AGNRA),float(AGNDec)
            
                # calculate separation in RA only
                dRA_agn = calculateImpactParameter(AGNRA,targetDec,targetRA,targetDec,galaxyDist)

                # calculate separation in Dec only
                dDec_agn = calculateImpactParameter(targetRA,AGNDec,targetRA,targetDec,galaxyDist)
                
            
                # add signs back into physical impact parameters
                if AGNRA < targetRA:
                    if RAeastLeft:
                        dRA_agn = -dRA_agn
                
                if AGNRA > targetRA:
                    if not RAeastLeft:
                        dRA_agn = -dRA_agn
                
                # 'edge' effects
                if AGNRA >= 359.0 and targetRA <= 1.0:
                    dRA_agn = -dRA_agn

                if AGNDec < targetDec:
                    dDec_agn = -dDec_agn
                
                plotPositionAGNRA.append(dRA_agn)
                plotPositionAGNDec.append(dDec_agn)
                
                print 'dRA_agn: ',dRA_agn
                print 'dDec_agn: ',dDec_agn
                print

                
        # loop through the returned galaxy environment data, making calculations and
        # populating lists as we go
        counter = 0
        for row in galaxyInfo:
            counter+=1
            vhel, galaxyRow = row
            AGNposition = eval(str(galaxyRow['degreesJ2000RA_DecAGN']))

            
            # crop off results that fall out of the 'velocityWindow' parameter
            # 'velocityWindow' = a cut in velocity space 
            if counter <= maxPlotObjects and float(vhel)-velocityWindow <= center and float(vhel)+velocityWindow >= center:
                separation = galaxyRow['impactParameter (kpc)']
                galaxyName = galaxyRow['galaxyName']
                galaxyPosition = eval(str(galaxyRow['degreesJ2000RA_DecGalaxy']))
                AGNposition = eval(str(galaxyRow['degreesJ2000RA_DecAGN']))
                galaxyDist = galaxyRow['distGalaxy (Mpc)']
                group = eval(str(galaxyRow['groups_dist_std (Mpc)']))
                major,minor = eval(str(galaxyRow['linDiameters (kpc)']))
                galaxyVcorr = galaxyRow['vcorrGalaxy (km/s)']
                galaxyVel = galaxyRow['radialVelocity (km/s)']
                morphology = galaxyRow['morphology']
                RC3Type = galaxyRow['RC3type']
                RC3PA = galaxyRow['RC3pa (deg)']
                RC3inc = galaxyRow['RC3inc (deg)']
                positionAngle = galaxyRow['positionAngle (deg)']
                inclination = galaxyRow['inclination (deg)']
                azimuth = galaxyRow['azimuth (deg)']
#                 include = galaxyRow['include']
                
                if not isNull(inclination):
                    inclination = round(eval(inclination),0)

                if not isNull(galaxyDist):
                    positions.append(galaxyPosition)
                    separations.append(float(separation))
                    
                    if not isNull(positionAngle):
                        pa = float(positionAngle)
                    elif not isNull(RC3PA):
                        pa = float(RC3PA)
                    else:
                        pa = 0
                        
                    #find plot placement
                    # calculate angular separations in ra and dec to determine positions on chart w.r.t. target AGN
                    gRA,gDec = float(galaxyPosition[0]),float(galaxyPosition[1])
                    agnRA,agnDec = float(AGNposition[0]),float(AGNposition[1])

                    
                    # calculate separation in RA only
#                     dRA = correlateSingle.calculateImpactParameter_slow(gRA,agnDec,agnRA,agnDec,galaxyDist)
                    dRA = calculateImpactParameter(gRA,agnDec,agnRA,agnDec,galaxyDist)

                    
                    # calculate separation in Dec only
#                     dDec = correlateSingle.calculateImpactParameter_slow(agnRA,gDec,agnRA,agnDec,galaxyDist)
                    dDec = calculateImpactParameter(agnRA,gDec,agnRA,agnDec,galaxyDist)
                    
                    # add signs back into physical impact parameters
#                     if gRA < agnRA:
#                         if not RAeastLeft:
# #                             dRA = -dRA
#                             dRA = dRA

                    if gRA < agnRA:
                        if RAeastLeft:
                            dRA = -dRA
                
                    if gRA > agnRA:
                        if not RAeastLeft:
                            dRA = -dRA

                    # 'edge' effects
                    if gRA >= 359.0 and agnRA <= 1.0:
                        dRA = -dRA
                        
                    if gDec < agnDec:
                        dDec = -dDec
                        
                        
                    # calculate size by finding radius
                    noSize = False
                    if not isNull(major):
                        averageSize = float(major)
                    elif not isNull(minor) and not isNull(inclination):
                        averageSize = float(minor) / math.cos(float(inclination) * math.pi/180)
                    elif not isNull(minor) and isNull(inclination) and isNull(major):
                        averageSize = float(minor)
                        inclination = 0
                    else:
                        averageSize = 2
                        inclination = 0
                        noSize = True
                                            
#                     if not isNumber(galaxyVcorr):
#                         galaxyVcorr = 0
#                         galaxyVel = 0

                    localType = 'x'
                    rc3Type = RC3Type.lower()
                    r = rc3Type[:10]
                    morphology = morphology.lower()
                    m = morphology[:10]
                    print 'm vs morphology: ',m, ' vs ',morphology


                    if not isNull(morphology):
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
                        
                        elif bfind(m,'ir') or bfind(m,'im') or bfind(m,'i ') or bfind(m,'ia'):
                            localType = 's'
                        
                        else:
                            localType = 'x'
                    
                    # try RC3 types
                    elif not isNull(rc3Type):
                        print 'Not null rc3 morphology: ',morphology
                        print
                        if bfind(r,'s'):
                            if not bfind(r,'s0'):
                                # straight spiral type
                                localType = 's'
                            else:
                                # lenticular or S0 type
                                localType = 'e'
                        elif bfind(r,'e') or bfind(r,'dwarf') or bfind(r,'pec'):
                            localType = 'e'
                        elif bfind(r,'len'):
                            localType = 'e'
                            
                        elif bfind(r,'ir') or bfind(r,'im') or bfind(r,'i ') or bfind(r,'ia'):
                            # irregular type
                            localType = 's'
                        else:
                            localType = 'x'
                            
                    # Latex format galaxyName
                    galaxyName = r'{0}'.format(galaxyName)
                    galaxyName = galaxyName.replace('_','\_')
                            
                    infoDict[galaxyName]=galaxyRow
                    
                    galaxyType = {'s':2,'e':1,'x':3}

                    
                    if noSize:
                        galaxyNames.append('*'+galaxyName)
                    else:
                        galaxyNames.append(galaxyName)
                    
                    print 'galaxyName - morph - localType: ',galaxyName, ' - ',morphology,' - ',localType
                    print
                        
                    plotPositionsRA.append(dRA)
                    plotPositionsDec.append(dDec)
                    plotSizes.append(averageSize)
                    plotVelocity.append(float(galaxyVel))
                    typeList.append(galaxyType[localType])
                    PA.append(pa)
                    inc.append(inclination)
                    plotRA.append(gRA)
                    plotDec.append(gDec)
                    plotAGNposition.append(AGNposition)
                    

        if len(galaxyNames) >=1:                        
            ##################################
            print 'starting first plot...'
            
            x = arange(len(galaxyNames))+1
            fig = figure(figsize=(9,7))
            ax = fig.add_subplot(111)
            width = 0.30
            
            # format the axes:
            #
            # x-axis
            majorLocator   = MultipleLocator(xAxisMajorTicks)
            majorFormatter = FormatStrFormatter(r'$\rm %d$')
            minorLocator   = MultipleLocator(xAxisMinorTicks)
            ax.xaxis.set_major_locator(majorLocator)
            ax.xaxis.set_major_formatter(majorFormatter)
            ax.xaxis.set_minor_locator(minorLocator)
        
            # y axis
            majorLocator   = MultipleLocator(yAxisMajorTicks)
            majorFormatter = FormatStrFormatter(r'$\rm %d$')
            minorLocator   = MultipleLocator(yAxisMinorTicks)
            ax.yaxis.set_major_locator(majorLocator)
            ax.yaxis.set_major_formatter(majorFormatter)
            ax.yaxis.set_minor_locator(minorLocator)            
            
#             ax.xaxis.set_tick_params(labelweight='normal')
#             ax.yaxis.set_tick_params(labelweight='normal')
            
            # scale sizes and velocities
            maxSize = 300
            minSize = 80
            largest = float(max(plotSizes))
            smallest = float(min(plotSizes))
            
            newSizes = []
        
            # multiply sizes of galaxies by 10
            for s in plotSizes:
                new = s*10
                newSizes.append(new)
            
            vmaxVal = velocityWindow
            vminVal = -velocityWindow

            # +/- 400 km/s around the center
            largestVelocity = velocityWindow
            smallestVelocity = -velocityWindow

            newVelocities = []
            # check if there's more than one
            for v in plotVelocity:
                # convert to delta-v = v_absorber - v_galaxy (neg = absorber is blueward of galaxy)
#                 velocity = center - v
                velocity = v - center

#                 newVelocity = ((float(velocity) - smallestVelocity)/(largestVelocity-smallestVelocity)) * (vmaxVal-0)+0
#                 newVelocities.append(newVelocity)

                newVelocities.append(velocity)

            
#             rounding = -1
#             step = round(((largestVelocity-smallestVelocity)/vmaxVal),rounding)
#             ticks = arange(int(round(smallestVelocity,-2)),int(round(largestVelocity,-2))+int(step),int(step))

            rounding = -1
            step = velocityStepSize
            ticks = arange(-velocityWindow,velocityWindow+step,int(step))
            
            norm = matplotlib.colors.Normalize(vmin = vminVal, vmax = vmaxVal)
            m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colmap)

            numberGalaxies = len(newVelocities)
            galaxyIndices = range(numberGalaxies)
            
            # make ellipses to indicate position angle and inclination at once
            if len(plotVelocity) !=0:
                mapping = dict((i,n) for i,n in zip(galaxyIndices,galaxyNames))
                
                for i in range(numberGalaxies):
                    if typeList[i] == 1 or typeList[i] == 2:
                        # typeList[i] == 2 is a spiral, 1 is E type
                        print 'i: ',i,galaxyIndices[i]
                        if RAeastLeft:
                            e = Ellipse(xy=(plotPositionsRA[i],plotPositionsDec[i]),\
                            width=newSizes[i] * math.cos(float(inc[i]) * math.pi/180)/2,\
                            height=float(newSizes[i])/2, angle=float(-PA[i]))
                            print 'PA = ',PA
                            print 'angle=float(-PA[i]) : ',float(-PA[i])
                            print
                        
                        else:
                            e = Ellipse(xy=(plotPositionsRA[i],plotPositionsDec[i]),\
                            width=newSizes[i] * math.cos(float(inc[i]) * math.pi/180)/2,\
                            height=float(newSizes[i])/2, angle=float(PA[i]))
                            
                        # no transparency
                        e.set_alpha(0.9)
                        
                        if typeList[i] == 2:
                            # spiral - edge color is black
                            
                            ax.add_artist(e)
                            e.set_facecolor(m.to_rgba(newVelocities[i]))
                            e.set_edgecolor('black')
                        
                        if typeList[i] == 1:
                            # elliptical - edge color is same as galaxy
                            
                            ax.add_artist(e)
                            e.set_facecolor(m.to_rgba(newVelocities[i]))
                            e.set_edgecolor(m.to_rgba(newVelocities[i]))
                        
                    else:
                        if RAeastLeft:
                            ax.scatter(plotPositionsRA[i],plotPositionsDec[i],\
                            s=newSizes[i]*4,c=newVelocities[i],vmin=vminVal,vmax=vmaxVal,\
                            marker=(2,1,float(PA[i])),lw=1,cmap=colmap)
                            
                        else:
                            ax.scatter(plotPositionsRA[i],plotPositionsDec[i],\
                            s=newSizes[i]*4,c=newVelocities[i],vmin=vminVal,vmax=vmaxVal\
                            ,marker=(2,1,float(-PA[i])),lw=1,cmap=colmap)

                        
                    # annotate with galaxy names if includeNameTags == True
                    if includeNameTags:
                        if galaxyNames[i].find('*')!=-1:
                            # this indicates no size data is available
#                             plt.annotate('*'+str(galaxyNames[i]),xy=(plotPositionsRA[i],\
#                             plotPositionsDec[i]),xytext=(xTagOffset,newSizes[i]/yTagOffset),textcoords='offset points',size=nameTagFont)
#                             newgalaxyNames = galaxyNames[i].replace('_','\_')
                            if bfind(galaxyNames[i],'_') and not bfind(galaxyNames[i],'\_'):
                                newgalaxyNames = galaxyNames[i].replace('_','\_')
                            else:
                                newgalaxyNames = galaxyNames[i]

                            plt.annotate('*'+str(newgalaxyNames),xy=(plotPositionsRA[i],\
                            plotPositionsDec[i]),xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=nameTagFont)

#                             plt.annotate('*'+str(galaxyNames[i]),xy=(plotPositionsRA[i],\
#                             plotPositionsDec[i]),xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=nameTagFont)

                        else:
#                             plt.annotate(galaxyNames[i],xy=(plotPositionsRA[i],plotPositionsDec[i]),\
#                             xytext=(xTagOffset,newSizes[i]/yTagOffset),textcoords='offset points',size=nameTagFont)
                            plt.annotate(galaxyNames[i],xy=(plotPositionsRA[i],plotPositionsDec[i]),\
                            xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=nameTagFont)
            
            plot1 = ax.scatter(plotPositionsRA,plotPositionsDec,s=0,c=newVelocities,vmin=vminVal,vmax=vmaxVal,marker='.',cmap=colmap)
            
            # if we chose to plot nearby AGN as well, plot them
            if includeAGN:
                for i in range(len(plotPositionAGNRA)):
                    plot_agn = ax.scatter(plotPositionAGNRA[i],plotPositionAGNDec[i],s=AGNsize,c='green',marker='*',lw=0.5)
                
                    if includeAGNnameTags:
                        if bfind(AGNnames[i],'_') and not bfind(AGNnames[i],'\_'):
                            newAGNname = AGNnames[i].replace('_','\_')
                        else:
                            newAGNname = AGNnames[i]
#                         plt.annotate(AGNnames[i],xy=(plotPositionAGNRA[i],plotPositionAGNDec[i]),\
#                         xytext=(xTagOffset,AGNsize/yTagOffset),textcoords='offset points',size=nameTagFont)
                        plt.annotate(newAGNname,xy=(plotPositionAGNRA[i],plotPositionAGNDec[i]),\
                        xytext=(xTagOffset,yTagOffset),textcoords='offset points',size=nameTagFont)


#########################################################################################
#########################################################################################
            # prepare all the stuff for writing. Deciding to write happens later

            fieldnames = ('targetName',\
            'center',\
            'galaxyName',\
            'environment',\
            'degreesJ2000RA_DecAGN',\
            'degreesJ2000RA_DecGalaxy',\
            'likelihood',\
            'likelihood_1.5',\
            'virialRadius',\
            'd^1.5',\
            'impactParameter (kpc)',\
            'redshiftDistances',\
            'vcorrGalaxy (km/s)',\
            'radialVelocity (km/s)',\
            'vel_diff',\
            'distGalaxy (Mpc)',\
            'AGN S/N',\
            'majorAxis (kpc)',\
            'minorAxis (kpc)',\
            'inclination (deg)',\
            'positionAngle (deg)',\
            'azimuth (deg)',\
            'RC3flag',\
            'RC3type',\
            'RC3inc (deg)',\
            'RC3pa (deg)',\
            'morphology',\
            'final_morphology',\
            'galaxyRedshift')
            
            virList = []
            customList = []
            environment = numberGalaxies
            
            for number in range(numberGalaxies):
                if mapping[number].find("*")!=-1:
                    infoRow = infoDict[mapping[number].strip('*')]
                else:
                    infoRow = infoDict[mapping[number]]
                    
                inc = infoRow['inclination (deg)']
                if not isNull(inc):
                    inc = round(float(inc),1)

                pa = infoRow['positionAngle (deg)']
                RC3pa = infoRow['RC3pa (deg)']
                if not isNull(pa):
                    pa = round(float(pa),1)
                if not isNull(RC3pa):
                    RC3pa = round(float(RC3pa),1)
                    
                az = infoRow['azimuth (deg)']
                if isNumber(az):
                    az = round(float(az),1)
                else:
                    az = 'x'
                    
                vcorr = infoRow['vcorrGalaxy (km/s)']
                vhel = infoRow['radialVelocity (km/s)']
                if not isNull(vcorr):
                    vcorr = round(float(vcorr),1)
                    vhel = round(float(vhel),1)
                    vel_dif = vhel - float(center)
                else:
                    vcorr = 'x'
                    vhel = 'x'
                    vel_dif = 'x'
                
                try:
                    major,minor = eval(infoRow['linDiameters (kpc)'])
                except Exception,e:
                    major,minor = infoRow['linDiameters (kpc)']
                
                
                if not isNull(major):
                    major = round(float(major),1)
                if not isNull(minor):
                    minor = round(float(minor),1)
    
                agn_sn = 'x'
                corrected_az = az
                redshiftDistance = 'x'
                impact = float(infoRow['impactParameter (kpc)'])
                
                if not isNull(major) and float(major) >= 0:
                    rVir = calculateVirialRadius(major)

                    # try this "sphere of influence" value instead
                    m15 = major**1.5

                    # first for the virial radius
                    likelihood = math.exp(-(impact/rVir)**2) * math.exp(-(vel_dif/200.)**2)
                    
                    if rVir>= impact:
                        likelihood = likelihood*2
                        
                    # then for the second 'virial like' m15 radius
                    likelihoodm15 = math.exp(-(impact/m15)**2) * math.exp(-(vel_dif/200.)**2)
                    
                    if m15>= impact:
                        likelihoodm15 = likelihoodm15*2
                        
                    # should be like 33% at 1R_v, and linear down after that.
                    # 66% at 0.5R_v, something quadratic, or logarithmic
                    # look up the "sphere of influence" for where the probability drops
                    # down to 10%
                    
                    # use M/L, surpress environment by M/L ratios
                    # also include velocity difference to make the virial radius 3D - 
                    # really a 3D impact parameter
                    
                else:
                    likelihood = nullFloat
                    likelihoodm15 = nullFloat
                    rVir = nullFloat
                    m15 = nullFloat
                
                objectInfoList = [targetName,\
                center,\
                infoRow['galaxyName'],\
                environment,\
                infoRow['degreesJ2000RA_DecAGN'],\
                infoRow['degreesJ2000RA_DecGalaxy'],\
                likelihood,\
                likelihoodm15,\
                rVir,\
                m15,\
                impact,\
                redshiftDistance,\
                vcorr,\
                vhel,\
                vel_dif,\
                infoRow['distGalaxy (Mpc)'],\
                agn_sn,\
                major,\
                minor,\
                inc,\
                pa,\
                az,\
                infoRow['RC3flag'],\
                infoRow['RC3type'],\
                RC3inc,\
                RC3pa,\
                infoRow['morphology'],\
                infoRow['morphology'],\
                infoRow['galaxyRedshift']]
                
                virList.append([likelihood,objectInfoList])
                customList.append([likelihoodm15,objectInfoList])
                
                
                
            # now sort the list of galaxies by likelihood:
            sorted_virList = sorted(virList,reverse=True)
            sorted_cusList = sorted(customList,reverse=True)
            
            print targetName,' - ',center,': ', sorted_virList
            print 'and : ', sorted_cusList
            
            # if there are enough galaxies, compare the best two : VIRIAL
            if len(sorted_virList) >=2:
                first_vir = sorted_virList[0][0]
                second_vir = sorted_virList[1][0]
            
                # the most likely galaxy must have rigor * the likelihood of the #2 galaxy
                if second_vir*rigor <= first_vir:
                    # larger than the min?
                    if first_vir >= l_min:
                        include_vir = True
                    else:
                        include_vir = False
                else:
                    include_vir = False
            
            
            # include lone galaxies if loner = TRUE, otherwise same l_min requirement
            elif len(sorted_virList) == 1:
                first_vir = sorted_virList[0][0]
                
                if loner:
                    # include by default
                    include_vir = True
                    
                elif first_vir >= l_min:
                    # include only if above l_min threshold
                    print 'first_vir >= l_min: ',first_vir
                    print 'targetName, center = ',targetName, center
                    include_vir = True
                else:
                    include_vir = False
                    print 'first_vir < l_min: ',first_vir
                    print 'targetName, center = ',targetName, center
            
            else:
                # if there are no galaxies, don't include
                include_vir = False
                

            # if there are enough galaxies, compare the best two : CUSTOM
            if len(sorted_cusList) >=2:
                first_cus = sorted_cusList[0][0]
                second_cus = sorted_cusList[1][0]
            
                # the most likely galaxy must have rigor * the likelihood of the #2 galaxy
                if second_cus*rigor <= first_cus:
                    # above the hard lower limit?
                    if first_cus >= l_min:
                        include_cus = True
                    else:
                        include_cus = False
                else:
                    include_cus = False                
            
            # include lone galaxies if loner = TRUE
            elif len(sorted_cusList) == 1:
                first_cus = sorted_cusList[0][0]
                
                if loner:
                    # include by default
                    include_cus = True
                    
                elif first_cus >= l_min:
                    # include only if above l_min threshold
                    print 'first_cus >= l_min: ',first_cus
                    print 'targetName, center = ',targetName, center
                    include_cus = True
                else:
                    # don't include
                    include_cus = False
                    print 'first_cus < l_min: ',first_cus
                    print 'targetName, center = ',targetName, center
            
            else:
                # if there are no galaxies, don't include
                include_cus = False
                
                
            # where to put it?
            if sortIntoFolders:
                if include_cus and include_vir:
                    include_folder = 'associated'
                elif include_cus or include_vir:
                    include_folder = 'ambiguous'
                    print 'include_cus, include_vir = ',include_cus, include_vir
                else:
                    include_folder = 'nonassociated'
            else:
                include_folder = ''
                
            # append this result to the master lists
            masterVirList.append(sorted_virList)
            masterCustomList.append(sorted_cusList)
                
            # include a simple position skyplot?
            if includeSkyPlot:
                # plot
                fig = figure(figsize=(12,10))
                ax = fig.add_subplot(111)
                plot1 = ax.scatter(plotRA,plotDec,marker='d')
                plot2 = ax.scatter(float(AGNposition[0]),float(AGNposition[1]),marker='*')
                ax.set_xlabel('RA')
                ax.set_ylabel('Dec')
                
                if RAeastLeft:
                    ax.invert_yaxis()
        
                savefig('{0}{1}/map2_{2}_{3}_simple.pdf'.format(saveDirectory,include_folder,targetName,center),format='pdf')

            # save the map plot tables
            if saveMapTables:
                writerOutFile = open('{0}{1}/map_{2}_{3}_table.csv'.format(saveDirectory,include_folder,targetName,center),'wt')

                writer = csv.DictWriter(writerOutFile, fieldnames=fieldnames)
                headers = dict((n,n) for n in fieldnames)
                writer.writerow(headers)
               
                for s in sorted_virList:
                    likelihood, rest = s
                    row = dict((f,o) for f,o in zip(fieldnames,rest))
                    writer.writerow(row)
    
                writerOutFile.close()
            
            else:
                for s in sorted_virList:
                    likelihood, rest = s
                    row = dict((f,o) for f,o in zip(fieldnames,rest))
                    print
                    print 'L = {0} : {1}, D = {2}, dv = {3}'.format(likelihood,row['galaxyName'],row['majorAxis (kpc)'],row['vel_diff'])
            
##########################################################################################
##########################################################################################

            # plot AGN target in center
            ax.scatter(0,0,s=200,c='black',marker='*')
    
#             cbar = plt.colorbar(plot1,ticks=range(0,21),cmap=colmap,orientation='vertical')
#             cbar = plt.colorbar(plot1,ticks=ticks,cmap=colmap,orientation='vertical')
            
            cbar = plt.colorbar(plot1,ticks=ticks,format=r'$\rm %d$',cmap=colmap,orientation='vertical')
            
#             cbar.ax.set_yticklabels([r'$\rm %d$' % i for i ticks])
#             cbar.ax.set_yticklabels(['%d' % i for i ticks])
#             ax.yaxis.set_ticklabels(['%.2f' % 0.1/100*i for i in np.arange(0,100,10)]) 
            cbar.set_label(r'$\rm \Delta v ~[km ~s^{-1}]$')
        
            ax.grid(b=None,which='major',axis='both')
            ax.set_ylim(-maxSep,maxSep)
            ax.set_xlim(-maxSep,maxSep)
            
            if RAeastLeft:
                ax.invert_xaxis()
            
            ax.set_xlabel(r'$\rm R.A. ~Separation ~[kpc]$')
            ax.set_ylabel(r'$\rm Dec. ~Separation ~[kpc]$')
            if includeTitle:
                if bfind(targetName,'_') and not bfind(targetName,'\_'):
                    newTargetName = targetName.replace('_','\\_')
                else:
                    newTargetName = targetName
#                 title("{0} centered velocity = {1} +/- {2} km/s".format(targetName,center,velocityWindow))
                title(r'$\rm {0} ~ centered ~ velocity = {1} +/- {2}~ km/s $'.format(targetName,center,velocityWindow))

            # now write it all to file, or display the finished figure
            if saveMaps:
                savefig('{0}{1}/map_{2}_{3}.pdf'.format(saveDirectory,include_folder,targetName,center),\
                bbox_inches='tight',format='pdf')
            else:
                show()
        
        else:
            print 'TARGET: {0} = NO GALAXIES!!'.format(i)
            
            
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
            
            
    def writeFullResults(listVir,listCustom):
        # write the full results to a new file,  specified by 'outputFile'
        # 
        # this new file will contain all results, both positive and negative
        '''
            list contains entries like [likelihood, [all]], where all is:
            
                ('targetName',\
                'center',\
                'galaxyName',\
                'environment',\
                'degreesJ2000RA_DecAGN',\
                'degreesJ2000RA_DecGalaxy',\
                'likelihood',\
                'likelihood_1.5'
                'virialRadius',\
                'd^1.5',\
                'impactParameter (kpc)',\
                'redshiftDistances',\
                'vcorrGalaxy (km/s)',\
                'vel_diff',\
                'distGalaxy (Mpc)',\
                'AGN S/N',\
                'majorAxis (kpc)',\
                'minorAxis (kpc)',\
                'inclination (deg)',\
                'positionAngle (deg)',\
                'azimuth (deg)',\
                'RC3flag',\
                'RC3type',\
                'RC3inc (deg)',\
                'RC3pa (deg)',\
                'morphology',\
                'galaxyRedshift')
                
            listVir is sorted by likelihood, listCustom is sorted by custom likelihood
                
        '''
        
        # open the origin file
        f = open(targetFile,'rU')
        reader = csv.DictReader(f)
        
        
        # fieldnames for the new outputFile
        fieldnames = ('targetName',\
        'center',\
        'galaxyName',\
        'environment',\
        'degreesJ2000RA_DecAGN',\
        'degreesJ2000RA_DecGalaxy',\
        'likelihood',\
        'likelihood_1.5',\
        'virialRadius',\
        'd^1.5',\
        'impactParameter (kpc)',\
        'redshiftDistances',\
        'vcorrGalaxy (km/s)',\
        'radialVelocity (km/s)',\
        'vel_diff',\
        'distGalaxy (Mpc)',\
        'AGN S/N',\
        'majorAxis (kpc)',\
        'minorAxis (kpc)',\
        'inclination (deg)',\
        'positionAngle (deg)',\
        'azimuth (deg)',\
        'RC3flag',\
        'RC3type',\
        'RC3inc (deg)',\
        'RC3pa (deg)',\
        'morphology',\
        'final_morphology',\
        'galaxyRedshift',\
        'AGNredshift',\
        'spectrumStatus',\
        'include',\
        'include_vir',\
        'include_custom',\
        'Lya_v',\
        'vlimits',\
        'Lya_W',\
        'Na',\
        'b',\
        'identified',\
        'comment')
        
        # open the output file
        writerOutFile = open(outputFile,'wt')
        writer = csv.DictWriter(writerOutFile, fieldnames=fieldnames)
        headers = dict((n,n) for n in fieldnames)
        writer.writerow(headers)
        
        count = -1
#         lenList = len(listVir)
        for l in reader:
            targetName = l[targetHeader]
            Lya_v = l[velocityHeader]
            AGNredshift = l[z_targetHeader]
            spectrumStatus = 'x'
            v_limits = l[v_limitsHeader]
#             Lya_v = l[Lya_vHeader]
            Lya_W = l[Lya_WHeader]
            Na = l[NaHeader]
            b = l[bHeader]
            identified = l[identifiedHeader]
            comment = l['comment']
#             degreesJ2000RA_DecAGN = (l['RAdeg_target'],l['DEdeg_target'])
            degreesJ2000RA_DecAGN = eval(l[AGN_coordsHeader])

            print 'targetName: ',targetName
            print 'Lya_v: ',Lya_v
            print 'listVir: ',listVir
            print 'count: ',count
        
            if isNumber(Lya_v):
                # this should then correspond to an element in 'list'
                count +=1
                Lya_v = float(Lya_v)
                
                # grab the corresponding data. list[count] is all the galaxies around a
                # particular target and line
                if count < len(listVir):
                    systemsV = listVir[count]
                    systemsC = listCustom[count]
                
                    # test to make sure this is the correct one. If there were no galaxies,
                    # the count will get off
                    v1 = systemsV[0][1]
                    v1_targetName = v1[0]
                    v1_center = v1[1]
                
                    print 'v1: ',v1
                    print 'v1_targetName: ',v1_targetName
                    print 'v1_center: ',v1_center
                
                    if v1_targetName == targetName and v1_center == Lya_v:
                        match = True
                    else:
                        match = False
                        count-=1
                
                    # now find the 1st, 2nd highest likelihood for both estimates
                    sysLen = len(systemsV)
                
                    finalEntry = []
                    
                else:
                    match = False
                
                if sysLen >=2 and match:
                    # systemsV[0] is the most likely as such: [likelihood, [all]]
                    # thus systemsV[1] is the second most likely
                    vir_1 = systemsV[0]
                    vir_1_like = vir_1[0]
                    vir_1_all = vir_1[1]
                    vir_1_galaxy = vir_1_all[2]
                    
                    vir_2 = systemsV[1]
                    vir_2_like = vir_2[0]
                    vir_2_all = vir_2[1]
                    vir_2_galaxy = vir_2_all[2]
                                    
                    # now for the second, 'custom' likelihood
                    cus_1 = systemsC[0]
                    cus_1_like = cus_1[0]
                    cus_1_all = cus_1[1]
                    cus_1_galaxy = cus_1_all[2]
                    
                    cus_2 = systemsC[1]
                    cus_2_like = cus_2[0]
                    cus_2_all = cus_2[1]
                    cus_2_galaxy = cus_2_all[2]
                    
                    # now decide to marked "include" for each likelihood estimate
                    #
                    # first: for virial radius based likelihood                    
                    if vir_2_like * rigor <= vir_1_like:
                        if vir_1_like >= l_min:
                            virInclude = True
                        else:
                            virInclude = False
                    else:
                        virInclude = False
                    
                    # second: for custom based likelihood                    
                    if cus_2_like * rigor <= cus_1_like:
                        if cus_1_like >= l_min:
                            cusInclude = True
                        else:
                            cusInclude = False                            
                    else:
                        cusInclude = False
                        
                    print
                    print
                    print 'virInclude: ',virInclude
                    print 'vir_1 : ',vir_1
                    print 'vir_2 : ',vir_2
                    
                    print
                    
                    print 'cus_Include: ',cusInclude
                    print 'cus_1: ',cus_1
                    print 'cus_2: ',cus_2
                             
                    
                    # decide on the final thing to include.
                    # first, check if both methods agree that a best galaxy is found
                    if cusInclude and virInclude:
                        
                        # now check if they find the SAME galaxy
                        if vir_1_galaxy == cus_1_galaxy:
                        
                            # in this case vir_1_all and cus_1_all are the same, so just 
                            # pick one to include
                            entry = vir_1_all
                            entry.append(AGNredshift)
                            entry.append(spectrumStatus)
                            entry.append('?')
                            entry.append(virInclude)
                            entry.append(cusInclude)
                            entry.append(Lya_v)
                            entry.append(v_limits)
                            entry.append(Lya_W)
                            entry.append(Na)
                            entry.append(b)
                            entry.append(identified)
                            entry.append(comment)
                            finalEntry = [entry]
                                                    
                        
                        # if they don't find the same galaxy...
                        else:
                            # enter two lines, first for the virial result, then for 
                            # the custom one later
                            entry = vir_1_all
                            entry.append(AGNredshift)
                            entry.append(spectrumStatus)
                            entry.append('?')
                            
                            # force True for virInclude
                            entry.append(True)
                            # force False for cusInclude
                            entry.append(False)
                            
                            entry.append(Lya_v)
                            entry.append(v_limits)
                            entry.append(Lya_W)
                            entry.append(Na)
                            entry.append(b)
                            entry.append(identified)
                            
                            # update the comment to include this note
                            comment1 = str(comment) + ' : VIRIAL VS CUSTOM RESULT MISMATCH'
                            entry.append(comment1)
                            
                            
                            # now enter a second line for the custom result
                            entry2 = cus_1_all
                            entry.append(AGNredshift)
                            entry.append(spectrumStatus)
                            entry2.append('?')
                            
                            # now force False for virInclude
                            entry2.append(False)
                            # and True for cusInclude
                            entry2.append(True)
                            
                            entry.append(Lya_v)
                            entry.append(v_limits)
                            entry.append(Lya_W)
                            entry.append(Na)
                            entry.append(b)
                            entry.append(identified)
                            
                            # include the same updated comment as above
                            entry2.append(comment1)
                            
                            finalEntry = [entry,entry2]
                            
                    # if the methods do not agree (i.e. one finds a best galaxy, and the
                    # other finds NO best galaxy)
                    else:
                        # in this case vir_1_all and cus_1_all will not be the same, so I
                        # need to pick the right one
                        if virInclude:
                            entry = vir_1_all
                        else:
                            entry = cus_1_all
                            
                        entry.append(AGNredshift)
                        entry.append(spectrumStatus)
                        entry.append('?')
                        entry.append(virInclude)
                        entry.append(cusInclude)
                        entry.append(Lya_v)
                        entry.append(v_limits)
                        entry.append(Lya_W)
                        entry.append(Na)
                        entry.append(b)
                        entry.append(identified)
                        
                        # update the comment to say the includes do not match
                        comment1 = str(comment) + " : VIRIAL AND CUSTOM INCLUDE MISMATCH"
                        entry.append(comment1)
                        
                        finalEntry = [entry]              
                            
                        
                
                # now if there is only one galaxy
                elif sysLen == 1 and match:
        
                    # is the likelihood big enough? (>= l_min)
                    vir_1 = systemsV[0]
                    vir_1_like = vir_1[0]
                    if float(vir_1_like) >= l_min:
                        virInclude = True
                    else:
                        virInclude = False
                
                    cus_1 = systemsC[0]
                    cus_1_like = cus_1[0]
                    if float(cus_1_like) >= l_min:
                        cusInclude = True
                    else:
                        cusInclude = False
                    
                    entry = vir_1[1]
                    entry.append(AGNredshift)
                    entry.append(spectrumStatus)
                    entry.append('?')
                    entry.append(virInclude)
                    entry.append(cusInclude)
                    entry.append(Lya_v)
                    entry.append(v_limits)
                    entry.append(Lya_W)
                    entry.append(Na)
                    entry.append(b)
                    entry.append(identified)
                
                    # include the same updated comment as above
                    comment1 = str(comment) + " : ONLY 1 GALAXY"
                    entry.append(comment1)
                
                    finalEntry = [entry]

                # no galaxies, mark it so in the comments
                else:
                    virInclude = False
                    cusInclude = False
                    
                    # also finalInclude is forced to FALSE - see below

                    comment1 = str(comment) + " : NO GALAXIES"
                    
                    finalEntry = [[l[targetHeader],\
                    Lya_v,\
                    'x',\
                    0,\
                    degreesJ2000RA_DecAGN,\
                    ('x','x'),\
                    -99,\
                    -99,\
                    -99,\
                    -99,\
                    'x',\
                    'x',\
                    'x',\
                    'x',\
                    'x',\
                    'x',\
                    'x',\
                    'x',\
                    'x',\
                    'x',\
                    'x',\
                    'x',\
                    'x',\
                    'x',\
                    'x',\
                    'x',\
                    'x',\
                    'x',\
                    'x',\
                    AGNredshift,\
                    spectrumStatus,\
                    False,\
                    virInclude,\
                    virInclude,\
                    Lya_v,\
                    v_limits,\
                    Lya_W,\
                    Na,\
                    b,\
                    identified,\
                    comment1]]
                    
                                        
                # write it all out to file
                for entry in finalEntry:
#                     print 'entry: ',entry

                    row = dict((f,o) for f,o in zip(fieldnames,entry))
                    writer.writerow(row)
            
            
        # close the files
        f.close()
        writerOutFile.close()
            
            
##########################################################################################
##########################################################################################
       
    # now actually do it
    if saveResults:
        print
        print
        print
        print
        print 'masterVirList: ',masterVirList
        print
        print
        print
        print
        print
        writeFullResults(masterVirList,masterCustomList)
    
    print "Finished."
    

if __name__=="__main__":
    main()
