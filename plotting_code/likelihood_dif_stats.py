#!/usr/bin/env python

'''
By David French (frenchd@astro.wisc.edu)

$Id:  pilot_stats.py, v 1.4 09/23/16

Print out all the relevant stats on the dataset for the pilot paper (01/04/2016)

v1.1: updates for v_hel velocity and probably something else? (2/22/16)

v1.2: updates for the new large galaxy sample (07/14/16) -> /plots4/

v1.3: added ability to limit results by environment variable (7/14/16)
    also same for likelihood values, including some stats about them also
    
v1.4: updated to LG_correlation_combined5_11_25cut_edit4.csv and plots5 (9/23/16)

'''

import sys
import os
import csv

from pylab import *
from scipy import stats
# import atpy
from math import *
from utilities import *
import getpass
import pickle

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



###########################################################################

    
def main():
    # assuming 'theFile' contains one name per line, read the file
    
    if getpass.getuser() == 'frenchd':

#         pickleFilename = '/Users/frenchd/Research/inclination/git_inclination/picklePilot_plusSALT_14.p'
#         gtPickleFilename = '/Users/frenchd/Research/inclination/git_inclination/pickleGT.p'
#         saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/plotting_code/figs'

#         gtPickleFilename = '/Users/frenchd/Research/inclination/git_inclination/pickleGT_filteredAll.p'
        gtPickleFilename = '/Users/frenchd/Research/GT_update2/pickleGT_filteredAll.p'

        saveDirectory = '/Users/frenchd/Research/inclination/git_inclination/plotting_code/figs/'
        
        isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/isolated8.p'
        L_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_isolated8.p'
        L_associated_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated8.p'
        L_associated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated8.p'
        L_nonassociated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated8.p'
        L_two_filename = '/Users/frenchd/Research/inclination/git_inclination/L_two8.p'
        L_three_plus_filename = '/Users/frenchd/Research/inclination/git_inclination/L_three_plus8.p'
        L_group_filename = '/Users/frenchd/Research/inclination/git_inclination/L_group8.p'
        L_summed_filename = '/Users/frenchd/Research/inclination/git_inclination/L_summed8.p'
        all_filename = '/Users/frenchd/Research/inclination/git_inclination/all8.p'


        # double 
        isolated_filename_double = '/Users/frenchd/Research/inclination/git_inclination/isolated8_double.p'
        L_isolated_filename_double = '/Users/frenchd/Research/inclination/git_inclination/L_isolated8_double.p'
        L_associated_isolated_filename_double = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated8_double.p'
        L_associated_filename_double = '/Users/frenchd/Research/inclination/git_inclination/L_associated8_double.p'
        L_nonassociated_filename_double = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated8_double.p'
        L_two_filename_double = '/Users/frenchd/Research/inclination/git_inclination/L_two8_double.p'
        L_three_plus_filename_double = '/Users/frenchd/Research/inclination/git_inclination/L_three_plus8_double.p'
        L_group_filename_double = '/Users/frenchd/Research/inclination/git_inclination/L_group8_double.p'
        L_summed_filename_double = '/Users/frenchd/Research/inclination/git_inclination/L_summed8_double.p'
        all_filename_double = '/Users/frenchd/Research/inclination/git_inclination/all8_double.p'


        # min001
        data_set = '_min001'
        isolated_filename_min001 = '/Users/frenchd/Research/inclination/git_inclination/isolated8{0}.p'.format(data_set)
        L_isolated_filename_min001 = '/Users/frenchd/Research/inclination/git_inclination/L_isolated8{0}.p'.format(data_set)
        L_associated_isolated_filename_min001 = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated8{0}.p'.format(data_set)
        L_associated_filename_min001 = '/Users/frenchd/Research/inclination/git_inclination/L_associated8{0}.p'.format(data_set)
        L_nonassociated_filename_min001 = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated8{0}.p'.format(data_set)
        L_two_filename_min001 = '/Users/frenchd/Research/inclination/git_inclination/L_two8{0}.p'.format(data_set)
        L_three_plus_filename_min001 = '/Users/frenchd/Research/inclination/git_inclination/L_three_plus8{0}.p'.format(data_set)
        L_group_filename_min001 = '/Users/frenchd/Research/inclination/git_inclination/L_group8{0}.p'.format(data_set)
        L_summed_filename_min001 = '/Users/frenchd/Research/inclination/git_inclination/L_summed8{0}.p'.format(data_set)
        all_filename_min001 = '/Users/frenchd/Research/inclination/git_inclination/all8{0}.p'.format(data_set)

        # min001, rigor6
        data_set = '_min001_rigor6'
        isolated_filename_min001_rigor6 = '/Users/frenchd/Research/inclination/git_inclination/isolated8{0}.p'.format(data_set)
        L_isolated_filename_min001_rigor6 = '/Users/frenchd/Research/inclination/git_inclination/L_isolated8{0}.p'.format(data_set)
        L_associated_isolated_filename_min001_rigor6 = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated8{0}.p'.format(data_set)
        L_associated_filename_min001_rigor6 = '/Users/frenchd/Research/inclination/git_inclination/L_associated8{0}.p'.format(data_set)
        L_nonassociated_filename_min001_rigor6 = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated8{0}.p'.format(data_set)
        L_two_filename_min001_rigor6 = '/Users/frenchd/Research/inclination/git_inclination/L_two8{0}.p'.format(data_set)
        L_three_plus_filename_min001_rigor6 = '/Users/frenchd/Research/inclination/git_inclination/L_three_plus8{0}.p'.format(data_set)
        L_group_filename_min001_rigor6 = '/Users/frenchd/Research/inclination/git_inclination/L_group8{0}.p'.format(data_set)
        L_summed_filename_min001_rigor6 = '/Users/frenchd/Research/inclination/git_inclination/L_summed8{0}.p'.format(data_set)
        all_filename_min001_rigor6 = '/Users/frenchd/Research/inclination/git_inclination/all8{0}.p'.format(data_set)
        
        # min001, rigor7
        data_set = '_min001_rigor7'
        isolated_filename_min001_rigor7 = '/Users/frenchd/Research/inclination/git_inclination/isolated8{0}.p'.format(data_set)
        L_isolated_filename_min001_rigor7 = '/Users/frenchd/Research/inclination/git_inclination/L_isolated8{0}.p'.format(data_set)
        L_associated_isolated_filename_min001_rigor7 = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated8{0}.p'.format(data_set)
        L_associated_filename_min001_rigor7 = '/Users/frenchd/Research/inclination/git_inclination/L_associated8{0}.p'.format(data_set)
        L_nonassociated_filename_min001_rigor7 = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated8{0}.p'.format(data_set)
        L_two_filename_min001_rigor7 = '/Users/frenchd/Research/inclination/git_inclination/L_two8{0}.p'.format(data_set)
        L_three_plus_filename_min001_rigor7 = '/Users/frenchd/Research/inclination/git_inclination/L_three_plus8{0}.p'.format(data_set)
        L_group_filename_min001_rigor7 = '/Users/frenchd/Research/inclination/git_inclination/L_group8{0}.p'.format(data_set)
        L_summed_filename_min001_rigor7 = '/Users/frenchd/Research/inclination/git_inclination/L_summed8{0}.p'.format(data_set)
        all_filename_min001_rigor7 = '/Users/frenchd/Research/inclination/git_inclination/all8{0}.p'.format(data_set)

        # min001, rigor8
        data_set = '_min001_rigor8'
        isolated_filename_min001_rigor8 = '/Users/frenchd/Research/inclination/git_inclination/isolated8{0}.p'.format(data_set)
        L_isolated_filename_min001_rigor8 = '/Users/frenchd/Research/inclination/git_inclination/L_isolated8{0}.p'.format(data_set)
        L_associated_isolated_filename_min001_rigor8 = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated8{0}.p'.format(data_set)
        L_associated_filename_min001_rigor8 = '/Users/frenchd/Research/inclination/git_inclination/L_associated8{0}.p'.format(data_set)
        L_nonassociated_filename_min001_rigor8 = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated8{0}.p'.format(data_set)
        L_two_filename_min001_rigor8 = '/Users/frenchd/Research/inclination/git_inclination/L_two8{0}.p'.format(data_set)
        L_three_plus_filename_min001_rigor8 = '/Users/frenchd/Research/inclination/git_inclination/L_three_plus8{0}.p'.format(data_set)
        L_group_filename_min001_rigor8 = '/Users/frenchd/Research/inclination/git_inclination/L_group8{0}.p'.format(data_set)
        L_summed_filename_min001_rigor8 = '/Users/frenchd/Research/inclination/git_inclination/L_summed8{0}.p'.format(data_set)
        all_filename_min001_rigor8 = '/Users/frenchd/Research/inclination/git_inclination/all8{0}.p'.format(data_set)

        # min001, cus
        data_set = '_min001_cus'
        isolated_filename_min001_cus = '/Users/frenchd/Research/inclination/git_inclination/isolated8{0}.p'.format(data_set)
        L_isolated_filename_min001_cus = '/Users/frenchd/Research/inclination/git_inclination/L_isolated8{0}.p'.format(data_set)
        L_associated_isolated_filename_min001_cus = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated8{0}.p'.format(data_set)
        L_associated_filename_min001_cus = '/Users/frenchd/Research/inclination/git_inclination/L_associated8{0}.p'.format(data_set)
        L_nonassociated_filename_min001_cus = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated8{0}.p'.format(data_set)
        L_two_filename_min001_cus = '/Users/frenchd/Research/inclination/git_inclination/L_two8{0}.p'.format(data_set)
        L_three_plus_filename_min001_cus = '/Users/frenchd/Research/inclination/git_inclination/L_three_plus8{0}.p'.format(data_set)
        L_group_filename_min001_cus = '/Users/frenchd/Research/inclination/git_inclination/L_group8{0}.p'.format(data_set)
        L_summed_filename_min001_cus = '/Users/frenchd/Research/inclination/git_inclination/L_summed8{0}.p'.format(data_set)
        all_filename_min001_cus = '/Users/frenchd/Research/inclination/git_inclination/all8{0}.p'.format(data_set)

        # min001, double
        data_set = '_min001_double'
        isolated_filename_min001_double = '/Users/frenchd/Research/inclination/git_inclination/isolated8{0}.p'.format(data_set)
        L_isolated_filename_min001_double = '/Users/frenchd/Research/inclination/git_inclination/L_isolated8{0}.p'.format(data_set)
        L_associated_isolated_filename_min001_double = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated8{0}.p'.format(data_set)
        L_associated_filename_min001_double = '/Users/frenchd/Research/inclination/git_inclination/L_associated8{0}.p'.format(data_set)
        L_nonassociated_filename_min001_double = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated8{0}.p'.format(data_set)
        L_two_filename_min001_double = '/Users/frenchd/Research/inclination/git_inclination/L_two8{0}.p'.format(data_set)
        L_three_plus_filename_min001_double = '/Users/frenchd/Research/inclination/git_inclination/L_three_plus8{0}.p'.format(data_set)
        L_group_filename_min001_double = '/Users/frenchd/Research/inclination/git_inclination/L_group8{0}.p'.format(data_set)
        L_summed_filename_min001_double = '/Users/frenchd/Research/inclination/git_inclination/L_summed8{0}.p'.format(data_set)
        all_filename_min001_double = '/Users/frenchd/Research/inclination/git_inclination/all8{0}.p'.format(data_set)

       # cus
        data_set = '_cus'
        isolated_filename_cus = '/Users/frenchd/Research/inclination/git_inclination/isolated8{0}.p'.format(data_set)
        L_isolated_filename_cus = '/Users/frenchd/Research/inclination/git_inclination/L_isolated8{0}.p'.format(data_set)
        L_associated_isolated_filename_cus = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated8{0}.p'.format(data_set)
        L_associated_filename_cus = '/Users/frenchd/Research/inclination/git_inclination/L_associated8{0}.p'.format(data_set)
        L_nonassociated_filename_cus = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated8{0}.p'.format(data_set)
        L_two_filename_cus = '/Users/frenchd/Research/inclination/git_inclination/L_two8{0}.p'.format(data_set)
        L_three_plus_filename_cus = '/Users/frenchd/Research/inclination/git_inclination/L_three_plus8{0}.p'.format(data_set)
        L_group_filename_cus = '/Users/frenchd/Research/inclination/git_inclination/L_group8{0}.p'.format(data_set)
        L_summed_filename_cus = '/Users/frenchd/Research/inclination/git_inclination/L_summed8{0}.p'.format(data_set)
        all_filename_cus = '/Users/frenchd/Research/inclination/git_inclination/all8{0}.p'.format(data_set)

       # min005_v150
        data_set = '_min005_v150'
        isolated_filename_min005_v150 = '/Users/frenchd/Research/inclination/git_inclination/isolated8{0}.p'.format(data_set)
        L_isolated_filename_min005_v150 = '/Users/frenchd/Research/inclination/git_inclination/L_isolated8{0}.p'.format(data_set)
        L_associated_isolated_filename_min005_v150 = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated8{0}.p'.format(data_set)
        L_associated_filename_min005_v150 = '/Users/frenchd/Research/inclination/git_inclination/L_associated8{0}.p'.format(data_set)
        L_nonassociated_filename_min005_v150 = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated8{0}.p'.format(data_set)
        L_two_filename_min005_v150 = '/Users/frenchd/Research/inclination/git_inclination/L_two8{0}.p'.format(data_set)
        L_three_plus_filename_min005_v150 = '/Users/frenchd/Research/inclination/git_inclination/L_three_plus8{0}.p'.format(data_set)
        L_group_filename_min005_v150 = '/Users/frenchd/Research/inclination/git_inclination/L_group8{0}.p'.format(data_set)
        L_summed_filename_min005_v150 = '/Users/frenchd/Research/inclination/git_inclination/L_summed8{0}.p'.format(data_set)
        all_filename_min005_v150 = '/Users/frenchd/Research/inclination/git_inclination/all8{0}.p'.format(data_set)

       # min005_v250
        data_set = '_min005_v250'
        isolated_filename_min005_v250 = '/Users/frenchd/Research/inclination/git_inclination/isolated8{0}.p'.format(data_set)
        L_isolated_filename_min005_v250 = '/Users/frenchd/Research/inclination/git_inclination/L_isolated8{0}.p'.format(data_set)
        L_associated_isolated_filename_min005_v250 = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated8{0}.p'.format(data_set)
        L_associated_filename_min005_v250 = '/Users/frenchd/Research/inclination/git_inclination/L_associated8{0}.p'.format(data_set)
        L_nonassociated_filename_min005_v250 = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated8{0}.p'.format(data_set)
        L_two_filename_min005_v250 = '/Users/frenchd/Research/inclination/git_inclination/L_two8{0}.p'.format(data_set)
        L_three_plus_filename_min005_v250 = '/Users/frenchd/Research/inclination/git_inclination/L_three_plus8{0}.p'.format(data_set)
        L_group_filename_min005_v250 = '/Users/frenchd/Research/inclination/git_inclination/L_group8{0}.p'.format(data_set)
        L_summed_filename_min005_v250 = '/Users/frenchd/Research/inclination/git_inclination/L_summed8{0}.p'.format(data_set)
        all_filename_min005_v250 = '/Users/frenchd/Research/inclination/git_inclination/all8{0}.p'.format(data_set)

    else:
        print 'Could not determine username. Exiting.'
        sys.exit()
    
    # pickle file for the whole galaxy table:
    gtPickleFile = open(gtPickleFilename,'rU')
    gtDict = pickle.load(gtPickleFile)
    gtPickleFile.close()
    
    
    # open all the pickle files
    isolated_file = open(isolated_filename,'r')
    L_isolated_file = open(L_isolated_filename,'r')
    L_associated_isolated_file = open(L_associated_isolated_filename,'r')
    L_associated_file = open(L_associated_filename,'r')
    L_nonassociated_file = open(L_nonassociated_filename,'r')
    L_two_file = open(L_two_filename,'r')
    L_three_plus_file = open(L_three_plus_filename,'r')
    L_group_file = open(L_group_filename,'r')

    # unload the data from them
    isolated = pickle.load(isolated_file)
    L_isolated = pickle.load(L_isolated_file)
    L_associated_isolated = pickle.load(L_associated_isolated_file)
    L_associated = pickle.load(L_associated_file)
    L_nonassociated = pickle.load(L_nonassociated_file)
    L_two = pickle.load(L_two_file)
    L_three_plus = pickle.load(L_three_plus_file)
    L_group = pickle.load(L_group_file)
    
    # close the files
    isolated_file.close()
    L_isolated_file.close()
    L_associated_isolated_file.close()
    L_associated_file.close()
    L_nonassociated_file.close()
    L_two_file.close()
    L_three_plus_file.close()
    L_group_file.close()
    
##########################################################################################
##########################################################################################
    # double

    # open all the pickle files
    isolated_file_double = open(isolated_filename_double,'r')
    L_isolated_file_double = open(L_isolated_filename_double,'r')
    L_associated_isolated_file_double = open(L_associated_isolated_filename_double,'r')
    L_associated_file_double = open(L_associated_filename_double,'r')
    L_nonassociated_file_double = open(L_nonassociated_filename_double,'r')
    L_two_file_double = open(L_two_filename_double,'r')
    L_three_plus_file_double = open(L_three_plus_filename_double,'r')
    L_group_file_double = open(L_group_filename_double,'r')

    # unload the data from them
    isolated_double = pickle.load(isolated_file_double)
    L_isolated_double = pickle.load(L_isolated_file_double)
    L_associated_isolated_double = pickle.load(L_associated_isolated_file_double)
    L_associated_double = pickle.load(L_associated_file_double)
    L_nonassociated_double = pickle.load(L_nonassociated_file_double)
    L_two_double = pickle.load(L_two_file_double)
    L_three_plus_double = pickle.load(L_three_plus_file_double)
    L_group_double = pickle.load(L_group_file_double)
    
    # close the files
    isolated_file_double.close()
    L_isolated_file_double.close()
    L_associated_isolated_file_double.close()
    L_associated_file_double.close()
    L_nonassociated_file_double.close()
    L_two_file_double.close()
    L_three_plus_file_double.close()
    L_group_file_double.close()
    
##########################################################################################
##########################################################################################


    # _min001
    # open all the pickle files
    isolated_file_min001 = open(isolated_filename_min001,'r')
    L_isolated_file_min001 = open(L_isolated_filename_min001,'r')
    L_associated_isolated_file_min001 = open(L_associated_isolated_filename_min001,'r')
    L_associated_file_min001 = open(L_associated_filename_min001,'r')
    L_nonassociated_file_min001 = open(L_nonassociated_filename_min001,'r')
    L_two_file_min001 = open(L_two_filename_min001,'r')
    L_three_plus_file_min001 = open(L_three_plus_filename_min001,'r')
    L_group_file_min001 = open(L_group_filename_min001,'r')

    # unload the data from them
    isolated_min001 = pickle.load(isolated_file_min001)
    L_isolated_min001 = pickle.load(L_isolated_file_min001)
    L_associated_isolated_min001 = pickle.load(L_associated_isolated_file_min001)
    L_associated_min001 = pickle.load(L_associated_file_min001)
    L_nonassociated_min001 = pickle.load(L_nonassociated_file_min001)
    L_two_min001 = pickle.load(L_two_file_min001)
    L_three_plus_min001 = pickle.load(L_three_plus_file_min001)
    L_group_min001 = pickle.load(L_group_file_min001)
    
    # close the files
    isolated_file_min001.close()
    L_isolated_file_min001.close()
    L_associated_isolated_file_min001.close()
    L_associated_file_min001.close()
    L_nonassociated_file_min001.close()
    L_two_file_min001.close()
    L_three_plus_file_min001.close()
    L_group_file_min001.close()
    
##########################################################################################
##########################################################################################
    # _min001_rigor6
    # open all the pickle files
    isolated_file_min001_rigor6 = open(isolated_filename_min001_rigor6,'r')
    L_isolated_file_min001_rigor6 = open(L_isolated_filename_min001_rigor6,'r')
    L_associated_isolated_file_min001_rigor6 = open(L_associated_isolated_filename_min001_rigor6,'r')
    L_associated_file_min001_rigor6 = open(L_associated_filename_min001_rigor6,'r')
    L_nonassociated_file_min001_rigor6 = open(L_nonassociated_filename_min001_rigor6,'r')
    L_two_file_min001_rigor6 = open(L_two_filename_min001_rigor6,'r')
    L_three_plus_file_min001_rigor6 = open(L_three_plus_filename_min001_rigor6,'r')
    L_group_file_min001_rigor6 = open(L_group_filename_min001_rigor6,'r')

    # unload the data from them
    isolated_min001_rigor6 = pickle.load(isolated_file_min001_rigor6)
    L_isolated_min001_rigor6 = pickle.load(L_isolated_file_min001_rigor6)
    L_associated_isolated_min001_rigor6 = pickle.load(L_associated_isolated_file_min001_rigor6)
    L_associated_min001_rigor6 = pickle.load(L_associated_file_min001_rigor6)
    L_nonassociated_min001_rigor6 = pickle.load(L_nonassociated_file_min001_rigor6)
    L_two_min001_rigor6 = pickle.load(L_two_file_min001_rigor6)
    L_three_plus_min001_rigor6 = pickle.load(L_three_plus_file_min001_rigor6)
    L_group_min001_rigor6 = pickle.load(L_group_file_min001_rigor6)
    
    # close the files
    isolated_file_min001_rigor6.close()
    L_isolated_file_min001_rigor6.close()
    L_associated_isolated_file_min001_rigor6.close()
    L_associated_file_min001_rigor6.close()
    L_nonassociated_file_min001_rigor6.close()
    L_two_file_min001_rigor6.close()
    L_three_plus_file_min001_rigor6.close()
    L_group_file_min001_rigor6.close()

##########################################################################################
##########################################################################################
    # _min001_rigor7
    # open all the pickle files
    isolated_file_min001_rigor7 = open(isolated_filename_min001_rigor7,'r')
    L_isolated_file_min001_rigor7 = open(L_isolated_filename_min001_rigor7,'r')
    L_associated_isolated_file_min001_rigor7 = open(L_associated_isolated_filename_min001_rigor7,'r')
    L_associated_file_min001_rigor7 = open(L_associated_filename_min001_rigor7,'r')
    L_nonassociated_file_min001_rigor7 = open(L_nonassociated_filename_min001_rigor7,'r')
    L_two_file_min001_rigor7 = open(L_two_filename_min001_rigor7,'r')
    L_three_plus_file_min001_rigor7 = open(L_three_plus_filename_min001_rigor7,'r')
    L_group_file_min001_rigor7 = open(L_group_filename_min001_rigor7,'r')

    # unload the data from them
    isolated_min001_rigor7 = pickle.load(isolated_file_min001_rigor7)
    L_isolated_min001_rigor7 = pickle.load(L_isolated_file_min001_rigor7)
    L_associated_isolated_min001_rigor7 = pickle.load(L_associated_isolated_file_min001_rigor7)
    L_associated_min001_rigor7 = pickle.load(L_associated_file_min001_rigor7)
    L_nonassociated_min001_rigor7 = pickle.load(L_nonassociated_file_min001_rigor7)
    L_two_min001_rigor7 = pickle.load(L_two_file_min001_rigor7)
    L_three_plus_min001_rigor7 = pickle.load(L_three_plus_file_min001_rigor7)
    L_group_min001_rigor7 = pickle.load(L_group_file_min001_rigor7)
    
    # close the files
    isolated_file_min001_rigor7.close()
    L_isolated_file_min001_rigor7.close()
    L_associated_isolated_file_min001_rigor7.close()
    L_associated_file_min001_rigor7.close()
    L_nonassociated_file_min001_rigor7.close()
    L_two_file_min001_rigor7.close()
    L_three_plus_file_min001_rigor7.close()
    L_group_file_min001_rigor7.close()
    
##########################################################################################
##########################################################################################
    # _min001_rigor8
    # open all the pickle files
    isolated_file_min001_rigor8 = open(isolated_filename_min001_rigor8,'r')
    L_isolated_file_min001_rigor8 = open(L_isolated_filename_min001_rigor8,'r')
    L_associated_isolated_file_min001_rigor8 = open(L_associated_isolated_filename_min001_rigor8,'r')
    L_associated_file_min001_rigor8 = open(L_associated_filename_min001_rigor8,'r')
    L_nonassociated_file_min001_rigor8 = open(L_nonassociated_filename_min001_rigor8,'r')
    L_two_file_min001_rigor8 = open(L_two_filename_min001_rigor8,'r')
    L_three_plus_file_min001_rigor8 = open(L_three_plus_filename_min001_rigor8,'r')
    L_group_file_min001_rigor8 = open(L_group_filename_min001_rigor8,'r')

    # unload the data from them
    isolated_min001_rigor8 = pickle.load(isolated_file_min001_rigor8)
    L_isolated_min001_rigor8 = pickle.load(L_isolated_file_min001_rigor8)
    L_associated_isolated_min001_rigor8 = pickle.load(L_associated_isolated_file_min001_rigor8)
    L_associated_min001_rigor8 = pickle.load(L_associated_file_min001_rigor8)
    L_nonassociated_min001_rigor8 = pickle.load(L_nonassociated_file_min001_rigor8)
    L_two_min001_rigor8 = pickle.load(L_two_file_min001_rigor8)
    L_three_plus_min001_rigor8 = pickle.load(L_three_plus_file_min001_rigor8)
    L_group_min001_rigor8 = pickle.load(L_group_file_min001_rigor8)
    
    # close the files
    isolated_file_min001_rigor8.close()
    L_isolated_file_min001_rigor8.close()
    L_associated_isolated_file_min001_rigor8.close()
    L_associated_file_min001_rigor8.close()
    L_nonassociated_file_min001_rigor8.close()
    L_two_file_min001_rigor8.close()
    L_three_plus_file_min001_rigor8.close()
    L_group_file_min001_rigor8.close()

##########################################################################################
##########################################################################################
    # _min001_cus
    # open all the pickle files
    isolated_file_min001_cus = open(isolated_filename_min001_cus,'r')
    L_isolated_file_min001_cus = open(L_isolated_filename_min001_cus,'r')
    L_associated_isolated_file_min001_cus = open(L_associated_isolated_filename_min001_cus,'r')
    L_associated_file_min001_cus = open(L_associated_filename_min001_cus,'r')
    L_nonassociated_file_min001_cus = open(L_nonassociated_filename_min001_cus,'r')
    L_two_file_min001_cus = open(L_two_filename_min001_cus,'r')
    L_three_plus_file_min001_cus = open(L_three_plus_filename_min001_cus,'r')
    L_group_file_min001_cus = open(L_group_filename_min001_cus,'r')

    # unload the data from them
    isolated_min001_cus = pickle.load(isolated_file_min001_cus)
    L_isolated_min001_cus = pickle.load(L_isolated_file_min001_cus)
    L_associated_isolated_min001_cus = pickle.load(L_associated_isolated_file_min001_cus)
    L_associated_min001_cus = pickle.load(L_associated_file_min001_cus)
    L_nonassociated_min001_cus = pickle.load(L_nonassociated_file_min001_cus)
    L_two_min001_cus = pickle.load(L_two_file_min001_cus)
    L_three_plus_min001_cus = pickle.load(L_three_plus_file_min001_cus)
    L_group_min001_cus = pickle.load(L_group_file_min001_cus)
    
    # close the files
    isolated_file_min001_cus.close()
    L_isolated_file_min001_cus.close()
    L_associated_isolated_file_min001_cus.close()
    L_associated_file_min001_cus.close()
    L_nonassociated_file_min001_cus.close()
    L_two_file_min001_cus.close()
    L_three_plus_file_min001_cus.close()
    L_group_file_min001_cus.close()
    
##########################################################################################
##########################################################################################
    # _min001_double
    # open all the pickle files
    isolated_file_min001_double = open(isolated_filename_min001_double,'r')
    L_isolated_file_min001_double = open(L_isolated_filename_min001_double,'r')
    L_associated_isolated_file_min001_double = open(L_associated_isolated_filename_min001_double,'r')
    L_associated_file_min001_double = open(L_associated_filename_min001_double,'r')
    L_nonassociated_file_min001_double = open(L_nonassociated_filename_min001_double,'r')
    L_two_file_min001_double = open(L_two_filename_min001_double,'r')
    L_three_plus_file_min001_double = open(L_three_plus_filename_min001_double,'r')
    L_group_file_min001_double = open(L_group_filename_min001_double,'r')

    # unload the data from them
    isolated_min001_double = pickle.load(isolated_file_min001_double)
    L_isolated_min001_double = pickle.load(L_isolated_file_min001_double)
    L_associated_isolated_min001_double = pickle.load(L_associated_isolated_file_min001_double)
    L_associated_min001_double = pickle.load(L_associated_file_min001_double)
    L_nonassociated_min001_double = pickle.load(L_nonassociated_file_min001_double)
    L_two_min001_double = pickle.load(L_two_file_min001_double)
    L_three_plus_min001_double = pickle.load(L_three_plus_file_min001_double)
    L_group_min001_double = pickle.load(L_group_file_min001_double)
    
    # close the files
    isolated_file_min001_double.close()
    L_isolated_file_min001_double.close()
    L_associated_isolated_file_min001_double.close()
    L_associated_file_min001_double.close()
    L_nonassociated_file_min001_double.close()
    L_two_file_min001_double.close()
    L_three_plus_file_min001_double.close()
    L_group_file_min001_double.close()
    
##########################################################################################
##########################################################################################
    # cus
    # open all the pickle files
    isolated_file_cus = open(isolated_filename_cus,'r')
    L_isolated_file_cus = open(L_isolated_filename_cus,'r')
    L_associated_isolated_file_cus = open(L_associated_isolated_filename_cus,'r')
    L_associated_file_cus = open(L_associated_filename_cus,'r')
    L_nonassociated_file_cus = open(L_nonassociated_filename_cus,'r')
    L_two_file_cus = open(L_two_filename_cus,'r')
    L_three_plus_file_cus = open(L_three_plus_filename_cus,'r')
    L_group_file_cus = open(L_group_filename_cus,'r')

    # unload the data from them
    isolated_cus = pickle.load(isolated_file_cus)
    L_isolated_cus = pickle.load(L_isolated_file_cus)
    L_associated_isolated_cus = pickle.load(L_associated_isolated_file_cus)
    L_associated_cus = pickle.load(L_associated_file_cus)
    L_nonassociated_cus = pickle.load(L_nonassociated_file_cus)
    L_two_cus = pickle.load(L_two_file_cus)
    L_three_plus_cus = pickle.load(L_three_plus_file_cus)
    L_group_cus = pickle.load(L_group_file_cus)
    
    # close the files
    isolated_file_cus.close()
    L_isolated_file_cus.close()
    L_associated_isolated_file_cus.close()
    L_associated_file_cus.close()
    L_nonassociated_file_cus.close()
    L_two_file_cus.close()
    L_three_plus_file_cus.close()
    L_group_file_cus.close()

##########################################################################################
##########################################################################################
    # min005_v150
    # open all the pickle files
    isolated_file_min005_v150 = open(isolated_filename_min005_v150,'r')
    L_isolated_file_min005_v150 = open(L_isolated_filename_min005_v150,'r')
    L_associated_isolated_file_min005_v150 = open(L_associated_isolated_filename_min005_v150,'r')
    L_associated_file_min005_v150 = open(L_associated_filename_min005_v150,'r')
    L_nonassociated_file_min005_v150 = open(L_nonassociated_filename_min005_v150,'r')
    L_two_file_min005_v150 = open(L_two_filename_min005_v150,'r')
    L_three_plus_file_min005_v150 = open(L_three_plus_filename_min005_v150,'r')
    L_group_file_min005_v150 = open(L_group_filename_min005_v150,'r')

    # unload the data from them
    isolated_min005_v150 = pickle.load(isolated_file_min005_v150)
    L_isolated_min005_v150 = pickle.load(L_isolated_file_min005_v150)
    L_associated_isolated_min005_v150 = pickle.load(L_associated_isolated_file_min005_v150)
    L_associated_min005_v150 = pickle.load(L_associated_file_min005_v150)
    L_nonassociated_min005_v150 = pickle.load(L_nonassociated_file_min005_v150)
    L_two_min005_v150 = pickle.load(L_two_file_min005_v150)
    L_three_plus_min005_v150 = pickle.load(L_three_plus_file_min005_v150)
    L_group_min005_v150 = pickle.load(L_group_file_min005_v150)
    
    # close the files
    isolated_file_min005_v150.close()
    L_isolated_file_min005_v150.close()
    L_associated_isolated_file_min005_v150.close()
    L_associated_file_min005_v150.close()
    L_nonassociated_file_min005_v150.close()
    L_two_file_min005_v150.close()
    L_three_plus_file_min005_v150.close()
    L_group_file_min005_v150.close()
    
##########################################################################################
##########################################################################################
    # min005_v250
    # open all the pickle files
    isolated_file_min005_v250 = open(isolated_filename_min005_v250,'r')
    L_isolated_file_min005_v250 = open(L_isolated_filename_min005_v250,'r')
    L_associated_isolated_file_min005_v250 = open(L_associated_isolated_filename_min005_v250,'r')
    L_associated_file_min005_v250 = open(L_associated_filename_min005_v250,'r')
    L_nonassociated_file_min005_v250 = open(L_nonassociated_filename_min005_v250,'r')
    L_two_file_min005_v250 = open(L_two_filename_min005_v250,'r')
    L_three_plus_file_min005_v250 = open(L_three_plus_filename_min005_v250,'r')
    L_group_file_min005_v250 = open(L_group_filename_min005_v250,'r')

    # unload the data from them
    isolated_min005_v250 = pickle.load(isolated_file_min005_v250)
    L_isolated_min005_v250 = pickle.load(L_isolated_file_min005_v250)
    L_associated_isolated_min005_v250 = pickle.load(L_associated_isolated_file_min005_v250)
    L_associated_min005_v250 = pickle.load(L_associated_file_min005_v250)
    L_nonassociated_min005_v250 = pickle.load(L_nonassociated_file_min005_v250)
    L_two_min005_v250 = pickle.load(L_two_file_min005_v250)
    L_three_plus_min005_v250 = pickle.load(L_three_plus_file_min005_v250)
    L_group_min005_v250 = pickle.load(L_group_file_min005_v250)
    
    # close the files
    isolated_file_min005_v250.close()
    L_isolated_file_min005_v250.close()
    L_associated_isolated_file_min005_v250.close()
    L_associated_file_min005_v250.close()
    L_nonassociated_file_min005_v250.close()
    L_two_file_min005_v250.close()
    L_three_plus_file_min005_v250.close()
    L_group_file_min005_v250.close()
    
    
    
############################################################################################################
    # standard
    
    # associated-isolated
    Lya_Ws_associated_isolated = L_associated_isolated['Lya_Ws']
    Nas_associated_isolated = L_associated_isolated['Nas']
    bs_associated_isolated = L_associated_isolated['bs']
    
    # associated
    Lya_Ws_associated = L_associated['Lya_Ws']
    Nas_associated = L_associated['Nas']
    bs_associated = L_associated['bs']

    # isolated
    Lya_Ws_isolated = isolated['Lya_Ws']
    Nas_isolated = isolated['Nas']
    bs_isolated = isolated['bs']
    
    # L-isolated
    Lya_Ws_L_isolated = L_isolated['Lya_Ws']
    Nas_L_isolated = L_isolated['Nas']
    bs_L_isolated = L_isolated['bs']
    
    # two_plus
    Lya_Ws_two_plus = np.array(list(L_two['Lya_Ws']) + list(L_three_plus['Lya_Ws']))
    Nas_two_plus = np.array(list(L_two['Nas']) + list(L_three_plus['Lya_Ws']))
    bs_two_plus = np.array(list(L_two['bs']) + list(L_three_plus['bs']))
    
    ############################################################
    # double
    
    # associated-isolated
    Lya_Ws_associated_isolated_double = L_associated_isolated_double['Lya_Ws']
    Nas_associated_isolated_double = L_associated_isolated_double['Nas']
    bs_associated_isolated_double = L_associated_isolated_double['bs']
    incs_associated_isolated_double = L_associated_isolated_double['adjustedIncs']
    azimuths_associated_isolated_double = L_associated_isolated_double['azimuths']
    
    azimuths_associated_isolated_double_min300 = []
    incs_associated_isolated_double_min300 = []
    
    for az, inc, W in zip(azimuths_associated_isolated_double, incs_associated_isolated_double,Lya_Ws_associated_isolated_double):
        W = float(W)
        az = float(az)
        inc = float(inc)
        if W >= 300:
            if az >= 0:
                azimuths_associated_isolated_double_min300.append(az)
            
            if inc >= 0:
                incs_associated_isolated_double_min300.append(inc)


    # associated
    Lya_Ws_associated_double = L_associated_double['Lya_Ws']
    Nas_associated_double = L_associated_double['Nas']
    bs_associated_double = L_associated_double['bs']
    incs_associated_double = L_associated_double['adjustedIncs']
    azimuths_associated_double = L_associated_double['azimuths']

    azimuths_associated_double_min300 = []
    incs_associated_double_min300 = []
    for az, inc, W in zip(azimuths_associated_double, incs_associated_double, Lya_Ws_associated_double):
        W = float(W)
        az = float(az)
        inc = float(inc)
        if W >= 300:
            if az >= 0:
                azimuths_associated_double_min300.append(az)
            
            if inc >= 0:
                incs_associated_double_min300.append(inc)


    # isolated
    Lya_Ws_isolated_double = isolated_double['Lya_Ws']
    Nas_isolated_double = isolated_double['Nas']
    bs_isolated_double = isolated_double['bs']
    
    # L-isolated
    Lya_Ws_L_isolated_double = L_isolated_double['Lya_Ws']
    Nas_L_isolated_double = L_isolated_double['Nas']
    bs_L_isolated_double = L_isolated_double['bs']
    
    # two_plus
    Lya_Ws_two_plus_double = np.array(list(L_two_double['Lya_Ws']) + list(L_three_plus_double['Lya_Ws']))
    Nas_two_plus_double = np.array(list(L_two_double['Nas']) + list(L_three_plus_double['Lya_Ws']))
    bs_two_plus_double = np.array(list(L_two_double['bs']) + list(L_three_plus_double['bs']))
    incs_two_plus_double = np.array(list(L_two_double['adjustedIncs']) + list(L_three_plus_double['adjustedIncs']))
    azimuths_two_plus_double = np.array(list(L_two_double['azimuths']) + list(L_three_plus_double['azimuths']))

    ############################################################
    # min001
    
    # associated-isolated
    Lya_Ws_associated_isolated_min001 = L_associated_isolated_min001['Lya_Ws']
    Nas_associated_isolated_min001 = L_associated_isolated_min001['Nas']
    bs_associated_isolated_min001 = L_associated_isolated_min001['bs']
    
    # associated
    Lya_Ws_associated_min001 = L_associated_min001['Lya_Ws']
    Nas_associated_min001 = L_associated_min001['Nas']
    bs_associated_min001 = L_associated_min001['bs']

    # isolated
    Lya_Ws_isolated_min001 = isolated_min001['Lya_Ws']
    Nas_isolated_min001 = isolated_min001['Nas']
    bs_isolated_min001 = isolated_min001['bs']
    
    # L-isolated
    Lya_Ws_L_isolated_min001 = L_isolated_min001['Lya_Ws']
    Nas_L_isolated_min001 = L_isolated_min001['Nas']
    bs_L_isolated_min001 = L_isolated_min001['bs']
    
    # two-plus
    Lya_Ws_two_plus_min001 = np.array(list(L_two_min001['Lya_Ws']) + list(L_three_plus_min001['Lya_Ws']))
    Nas_two_plus_min001 = np.array(list(L_two_min001['Nas']) + list(L_three_plus_min001['Lya_Ws']))
    bs_two_plus_min001 = np.array(list(L_two_min001['bs']) + list(L_three_plus_min001['bs']))

    ############################################################
    # min001, rigor6
    
    # associated-isolated
    Lya_Ws_associated_isolated_min001_rigor6 = L_associated_isolated_min001_rigor6['Lya_Ws']
    Nas_associated_isolated_min001_rigor6 = L_associated_isolated_min001_rigor6['Nas']
    bs_associated_isolated_min001_rigor6 = L_associated_isolated_min001_rigor6['bs']
    
    # associated
    Lya_Ws_associated_min001_rigor6 = L_associated_min001_rigor6['Lya_Ws']
    Nas_associated_min001_rigor6 = L_associated_min001_rigor6['Nas']
    bs_associated_min001_rigor6 = L_associated_min001_rigor6['bs']

    # isolated
    Lya_Ws_isolated_min001_rigor6 = isolated_min001_rigor6['Lya_Ws']
    Nas_isolated_min001_rigor6 = isolated_min001_rigor6['Nas']
    bs_isolated_min001_rigor6 = isolated_min001_rigor6['bs']
    
    # L-isolated
    Lya_Ws_L_isolated_min001_rigor6 = L_isolated_min001_rigor6['Lya_Ws']
    Nas_L_isolated_min001_rigor6 = L_isolated_min001_rigor6['Nas']
    bs_L_isolated_min001_rigor6 = L_isolated_min001_rigor6['bs']
    
    # two-plus
    Lya_Ws_two_plus_min001_rigor6 = np.array(list(L_two_min001_rigor6['Lya_Ws']) + list(L_three_plus_min001_rigor6['Lya_Ws']))
    Nas_two_plus_min001_rigor6 = np.array(list(L_two_min001_rigor6['Nas']) + list(L_three_plus_min001_rigor6['Lya_Ws']))
    bs_two_plus_min001_rigor6 = np.array(list(L_two_min001_rigor6['bs']) + list(L_three_plus_min001_rigor6['bs']))

    ############################################################
    # min001, rigor7
    
    # associated-isolated
    Lya_Ws_associated_isolated_min001_rigor7 = L_associated_isolated_min001_rigor7['Lya_Ws']
    Nas_associated_isolated_min001_rigor7 = L_associated_isolated_min001_rigor7['Nas']
    bs_associated_isolated_min001_rigor7 = L_associated_isolated_min001_rigor7['bs']
    
    # associated
    Lya_Ws_associated_min001_rigor7 = L_associated_min001_rigor7['Lya_Ws']
    Nas_associated_min001_rigor7 = L_associated_min001_rigor7['Nas']
    bs_associated_min001_rigor7 = L_associated_min001_rigor7['bs']

    # isolated
    Lya_Ws_isolated_min001_rigor7 = isolated_min001_rigor7['Lya_Ws']
    Nas_isolated_min001_rigor7 = isolated_min001_rigor7['Nas']
    bs_isolated_min001_rigor7 = isolated_min001_rigor7['bs']
    
    # L-isolated
    Lya_Ws_L_isolated_min001_rigor7 = L_isolated_min001_rigor7['Lya_Ws']
    Nas_L_isolated_min001_rigor7 = L_isolated_min001_rigor7['Nas']
    bs_L_isolated_min001_rigor7 = L_isolated_min001_rigor7['bs']
    
    # two-plus
    Lya_Ws_two_plus_min001_rigor7 = np.array(list(L_two_min001_rigor7['Lya_Ws']) + list(L_three_plus_min001_rigor7['Lya_Ws']))
    Nas_two_plus_min001_rigor7 = np.array(list(L_two_min001_rigor7['Nas']) + list(L_three_plus_min001_rigor7['Lya_Ws']))
    bs_two_plus_min001_rigor7 = np.array(list(L_two_min001_rigor7['bs']) + list(L_three_plus_min001_rigor7['bs']))

    ############################################################
    # min001, rigor8
    
    # associated-isolated
    Lya_Ws_associated_isolated_min001_rigor8 = L_associated_isolated_min001_rigor8['Lya_Ws']
    Nas_associated_isolated_min001_rigor8 = L_associated_isolated_min001_rigor8['Nas']
    bs_associated_isolated_min001_rigor8 = L_associated_isolated_min001_rigor8['bs']
    
    # associated
    Lya_Ws_associated_min001_rigor8 = L_associated_min001_rigor8['Lya_Ws']
    Nas_associated_min001_rigor8 = L_associated_min001_rigor8['Nas']
    bs_associated_min001_rigor8 = L_associated_min001_rigor8['bs']

    # isolated
    Lya_Ws_isolated_min001_rigor8 = isolated_min001_rigor8['Lya_Ws']
    Nas_isolated_min001_rigor8 = isolated_min001_rigor8['Nas']
    bs_isolated_min001_rigor8 = isolated_min001_rigor8['bs']
    
    # L-isolated
    Lya_Ws_L_isolated_min001_rigor8 = L_isolated_min001_rigor8['Lya_Ws']
    Nas_L_isolated_min001_rigor8 = L_isolated_min001_rigor8['Nas']
    bs_L_isolated_min001_rigor8 = L_isolated_min001_rigor8['bs']
    
    # two-plus
    Lya_Ws_two_plus_min001_rigor8 = np.array(list(L_two_min001_rigor8['Lya_Ws']) + list(L_three_plus_min001_rigor8['Lya_Ws']))
    Nas_two_plus_min001_rigor8 = np.array(list(L_two_min001_rigor8['Nas']) + list(L_three_plus_min001_rigor8['Lya_Ws']))
    bs_two_plus_min001_rigor8 = np.array(list(L_two_min001_rigor8['bs']) + list(L_three_plus_min001_rigor8['bs']))

    ############################################################
    # min001, cus
    
    # associated-isolated
    Lya_Ws_associated_isolated_min001_cus = L_associated_isolated_min001_cus['Lya_Ws']
    Nas_associated_isolated_min001_cus = L_associated_isolated_min001_cus['Nas']
    bs_associated_isolated_min001_cus = L_associated_isolated_min001_cus['bs']
    
    # associated
    Lya_Ws_associated_min001_cus = L_associated_min001_cus['Lya_Ws']
    Nas_associated_min001_cus = L_associated_min001_cus['Nas']
    bs_associated_min001_cus = L_associated_min001_cus['bs']

    # isolated
    Lya_Ws_isolated_min001_cus = isolated_min001_cus['Lya_Ws']
    Nas_isolated_min001_cus = isolated_min001_cus['Nas']
    bs_isolated_min001_cus = isolated_min001_cus['bs']
    
    # L-isolated
    Lya_Ws_L_isolated_min001_cus = L_isolated_min001_cus['Lya_Ws']
    Nas_L_isolated_min001_cus = L_isolated_min001_cus['Nas']
    bs_L_isolated_min001_cus = L_isolated_min001_cus['bs']
    
    # two-plus
    Lya_Ws_two_plus_min001_cus = np.array(list(L_two_min001_cus['Lya_Ws']) + list(L_three_plus_min001_cus['Lya_Ws']))
    Nas_two_plus_min001_cus = np.array(list(L_two_min001_cus['Nas']) + list(L_three_plus_min001_cus['Lya_Ws']))
    bs_two_plus_min001_cus = np.array(list(L_two_min001_cus['bs']) + list(L_three_plus_min001_cus['bs']))

    ############################################################
    # min001, double
    
    # associated-isolated
    Lya_Ws_associated_isolated_min001_double = L_associated_isolated_min001_double['Lya_Ws']
    Nas_associated_isolated_min001_double = L_associated_isolated_min001_double['Nas']
    bs_associated_isolated_min001_double = L_associated_isolated_min001_double['bs']
    
    # associated
    Lya_Ws_associated_min001_double = L_associated_min001_double['Lya_Ws']
    Nas_associated_min001_double = L_associated_min001_double['Nas']
    bs_associated_min001_double = L_associated_min001_double['bs']

    # isolated
    Lya_Ws_isolated_min001_double = isolated_min001_double['Lya_Ws']
    Nas_isolated_min001_double = isolated_min001_double['Nas']
    bs_isolated_min001_double = isolated_min001_double['bs']
    
    # L-isolated
    Lya_Ws_L_isolated_min001_double = L_isolated_min001_double['Lya_Ws']
    Nas_L_isolated_min001_double = L_isolated_min001_double['Nas']
    bs_L_isolated_min001_double = L_isolated_min001_double['bs']
    
    # two-plus
    Lya_Ws_two_plus_min001_double = np.array(list(L_two_min001_double['Lya_Ws']) + list(L_three_plus_min001_double['Lya_Ws']))
    Nas_two_plus_min001_double = np.array(list(L_two_min001_double['Nas']) + list(L_three_plus_min001_double['Lya_Ws']))
    bs_two_plus_min001_double = np.array(list(L_two_min001_double['bs']) + list(L_three_plus_min001_double['bs']))

    ############################################################
    # cus
    
    # associated-isolated
    Lya_Ws_associated_isolated_cus = L_associated_isolated_cus['Lya_Ws']
    Nas_associated_isolated_cus = L_associated_isolated_cus['Nas']
    bs_associated_isolated_cus = L_associated_isolated_cus['bs']
    
    # associated
    Lya_Ws_associated_cus = L_associated_cus['Lya_Ws']
    Nas_associated_cus = L_associated_cus['Nas']
    bs_associated_cus = L_associated_cus['bs']

    # isolated
    Lya_Ws_isolated_cus = isolated_cus['Lya_Ws']
    Nas_isolated_cus = isolated_cus['Nas']
    bs_isolated_cus = isolated_cus['bs']
    
    # L-isolated
    Lya_Ws_L_isolated_cus = L_isolated_cus['Lya_Ws']
    Nas_L_isolated_cus = L_isolated_cus['Nas']
    bs_L_isolated_cus = L_isolated_cus['bs']
    
    # two-plus
    Lya_Ws_two_plus_cus = np.array(list(L_two_cus['Lya_Ws']) + list(L_three_plus_cus['Lya_Ws']))
    Nas_two_plus_cus = np.array(list(L_two_cus['Nas']) + list(L_three_plus_cus['Lya_Ws']))
    bs_two_plus_cus = np.array(list(L_two_cus['bs']) + list(L_three_plus_cus['bs']))

    ############################################################
    # min005, v150
    
    # associated-isolated
    Lya_Ws_associated_isolated_min005_v150 = L_associated_isolated_min005_v150['Lya_Ws']
    Nas_associated_isolated_min005_v150 = L_associated_isolated_min005_v150['Nas']
    bs_associated_isolated_min005_v150 = L_associated_isolated_min005_v150['bs']
    
    # associated
    Lya_Ws_associated_min005_v150 = L_associated_min005_v150['Lya_Ws']
    Nas_associated_min005_v150 = L_associated_min005_v150['Nas']
    bs_associated_min005_v150 = L_associated_min005_v150['bs']

    # isolated
    Lya_Ws_isolated_min005_v150 = isolated_min005_v150['Lya_Ws']
    Nas_isolated_min005_v150 = isolated_min005_v150['Nas']
    bs_isolated_min005_v150 = isolated_min005_v150['bs']
    
    # L-isolated
    Lya_Ws_L_isolated_min005_v150 = L_isolated_min005_v150['Lya_Ws']
    Nas_L_isolated_min005_v150 = L_isolated_min005_v150['Nas']
    bs_L_isolated_min005_v150 = L_isolated_min005_v150['bs']
    
    # two-plus
    Lya_Ws_two_plus_min005_v150 = np.array(list(L_two_min005_v150['Lya_Ws']) + list(L_three_plus_min005_v150['Lya_Ws']))
    Nas_two_plus_min005_v150 = np.array(list(L_two_min005_v150['Nas']) + list(L_three_plus_min005_v150['Lya_Ws']))
    bs_two_plus_min005_v150 = np.array(list(L_two_min005_v150['bs']) + list(L_three_plus_min005_v150['bs']))

    ############################################################
    # min005, v250
    
    # associated-isolated
    Lya_Ws_associated_isolated_min005_v250 = L_associated_isolated_min005_v250['Lya_Ws']
    Nas_associated_isolated_min005_v250 = L_associated_isolated_min005_v250['Nas']
    bs_associated_isolated_min005_v250 = L_associated_isolated_min005_v250['bs']
    
    # associated
    Lya_Ws_associated_min005_v250 = L_associated_min005_v250['Lya_Ws']
    Nas_associated_min005_v250 = L_associated_min005_v250['Nas']
    bs_associated_min005_v250 = L_associated_min005_v250['bs']

    # isolated
    Lya_Ws_isolated_min005_v250 = isolated_min005_v250['Lya_Ws']
    Nas_isolated_min005_v250 = isolated_min005_v250['Nas']
    bs_isolated_min005_v250 = isolated_min005_v250['bs']
    
    # L-isolated
    Lya_Ws_L_isolated_min005_v250 = L_isolated_min005_v250['Lya_Ws']
    Nas_L_isolated_min005_v250 = L_isolated_min005_v250['Nas']
    bs_L_isolated_min005_v250 = L_isolated_min005_v250['bs']
    
    # two-plus
    Lya_Ws_two_plus_min005_v250 = np.array(list(L_two_min005_v250['Lya_Ws']) + list(L_three_plus_min005_v250['Lya_Ws']))
    Nas_two_plus_min005_v250 = np.array(list(L_two_min005_v250['Nas']) + list(L_three_plus_min005_v250['Lya_Ws']))
    bs_two_plus_min005_v250 = np.array(list(L_two_min005_v250['bs']) + list(L_three_plus_min005_v250['bs']))


    majorAxisL = gtDict['majorAxis']
    incL = gtDict['inc']
    adjustedIncL = gtDict['adjustedInc']
    paL = gtDict['PA']
    BmagL = gtDict['Bmag']
    RID_medianL = gtDict['RID_median']
    RID_meanL = gtDict['RID_mean']
    RID_stdL = gtDict['RID_std']
    VhelL = gtDict['Vhel']
    RAdegL = gtDict['RAdeg']
    DEdegL = gtDict['DEdeg']
    NameL= gtDict['Name']
    
    allPA = paL
    allInclinations = []
    allAdjustedIncs = []
    allCosInclinations = []

#     print 'type: ',type(incL)
    for i in incL:
        if i != -99:
            i = float(i)
            allInclinations.append(i)
            
            i2 = pi/180. * i
            cosi2 = cos(i)
            allCosInclinations.append(cosi2)
            
    allCosFancyCosInclinations = []
    for i in adjustedIncL:
        if str(i) != '-99':
            i = float(i)

            allAdjustedIncs.append(i)
            
            i2 = pi/180. * i
            cosi2 = cos(i)
            allCosFancyCosInclinations.append(cosi2)
            
    allDiameter = majorAxisL

    print 'finished with this shit'
    print 'len(allAdjustedIncs): ',len(allAdjustedIncs)
    print


########################################################################################
#########################################################################################
    
            
##########################################################################################
##########################################################################################
    
    # Standard stats
    print '---------------------------- Standard ----------------------------'
    print
    print 'stats.stats(Lya_Ws_isolated): ',stats.describe(Lya_Ws_isolated)
    print 'stats.stats(Lya_Ws_L_isolated): ',stats.describe(Lya_Ws_L_isolated)
    print 'stats.stats(Lya_Ws_L_associated_isolated): ',stats.describe(Lya_Ws_associated_isolated)
    print 'stats.stats(Lya_Ws_L_associated): ',stats.describe(Lya_Ws_associated)
    print 'stats.stats(Lya_Ws_L_two_plus): ',stats.describe(Lya_Ws_two_plus)
    print
    print
    # associated vs associated_isolated
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated, Lya_Ws_associated)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated, Lya_Ws_associated])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated, Lya_Ws_associated)

    print 'KS Lya_Ws: L_associated_isolated vs L_associated: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_associated: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_associated: ', p_val
    print
    print
    
    # associated_isolated vs two_plus
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated, Lya_Ws_two_plus)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated, Lya_Ws_two_plus])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated, Lya_Ws_two_plus)

    print 'KS Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_two_plus: ', p_val
    print
    
    # b-parameters
    ans1 = stats.ks_2samp(bs_associated_isolated, bs_two_plus)
    ans1a = stats.anderson_ksamp([bs_associated_isolated, bs_two_plus])
    z_stat, p_val = stats.ranksums(bs_associated_isolated, bs_two_plus)
    print 'KS bs: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD bs: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum bs: L_associated_isolated vs L_two_plus: ', p_val
    print

###############################################################################

    
    # double stats
    print '---------------------------- Double ----------------------------'
    print
    print 'len(Lya_Ws_isolated_double): ',len(Lya_Ws_isolated_double)
    print 'len(Lya_Ws_L_isolated_double): ',len(Lya_Ws_L_isolated_double)
    print 'len(Lya_Ws_associated_isolated_double): ',len(Lya_Ws_associated_isolated_double)
    print 'len(Lya_Ws_associated_double): ',len(Lya_Ws_associated_double)
    print 'len(Lya_Ws_two_plus_double): ',len(Lya_Ws_two_plus_double)
    print
    
    
    
    print
    print 'stats.stats(Lya_Ws_isolated): ',stats.describe(Lya_Ws_isolated_double)
    print 'stats.stats(Lya_Ws_L_isolated): ',stats.describe(Lya_Ws_L_isolated_double)
    print 'stats.stats(Lya_Ws_L_associated_isolated): ',stats.describe(Lya_Ws_associated_isolated_double)
    print 'stats.stats(Lya_Ws_L_associated): ',stats.describe(Lya_Ws_associated_double)
    print 'stats.stats(Lya_Ws_L_two_plus): ',stats.describe(Lya_Ws_two_plus_double)
    print
    print
    # associated vs associated_isolated
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_double, Lya_Ws_associated_double)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_double, Lya_Ws_associated_double])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_double, Lya_Ws_associated_double)

    print 'KS Lya_Ws: L_associated_isolated vs L_associated: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_associated: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_associated: ', p_val
    print
    print
    
    # associated_isolated vs two_plus
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_double, Lya_Ws_two_plus_double)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_double, Lya_Ws_two_plus_double])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_double, Lya_Ws_two_plus_double)

    print 'KS Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_two_plus: ', p_val
    print
    
    # associated_isolated vs L_isolated
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_double, Lya_Ws_isolated_double)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_double, Lya_Ws_isolated_double])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_double, Lya_Ws_isolated_double)

    print 'KS Lya_Ws: L_associated_isolated vs Lya_Ws_isolated_double: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs Lya_Ws_isolated_double: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs Lya_Ws_isolated_double: ', p_val
    print
    print

    # associated vs L_isolated
    ans1 = stats.ks_2samp(Lya_Ws_associated_double, Lya_Ws_isolated_double)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_double, Lya_Ws_isolated_double])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_double, Lya_Ws_isolated_double)

    print 'KS Lya_Ws: Lya_Ws_associated_double vs Lya_Ws_isolated_double: ',ans1
    print 'AD Lya_Ws: Lya_Ws_associated_double vs Lya_Ws_isolated_double: ',ans1a
    print 'Ranksum Lya_Ws: Lya_Ws_associated_double vs Lya_Ws_isolated_double: ', p_val
    print
    print
    
    # associated vs L_isolated
    ans1 = stats.ks_2samp(Lya_Ws_two_plus_double, Lya_Ws_isolated_double)
    ans1a = stats.anderson_ksamp([Lya_Ws_two_plus_double, Lya_Ws_isolated_double])
    z_stat, p_val = stats.ranksums(Lya_Ws_two_plus_double, Lya_Ws_isolated_double)

    print 'KS Lya_Ws: Lya_Ws_two_plus_double vs Lya_Ws_isolated_double: ',ans1
    print 'AD Lya_Ws: Lya_Ws_two_plus_double vs Lya_Ws_isolated_double: ',ans1a
    print 'Ranksum Lya_Ws: Lya_Ws_two_plus_double vs Lya_Ws_isolated_double: ', p_val
    print
    print
    
    
    
    # b-parameters
    ans1 = stats.ks_2samp(bs_associated_isolated_double, bs_two_plus_double)
    ans1a = stats.anderson_ksamp([bs_associated_isolated_double, bs_two_plus_double])
    z_stat, p_val = stats.ranksums(bs_associated_isolated_double, bs_two_plus_double)
    print 'KS bs: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD bs: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum bs: L_associated_isolated vs L_two_plus: ', p_val
    print
    
    # inclinations
    ans1 = stats.ks_2samp(incs_associated_isolated_double, allAdjustedIncs)
    ans1a = stats.anderson_ksamp([incs_associated_isolated_double, allAdjustedIncs])
    z_stat, p_val = stats.ranksums(incs_associated_isolated_double, allAdjustedIncs)
    print 'KS adjustedIncs: L_associated_isolated vs all: ',ans1
    print 'AD adjustedIncs: L_associated_isolated vs all: ',ans1a
    print 'Ranksum adjustedIncs: L_associated_isolated vs all: ', p_val
    print
    
    ans1 = stats.ks_2samp(incs_associated_double, allAdjustedIncs)
    ans1a = stats.anderson_ksamp([incs_associated_double, allAdjustedIncs])
    z_stat, p_val = stats.ranksums(incs_associated_double, allAdjustedIncs)
    print 'KS adjustedIncs: L_associated vs all: ',ans1
    print 'AD adjustedIncs: L_associated vs all: ',ans1a
    print 'Ranksum adjustedIncs: L_associated vs all: ', p_val
    print
    
    ans1 = stats.ks_2samp(incs_two_plus_double, allAdjustedIncs)
    ans1a = stats.anderson_ksamp([incs_two_plus_double, allAdjustedIncs])
    z_stat, p_val = stats.ranksums(incs_two_plus_double, allAdjustedIncs)
    print 'KS adjustedIncs: L_two+ vs all: ',ans1
    print 'AD adjustedIncs: L_two+ vs all: ',ans1a
    print 'Ranksum adjustedIncs: L_two+ vs all: ', p_val
    print

    all_assoc_incs = np.array(list(incs_associated_isolated_double) + list(incs_associated_double))
    all_assoc_incs_min300 = np.array(list(incs_associated_isolated_double_min300) + list(incs_associated_double_min300))
    
    ans1 = stats.ks_2samp(all_assoc_incs, allAdjustedIncs)
    ans1a = stats.anderson_ksamp([all_assoc_incs, allAdjustedIncs])
    z_stat, p_val = stats.ranksums(all_assoc_incs, allAdjustedIncs)
    print 'KS adjustedIncs: all_assoc_incs vs all: ',ans1
    print 'AD adjustedIncs: all_assoc_incs vs all: ',ans1a
    print 'Ranksum adjustedIncs: all_assoc_incs vs all: ', p_val
    print

    ans1 = stats.ks_2samp(all_assoc_incs, incs_two_plus_double)
    ans1a = stats.anderson_ksamp([all_assoc_incs, incs_two_plus_double])
    z_stat, p_val = stats.ranksums(all_assoc_incs, incs_two_plus_double)
    print 'KS adjustedIncs: all_assoc_incs vs incs_two_plus_double: ',ans1
    print 'AD adjustedIncs: all_assoc_incs vs incs_two_plus_double: ',ans1a
    print 'Ranksum adjustedIncs: all_assoc_incs vs incs_two_plus_double: ', p_val
    print
    ans1 = stats.ks_2samp(all_assoc_incs_min300, allAdjustedIncs)
    ans1a = stats.anderson_ksamp([all_assoc_incs_min300, allAdjustedIncs])
    z_stat, p_val = stats.ranksums(all_assoc_incs_min300, allAdjustedIncs)
    print 'KS adjustedIncs: all_aall_assoc_incs_min300ssoc_incs vs all: ',ans1
    print 'AD adjustedIncs: all_assoc_incs_min300 vs all: ',ans1a
    print 'Ranksum adjustedIncs: all_assoc_incs_min300 vs all: ', p_val
    print

    ans1 = stats.ks_2samp(all_assoc_incs_min300, incs_two_plus_double)
    ans1a = stats.anderson_ksamp([all_assoc_incs_min300, incs_two_plus_double])
    z_stat, p_val = stats.ranksums(all_assoc_incs_min300, incs_two_plus_double)
    print 'KS adjustedIncs: all_assoc_incs_min300 vs incs_two_plus_double: ',ans1
    print 'AD adjustedIncs: all_assoc_incs_min300 vs incs_two_plus_double: ',ans1a
    print 'Ranksum adjustedIncs: all_assoc_incs_min300 vs incs_two_plus_double: ', p_val
    print


    # azimuths
    number_az = 216
    
    flat_az = np.ones(9)*24

    all_assoc_azimuths = np.array(list(azimuths_associated_isolated_double) + list(azimuths_associated_double))
    all_assoc_azimuths_min300 = np.array(list(azimuths_associated_isolated_double_min300) + list(azimuths_associated_double_min300))

    ans1 = stats.ks_2samp(all_assoc_azimuths, azimuths_two_plus_double)
    ans1a = stats.anderson_ksamp([all_assoc_azimuths, azimuths_two_plus_double])
    z_stat, p_val = stats.ranksums(all_assoc_azimuths, azimuths_two_plus_double)
    print 'KS Azimuths: all_assoc_azimuths vs azimuths_two_plus_double: ',ans1
    print 'AD Azimuths: all_assoc_azimuths vs azimuths_two_plus_double: ',ans1a
    print 'Ranksum Azimuths: all_assoc_azimuths vs azimuths_two_plus_double: ', p_val
    print
    
    ans1 = stats.ks_2samp(all_assoc_azimuths, flat_az)
    ans1a = stats.anderson_ksamp([all_assoc_azimuths, flat_az])
    z_stat, p_val = stats.ranksums(all_assoc_azimuths, flat_az)
    print 'KS Azimuths: all_assoc_azimuths vs flat_az: ',ans1
    print 'AD Azimuths: all_assoc_azimuths vs flat_az: ',ans1a
    print 'Ranksum Azimuths: all_assoc_azimuths vs flat_az: ', p_val
    print
    
    ans1 = stats.ks_2samp(azimuths_two_plus_double, flat_az)
    ans1a = stats.anderson_ksamp([azimuths_two_plus_double, flat_az])
    z_stat, p_val = stats.ranksums(azimuths_two_plus_double, flat_az)
    print 'KS Azimuths: azimuths_two_plus_double vs flat_az: ',ans1
    print 'AD Azimuths: azimuths_two_plus_double vs flat_az: ',ans1a
    print 'Ranksum Azimuths: azimuths_two_plus_double vs flat_az: ', p_val
    print
    
    ans1 = stats.ks_2samp(all_assoc_azimuths_min300, azimuths_two_plus_double)
    ans1a = stats.anderson_ksamp([all_assoc_azimuths_min300, azimuths_two_plus_double])
    z_stat, p_val = stats.ranksums(all_assoc_azimuths_min300, azimuths_two_plus_double)
    print 'KS Azimuths: all_assoc_azimuths_min300 vs azimuths_two_plus_double: ',ans1
    print 'AD Azimuths: all_assoc_azimuths_min300 vs azimuths_two_plus_double: ',ans1a
    print 'Ranksum Azimuths: all_assoc_azimuths_min300 vs azimuths_two_plus_double: ', p_val
    print
    
    ans1 = stats.ks_2samp(all_assoc_azimuths_min300, flat_az)
    ans1a = stats.anderson_ksamp([all_assoc_azimuths_min300, flat_az])
    z_stat, p_val = stats.ranksums(all_assoc_azimuths_min300, flat_az)
    print 'KS Azimuths: all_assoc_azimuths_min300 vs flat_az: ',ans1
    print 'AD Azimuths: all_assoc_azimuths_min300 vs flat_az: ',ans1a
    print 'Ranksum Azimuths: all_assoc_azimuths_min300 vs flat_az: ', p_val
    print
    
    ans1 = stats.ks_2samp(all_assoc_azimuths_min300, all_assoc_azimuths)
    ans1a = stats.anderson_ksamp([all_assoc_azimuths_min300, all_assoc_azimuths])
    z_stat, p_val = stats.ranksums(azimuths_two_plus_double, all_assoc_azimuths)
    print 'KS Azimuths: all_assoc_azimuths_min300 vs all_assoc_azimuths: ',ans1
    print 'AD Azimuths: all_assoc_azimuths_min300 vs all_assoc_azimuths: ',ans1a
    print 'Ranksum Azimuths: all_assoc_azimuths_min300 vs all_assoc_azimuths: ', p_val
    print
    

###############################################################################

    # min001 stats
    print '---------------------------- min001 ----------------------------'
    print
    print 'stats.describe(Lya_Ws_isolated): ',stats.describe(Lya_Ws_isolated_min001)
    print 'stats.describe(Lya_Ws_L_isolated): ',stats.describe(Lya_Ws_L_isolated_min001)
    print 'stats.describe(Lya_Ws_L_associated_isolated): ',stats.describe(Lya_Ws_associated_isolated_min001)
    print 'stats.describe(Lya_Ws_L_associated): ',stats.describe(Lya_Ws_associated_min001)
    print 'stats.describe(Lya_Ws_L_two_plus): ',stats.describe(Lya_Ws_two_plus_min001)
    print
    print
    # associated vs associated_isolated
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_min001, Lya_Ws_associated_min001)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_min001, Lya_Ws_associated_min001])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_min001, Lya_Ws_associated_min001)

    print 'KS Lya_Ws: L_associated_isolated vs L_associated: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_associated: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_associated: ', p_val
    print
    print
    
    # associated_isolated vs two_plus
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_min001, Lya_Ws_two_plus_min001)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_min001, Lya_Ws_two_plus_min001])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_min001, Lya_Ws_two_plus_min001)

    print 'KS Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_two_plus: ', p_val
    print

    # b-parameters
    ans1 = stats.ks_2samp(bs_associated_isolated_min001, bs_two_plus_min001)
    ans1a = stats.anderson_ksamp([bs_associated_isolated_min001, bs_two_plus_min001])
    z_stat, p_val = stats.ranksums(bs_associated_isolated_min001, bs_two_plus_min001)

    print 'KS bs: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD bs: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum bs: L_associated_isolated vs L_two_plus: ', p_val
    print

###############################################################################

    # min001, rigor6 stats
    print '---------------------------- min001, rigor6 ----------------------------'
    print
    print 'stats.describe(Lya_Ws_isolated): ',stats.describe(Lya_Ws_isolated_min001_rigor6)
    print 'stats.describe(Lya_Ws_L_isolated): ',stats.describe(Lya_Ws_L_isolated_min001_rigor6)
    print 'stats.describe(Lya_Ws_L_associated_isolated): ',stats.describe(Lya_Ws_associated_isolated_min001_rigor6)
    print 'stats.describe(Lya_Ws_L_associated): ',stats.describe(Lya_Ws_associated_min001_rigor6)
    print 'stats.describe(Lya_Ws_L_two_plus): ',stats.describe(Lya_Ws_two_plus_min001_rigor6)
    print
    print
    # associated vs associated_isolated
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_min001_rigor6, Lya_Ws_associated_min001_rigor6)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_min001_rigor6, Lya_Ws_associated_min001_rigor6])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_min001_rigor6, Lya_Ws_associated_min001_rigor6)

    print 'KS Lya_Ws: L_associated_isolated vs L_associated: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_associated: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_associated: ', p_val
    print
    print
    
    # associated_isolated vs two_plus
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_min001_rigor6, Lya_Ws_two_plus_min001_rigor6)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_min001_rigor6, Lya_Ws_two_plus_min001_rigor6])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_min001_rigor6, Lya_Ws_two_plus_min001_rigor6)

    print 'KS Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_two_plus: ', p_val
    print
    
    # b-parameters
    ans1 = stats.ks_2samp(bs_associated_isolated_min001_rigor6, bs_two_plus_min001_rigor6)
    ans1a = stats.anderson_ksamp([bs_associated_isolated_min001_rigor6, bs_two_plus_min001_rigor6])
    z_stat, p_val = stats.ranksums(bs_associated_isolated_min001_rigor6, bs_two_plus_min001_rigor6)

    print 'KS bs: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD bs: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum bs: L_associated_isolated vs L_two_plus: ', p_val
    print

###############################################################################

    # min001, rigor7 stats
    print '---------------------------- min001, rigor7 ----------------------------'
    print
    print 'stats.describe(Lya_Ws_isolated): ',stats.describe(Lya_Ws_isolated_min001_rigor7)
    print 'stats.describe(Lya_Ws_L_isolated): ',stats.describe(Lya_Ws_L_isolated_min001_rigor7)
    print 'stats.describe(Lya_Ws_L_associated_isolated): ',stats.describe(Lya_Ws_associated_isolated_min001_rigor7)
    print 'stats.describe(Lya_Ws_L_associated): ',stats.describe(Lya_Ws_associated_min001_rigor7)
    print 'stats.describe(Lya_Ws_L_two_plus): ',stats.describe(Lya_Ws_two_plus_min001_rigor7)
    print
    print
    # associated vs associated_isolated
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_min001_rigor7, Lya_Ws_associated_min001_rigor7)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_min001_rigor7, Lya_Ws_associated_min001_rigor7])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_min001_rigor7, Lya_Ws_associated_min001_rigor7)

    print 'KS Lya_Ws: L_associated_isolated vs L_associated: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_associated: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_associated: ', p_val
    print
    print
    
    # associated_isolated vs two_plus
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_min001_rigor7, Lya_Ws_two_plus_min001_rigor7)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_min001_rigor7, Lya_Ws_two_plus_min001_rigor7])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_min001_rigor7, Lya_Ws_two_plus_min001_rigor7)

    print 'KS Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_two_plus: ', p_val
    print

    # b-parameters
    ans1 = stats.ks_2samp(bs_associated_isolated_min001_rigor7, bs_two_plus_min001_rigor7)
    ans1a = stats.anderson_ksamp([bs_associated_isolated_min001_rigor7, bs_two_plus_min001_rigor7])
    z_stat, p_val = stats.ranksums(bs_associated_isolated_min001_rigor7, bs_two_plus_min001_rigor7)

    print 'KS bs: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD bs: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum bs: L_associated_isolated vs L_two_plus: ', p_val
    print

###############################################################################

    # min001, rigor8 stats
    print '---------------------------- min001, rigor8 ----------------------------'
    print
    print 'stats.describe(Lya_Ws_isolated): ',stats.describe(Lya_Ws_isolated_min001_rigor8)
    print 'stats.describe(Lya_Ws_L_isolated): ',stats.describe(Lya_Ws_L_isolated_min001_rigor8)
    print 'stats.describe(Lya_Ws_L_associated_isolated): ',stats.describe(Lya_Ws_associated_isolated_min001_rigor8)
    print 'stats.describe(Lya_Ws_L_associated): ',stats.describe(Lya_Ws_associated_min001_rigor8)
    print 'stats.describe(Lya_Ws_L_two_plus): ',stats.describe(Lya_Ws_two_plus_min001_rigor8)
    print
    print
    # associated vs associated_isolated
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_min001_rigor8, Lya_Ws_associated_min001_rigor8)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_min001_rigor8, Lya_Ws_associated_min001_rigor8])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_min001_rigor8, Lya_Ws_associated_min001_rigor8)

    print 'KS Lya_Ws: L_associated_isolated vs L_associated: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_associated: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_associated: ', p_val
    print
    print
    
    # associated_isolated vs two_plus
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_min001_rigor8, Lya_Ws_two_plus_min001_rigor8)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_min001_rigor8, Lya_Ws_two_plus_min001_rigor8])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_min001_rigor8, Lya_Ws_two_plus_min001_rigor8)

    print 'KS Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_two_plus: ', p_val
    print
    
    # b-parameters
    ans1 = stats.ks_2samp(bs_associated_isolated_min001_rigor8, bs_two_plus_min001_rigor8)
    ans1a = stats.anderson_ksamp([bs_associated_isolated_min001_rigor8, bs_two_plus_min001_rigor8])
    z_stat, p_val = stats.ranksums(bs_associated_isolated_min001_rigor8, bs_two_plus_min001_rigor8)

    print 'KS bs: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD bs: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum bs: L_associated_isolated vs L_two_plus: ', p_val
    print

###############################################################################

    # min001, cus stats
    print '---------------------------- min001, cus ----------------------------'
    print
    print 'stats.describe(Lya_Ws_isolated): ',stats.describe(Lya_Ws_isolated_min001_cus)
    print 'stats.describe(Lya_Ws_L_isolated): ',stats.describe(Lya_Ws_L_isolated_min001_cus)
    print 'stats.describe(Lya_Ws_L_associated_isolated): ',stats.describe(Lya_Ws_associated_isolated_min001_cus)
    print 'stats.describe(Lya_Ws_L_associated): ',stats.describe(Lya_Ws_associated_min001_cus)
    print 'stats.describe(Lya_Ws_L_two_plus): ',stats.describe(Lya_Ws_two_plus_min001_cus)
    print
    print
    # associated vs associated_isolated
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_min001_cus, Lya_Ws_associated_min001_cus)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_min001_cus, Lya_Ws_associated_min001_cus])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_min001_cus, Lya_Ws_associated_min001_cus)

    print 'KS Lya_Ws: L_associated_isolated vs L_associated: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_associated: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_associated: ', p_val
    print
    print
    
    # associated_isolated vs two_plus
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_min001_cus, Lya_Ws_two_plus_min001_cus)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_min001_cus, Lya_Ws_two_plus_min001_cus])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_min001_cus, Lya_Ws_two_plus_min001_cus)

    print 'KS Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_two_plus: ', p_val
    print
    
    # b-parameters
    ans1 = stats.ks_2samp(bs_associated_isolated_min001_cus, bs_two_plus_min001_cus)
    ans1a = stats.anderson_ksamp([bs_associated_isolated_min001_cus, bs_two_plus_min001_cus])
    z_stat, p_val = stats.ranksums(bs_associated_isolated_min001_cus, bs_two_plus_min001_cus)

    print 'KS bs: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD bs: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum bs: L_associated_isolated vs L_two_plus: ', p_val
    print

###############################################################################

    # min001, double stats
    print '---------------------------- min001, double ----------------------------'
    print
    print 'stats.describe(Lya_Ws_isolated): ',stats.describe(Lya_Ws_isolated_min001_double)
    print 'stats.describe(Lya_Ws_L_isolated): ',stats.describe(Lya_Ws_L_isolated_min001_double)
    print 'stats.describe(Lya_Ws_L_associated_isolated): ',stats.describe(Lya_Ws_associated_isolated_min001_double)
    print 'stats.describe(Lya_Ws_L_associated): ',stats.describe(Lya_Ws_associated_min001_double)
    print 'stats.describe(Lya_Ws_L_two_plus): ',stats.describe(Lya_Ws_two_plus_min001_double)
    print
    print
    # associated vs associated_isolated
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_min001_double, Lya_Ws_associated_min001_double)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_min001_double, Lya_Ws_associated_min001_double])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_min001_double, Lya_Ws_associated_min001_double)

    print 'KS Lya_Ws: L_associated_isolated vs L_associated: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_associated: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_associated: ', p_val
    print
    print
    
    # associated_isolated vs two_plus
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_min001_double, Lya_Ws_two_plus_min001_double)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_min001_double, Lya_Ws_two_plus_min001_double])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_min001_double, Lya_Ws_two_plus_min001_double)

    print 'KS Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_two_plus: ', p_val
    print
    
    # b-parameters
    ans1 = stats.ks_2samp(bs_associated_isolated_min001_double, bs_two_plus_min001_double)
    ans1a = stats.anderson_ksamp([bs_associated_isolated_min001_double, bs_two_plus_min001_double])
    z_stat, p_val = stats.ranksums(bs_associated_isolated_min001_double, bs_two_plus_min001_double)

    print 'KS bs: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD bs: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum bs: L_associated_isolated vs L_two_plus: ', p_val
    print

###############################################################################

    # min001, v150 stats
    print '---------------------------- min005, v150 ----------------------------'
    print
    print 'stats.describe(Lya_Ws_isolated): ',stats.describe(Lya_Ws_isolated_min005_v150)
    print 'stats.describe(Lya_Ws_L_isolated): ',stats.describe(Lya_Ws_L_isolated_min005_v150)
    print 'stats.describe(Lya_Ws_L_associated_isolated): ',stats.describe(Lya_Ws_associated_isolated_min005_v150)
    print 'stats.describe(Lya_Ws_L_associated): ',stats.describe(Lya_Ws_associated_min005_v150)
    print 'stats.describe(Lya_Ws_L_two_plus): ',stats.describe(Lya_Ws_two_plus_min005_v150)
    print
    print
    # associated vs associated_isolated
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_min005_v150, Lya_Ws_associated_min005_v150)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_min005_v150, Lya_Ws_associated_min005_v150])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_min005_v150, Lya_Ws_associated_min005_v150)

    print 'KS Lya_Ws: L_associated_isolated vs L_associated: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_associated: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_associated: ', p_val
    print
    print
    
    # associated_isolated vs two_plus
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_min005_v150, Lya_Ws_two_plus_min005_v150)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_min005_v150, Lya_Ws_two_plus_min005_v150])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_min005_v150, Lya_Ws_two_plus_min005_v150)

    print 'KS Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_two_plus: ', p_val
    print
    
    # b-parameters
    ans1 = stats.ks_2samp(bs_associated_isolated_min005_v150, bs_two_plus_min005_v150)
    ans1a = stats.anderson_ksamp([bs_associated_isolated_min005_v150, bs_two_plus_min005_v150])
    z_stat, p_val = stats.ranksums(bs_associated_isolated_min005_v150, bs_two_plus_min005_v150)

    print 'KS bs: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD bs: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum bs: L_associated_isolated vs L_two_plus: ', p_val
    print

###############################################################################

    # min001, v250 stats
    print '---------------------------- min005, v250 ----------------------------'
    print
    print 'stats.describe(Lya_Ws_isolated): ',stats.describe(Lya_Ws_isolated_min005_v250)
    print 'stats.describe(Lya_Ws_L_isolated): ',stats.describe(Lya_Ws_L_isolated_min005_v250)
    print 'stats.describe(Lya_Ws_L_associated_isolated): ',stats.describe(Lya_Ws_associated_isolated_min005_v250)
    print 'stats.describe(Lya_Ws_L_associated): ',stats.describe(Lya_Ws_associated_min005_v250)
    print 'stats.describe(Lya_Ws_L_two_plus): ',stats.describe(Lya_Ws_two_plus_min005_v250)
    print
    print
    # associated vs associated_isolated
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_min005_v250, Lya_Ws_associated_min005_v250)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_min005_v250, Lya_Ws_associated_min005_v250])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_min005_v250, Lya_Ws_associated_min005_v250)

    print 'KS Lya_Ws: L_associated_isolated vs L_associated: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_associated: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_associated: ', p_val
    print
    print
    
    # associated_isolated vs two_plus
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_min005_v250, Lya_Ws_two_plus_min005_v250)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_min005_v250, Lya_Ws_two_plus_min005_v250])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_min005_v250, Lya_Ws_two_plus_min005_v250)

    print 'KS Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_two_plus: ', p_val
    print
    
    # b-parameters
    ans1 = stats.ks_2samp(bs_associated_isolated_min005_v250, bs_two_plus_min005_v250)
    ans1a = stats.anderson_ksamp([bs_associated_isolated_min005_v250, bs_two_plus_min005_v250])
    z_stat, p_val = stats.ranksums(bs_associated_isolated_min005_v250, bs_two_plus_min005_v250)

    print 'KS bs: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD bs: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum bs: L_associated_isolated vs L_two_plus: ', p_val
    print
    
###############################################################################

    # cus stats
    print '---------------------------- cus ----------------------------'
    print
    print 'stats.describe(Lya_Ws_isolated): ',stats.describe(Lya_Ws_isolated_cus)
    print 'stats.describe(Lya_Ws_L_isolated): ',stats.describe(Lya_Ws_L_isolated_cus)
    print 'stats.describe(Lya_Ws_L_associated_isolated): ',stats.describe(Lya_Ws_associated_isolated_cus)
    print 'stats.describe(Lya_Ws_L_associated): ',stats.describe(Lya_Ws_associated_cus)
    print 'stats.describe(Lya_Ws_L_two_plus): ',stats.describe(Lya_Ws_two_plus_cus)
    print
    print
    # associated vs associated_isolated
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_cus, Lya_Ws_associated_cus)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_cus, Lya_Ws_associated_cus])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_cus, Lya_Ws_associated_cus)

    print 'KS Lya_Ws: L_associated_isolated vs L_associated: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_associated: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_associated: ', p_val
    print
    print
    
    # associated_isolated vs two_plus
    ans1 = stats.ks_2samp(Lya_Ws_associated_isolated_cus, Lya_Ws_two_plus_cus)
    ans1a = stats.anderson_ksamp([Lya_Ws_associated_isolated_cus, Lya_Ws_two_plus_cus])
    z_stat, p_val = stats.ranksums(Lya_Ws_associated_isolated_cus, Lya_Ws_two_plus_cus)

    print 'KS Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD Lya_Ws: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum Lya_Ws: L_associated_isolated vs L_two_plus: ', p_val
    print
    
    # b-parameters
    ans1 = stats.ks_2samp(bs_associated_isolated_cus, bs_two_plus_cus)
    ans1a = stats.anderson_ksamp([bs_associated_isolated_cus, bs_two_plus_cus])
    z_stat, p_val = stats.ranksums(bs_associated_isolated_cus, bs_two_plus_cus)

    print 'KS bs: L_associated_isolated vs L_two_plus: ',ans1
    print 'AD bs: L_associated_isolated vs L_two_plus: ',ans1a
    print 'Ranksum bs: L_associated_isolated vs L_two_plus: ', p_val
    print
    

###############################################################################
###############################################################################
###############################################################################


if __name__=="__main__":
    # do the work
    main()
    