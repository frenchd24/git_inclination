#!/usr/bin/env python



import sys
import csv
import math
from pylab import *
import correlateSingle11 as correlateSingle
from utilities import *
import getpass
import pickle

maxSep = 500.
agnSeparation = 8000.
minVcorr = 450.
minSize = 0.
max_deltav = 400.
min_likelihood = 0.01
rigor = 5
target = 'MRK290'
                
def do_a_thing(a):
    a*5
    print 'yes'


# correlation = correlateSingle.correlateTarget(target, maxSep, agnSeparation, minVcorr, minSize, slow=False, searchAll=False, returnDict=True)

# print 'correlation: ',correlation

# for c in correlation:
#     print 'c: ',c, correlation[c]['Vhel']
#     
# 
# print
# do_a_thing(5)


isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/isolated8_double.p'
L_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_isolated8_double.p'
L_associated_isolated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated_isolated8_double.p'
L_associated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_associated8_double.p'
L_nonassociated_filename = '/Users/frenchd/Research/inclination/git_inclination/L_nonassociated8_double.p'
L_two_filename = '/Users/frenchd/Research/inclination/git_inclination/L_two8_double.p'
L_three_plus_filename = '/Users/frenchd/Research/inclination/git_inclination/L_three_plus8_double.p'
L_group_filename = '/Users/frenchd/Research/inclination/git_inclination/L_group8_double.p'
all_filename = '/Users/frenchd/Research/inclination/git_inclination/all8_double.p'


morph_dict_filename = '/Users/frenchd/Research/inclination/git_inclination/morph_dict2.p'
new_morph_dict_filename = '/Users/frenchd/Research/inclination/git_inclination/morph_dict3.p'


isolated_file = open(isolated_filename,'r')
L_isolated_file = open(L_isolated_filename,'r')
L_associated_isolated_file = open(L_associated_isolated_filename,'r')
L_associated_file = open(L_associated_filename,'r')
L_nonassociated_file = open(L_nonassociated_filename,'r')
L_two_file = open(L_two_filename,'r')
L_three_plus_file = open(L_three_plus_filename,'r')
L_group_file = open(L_group_filename,'r')
all_file = open(all_filename,'r')

isolated_pickle = pickle.load(isolated_file)
L_isolated_pickle = pickle.load(L_isolated_file)
L_associated_isolated_pickle = pickle.load(L_associated_isolated_file)
L_associated_pickle = pickle.load(L_associated_file)
L_nonassociated_pickle = pickle.load(L_nonassociated_file)
L_two_pickle = pickle.load(L_two_file)
L_three_plus_pickle = pickle.load(L_three_plus_file)
L_group_pickle = pickle.load(L_group_file)
all_pickle = pickle.load(all_file)

new_morph_dict_file = open(new_morph_dict_filename,'wt')

morph_dict_file = open(morph_dict_filename,'rU')
morph_dict = pickle.load(morph_dict_file)
morph_dict_file.close()


# print "isolated_pickle['Names']: ", isolated_pickle['Names']
# print "isolated_pickle['targets']: ", isolated_pickle['targets']
# print "isolated_pickle['Lya_vs']: ", isolated_pickle['Lya_vs']
# print
# print
# print
# print "L_isolated_pickle['Names']: ",L_isolated_pickle['Names']
# print "L_isolated_pickle['targets']: ",L_isolated_pickle['targets']
# print "L_isolated_pickle['Lya_vs']: ",L_isolated_pickle['Lya_vs']
# print
# print
# print
# print
# print 'L_associated_isolated_pickle: ',L_associated_isolated_pickle
# print "L_associated_isolated_pickle['Names']: ",L_associated_isolated_pickle['Names']
# print "L_associated_isolated_pickle['targets']: ",L_associated_isolated_pickle['targets']
# print "L_associated_isolated_pickle['Lya_vs']: ",L_associated_isolated_pickle['Lya_vs']
# print
# print
# print
# print

# Names = L_associated_isolated_pickle['Names']
# targets = L_associated_isolated_pickle['targets']
# Lya_vs = L_associated_isolated_pickle['Lya_vs']
# ls = L_associated_isolated_pickle['ls']


# Names = L_associated_pickle['Names']
# targets = L_associated_pickle['targets']
# Lya_vs = L_associated_pickle['Lya_vs']
# ls = L_associated_pickle['ls']


# Names = L_associated_pickle['Names']
# targets = L_associated_pickle['targets']
# Lya_vs = L_associated_pickle['Lya_vs']
# ls = L_associated_pickle['ls']


# Names = L_group_pickle['Names']
# targets = L_group_pickle['targets']
# Lya_vs = L_group_pickle['Lya_vs']
# ls = L_group_pickle['ls']

all_Lya_Ws = all_pickle['Lya_Ws']

all_Lya_Ws.sort(reverse=True)



print "Max(all['Lya_Ws']) = ",max(all_Lya_Ws)
print
print "highest 10: ", all_Lya_Ws[:10]
print
print 'len(all_Lya_Ws): ',len(all_Lya_Ws)
print
print "len(isolated_pickle): ", len(isolated_pickle["Lya_Ws"])
print
print "len(L_isolated_pickle): ", len(L_isolated_pickle["Lya_Ws"])
print
print "len(L_associated_isolated_pickle): ", len(L_associated_isolated_pickle["Lya_Ws"])
print
print "len(L_associated_pickle): ", len(L_associated_pickle["Lya_Ws"])
print
print "len(L_nonassociated_pickle): ", len(L_nonassociated_pickle["Lya_Ws"])
print
print "len(L_two_pickle): ", len(L_two_pickle["Lya_Ws"])
print
print "len(L_three_plus_pickle): ", len(L_three_plus_pickle["Lya_Ws"])
print
print "len(L_group_pickle): ", len(L_group_pickle["Lya_Ws"])
print

target_dict = {}
for i in all_pickle['targets']:
    if not target_dict.has_key(i):
        target_dict[i] = 1
        
old_morph_dict = {}
new_morph_dict = {}

print
print 'len(targets_dict): ',len(target_dict)

print
print
print '----- L_associated ----'
for morph, target, name in zip(L_associated_pickle['MTypes'], L_associated_pickle['targets'], L_associated_pickle['Names']):
    if not old_morph_dict.has_key(name):
        old_morph_dict[name] = morph
    

print '----- L_associated_isolated ----'
for i, target, name in zip(L_associated_isolated_pickle['MTypes'], L_associated_isolated_pickle['targets'], L_associated_isolated_pickle['Names']):
    if not old_morph_dict.has_key(name):
        old_morph_dict[name] = morph




for name in old_morph_dict:
    if not morph_dict.has_key(name):
        print '{0} - {1}'.format(name, old_morph_dict[name])
    
        ans = 'x'
        while ans != 's0' and ans != 'e' and ans != 'i' and ans != 'sa' and ans != 'sb' and ans != 'q':
            ans = raw_input("Enter type: 's0' = 'S0', 'e' = 'E', 'i' = 'Im', 'sa' = 'SA', 'sb' = 'SB' : ")
        
        if ans == 'q':
            sys.exit()

        else:
            print 'Saving {0} = {1}'.format(name,ans)
            print
            new_morph_dict[name] = ans

    else:
        new_morph_dict[name] = morph_dict[name]
    
    

# for name, target, lya_v, l in zip(Names, targets, Lya_vs, ls):
#     if target == 'PG0003+158':
#         print '{0} : {1} - {2} : {3}'.format(name, target, lya_v, l)
#     
#     
# print
# print 'Total number: ',len(Names)
# print
# print 'Done'


isolated_file.close()
L_isolated_file.close()
L_associated_isolated_file.close()
L_associated_file.close()
L_nonassociated_file.close()
L_two_file.close()
L_three_plus_file.close()
L_group_file.close()


pickle.dump(new_morph_dict, new_morph_dict_file)
new_morph_dict_file.close()

print 'done'