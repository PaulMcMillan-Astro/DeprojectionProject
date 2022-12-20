#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 15:20:40 2021

@author: daniel
"""
import os
os.chdir('/home/daniel/DeprojectionProject')
from Profile_v1_0 import *
import time
from astropy.io.ascii import read as tableread
from astropy.table import QTable
from astropy.stats import bootstrap
from datetime import date
import builtins
import string
import inspect
from termcolor import colored
from Deproject_test import sanity_check
import matplotlib
#from memory_profiler import LogFile

matplotlib.use('TkAgg')
from Deproject_plots import plot_fv, plot_L_and_dL, DM_plt_prefs
DM_plt_prefs()
"""Script that when called by Python in terminal will perform the computations needed
to run a maximisation scheme. It will also allow you to plot the results and the
change in L & dL over all iterations.

If provided the file 'vars.ini' as an argument such that the terminal command is
'python Deproject_exe.py -i vars.ini' it will read the input variables automatically"""

def process_input(file_): #Function that reads the input file if any
    h_file = []
    input_ = open(file_, 'r')
    for line in input_:
        if line.find("#") != -1:
            continue
        elif line.find("\n") == 0:
            continue
        else:
            h_file.append(line.split('\t'))
    return h_file

# Read the input file
guessfile = 'profile_vars.ini' #The vars.ini file
vars_ = process_input(guessfile)
N = int(vars_[0][0])
n = np.array(vars_[1][0].split(',')[:],dtype=int)
vmin = np.array(vars_[2][0].split(',')[:],dtype=float)
dv = np.array(vars_[3][0].split(',')[:],dtype=float)
use_guess = bool(int(vars_[4][0]))
non_iso = bool(int(vars_[5][0]))
v_guess = np.array(vars_[6][0].split(',')[:],dtype=float)
disp_guess = np.array(vars_[7][0].split(',')[:],dtype=float)
alpha = float(vars_[8][0])
datafile = vars_[9][0].rstrip('\n')
logging = bool(int(vars_[10][0]))

os.chdir("DATA/")

data_raw = QTable.read(datafile, format='votable')
dist = data_raw['dist']
RA = data_raw['ra']
DEC = data_raw['dec']
plx = data_raw['parallax']
pm_RA = data_raw['pmra']
pm_DEC = data_raw['pmdec']


sample_icrs = coord.ICRS(ra = RA, dec = DEC, pm_ra_cosdec = pm_RA, pm_dec = pm_DEC, distance=coord.Distance(parallax=plx))

sample = sample_icrs.transform_to(coord.Galactic)

pvals, rhatvals = calc_p_rhat(sample)
pvals = pvals.value

ti = time.time()
builtins.ti = ti

mxl, fmin_it = multigrid_max_L(alpha, pvals, rhatvals, vmin, dv, n ,v0_guess=v_guess, disp_guess=disp_guess, noniso=non_iso)
    
tf = time.time()
endtime = (tf - builtins.ti)/60
print("\nThe run took: ", endtime, 'mins')