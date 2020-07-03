from Deproject_v1_0_testing import *
from sys import argv
import os
import time
from astropy.io.ascii import read as tableread

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

ti_a = time.time()

guessfile = "vars.ini"

# READ VARS.INI
#guessfile = argv[1] #The vars.ini file
vars_ = process_input(guessfile)
N = int(vars_[0][0])
n = np.array(vars_[1][0].split(',')[:],dtype=int)
vmin = np.array(vars_[2][0].split(',')[:],dtype=float)
dv = np.array(vars_[3][0].split(',')[:],dtype=float)
v_guess = np.array(vars_[4][0].split(',')[:],dtype=float)
disp_guess = np.array(vars_[5][0].split(',')[:],dtype=float)
alpha = float(vars_[6][0])
datafile = vars_[7][0].rstrip('\n')

# EXTRACT SAMPLE DATA
try:
        os.chdir("DATA/")
except FileNotFoundError:
    print('Edit your desired path in the script')
data_raw = tableread(str(datafile))

dist = 1000/data_raw['parallax']*u.pc        
RA = (data_raw['ra']*u.degree)
DEC = (data_raw['dec']*u.degree)
plx = (data_raw['parallax']*u.mas)
pm_RA = (data_raw['pmra']*u.mas/u.yr)
pm_DEC = (data_raw['pmdec']*u.mas/u.yr)
    
    
    
    
#SET UP DATA

sample_icrs = coord.ICRS(ra = RA, dec = DEC, pm_ra_cosdec = pm_RA, pm_dec = pm_DEC,distance=dist)

sample = sample_icrs.transform_to(coord.Galactic)
   
# CALC THINGS 
pvals, rhatvals = calc_p_rhat(sample)

max_L(alpha, pvals, rhatvals, vmin, dv, 1, n,v0_guess=v_guess, disp_guess=disp_guess)
max_L(alpha, pvals, rhatvals, vmin, dv, 2, n,v0_guess=v_guess, disp_guess=disp_guess)

tf_a = time.time()
print("The run took", tf_a - ti_a, 's: a')