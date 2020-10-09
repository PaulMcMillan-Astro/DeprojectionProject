import os
os.chdir('/home/mikkola/Documents/DeprojectionProject')
from Deproject_v1_0 import *
import time
from astropy.io.ascii import read as tableread
from datetime import date
import builtins
import string
import inspect
from termcolor import colored
from Deproject_test import sanity_check
import matplotlib
matplotlib.use('TkAgg')
from Deproject_plots import plot_fv, plot_L_and_dL, DM_plt_prefs
DM_plt_prefs()
"""Script that when called by Python in terminal will perform the computations needed
to run a maximisation scheme. It will also allow you to plot the results and the 
change in L & dL over all iterations. 

If provided the file 'vars.ini' as an argument such that the terminal command is
'python Deproject_exe.py -i vars.ini' it will read the input variables automatically"""

import argparse

parser = argparse.ArgumentParser()
parser._action_groups.pop()

required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

required.add_argument('-i',
                    default='vars.ini',
                    help='datafile, typically vars.ini',
                    metavar='vars.ini',
                    required=True)

optional.add_argument("-a",
                    default='False',
                    const='False',
                    nargs='?',
                    choices=['True', 'False'],
                    help='Automatic plotting [default : False]',
                    metavar='bool')
optional.add_argument("-f",
                    default='multigrid_max_L',
                    const='multigrid_max_L',
                    nargs='?',
                    choices=['max_L', 'multigrid_max_L'],
                    help='Select PMLE scheme [default : multigrid_max_L]',
                    metavar='func')
                
args = parser.parse_args()

if args.f == 'max_L':
    from Deproject_v1_0 import max_L as max_func
elif args.f == 'multigrid_max_L':
    from Deproject_v1_0 import multigrid_max_L as max_func

##########You will need to change the directory path to your data#################
ti_a = time.time()

# Function that creates a non-existing folder to store output data
def make_folder():
    '''Function that writes a folder with name YYYY-MM-DDx where x 
    is the first letter that is a non-existing directory'''
    alp = string.ascii_lowercase
    date_str = str(date.today())
    
    list_str = []
    for letter in alp:
        list_str.append(date_str + letter)

    for letter1 in alp:
        for letter2 in alp:
            list_str.append(date_str + letter1 + letter2)

    list_str = np.array(list_str)


    os_dirs = os.popen('ls RUNS | grep "%s"' % date_str).read().split("\n")[:-1]
    os_dirs.sort(key=lambda x: len(x))
    os_dirs = np.array(os_dirs)

    existing_dirs = np.zeros(len(list_str),dtype='<U12')
    existing_dirs[:len(os_dirs)] = os_dirs

    folder_name=list_str[list_str != existing_dirs][0]
    os.system('mkdir RUNS/' + folder_name)
                
    return folder_name


##########You will need to change the directory path to your data#################
ti_a = time.time()

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
guessfile = args.i #The vars.ini file
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
print("Sample has " + str(len(dist)) + " stars\n")



sample_icrs = coord.ICRS(ra = RA, dec = DEC, pm_ra_cosdec = pm_RA, pm_dec = pm_DEC,distance=dist)

sample = sample_icrs.transform_to(coord.Galactic)


pvals, rhatvals = calc_p_rhat(sample) 
    
mxl, fmin_it, mxls, ns = max_func(alpha, pvals, rhatvals, vmin, dv, n ,v0_guess=v_guess, disp_guess=disp_guess, noniso=non_iso)
tf_a = time.time()
endtime = (tf_a - ti_a)/60
print("\nThe run took: ", endtime, 'mins')


def getXY(dv, vmin, n, nmax):
    dx, dy = dv[0], dv[1]
    vxmin, vymin = vmin[0], vmin[1]
    nmaxx, nmaxy = nmax[0], nmax[1]
    nx, ny = n[0], n[1]
    
    vxmax, vymax = vxmin+nmaxx*dx,vymin+nmaxy*dy

    x0, y0 = vxmin, vymin
    x1, y1 = vxmax, vymax

    xbins = np.linspace(x0,x1,nx+1)
    ybins = np.linspace(y0,y1,ny+1)

    [X,Y] = np.meshgrid(xbins,ybins);
    return X,Y

f,ax = plt.subplots(len(ns),1,sharex=True,sharey=True,figsize=(20/len(ns),20),frameon=False)
for i, mxl in enumerate(mxls):

    X,Y = getXY(dv, vmin, ns[i], ns[-1])
    twodfv = np.sum(np.exp(mxl.reshape(ns[i]))+1, axis=2)
    ax[i].pcolormesh(X,Y,twodfv.T, cmap = plt.cm.get_cmap('bone_r'))
    ax[i].set_title('Deprojected fv')


plt.show()
