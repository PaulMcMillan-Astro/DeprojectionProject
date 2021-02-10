import os
os.chdir('/home/daniel/DeprojectionProject')
from Deproject_v1_0 import *
from Deproject_oa_scripts import *
import time
from astropy.io.ascii import read as tableread
from astropy.table import QTable
from datetime import date
import builtins
import string
import inspect
from termcolor import colored
from Deproject_test import sanity_check
import matplotlib
from memory_profiler import LogFile
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
                    default='0',
                    const='0',
                    nargs='?',
                    choices=['1', '0'],
                    help='Automatic plotting [default : 0]',
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
    #from Deproject_v1_0 import max_L as max_func
    from Deproject_oa_scripts import max_L_oa as max_func
elif args.f == 'multigrid_max_L':
    # from Deproject_v1_0 import multigrid_max_L as max_func
    from Deproject_oa_scripts import multigrid_max_L as max_func

    
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

os.chdir("DATA/")

data = Table.read(datafile, format='votable')

pvals = data['pval']
rhatvals = data['rhat']


if use_guess:
    mxl, fmin_it = max_func(alpha, pvals, rhatvals, vmin, dv, n ,v0_guess=v_guess, disp_guess=disp_guess, noniso=non_iso)
elif not use_guess:
    mxl, fmin_it = max_func(alpha, pvals, rhatvals, vmin, dv, n, noniso=non_iso)
tf_a = time.time()
endtime = (tf_a - ti_a)/60
print("\nThe run took: ", endtime, 'mins')

n = builtins.n
dv = builtins.dv

if logging:
    # Create a folder for the run and save mxl data
    os.chdir('/home/daniel/DeprojectionProject')
    folder = make_folder()
    np.save('RUNS/' + folder + '/mxl_data',mxl)
    
    # Save output
    # RUNS folder identifier
    with open('logs/log_dir_identifier.txt', 'a') as logfile:
        if folder[10] == 'a' and len(folder) == 11:
            mark = '='
            logfile.write('\n' + mark*120 + '\n')
            logfile.write(mark*55 + folder[:10] + mark*55 + '\n')
            logfile.write(mark*120 + '\n')
            
        logfile.write('\nFolder name : ' + folder + '\n')
        logfile.write('Datafile    : ' + datafile + '\n')
        logfile.write('fmin its    : ' + str(fmin_it) + '\n')
        logfile.write('Time needed : ' + str(endtime/60) + ' mins\n')
        
    # Logfile in RUNS folder 
    with open('RUNS/' + folder + '/log.txt', 'a') as logfile:
        logfile.write('Datafile    : ' + datafile + '\n')
        logfile.write('fmin its    : ' + str(fmin_it) + '\n')
        logfile.write('Time needed : ' + str(endtime/60) + ' mins\n')
        logfile.write('Labels      : Nbins[1x3], vmin[1x3], bin size, use_guess, noniso, mu_guess, sigma_guess, alpha\n')    
        value_string=str((str(list(n)).replace(",",":").replace(":",""),str(list(vmin)).replace(",",":").replace(":","")
                          ,str(list(dv)).replace(",",":").replace(":",""),use_guess,non_iso,str(list(v_guess)).replace(",",":").replace(":","")
                          , str(list(disp_guess)).replace(",",":").replace(":",""), alpha)).replace("'","")[1:-1]
        
        logfile.write("Values      : " + value_string + '\n')
        
    # MORE IS SAVE BY SANITY_CHECK()
else:
    folder = ''
        
builtins.autoplot = bool(int(args.a))
if not autoplot:
    sane = input('Do you want to perform a sanity check [y/n]? ')
    while sane != 'y' and sane != 'n':   
        sane = input('Incorrect entry, try again [y/n]! ')
    
    if sane == 'y':
        
        from Deproject_test import sanity_check
        
        sanity_check(pvals,rhatvals,mxl,vmin,dv,n,logging)
    
    elif sane == 'n':
        
        print('Suit yourself')
        
    shouldiplot = input('Do you want to plot your results [y/n]? ')
    while shouldiplot != 'y' and shouldiplot != 'n':   
        shouldiplot = input('Incorrect entry, try again [y/n]! ')
    
    if shouldiplot == 'y':
        from Deproject_plots import plot_fv,plot_L_and_dL
        
        plot_fv(mxl,input('What plane should I project onto? '),vmin,dv,n,folder,logging)
        
        s=0
        
        while True:
            if s==2:
                break
            plotagain = input('Do you want to plot another plane [y/n]? ')
            if plotagain == 'y':
                plot_fv(mxl,input('What plane should I project onto [xy/yz/xz]? '),vmin,dv,n,folder,logging)
                s+=1
                continue
            else:
                break
    
    shouldiplotL = input('Do you want to plot the change in L and dL during the maximisation [y/n]? ')
    while shouldiplotL != 'y' and shouldiplotL != 'n':   
        shouldiplotL = input('Incorrect entry, try again [y/n]! ')
    
    if shouldiplotL=='y':
        plot_L_and_dL(folder,logging)
        
else:
    #sanity_check(pvals,rhatvals,mxl,vmin,dv,n,logging,folder)
    plot_fv(mxl,'xy',vmin,dv,n,folder,logging)  
    plot_fv(mxl,'yz',vmin,dv,n,folder,logging)  
    plot_fv(mxl,'xz',vmin,dv,n,folder,logging)    
    plot_L_and_dL(folder, logging)

print('My work here is done')
