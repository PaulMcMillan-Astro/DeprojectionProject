import os
os.chdir('/home/daniel/DeprojectionProject')
from Deproject_v1_0 import *
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

#@profile
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

optional.add_argument("-b",
                    default='0',
                    const='0',
                    nargs='?',
                    choices=['1', '0'],
                    help='Whether or not to resample the parallax and proper motions',
                    metavar='bool')

optional.add_argument("-bs",
                    default='0',
                    const='0',
                    nargs='?',
                    choices=['1', '0'],
                    help='Whether or not to bootstrap pvals and rhatvals',
                    metavar='bool')

args = parser.parse_args()

if args.f == 'max_L':
    from Deproject_v1_0 import max_L as max_func
elif args.f == 'multigrid_max_L':
    from Deproject_v1_0 import multigrid_max_L as max_func


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

def resample(plx, pmra, pmdec, plx_err, pmra_err, pmdec_err, plx_pmra_corr, plx_pmdec_corr, pmra_pmdec_corr):
    '''Use the uncertainties and correlations between plx, pmra, pmdec to resample the values.'''
    rng = np.random.default_rng()

    plx_err = plx_err.value
    pmra_err = pmra_err.value
    pmdec_err = pmdec_err.value

    for i in range(len(plx)):
        mean = [plx[i].value, pmra[i].value, pmdec[i].value]
        cov = np.array([
                       [plx_err[i]**2,
                        plx_pmra_corr[i]*plx_err[i]*pmra_err[i],
                        plx_pmdec_corr[i]*plx_err[i]*pmdec_err[i]],
                       [plx_pmra_corr[i]*plx_err[i]*pmra_err[i],
                        pmra_err[i]**2,
                        pmra_pmdec_corr[i]*pmra_err[i]*pmdec_err[i]],
                       [plx_pmdec_corr[i]*plx_err[i]*pmdec_err[i],
                        pmra_pmdec_corr[i]*pmra_err[i]*pmdec_err[i],
                        pmdec_err[i]**2]
                       ])

        plx[i], pmra[i], pmdec[i] = rng.multivariate_normal(mean, cov)*np.array([plx.unit, pmra.unit, pmdec.unit])
    return plx, pmra, pmdec

def bootstrap_sample(pvals, rhatvals):
    pvals, rhatvals = bootstrap(np.dstack([pvals, rhatvals]), bootnum=1).transpose(3, 0, 1, 2)
    return np.squeeze(pvals), np.squeeze(rhatvals)

##########You will need to change the directory path to your data#################
ti = time.time()
builtins.ti = ti
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
polar = bool(int(vars_[4][0]))
use_guess = bool(int(vars_[5][0]))
non_iso = bool(int(vars_[6][0]))
v_guess = np.array(vars_[7][0].split(',')[:],dtype=float)
disp_guess = np.array(vars_[8][0].split(',')[:],dtype=float)
alpha = float(vars_[9][0])
datafile = vars_[10][0].rstrip('\n')
logging = bool(int(vars_[11][0]))

os.chdir("DATA/")

if datafile[-4:] == '.vot':
    data_raw = QTable.read(datafile, format='votable')
    RA = data_raw['ra']
    DEC = data_raw['dec']
    plx = data_raw['parallax']
    pm_RA = data_raw['pmra']
    pm_DEC = data_raw['pmdec']
    if bool(int(args.b)):
        plx_err = data_raw['parallax_error']
        pm_RA_err = data_raw['pmra_error']
        pm_DEC_err = data_raw['pmdec_error']
        plx_pmra_corr = data_raw['parallax_pmra_corr']
        plx_pmdec_corr = data_raw['parallax_pmdec_corr']
        pmra_pmdec_corr = data_raw['pmra_pmdec_corr']

        plx, pm_RA, pm_DEC = resample(plx, pm_RA, pm_DEC, plx_err, pm_RA_err, pm_DEC_err, plx_pmra_corr, plx_pmdec_corr, pmra_pmdec_corr)
else:
    data_raw = tableread(datafile)
    dist = 1000/data_raw['parallax']*u.pc
    RA = (data_raw['ra']*u.degree)
    DEC = (data_raw['dec']*u.degree)
    plx = (data_raw['parallax']*u.mas)
    pm_RA = (data_raw['pmra']*u.mas/u.yr)
    pm_DEC = (data_raw['pmdec']*u.mas/u.yr)

print(f'Sample has {len(data_raw)} stars\n')

sample = coord.SkyCoord(ra=RA,
                        dec=DEC,
                        pm_ra_cosdec=pm_RA,
                        pm_dec=pm_DEC, 
                        distance=coord.Distance(parallax=plx),
                        radial_velocity=np.zeros(len(plx))*u.km/u.s, # This is necessary to extract pvals from astropy
                        frame='icrs').galactic

pvals = sample.velocity.d_xyz.T.unmasked.value
rhatvals = sample.spherical.unit_vectors()['distance'].xyz.T.unmasked.value

if polar:
    pvals, rhatvals = make_polar(sample, pvals, rhatvals)

#pvals, rhatvals = calc_p_rhat(sample, polar=polar)
#pvals = pvals.unmasked.value

if bool(int(args.bs)):
    pvals, rhatvals = bootstrap_sample(pvals, rhatvals)

if use_guess:
    mxl, fmin_it = max_func(alpha, pvals, rhatvals, vmin, dv, n ,v0_guess=v_guess, disp_guess=disp_guess, noniso=non_iso, polar=polar)
elif not use_guess:
    mxl, fmin_it = max_func(alpha, pvals, rhatvals, vmin, dv, n, noniso=non_iso, polar=polar)
tf = time.time()
endtime = (tf - builtins.ti)/60
print("\nThe run took: ", endtime, 'mins')

n = builtins.n
dv = builtins.dv

if logging:
    # Create a folder for the run and save mxl data
    os.chdir('/home/daniel/DeprojectionProject')
    folder = make_folder()
    np.save('RUNS/' + folder + '/mxl_data', mxl)

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
        logfile.write('Time needed : ' + str(endtime/60) + ' hrs\n')

    # Logfile in RUNS folder
    with open('RUNS/' + folder + '/log.txt', 'a') as logfile:
        logfile.write('Datafile    : ' + datafile + '\n')
        logfile.write('fmin its    : ' + str(fmin_it) + '\n')
        logfile.write('Time needed : ' + str(endtime/60) + ' hrs\n')
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

        sanity_check(pvals,rhatvals,list(mxl.values())[-1],vmin,dv,n,logging)

    elif sane == 'n':

        print('Suit yourself')

    shouldiplot = input('Do you want to plot your results [y/n]? ')
    while shouldiplot != 'y' and shouldiplot != 'n':
        shouldiplot = input('Incorrect entry, try again [y/n]! ')

    if shouldiplot == 'y':
        from Deproject_plots import plot_fv,plot_L_and_dL

        plot_fv(list(mxl.values())[-1],input('What plane should I project onto? '),vmin,dv,n,folder,logging)

        s=0

        while True:
            if s==2:
                break
            plotagain = input('Do you want to plot another plane [y/n]? ')
            if plotagain == 'y':
                plot_fv(list(mxl.values())[-1],input('What plane should I project onto [xy/yz/xz]? '),vmin,dv,n,folder,logging)
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
    if args.f == 'multigrid_max_L':
        #sanity_check(pvals,rhatvals,mxl,vmin,dv,n,logging,folder)
        plot_fv(list(mxl.values())[-1],'xy',vmin,dv,n,folder,logging, polar)
        plot_fv(list(mxl.values())[-1],'yz',vmin,dv,n,folder,logging, polar)
        plot_fv(list(mxl.values())[-1],'xz',vmin,dv,n,folder,logging, polar)
    else:
        plot_fv(mxl,'xy',vmin,dv,n,folder,logging, polar)
        plot_fv(mxl,'yz',vmin,dv,n,folder,logging, polar)
        plot_fv(mxl,'xz',vmin,dv,n,folder,logging, polar)
    plot_L_and_dL(folder, logging)

print('My work here is done')
