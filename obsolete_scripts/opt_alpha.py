import os
os.chdir('/home/daniel/DeprojectionProject')
from Deproject_v1_0 import *
from Deproject_oa_scripts import *
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
import argparse

matplotlib.use('TkAgg')
from Deproject_plots import plot_fv, plot_L_and_dL, DM_plt_prefs
DM_plt_prefs()


"""This script will automatically initiate an optimisation of the alpha
smoothing parameter. Given an input file, i.e. 'alpha_vars.ini', it will 
load all parameters for the dimensions of the box which represents v-space.

If an input file is not provided, one can input these values manually via
inputs."""

parser = argparse.ArgumentParser()
parser._action_groups.pop()

required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

required.add_argument('-i',
                    default='vars.ini',
                    help='datafile, typically vars.ini',
                    metavar='vars.ini',
                    required=True)

args = parser.parse_args()

date_str = str(date.today())
start_time = time.time()
def alpha_script(stdscr):
    stdscr.addstr(2, 1, 'Started runing opt_alpha.py', curses.color_pair(0) | curses.A_BOLD)
    stdscr.refresh()
    
    def process_input(file_):
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
    
    
    
    try:
        argv[1]
    except (IndexError,NameError):
        raise NameError("Couldn't access the input file. Are you sure you wrote it correctly?")
    else:
        # Read the input file
        guessfile = args.i
        vars_ = process_input(guessfile)
        N = int(vars_[0][0])
        M = int(vars_[1][0])
        n = np.array(vars_[2][0].split(',')[:],dtype=int)
        vmin = np.array(vars_[3][0].split(',')[:],dtype=float)
        dv = np.array(vars_[4][0].split(',')[:],dtype=float)
        alpha0 = float(vars_[5][0])
        non_iso = bool(int(vars_[6][0]))
        opt_tol = float(vars_[7][0])
        mise_tol = float(vars_[8][0])
        datafile = vars_[9][0].rstrip('\n')

        os.chdir("DATA/")
    
    if datafile[-4:] == '.vot':
        data_raw = QTable.read(datafile, format='votable')
        dist = data_raw['dist']
        RA = data_raw['ra']
        DEC = data_raw['dec']
        plx = data_raw['parallax']
        pm_RA = data_raw['pmra']
        pm_DEC = data_raw['pmdec']
   else:
        data_raw = tableread(datafile)
        dist = 1000/data_raw['parallax']*u.pc
        RA = (data_raw['ra']*u.degree)
        DEC = (data_raw['dec']*u.degree)
        plx = (data_raw['parallax']*u.mas)
        pm_RA = (data_raw['pmra']*u.mas/u.yr)
        pm_DEC = (data_raw['pmdec']*u.mas/u.yr)

    sample_icrs = coord.ICRS(ra = RA, dec = DEC, pm_ra_cosdec = pm_RA, pm_dec = pm_DEC, distance=coord.Distance(parallax=plx))

    sample = sample_icrs.transform_to(coord.Galactic)

    pvals, rhatvals = calc_p_rhat(sample)
    pvals = pvals.value

    # Running the optimisation           
    stdscr.addstr(5, 1, 'Started 1st gss iteration...', curses.color_pair(0) | curses.A_BOLD)
    stdscr.refresh()
    alpha_fin = opt_alpha_gss(stdscr, alpha0, M, N, sample, vmin, dv, n, opt_tol=opt_tol, mise_tol=mise_tol, noniso=non_iso)
        

    endtime = time.time() - start_time

    with open('logs/log_alpha-opt.txt', 'a') as logfile:
        logfile.write('\n____________________________' + date_str +'______________________________________\n')
        logfile.write('Datafile    : ' + datafile + '\n')
        logfile.write('Time needed : ' + str(endtime/60) + ' mins ('+ argv[2] + ')\n')
        logfile.write('Final alpha : ' + str(alpha_fin) + '\n')
        logfile.write('Labels      : N, M, Nbins[1x3], vmin[1x3], bin size, alpha0, noniso\n')
        
        value_string=str((N,M
                          ,str(list(n)).replace(",",":").replace(":",""),str(list(vmin)).replace(",",":").replace(":",""),str(list(dv)).replace(",",":").replace(":","")
                          ,alpha0,non_iso)).replace("'","")[1:-1]
        
        logfile.write("Values      : " + value_string + '\n')
    return
stdscr = curses.initscr()
debugger()
curses.wrapper(alpha_script)

print(str('The run took %s hours' % np.around(endtime/3600,decimals=2)))
