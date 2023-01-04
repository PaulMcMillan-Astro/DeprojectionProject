import time
from Deproject_v1_0 import *
import time
from astropy.io.ascii import read as tableread
import builtins
import inspect
from termcolor import colored
from input_manipulation import *
# from Deproject_test import sanity_check
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


ti = time.time()
builtins.ti = ti


dr = DataReader(filename = 'INIT/example.ini', resample=int(args.b))
pvals, rhatvals = dr.create_sample()

if dr.polar:
    pvals, rhatvals = dr.make_polar(pvals, rhatvals)
if bool(int(args.bs)):
    pvals, rhatvals = dr.bootstrap_sample(pvals, rhatvals)

if dr.use_guess:
    mxl, fmin_it = max_func(dr.alpha, pvals, rhatvals, dr.vmin, dr.dv, dr.n, v0_guess=dr.v_guess, disp_guess=dr.disp_guess, noniso=dr.non_iso, polar=dr.polar)
elif not dr.use_guess:
    mxl, fmin_it = max_func(dr.alpha, pvals, rhatvals, vdr.vmin, dr.dv, dr.n, noniso=dr.non_iso, polar=dr.polar)

tf = time.time()
endtime = (tf - builtins.ti)/60
print("\nThe run took: ", endtime, 'mins')

n = builtins.n
dv = builtins.dv

if logging:
    # Create a folder for the run and save mxl data
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
        plot_fv(list(mxl.values())[-1],'xy',vmin,dv,n,folder,logging, polar)
        plot_fv(list(mxl.values())[-1],'yz',vmin,dv,n,folder,logging, polar)
        plot_fv(list(mxl.values())[-1],'xz',vmin,dv,n,folder,logging, polar)
    else:
        plot_fv(mxl,'xy',vmin,dv,n,folder,logging, polar)
        plot_fv(mxl,'yz',vmin,dv,n,folder,logging, polar)
        plot_fv(mxl,'xz',vmin,dv,n,folder,logging, polar)
    plot_L_and_dL(folder, logging)

print('My work here is done')
