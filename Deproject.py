# Script that when called by Python in terminal will perform the computations needed to run a maximisation scheme. It will also allow you to plot the results and the change in L & dL over all iterations.

# If provided the file 'vars.ini' as an argument such that the terminal command is 'python Deproject_exe.py -i vars.ini' it will read the input variables automatically

import time
import builtins
from termcolor import colored
from input_manipulation import DataReader
from likelihood_estimation import Deprojection
from log_and_plot import Logger
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

ti = time.time()
builtins.ti = ti

dr = DataReader(inifile=args.i, resample=int(args.b))
pvals, rhatvals = dr.create_sample()

if dr.polar:
    pvals, rhatvals = dr.make_polar(pvals, rhatvals)
if bool(int(args.bs)):
    pvals, rhatvals = dr.bootstrap_sample(pvals, rhatvals)

if dr.use_guess:
    dep = Deprojection(dr.alpha, pvals, rhatvals, dr.vmin, dr.dv, dr.n, v0_guess=dr.v_guess, disp_guess=dr.disp_guess, noniso=dr.non_iso, polar=dr.polar)
elif not dr.use_guess:
    dep = Deprojection(dr.alpha, pvals, rhatvals, dr.vmin, dr.dv, dr.n, noniso=dr.non_iso, polar=dr.polar)

if args.f == 'max_L':
    mxl, fmin_it = dep.max_L()
elif args.f == 'multigrid_max_L':
    mxl, fmin_it = dep.multigrid_max_L()

tf = time.time()
endtime = (tf - builtins.ti)/60
print("\nThe run took: ", endtime, 'mins')

n = builtins.n
dv = builtins.dv

lp = Logger()
lp.make_folder()
lp.log_mle_run(mxl, fmin_it, endtime, n, dr.vmin, dv, dr.polar, dr.use_guess, dr.non_iso, dr.v_guess, dr.disp_guess, dr.alpha, dr.datafile)
lp.plot_and_save()
lp.plot_L_and_dL()

print('My work here is done')
