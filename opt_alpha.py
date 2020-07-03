from Deproject_v1_0 import *
from sys import argv
import time 
import os
from decimal import Decimal
from astropy.io.ascii import read as tableread
import pandas as pd
from astropy.io.votable import parse_single_table
os.chdir('/home/mikkola/Documents/DeprojectionProject')
from datetime import date
from alpha_debugger import *
import builtins


"""This script will automatically initiate an optimisation of the alpha
smoothing parameter. Given an input file, i.e. 'alpha_vars.ini', it will 
load all parameters for the dimensions of the box which represents v-space.

If an input file is not provided, one can input these values manually via
inputs."""
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
        guessfile = argv[1] #Here we read the input file and all parameters
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
    
    if datafile=="Distances_PJM2017.csv":
        try:
            data_raw
        except NameError:
            
            try:
                os.chdir("DATA/")
            except FileNotFoundError:
                print('Edit your desired path in the script')
            data_raw = Table.read(str(datafile))
            pass
        flagset = 'flag_any'
        data_raw = filt_data(data_raw,flagset)
    
        dist = data_raw['distance']*u.pc
        filt = dist < 200*u.pc
        dist = dist[filt]
        
        RA = (data_raw['RAdeg']*u.degree)[filt]
        DEC = (data_raw['DEdeg']*u.degree)[filt]
        plx = (data_raw['parallax']*u.mas)[filt]
        pm_RA = (data_raw['pmRA_TGAS']*u.mas/u.yr)[filt]
        pm_DEC = (data_raw['pmDE_TGAS']*u.mas/u.yr)[filt]
    
    
        filt_data(data_raw,flagset)
    
    
    
        #We now just set up our data and are ready to apply cuts to it
    
        sample_icrs = coord.ICRS(ra = RA, dec = DEC, pm_ra_cosdec = pm_RA, pm_dec = pm_DEC,distance=dist)        
        sample = sample_icrs.transform_to(coord.Galactic)
            
    elif datafile[-5:] == 'table':
        try:
            data_raw
        except NameError:
            
            try:
                os.chdir("DATA/")
            except FileNotFoundError:
                print('Edit your desired path in the script')
            data_raw = tableread(str(datafile))
            pass
        
    
        dist = 1000/data_raw['parallax']*u.pc        
        RA = (data_raw['ra']*u.degree)
        DEC = (data_raw['dec']*u.degree)
        plx = (data_raw['parallax']*u.mas)
        pm_RA = (data_raw['pmra']*u.mas/u.yr)
        pm_DEC = (data_raw['pmdec']*u.mas/u.yr)
    
    
    
        sample_icrs = coord.ICRS(ra = RA, dec = DEC, pm_ra_cosdec = pm_RA, pm_dec = pm_DEC,distance=dist)
        sample = sample_icrs.transform_to(coord.Galactic)
            
    elif datafile == 'dehnen_binney_sample1.txt':
        table = parse_single_table("hipparcos.vot").to_table()
        df1 = table.to_pandas()
        
        DB = np.loadtxt("dehnen_binney_sample1.txt",dtype=int)
        df2 = pd.DataFrame(DB,columns=["HIP"])
        
        df = pd.merge(df1, df2, how='inner', on='HIP')
        
        HIP      = df["HIP"].values
        RA       = (df["RAICRS"].values*u.degree)
        DEC      = (df["DEICRS"].values*u.degree)
        plx      = (df["Plx"].values*u.mas)
        pm_RA    = (df["pmRA"].values*u.mas/u.yr)
        pm_DEC   = (df["pmDE"].values*u.mas/u.yr)
        e_RA     = (df["e_RAICRS"].values*u.degree)
        e_DEC    = (df["e_DEICRS"].values*u.degree)
        e_plx    = (df["e_Plx"].values*u.mas)
        e_pm_RA  = (df["e_pmRA"].values*u.mas/u.yr)
        e_pm_DEC = (df["e_pmDE"].values*u.mas/u.yr)    
        
        sample_icrs = coord.ICRS(ra = RA, dec = DEC, pm_ra_cosdec = pm_RA, pm_dec = pm_DEC,distance=dist)
        sample = sample_icrs.transform_to(coord.Galactic)
    os.chdir("/home/mikkola/Documents/DeprojectionProject")
    # Running the optimisation           
    if argv[2] == 'tenstep':
        stdscr.addstr(5, 1, 'Started 1st tenstep iteration...', curses.color_pair(0) | curses.A_BOLD)
        stdscr.refresh()
        alpha_fin = opt_alpha(stdscr, alpha0, M, N, sample, vmin, dv, n, opt_tol=opt_tol, mise_tol=mise_tol, noniso=non_iso)
    elif argv[2] == 'ternary':
        stdscr.addstr(5, 1, 'Started 1st ternary iteration...', curses.color_pair(0) | curses.A_BOLD)
        stdscr.refresh()
        alpha_fin = opt_alpha_ternary(stdscr, alpha0, M, N, sample, vmin, dv, n, opt_tol=opt_tol, mise_tol=mise_tol, noniso=non_iso)
    elif argv[2] == 'gss':
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