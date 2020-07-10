from sys import argv
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
from Deproject_plots import plot_fv,plot_L
"""Script that when called by Python in terminal will perform the computations needed
to run a maximisation scheme. It will also allow you to plot the results and the 
change in L over all iterations. 

If provided the file 'vars.ini' as an argument such that the terminal command is
'python Deproject_exe.py vars.ini' it will read the input variables automatically"""
try: 
    in_str     = argv[0] + ' ' + argv[1] + ' ' + argv[2]
except IndexError:
    in_str     = argv[0] + ' ' + argv[1]
    in_allowed = 'Deproject_exe.py vars.ini'
else:
    in_allowed = 'Deproject_exe.py vars.ini autoplot'

if in_str != in_allowed:
    raise NameError(colored('\nThe function can only be called as either:','red',attrs=['bold']) +
                    colored('\nDeproject_exe.py vars.ini', 'red') + 
                    colored('\nDeproject_exe.py vars.ini autoplot', 'red'))


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
        

#Depending on which type of sample we want to use, we either read a fits file
#automatically or provide inputs manually
  
try:
    argv[1]
except (IndexError,NameError):
    raise NameError("Couldn't access the input file. Are you sure you wrote it correctly?")
else:
    guessfile = argv[1] #The vars.ini file
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
if datafile=="Distances_PJM2017.csv":
    try:
        os.chdir("DATA/")
    except FileNotFoundError:
        print('Edit your desired path in the script')
    data_raw = Table.read(str(datafile))
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
    print("Sample has " + str(len(dist)) + " stars")



    sample_icrs = coord.ICRS(ra = RA, dec = DEC, pm_ra_cosdec = pm_RA, pm_dec = pm_DEC,distance=dist)

    sample = sample_icrs.transform_to(coord.Galactic)


pvals, rhatvals = calc_p_rhat(sample) 
    
if use_guess:
    mxl, phi_all, fmin_it = max_L(alpha, pvals, rhatvals, vmin, dv, n,v0_guess=v_guess, disp_guess=disp_guess, noniso=non_iso)
elif not use_guess:
    mxl, phi_all, fmin_it = max_L(alpha, pvals, rhatvals, vmin, dv, n, noniso=non_iso)
tf_a = time.time()
endtime = (tf_a - ti_a)/60
print("\nThe run took", endtime, 'm:')


# if logging:
#     # Create a folder for the run and save mxl data
#     os.chdir('/home/mikkola/Documents/DeprojectionProject')
#     folder = make_folder()
#     np.save('RUNS/' + folder + '/mxl_data',mxl)
    
#     # Save output
#     # RUNS folder identifier
#     with open('logs/log_dir_identifier.txt', 'a') as logfile:
#         if folder[10] == 'a' and len(folder) == 11:
#             mark = '='
#             logfile.write('\n' + mark*120 + '\n')
#             logfile.write(mark*55 + folder[:10] + mark*55 + '\n')
#             logfile.write(mark*120 + '\n')
            
#         logfile.write('\nFolder name : ' + folder + '\n')
#         logfile.write('Datafile    : ' + datafile + '\n')
#         logfile.write('fmin its    : ' + str(fmin_it) + '\n')
#         logfile.write('Time needed : ' + str(endtime/60) + ' mins\n')
        
#     # Logfile in RUNS folder 
#     with open('RUNS/' + folder + '/log.txt', 'a') as logfile:
#         logfile.write('Datafile    : ' + datafile + '\n')
#         logfile.write('fmin its    :' + str(fmin_it) + '\n')
#         logfile.write('Time needed : ' + str(endtime/60) + ' mins\n')
#         logfile.write('Labels      : Nbins[1x3], vmin[1x3], bin size, use_guess, noniso, mu_guess, sigma_guess, alpha\n')    
#         value_string=str((str(list(n)).replace(",",":").replace(":",""),str(list(vmin)).replace(",",":").replace(":","")
#                           ,str(list(dv)).replace(",",":").replace(":",""),use_guess,non_iso,str(list(v_guess)).replace(",",":").replace(":","")
#                           , str(list(disp_guess)).replace(",",":").replace(":",""), alpha)).replace("'","")[1:-1]
        
#         logfile.write("Values      : " + value_string + '\n')
        
#     # MORE IS SAVE BY SANITY_CHECK()
        
# try:
#     builtins.autoplot = argv[2]
# except IndexError:
#     argv.append(0)
    
# if argv[2] != 'autoplot':
#     sane = input('Do you want to perform a sanity check [y/n]? ')
#     while sane != 'y' and sane != 'n':   
#         sane = input('Incorrect entry, try again [y/n]! ')
    
#     if sane == 'y':
        
#         from Deproject_test import sanity_check
        
#         sanity_check(pvals,rhatvals,mxl,vmin,dv,n,logging,folder)
    
#     elif sane == 'n':
        
#         print('Suit yourself')
        
#     shouldiplot = input('Do you want to plot your results [y/n]? ')
#     while shouldiplot != 'y' and shouldiplot != 'n':   
#         shouldiplot = input('Incorrect entry, try again [y/n]! ')
    
#     if shouldiplot == 'y':
#         from Deproject_plots import plot_fv,plot_L
        
#         plot_fv(mxl,input('What plane should I project onto? '),vmin,dv,n,logging,folder)
        
#         s=0
        
#         while True:
#             if s==2:
#                 break
#             plotagain = input('Do you want to plot another plane [y/n]? ')
#             if plotagain == 'y':
#                 plot_fv(mxl,input('What plane should I project onto [xy/yz/xz]? '),vmin,dv,n,logging,folder)
#                 s+=1
#                 continue
#             else:
#                 break
    
#     shouldiplotL = input('Do you want to plot the change in L during the maximisation [y/n]? ')
#     while shouldiplotL != 'y' and shouldiplotL != 'n':   
#         shouldiplotL = input('Incorrect entry, try again [y/n]! ')
    
#     if shouldiplotL=='y':
#         plot_L(phi_all,pvals,rhatvals,vmin,dv,n,alpha,logging,folder)
        
# else:
#     sanity_check(pvals,rhatvals,mxl,vmin,dv,n,logging,folder)
#     plot_fv(mxl,'xy',vmin,dv,n,logging,folder)
#     plot_fv(mxl,'yz',vmin,dv,n,logging,folder)
#     plot_fv(mxl,'xz',vmin,dv,n,logging,folder)

# #Should add a way of saving the mxl data        

# print('My work here is done')
