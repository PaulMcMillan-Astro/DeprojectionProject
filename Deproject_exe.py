from Deproject_v1_0 import *
from sys import argv
import os
import time

"""Script that when called by Python in terminal will perform the computations needed
to run a maximisation scheme. It will also allow you to plot the results and the 
change in L over all iterations. 

If provided the file 'vars.ini' as an argument such that the terminal command is
'python Deproject_exe.py vars.ini' it will read the input variables automatically"""

##########You will need to change the directory path to your data#################

start_time = time.time()

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
        
#Next we check if there's an input file
#If there is, we proceed, otherwise one can enter pseudo or data
#to manually input the parameters necessary

while True:
    data_list = ['data','pseudo']
    try:
        argv[1]
    except (IndexError,NameError):
        whatdata = input('Enter \"data\" for measurements or \"pseudo\" for Gaussian pseudodata: ')
        if any(whatdata in data_list for i in data_list):
            break
        else:
            raise ValueError('Not a valid input')
    except ValueError:
        print('Please enter \"data\" or pseudo')
        continue
    else:
        whatdata = 'pseudo'
        break

#Depending on which type of sample we want to use, we either read a fits file
#automatically or provide inputs manually

if whatdata == 'data':
    
    data_raw = Table.read(input('Enter name of data file: '))
    
elif whatdata == 'pseudo':
    
    while True:    
        try:
            argv[1]
        except (IndexError,NameError):
            try:
                N = int(input('Enter number of stars: ')) #Number of stars we want to use in our sample
                n_str = input('Enter # of bins in vx, vy, vz: ').split(',')
                n = np.array([int(i) for i in n_str])
                dv_str = input('Enter length of bins in vx, vy, vz: ').split(',')
                dv = np.array([int(i) for i in dv_str])
                vguess_str = input('Enter guess for mux, muy, muz: ').split(',')
                v_guess = np.array([int(i) for i in vguess_str])
                dispguess_str = input('Enter guess for sigmax, sigmay, sigmaz: ').split(',')
                disp_guess = np.array([int(i) for i in dispguess_str])
                alpha = float(input('Enter value of alpha: '))
                vmin = np.array([-200,-200,-200])
                if any(len(i)!=3 for i in [v0,v_disp,v_guess,disp_guess,n]):
                    raise ValueError
            except ValueError:
                print('Not a valid format. Try again')
            else:
                break
        else:
            guessfile = argv[1] #The vars.ini file
            vars_ = process_input(guessfile)
            N = int(vars_[0][0])
            n = np.array(vars_[1][0].split(',')[:],dtype=int)
            vmin = np.array(vars_[2][0].split(',')[:],dtype=float)
            dv = np.array(vars_[3][0].split(',')[:],dtype=float)
            v_guess = np.array(vars_[4][0].split(',')[:],dtype=float)
            disp_guess = np.array(vars_[5][0].split(',')[:],dtype=float)
            alpha = float(vars_[6][0])
            try:
                datafile = vars_[7][0].rstrip('\n')
            except IndexError:
                datafile = []
                pass
            break
        
    if any([len(datafile) != 0,whatdata == 'data']):
        
        try:
            data_raw
        except NameError:
            
            try:
                os.chdir("/home/jooehn/Documents/Summer project/GDR2 data")
            except FileNotFoundError:
                print('Edit your desired path in the script')
            data_raw = Table.read(str(datafile))
            pass
        
        try:
            RA = data_raw['ra']
            DEC = data_raw['dec']
            pm_RA = data_raw['pmra']
            pm_DEC = data_raw['pmdec']
            parallax_raw = data_raw['parallax'].to(u.mas)
        
            G_mean_raw = data_raw['phot_g_mean_mag'].to(u.mag)
            BP_RP_raw = data_raw['bp_rp'].to(u.mag)
            E_BP_RP = data_raw['e_bp_min_rp_val']
            phot_bp_mean_mag = data_raw['phot_bp_mean_mag']
            phot_rp_mean_mag = data_raw['phot_rp_mean_mag']
        
        except KeyError:    
            RA = data_raw['RA']*u.deg
            DEC = data_raw['DEC']*u.deg
            pm_RA = data_raw['PMRA']*u.mas/u.yr
            pm_DEC = data_raw['PMDEC']*u.mas/u.yr
            parallax_raw = data_raw['PARALLAX']*u.mas
        
            G_mean_raw = data_raw['PHOT_G_MEAN_MAG']*u.mag
            BP_RP_raw = data_raw['BP_RP']*u.mag
            E_BP_RP = data_raw['E_BP_MIN_RP_VAL']
            phot_bp_mean_mag = data_raw['PHOT_BP_MEAN_MAG']
            phot_rp_mean_mag = data_raw['PHOT_RP_MEAN_MAG']
            pass
        
        #Properties used for data selection
        
        try:
            parallax_over_error = data_raw['parallax_over_error']
            phot_g_mean_flux_over_error = data_raw['phot_g_mean_flux_over_error']
            phot_rp_mean_flux_over_error = data_raw['phot_rp_mean_flux_over_error']
            phot_bp_mean_flux_over_error = data_raw['phot_bp_mean_flux_over_error']
            phot_bp_rp_excess_factor = data_raw['phot_bp_rp_excess_factor']
            visibility_periods_used = data_raw['visibility_periods_used']
            astrometric_chi2_al = data_raw['astrometric_chi2_al']
            astrometric_n_good_obs_al = data_raw['astrometric_n_good_obs_al']
        except KeyError:
            parallax_over_error = data_raw['PARALLAX_OVER_ERROR']
            phot_g_mean_flux_over_error = data_raw['PHOT_G_MEAN_FLUX_OVER_ERROR']
            phot_rp_mean_flux_over_error = data_raw['PHOT_RP_MEAN_FLUX_OVER_ERROR']
            phot_bp_mean_flux_over_error = data_raw['PHOT_BP_MEAN_FLUX_OVER_ERROR']
            phot_bp_rp_excess_factor = data_raw['PHOT_BP_RP_EXCESS_FACTOR']
            visibility_periods_used = data_raw['VISIBILITY_PERIODS_USED']
            astrometric_chi2_al = data_raw['ASTROMETRIC_CHI2_AL']
            astrometric_n_good_obs_al = data_raw['ASTROMETRIC_N_GOOD_OBS_AL']
        
        parallax_raw = parallax_raw.to(u.mas) #This is done to avoid having to deal with inconsitencies in units of the datasets
        G_mean = G_mean_raw.to(u.mag)
        BP_RP = BP_RP_raw.to(u.mag)
        
        dist = parallax_raw.to(u.kpc,equivalencies=u.parallax())
        
        sample_icrs = coord.ICRS(ra = RA, dec = DEC, pm_ra_cosdec = pm_RA, pm_dec = pm_DEC,distance=dist)
        
        sample_raw = sample_icrs.transform_to(coord.Galactic)
            
        chinu_raw = np.stack((np.exp(-0.4*(G_mean.value-19.5)), np.ones((G_mean.shape))),axis = 0)
        
        chinu = chinu_raw.reshape((2,len(sample_raw)))
    
        """Next we perform the standard cuts which are provided in Appendix B of the Observational HR-Diagram paper of
            Gaia Collaboration et al (2018)"""
    
        standard_cuts = [parallax_over_error<10,visibility_periods_used<8,phot_g_mean_flux_over_error<50,\
                   phot_rp_mean_flux_over_error<20,phot_bp_mean_flux_over_error<20,\
                   phot_bp_rp_excess_factor>1.3+0.06*(phot_bp_mean_mag-phot_rp_mean_mag)**2,\
                   phot_bp_rp_excess_factor<1.0+0.015*(phot_bp_mean_mag-phot_rp_mean_mag)**2,\
                   astrometric_chi2_al/(astrometric_n_good_obs_al-5)>(1.44*np.amax(chinu,axis=0)).reshape((len(sample_raw)))]

        cut_list = standard_cuts
        
        cut_list_rs = [np.reshape(i,(len(sample_raw))) for i in cut_list]
        bad_idx_list = [np.argwhere(i) for i in cut_list_rs]
        bad_idx_arr = np.concatenate(bad_idx_list)
        
        bad_idx_arr1 = bad_idx_arr.reshape((len(bad_idx_arr)))
    
        nan_idx = np.concatenate([np.ravel(np.argwhere(np.isnan(G_mean_raw))),np.ravel(np.argwhere(np.isnan(BP_RP_raw)))]) 
        nanvals = np.unique(nan_idx)
    
        bad_idx_arr2 = np.concatenate([bad_idx_arr1,nanvals])
    
        bad_idx = np.unique(bad_idx_arr2)
    
        mask = np.full((sample_raw.shape),True)
        mask[bad_idx] = False
    
        sample = sample_raw[mask]
        G_mean_new = G_mean[mask]
        BP_RP_new = BP_RP[mask]
    else:
        sample = model_sample(N)
    
pvals, rhatvals = calc_p_rhat(sample)

mxl, phi_all = max_L(alpha, pvals, rhatvals, vmin, dv, n,v0_guess=v_guess, disp_guess=disp_guess)

print("The run took", time.time() - start_time, 's')

sane = input('Do you want to perform a sanity check [y/n]? ')

if sane == 'y':
    
    from Deproject_test import sanity_check
    
    print(sanity_check(pvals,rhatvals,mxl,vmin,dv,n))

elif sane == 'n':
    
    print('Suit yourself')
    
shouldiplot = input('Do you want to plot your results?[y/n] ')

if shouldiplot == 'y':
    from Deproject_plots import *
    
    plot_fv(mxl,input('What plane should I project onto? '),vmin,dv,n)
    
    s=0
    
    while True:
        if s==2:
            break
        plotagain = input('Do you want to plot another plane?[y/n] ')
        if plotagain == 'y':
            plot_fv(mxl,input('What plane should I project onto? '),vmin,dv,n)
            s+=1
            continue
        else:
            break
    
    shouldiplotL = input('Do you want to plot the change in L during the maximisation?[y/n] ')
    if shouldiplotL=='y':
        plot_L(phi_all,pvals,rhatvals,vmin,dv,n,alpha)
        

print('My work here is done')