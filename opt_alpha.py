from Deproject_v1_0 import *
from sys import argv
import time 
import os

"""This script will automatically initiate an optimisation of the alpha
smoothing parameter. Given an input file, i.e. 'alpha_vars.ini', it will 
load all parameters for the dimensions of the box which represents v-space.

If an input file is not provided, one can input these values manually via
inputs."""

start_time = time.time()

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


while True:    
    try:
        argv[1]
    except (IndexError,NameError):
        try:
            N = int(input('Enter number of stars: ')) #Number of stars we want to use in our sample
            n_str = input('Enter # of bins in vx, vy, vz: ').split(',')
            n = np.array([int(i) for i in n_str]) #Dimensions of box
            dv_str = input('Enter length of bins in vx, vy, vz: ').split(',')
            dv = np.array([int(i) for i in dv_str]) #The box widths in x,y,z
            alpha0 = float(input('Enter value of alpha: ')) #Initial guess of alpha
            vmin = np.array([-200,-200,-200]) #Anchor of the box
            if any(len(i)!=3 for i in [vmin,dv,n]):
                raise ValueError
        except ValueError:
            print('Not a valid format. Try again')
        else:
            break
    else:
        guessfile = argv[1]
        vars_ = process_input(guessfile)
        N = int(vars_[0][0])
        M = int(vars_[1][0])
        n = np.array(vars_[2][0].split(',')[:],dtype=int)
        vmin = np.array(vars_[3][0].split(',')[:],dtype=float)
        dv = np.array(vars_[4][0].split(',')[:],dtype=float)
        alpha0 = float(vars_[5][0])
        datafile = vars_[6][0].rstrip('\n')
        break
    
if len(datafile)!=0:

    try:
        os.chdir("/home/jooehn/Documents/Summer project/GDR2 data")
    except FileNotFoundError:
        print('Edit your desired path in the python file')
    data_raw = Table.read(str(datafile))
    
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

opt_alpha(alpha0,M,N,sample,vmin,dv,n,tol=0.01)

print("The run took", time.time() - start_time, 's')