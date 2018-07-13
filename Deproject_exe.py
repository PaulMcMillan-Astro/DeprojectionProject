from Deproject_v0 import *
from sys import argv

import time
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

if whatdata == 'data':
    data_raw = ascii.read(input('Enter name of data file: '), format='fast_csv')
    
    RA = data_raw['RAdeg']*u.degree
    DEC = data_raw['DEdeg']*u.degree
    pm_RA = data_raw['pmRA_TGAS']*u.mas/u.yr
    pm_DEC = data_raw['pmDE_TGAS']*u.mas/u.yr
    parallax = data_raw['parallax']*u.mas
    
    dist = parallax.to(u.kpc,equivalencies=u.parallax())
    
    near_stars = np.where(dist.value<0.1)
    
    sample_raw = coord.ICRS(ra = RA, dec = DEC, pm_ra_cosdec = pm_RA, pm_dec = pm_DEC,distance=dist)
    
    sample = sample_raw[near_stars]
    
    sample = sample.transform_to(coord.Galactic)
    
elif whatdata == 'pseudo':
    
    while True:    
        try:
            argv[1]
        except (IndexError,NameError):
            try:
                N = int(input('Enter number of stars: ')) #Number of stars we want to use in our sample
                v0_str = input('Enter values for mux, muy, muz: ').split(',')
                v0 = np.array([int(i) for i in v0_str])
                vdisp_str = input('Enter values for sigmax, sigmay, sigmaz: ').split(',')
                v_disp = np.array([int(i) for i in vdisp_str])
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
            guessfile = argv[1]
            vars_ = process_input(guessfile)
            N = int(vars_[0][0])
            v0 = np.array(vars_[1][0].split(',')[:],dtype=float)
            v_disp = np.array(vars_[2][0].split(',')[:],dtype=float)
            n = np.array(vars_[3][0].split(',')[:],dtype=int)
            vmin = np.array(vars_[4][0].split(',')[:],dtype=float)
            dv = np.array(vars_[5][0].split(',')[:],dtype=float)
            v_guess = np.array(vars_[6][0].split(',')[:],dtype=float)
            disp_guess = np.array(vars_[7][0].split(',')[:],dtype=float)
            alpha = float(vars_[8][0])
            break

    sample = model_sample(N,v0,v_disp)
    
#Oort constant values from Bovy (2018)
A = (15.3*(u.km/(u.s*u.kpc))).to(1/u.yr)
B = (-11.9*(u.km/(u.s*u.kpc))).to(1/u.yr)


bvals = sample.b.to(u.deg)
lvals = sample.l.to(u.deg)

mul_obs = sample.pm_l.to(1/u.yr,equivalencies = u.dimensionless_angles())
mub_obs = sample.pm_b.to(1/u.yr,equivalencies = u.dimensionless_angles())

"""Computation of the relevant quantities

    l,b: Galactic coordinates
    s: the distance obtained by inverting the parallax
    mul, mub: proper motion in l and b
    pvals: Tangential velocities obtained from eq. 2 in DB98
    rhatvals: The unit vector of each star
    vmin: Vector containing the minimum velocities in v-space
    n: The number of cells we want in each dimension of our v-space box
    dv: Step sizes for each dimension"""

b = np.deg2rad(bvals).value # just a test
l = np.deg2rad(lvals).value
cosl = np.cos(l)
cosb = np.cos(b)
sinl = np.sin(l)
sinb = np.sin(b)
s = sample.distance
    
mul = mul_obs - A*np.cos(2*l)-B
mub = mub_obs + A*np.sin(2*l)*cosb*sinb

pvals = s*np.array([-sinl*cosb*mul - cosl*sinb*mub,
                 cosl*cosb*mul - sinl*sinb*mub,
                 cosb*mub])/u.yr
    
rhatvals = np.array([cosb*cosl, cosb*sinl, sinb]).T
pvals = pvals.to(u.km/u.s).value.T

mxl, phi_all = max_L(v_guess, disp_guess,alpha, pvals, rhatvals, vmin, dv, n)

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
    plot_L(phi_all,pvals,rhatvals,vmin,dv,n,alpha)
    
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