from Deproject_test import * 
from sys import argv 

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
            v0_str = input('Enter values for mux, muy, muz: ').split(',')
            v0 = np.array([int(i) for i in v0_str])
            vdisp_str = input('Enter values for sigmax, sigmay, sigmaz: ').split(',')
            v_disp = np.array([int(i) for i in vdisp_str])
            n_str = input('Enter # of bins in vx, vy, vz: ').split(',')
            n = np.array([int(i) for i in n_str])
            dv_str = input('Enter length of bins in vx, vy, vz: ').split(',')
            dv = np.array([int(i) for i in dv_str])
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
        alpha = float(vars_[6][0])
        break

sample = model_sample(N,v0,v_disp)

#Oort constant values from Bovy (2018)
A = (15.3*(u.km/(u.s*u.kpc))).to(1/u.yr)
B = (-11.9*(u.km/(u.s*u.kpc))).to(1/u.yr)


bvals = sample.b.to(u.deg)
lvals = sample.l.to(u.deg)

mul_obs = sample.pm_l.to(1/u.yr,equivalencies = u.dimensionless_angles())
mub_obs = sample.pm_b.to(1/u.yr,equivalencies = u.dimensionless_angles())

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

test_L('xy', alpha, v0, v_disp, pvals, rhatvals, vmin, dv, n)