import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import astropy.units as u
import astropy.coordinates as coord
import scipy.stats as st
import os
import inspect
import string
import time
import sys
import psutil
import builtins
from termcolor import colored
from astropy.table import Table
from astropy.io import ascii
from scipy import sparse as scisp
from scipy.optimize import fmin_cg
from scipy.interpolate import interpn
from scipy.ndimage import zoom
from datetime import date
from decimal import Decimal
from alpha_debugger import *
from IPython.display import clear_output
from types import ModuleType, FunctionType
from gc import get_referents

def mod_dict(vdict):
    for var in vdict.keys():
        vdict[var] = getsize(vdict[var])
    return vdict

def getsize(obj):
    BLACKLIST = type, ModuleType, FunctionType
    """sum size of object & members."""
    if isinstance(obj, BLACKLIST):
        return 0
    seen_ids = set()
    size = 0
    objects = [obj]
    while objects:
        need_referents = []
        for obj in objects:
            if not isinstance(obj, BLACKLIST) and id(obj) not in seen_ids:
                seen_ids.add(id(obj))
                size += sys.getsizeof(obj)
                need_referents.append(obj)
        objects = get_referents(*need_referents)
    return size


#@profile
def multigrid_steps(n):
    '''This function determines how many multigrids steps can be performed and returns the box size
    in each step. In addition, if the starting box dimensions are not divisible by 2**steps the final box dimensions
    are changed to accomodate this.'''
    step = 5
    while any(np.round(n/(2**step)) < 10):
        step -= 1
    box_steps = (2**np.linspace(0,step,step+1).reshape(step+1,1) * np.round(n/(2**step))).astype(int)

    if not all(box_steps[-1] == n):
        print('To oct-split box the recommended %s times, dimensions were changed from (%s, %s, %s) to (%s, %s, %s)\n'
              % (len(box_steps),n[0],n[1],n[2],box_steps[-1,0],box_steps[-1,1],box_steps[-1,2]))

    return box_steps.astype(int)

#@profile

builtins.variables = dict()
def update_variables(loc, glob):
    loc_vars = list(loc.keys())
    glob_vars = list(glob.keys())
    all_vars = list(set(loc_vars) | set(glob_vars))
    
    for var in all_vars:
        try:
            builtins.variables[var] = sys.getsizeof(glob[var])/1e6
        except KeyError:
            builtins.variables[var] = sys.getsizeof(loc[var])/1e6
            
    variable_list = sorted([(name, value) for (name, value) in builtins.variables.items()],
                          key=lambda x: x[1],
                          reverse=True)
    os.system('clear')
    print('-'*60)
    [print(var, ':', size, 'MiB') for var, size in variable_list[:3]]
    return
    

def write_variables():
    cwd = os.getcwd()
    os.chdir('/home/daniel/DeprojectionProject')
    variables = list(builtins.lv.items())
    l_vars = []
    for var, obj in variables:
        l_vars += [(var.ljust(30), obj/1e9)]
    l_vars = sorted(l_vars, key=lambda mem: mem[1], reverse=True)
    with open('variable_memory.log', 'w+') as vlog:
        vlog.write('Variable name'.ljust(30) + 'Size'.ljust(15) + 'GB\n')
        for line in l_vars:
            vlog.write(line[0] + str(line[1]).ljust(15) + 'GB\n')

    os.chdir(cwd)
    return


#@profile
def zoomed_mxl(mxl):
    '''Takes a mxl result and upscales it with interpolation'''
    phi0_guess = zoom(mxl, zoom=2, order=3) - np.log(8)

    return phi0_guess

#@profile
def callback_mg(x):
    i = builtins.grid_step
    builtins.L[i].append(builtins.current_L)
    builtins.gradL[i].append(builtins.current_gradL)

#@profile
def callback(x):
    builtins.L.append(builtins.current_L)
    builtins.gradL.append(builtins.current_gradL)

#@profile
def make_treefile():
    '''Function that creates a .txt file with name YYYY-MM-DDx where x is the
    first letter that is a non-existing file.
    The file will contain the log(alpha) and MISE values used in each
    opt iteration. This can be used to recreate the tree search path.'''

    alp = string.ascii_lowercase
    date_str = str(date.today())

    list_str = []
    for letter in alp:
        list_str.append(date_str + letter)

    for letter1 in alp:
        for letter2 in alp:
            list_str.append(date_str + letter1 + letter2)

    list_str = np.array(list_str)

    os_dirs = os.popen('ls logs | grep "%s"' % date_str).read().split("\n")[:-1]
    os_dirs.sort(key=lambda x: len(x))
    os_dirs = np.array(os_dirs)

    existing_dirs = np.zeros(len(list_str), dtype='<U12')
    existing_dirs[:len(os_dirs)] = os_dirs

    file_name = 'logs/tree_' + list_str[list_str != existing_dirs][0] + '.txt'


    return file_name

#
#@profile
def calc_p_rhat(sample):
    """Function that computes the values of rhat and the tangential velocity for a given sample

       Computation of the relevant quantities

       l,b: Galactic coordinates
       s: the distance obtained by inverting the parallax
       mul, mub: proper motion in l and b
       pvals: Tangential velocities obtained from eq. 2 in DB98
       rhatvals: The unit vector of each star
       vmin: Vector containing the minimum velocities in v-space
       n: The number of cells we want in each dimension of our v-space box
       dv: Step sizes for each dimension"""
    # Oort constant values from Bovy (2018)
    A = (15.3 * (u.km / (u.s * u.kpc)))
    B = (-11.9 * (u.km / (u.s * u.kpc)))

    if isinstance(sample, Table):
        b = sample['b']
        l = sample['l']
        mul_obs = sample['pml_cosb'].to(u.km / (u.s * u.kpc), equivalencies=u.dimensionless_angles())
        mub_obs = sample['pmb'].to(u.km / (u.s * u.kpc), equivalencies=u.dimensionless_angles())
        p = sample['parallax'].to(u.kpc, equivalencies=u.parallax())

    elif isinstance(sample, coord.builtin_frames.galactic.Galactic):
        b = sample.b
        l = sample.l
        mul_obs = sample.pm_l_cosb.to(u.km / (u.s * u.kpc), equivalencies=u.dimensionless_angles())
        mub_obs = sample.pm_b.to(u.km / (u.s * u.kpc), equivalencies=u.dimensionless_angles())
        p = sample.distance.to(u.kpc)

    cosl = np.cos(l)
    cosb = np.cos(b)
    sinl = np.sin(l)
    sinb = np.sin(b)

    mul = mul_obs - A * np.cos(2 * l) - B
    mub = mub_obs + A * np.sin(2 * l) * cosb * sinb

    pvals = p * np.vstack((-sinl*mul - cosl*sinb*mub,
                           cosl*mul - sinl*sinb*mub,
                           cosb*mub))
    pvals = pvals.T

    rhatvals = np.array([cosb * cosl, cosb * sinl, sinb]).T

    return pvals, rhatvals


#
#@profile
def model_sample(N):
    """Generates a simple model solar neighbourhood star sample in a Galactic frame of reference assuming a
    velocity distribution that is a sum of three predefined Gaussians.

    If one wishes to use other mean velocities, change these in mu0, mu1 and mu2, while the dispersions are changed
    in disp0, disp1 and disp2. The relative weights, w1, w2, w3 determine how many of the stars that belong to
    each Gaussian. As for now, they're just set to be ~1/3.

    The stars have distances from the Sun that are in the range [10,100] pc.

    Takes the following arguments:

    N: Number of stars in the sample"""

    xmax, ymax, zmax = np.array([100,100,100])/np.sqrt(3)
    xmin, ymin, zmin = -xmax,-ymax,-zmax
    #w0 = 1.
    #w1 = 0.
    w0 = w1 = 0.33
    w2 = 1-(w0+w1)

    #mu0 = np.array([0,0,0])
    mu0 = np.array([30,30,30])
    mu1 = np.array([-20,-20,-20])
    mu2 = np.array([15,-15,15])

    disp0 = np.array([25,23,27])
    disp1 = np.array([13,17,15])
    disp2 = np.array([9,14,12])

    w_sample = np.random.random_sample(N)

    psx = (np.random.rand(N) * (xmax - xmin) + xmin) * u.pc
    psy = (np.random.rand(N) * (ymax - ymin) + ymin) * u.pc
    psz = (np.random.rand(N) * (zmax - zmin) + zmin) * u.pc

    from_g0 = np.where(w_sample < w0)  # We get the indices of the stars that belong to the first Gaussian
    from_g1 = np.where((w0 <= w_sample) & (w_sample < (w0 + w1)))
    from_g2 = np.where(w_sample > (1 - w2))

    scale0 = np.random.randn(len(from_g0[0]), 3)
    scale1 = np.random.randn(len(from_g1[0]), 3)
    scale2 = np.random.randn(len(from_g2[0]), 3)

    psvels = np.zeros((N, 3))

    psvels[
        from_g0] = mu0 + scale0 * disp0  # We exchange our empty velocity values with the ones obtained from each Gaussian
    psvels[from_g1] = mu1 + scale1 * disp1
    psvels[from_g2] = mu2 + scale2 * disp2

    psvx, psvy, psvz = psvels.T
    # We use Astropy's coord class which makes it easy to keep track of units and conversions

    psample = coord.Galactic(u=psx, v=psy, w=psz, U=psvx * (u.km / u.s),
                             V=psvy * (u.km / u.s), W=psvz * (u.km / u.s),
                             representation_type=coord.CartesianRepresentation,
                             differential_type=coord.CartesianDifferential)

    psample.set_representation_cls(coord.SphericalRepresentation, coord.SphericalCosLatDifferential)


    return psample, psvx, psvy, psvz


#
#@profile
def calc_K(pk, rhat, vmin, dv, n):
    '''Calculate the values of K simultaneously for all bins for a given star with tangential velocity
    vector pk and a unit vector rhat'''

    vxmin, vymin, vzmin = vmin
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    pkx, pky, pkz = pk
    rhatx, rhaty, rhatz = rhat

    vxmax, vymax, vzmax = vxmin + nx * dvx, vymin + ny * dvy, vzmin + nz * dvz

    """We now solve the line equation v = pk + vr*rhat, where v are the intersections of the tangential velocity line
    and the boundaries of each bin"""

    vx_bins = np.linspace(vxmin, vxmax, nx + 1)
    vy_bins = np.linspace(vymin, vymax, ny + 1)
    vz_bins = np.linspace(vzmin, vzmax, nz + 1)

    vrx = (vx_bins - pkx) / rhatx
    vry = (vy_bins - pky) / rhaty
    vrz = (vz_bins - pkz) / rhatz

    """After solving the line equation for each dim we remove the vr values which solve the equation
    for bins outside of our specified box."""

    vrmax = min(max(vrx), max(vry), max(vrz))
    vrmin = max(min(vrx), min(vry), min(vrz))

    vrx = vrx[(vrx <= vrmax) & (vrx >= vrmin)]
    vry = vry[(vry <= vrmax) & (vry >= vrmin)]
    vrz = vrz[(vrz <= vrmax) & (vrz >= vrmin)]

    vr = np.concatenate((vrx, vry, vrz))
    vr.sort()  # We obtain an ordered list with all vr values for intersections between entry and exit points

    K = np.zeros((nx, ny, nz))

    if len(vr) == 0:

        return K

    vr_prime = (vr[:-1] + vr[1:]) / 2
    pks = np.ones((len(vr_prime), len(pk))) * pk[np.newaxis]
    rhats = np.ones((len(vr_prime), len(rhat))) * rhat[np.newaxis]
    vmins = np.ones((len(vr_prime), len(vmin))) * vmin[np.newaxis]
    vr_primestack = np.ones((len(vr_prime), 3))*vr_prime[:, np.newaxis]

    """We now solve the line equation again for values in the middle of each bin with a line segment in it.
    This gives us the coordinates for each relevant bin, given in line_bins.
    Finally we compute the length of each segment and add said value to the relevant box in our K-space."""
    
    v_prime = pks + vr_primestack * rhats

    line_bins = np.floor((v_prime - vmins) / dv)
    line_bins = line_bins.astype(int)

    line_len = vr[1:] - vr[:-1]  # Gets the length of each line
    non_zero = np.nonzero(line_len)
    line_len = line_len[non_zero]  # Removes duplicate intersections
    line_bins = line_bins[non_zero]

    K[line_bins[:, 0], line_bins[:, 1], line_bins[:, 2]] = line_len / (dvx * dvy * dvz)

    return K


#
#@profile
def calc_sigma2(pvals, rhat, give_vmean=False, noniso=False):
    """Function that applies equation 12 of DB98 for a set of stars from their tangential velocities and unit vectors.
    Returns the velocity dispersion tensor.

    pvals: array of the pk vectors for all N stars in our sample. Should have dims (N,3)
    rhat: array of N unit vectors, one for each star
    give_vmean: if True, returns the computed mean velocity vector value for the given sample
    noniso: if True, we no longer assume that the sample is isotropic"""

    # pmean = pvals.mean(axis=0)
    pmean = pvals.mean(axis=-2)

    # Fast way for outer product eq (4) of DB98
    rhat_outer = np.expand_dims(rhat, axis=-1) * np.expand_dims(rhat, axis=-2)

    # Create N_star identity matrices
    A0 = np.identity(3) * np.expand_dims(np.ones(rhat.shape[:-1]), axis=(-1, -2))
    A = A0 - rhat_outer
    del rhat_outer


    A_mean_inv = np.linalg.inv(A.mean(axis=-3))
    v_mean = np.einsum('...ij,...j->...i', A_mean_inv, pmean)
    del A_mean_inv

    # Computes p' from equation (6) in DB98
    pp = pvals - np.einsum('...ikj,...j->...ik', A, v_mean)

    pp2mean = np.mean(np.square(pp), axis=-2)
    if noniso:
        # In our main method, we rely on built-in tensor algebra functions
        # that perform these computations for us. We set up B by computing the
        # mean of the outer product of the peculiar tangential velocities
        # p' with itself.

        B = np.mean(np.expand_dims(pp, axis=-1) * np.expand_dims(pp, axis=-2), axis=-3)
        del pp
        # Next, we use the tensordot function of numpy to obtain T for each star.
        # We resort to a simple list comprehension operation to obtain these
        # values

        # We could alternatively use np.einsum here
        T = np.asarray(np.einsum('...ij,...kl->...ijlk', A, A)).mean(axis=-5)
        del A

        # With the non-singular A tensor at hand we can easily solve for D
        # and obtain the velocity dispersion tensor

        if T.ndim == 5:

            D = np.array([np.linalg.tensorsolve(t, b, (0, 2)) for t, b in zip(T, B)])
            del T, B
        else:

            D = np.linalg.tensorsolve(T, B, (0, 2))
            del T, B
        sigma2 = np.diagonal(D, axis1=-2, axis2=-1)

    else:
        del pp

        B = np.array([[9, -1, -1], [-1, 9, -1], [-1, -1, 9]])

        sigma2 = (3 / 14) * np.einsum('ij,...j->...i', B, pp2mean)  # The velocity dispersion tensor

    if give_vmean == True:

        return sigma2, v_mean

    else:

        return sigma2

#
#@profile
def sec_der(phi, sigma2, dv):
    '''This function calculates eq. 29 (and in turn eq. 30) of WD98. Each cell m has a contribution
    of -2*phi_m, phi_(m+e_i), and phi_(m-e_i) from itself, its neighbor one step up in the i
    direction and one step down in the i direction respectively. The direction i loops over
    the x, y, and z directions.

    The contributions are added onto a zero array of the same shape as phi and only when
    both neighbours are available in the given dimension. e.g. there is no x term when
    treating on [1, :, :] or [:, :, -1]

    Inputs are:

    phi - log of the distribution [nx,ny,nz]

    sigma2 - squared dispersion in each direction [array of len 3]

    dv - cell widths [array of len 3]

    ----

    Output is:

    phi_arr - The contributions to each cell, equiv. to eq 29 in WD98'''

    phi_arr = np.zeros(phi.shape)

    phi_arr[1:-1,:,:] += (sigma2[0] / (dv[0]*dv[0])) * (-2*phi[1:-1, :, :] + phi[2:, :, :] + phi[:-2, :, :])
    phi_arr[:,1:-1,:] += (sigma2[1] / (dv[1]*dv[1])) * (-2*phi[:, 1:-1, :] + phi[:, 2:, :] + phi[:, :-2, :])
    phi_arr[:,:,1:-1] += (sigma2[2] / (dv[2]*dv[2])) * (-2*phi[:, :, 1:-1] + phi[:, :, 2:] + phi[:, :, :-2])


    return phi_arr


#
#@profile
def grad_sec_der(phi, sigma2, dv):
    '''Here we calculate the equivalent factor to sec_der for the gradients third term'''

    phi_arr = np.zeros(phi.shape)

    phi_arr_loc = sec_der(phi, sigma2, dv) # This gives us the matrix A_m for all m = (i,j,k) cells

    # The x contribution
    phi_arr[:-2,:,:]  += 2 * (phi_arr_loc[1:-1,:,:]) * sigma2[0]/(dv[0]*dv[0]) # Adds A_(m-1) contribution
    phi_arr[2:,:,:]   += 2 * (phi_arr_loc[1:-1,:,:]) * sigma2[0]/(dv[0]*dv[0]) # Adds A_(m+1) contribution
    phi_arr[1:-1,:,:] -= 4 * (phi_arr_loc[1:-1,:,:]) * sigma2[0]/(dv[0]*dv[0]) # Adds A_m contribution

    # The y contribution
    phi_arr[:,:-2,:]  += 2 * (phi_arr_loc[:,1:-1,:]) * sigma2[1]/(dv[1]*dv[1]) # Adds A_(m-1) contribution
    phi_arr[:,2:,:]   += 2 * (phi_arr_loc[:,1:-1,:]) * sigma2[1]/(dv[1]*dv[1]) # Adds A_(m+1) contribution
    phi_arr[:,1:-1,:] -= 4 * (phi_arr_loc[:,1:-1,:]) * sigma2[1]/(dv[1]*dv[1]) # Adds A_m contribution

    # The z contribution
    phi_arr[:,:,:-2]  += 2 * (phi_arr_loc[:,:,1:-1]) * sigma2[2]/(dv[2]*dv[2]) # Adds A_(m-1) contribution
    phi_arr[:,:,2:]   += 2 * (phi_arr_loc[:,:,1:-1]) * sigma2[2]/(dv[2]*dv[2]) # Adds A_(m+1) contribution
    phi_arr[:,:,1:-1] -= 4 * (phi_arr_loc[:,:,1:-1]) * sigma2[2]/(dv[2]*dv[2]) # Adds A_m contribution

    return phi_arr

#
#@profile
def phi_guess(v0, disp0, vmin, dv, n):
    """Provides an initial guess of the phi values in each bin given an assumed distribution f(v).
    For now only allows for a sum of two Gaussians type guess given arrays with mean velocities and
    dispersions for each dimension.

    v0: A vector containing all the mean velocity values for the Gaussian distribution

    disp0: A vector with the velocity dispersions for the Gaussian in x, y, z

    vmin: The anchor point of our velocity space box

    """

    vxmin, vymin, vzmin = vmin
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    v0x, v0y, v0z = v0
    dispA = disp0
    dispB = disp0*2

    vxmax, vymax, vzmax = vxmin + nx * dvx, vymin + ny * dvy, vzmin + nz * dvz

    vx_bins = np.linspace(vxmin, vxmax, nx+1)
    vy_bins = np.linspace(vymin, vymax, ny+1)
    vz_bins = np.linspace(vzmin, vzmax, nz+1)

    vxc = (vx_bins[1:] + vx_bins[:-1]) / 2
    vyc = (vy_bins[1:] + vy_bins[:-1]) / 2
    vzc = (vz_bins[1:] + vz_bins[:-1]) / 2

    """Given the velocities of each bin we compute the 3D Gaussian value."""

    pos = np.stack(np.meshgrid(vxc,vyc,vzc),axis=3)

    fA = st.multivariate_normal(mean=v0, cov=np.diag(dispA**2))
    fB = st.multivariate_normal(mean=v0, cov=np.diag(dispB**2))

    phi = np.log((fA.pdf(pos) + fB.pdf(pos))/2)

    return phi


#
#@profile
def get_neg_L(phi, *args):
    """The function that we wish to optimise. Corresponds to eq. 31 in D98.

    N: Number of stars in our sample

    Kvals: Array of dimensions (N,nx,ny,nz) containing the K-values for each star in our sample.

    alpha: Smoothing parameter that can be found using the function opt_alpha

    We use sparse matrices in our computation of L_tilde because our K- arrays have low
    density. Therefore we use the scipy.sparse package to convert our arrays to sparse arrays.
    Using coo-matrices is faster when building the matrix, but the csc-matrices are faster for
    arithmetic operations

    phi_unr regains the original shape of our phi guess and is used to compute the
    third term of L_tilde."""

    Kvals, N, alpha, dv, n, sigma2 = args

    exphir = np.exp(phi)
    exphi_csc = scisp.coo_matrix(exphir).tocsc()
    Kphi = Kvals.multiply(exphi_csc).sum(axis=1) # Order all Kphi values in 1D arrays and compute the sum of exp(phi)*K(k|l) for each star
    Kphi_sum_tot = np.log(Kphi[Kphi != 0]).sum() # To make sure we don't get infinities and .sum() gives the double sum in the first term

    phi_unr = np.reshape(phi, n)
    phixhi_sum = (sec_der(phi_unr, sigma2, dv) ** 2).sum()

    t1 = Kphi_sum_tot / N
    t2 = exphi_csc.sum()
    t3 = ((alpha * dv[0] * dv[1] * dv[2]) / 2) * phixhi_sum

    L_tilde = t1 - t2 - t3 # eq. 31 in DB98

    neg_L = -1 * L_tilde  # Since we want to maximize L_tilde, we should minimize -L_tilde
    builtins.current_L = L_tilde

    return neg_L


#
#@profile
def get_grad_neg_L(phi, *args):
    """In this function we compute the gradient of L. We compute the derivative for each cell and return a
    1D array of length (nx*ny*nz).

    args: see get_L

    """
    Kvals, N, alpha, dv, n, sigma2 = args


    exphir = np.exp(phi)
    exphi_csc = scisp.coo_matrix(exphir).tocsc()

    Kphi = exphi_csc.multiply(Kvals)
    Kphi_sum = Kphi.sum(axis=1)
    Kphi_sum[Kphi_sum.nonzero()] = 1 / Kphi_sum[Kphi_sum.nonzero()] # We compute the sum of exp(phi)*K(k|l) for each star
    K_term0 = Kphi.multiply(Kphi_sum)
    K_term = K_term0.sum(axis=0) # The final array with the first term for each cell

    phi_unr = np.reshape(phi,n)
    dphixhi = grad_sec_der(phi_unr, sigma2, dv)

    t1 = K_term/N
    t2 = exphi_csc
    t3 = ((alpha * dv[0] * dv[1] * dv[2]) / 2) * dphixhi.ravel()

    grad_L = np.asarray(t1 - t2 - t3).reshape(len(phi), )

    neg_grad_L = -1 * grad_L
    builtins.current_gradL = np.linalg.norm(neg_grad_L)

    return neg_grad_L


#
#@profile
def fmin_cg_output(fopt, fcalls, gcalls, flag, fmin_it):
    l1 = ['Optimization terminated sucessfully.\n',
          colored('Warning','red',attrs=['bold'])+': Maximum number of iterations has been exceeded.\n',
          colored('Warning','red',attrs=['bold'])+': Desired error not necessarily achieved due to precision loss.\n',
          colored('Warning','red',attrs=['bold'])+': NaN result encountered.\n']
    l2 = ('         Current function value : %f\n' % fopt)
    l3 = ('         Iterations             : %s\n' % fmin_it)
    l4 = ('         Function evaluations   : %s\n' % fcalls)
    l5 = ('         Gradient evaluations   : %s\n' % gcalls)
    print(l1[flag] + l2 + l3 + l4 + l5)

    return


#@profile
def max_L(alpha, pvals, rhatvals, vmin, dv, n, phi0_guess=[], v0_guess=[], disp_guess=[], noniso=False, printing=True):
    """Function that employs scipy.optimize.fmin_cg to maximise the function get_neg_L().
    It takes guesses of the distribution (currently only supports Gaussian guesses) and the relevant data from the
    star sample for which the velocity distribution is to be estimated."""

    builtins.L     = []
    builtins.gradL = []

    dvx, dvy, dvz = dv
    nx, ny, nz = n
    N = len(pvals)

    sigma2, vmean = calc_sigma2(pvals, rhatvals, True, noniso=noniso)
    if np.size(v0_guess) == 0:
        v0_guess = vmean
    if np.size(disp_guess) == 0:
        sigma = np.sqrt(sigma2)
        disp_guess = sigma
    if np.size(phi0_guess) == 0:
        phi0 = phi_guess(v0_guess, disp_guess, vmin, dv,n)
    else:
        phi0 = phi0_guess

    Kvals_args = (pvals, rhatvals, vmin, dv, n, N, printing)
    Kvals = Kvals_function_selector(Kvals_args)

    args = (Kvals, N, alpha, dv, n, sigma2)
    phi0r = np.ravel(phi0)  # fmin_cg only takes one-dimensional inputs for the initial guess

    print('Started fmin_cg... ',end='')

    builtins.L.append(-1*get_neg_L(phi0r,Kvals, N, alpha, dv, n, sigma2))
    builtins.gradL.append(np.linalg.norm(get_grad_neg_L(phi0r,Kvals, N, alpha, dv, n, sigma2)))

    mxl, fopt, fcalls, gcalls, flag, phi_all = fmin_cg(get_neg_L, phi0r, fprime=get_grad_neg_L, gtol=1e-6, args=args, retall=True, disp=False, full_output=True, callback=callback)
    print(colored('Finished!','green',attrs=['bold','underline']))

    fmin_it = np.shape(phi_all)[0] - 1
    fmin_cg_output(fopt, fcalls, gcalls, flag, fmin_it)

    builtins.n = n
    builtins.dv = dv

    mxlnew = mxl.reshape(n)

    return mxlnew, fmin_it


#@profile
def multigrid_max_L(alpha, pvals, rhatvals, vmin, dv, n, phi0_guess=[], v0_guess=[], disp_guess=[], noniso=False, printing=False):
    """Function that employs scipy.optimize.fmin_cg to maximise the function get_neg_L().
    It takes guesses of the distribution (currently only supports Gaussian guesses) and the relevant data from the
    star sample for which the velocity distribution is to be estimated.

    This function differs from max_L() in that it starts with a crude box (nx,ny,nz) and itertively increases size
    to reach the final box size. This will significantly improve runtime"""

    N = len(pvals)
    vmax = vmin + dv*n

    box_steps = multigrid_steps(n)

    n = box_steps[0]
    dv = (vmax-vmin)/n

    sigma2, vmean = calc_sigma2(pvals, rhatvals, True, noniso=noniso)
    if np.size(v0_guess) == 0:
        v0_guess = vmean
    if np.size(disp_guess) == 0:
        sigma = np.sqrt(sigma2)
        disp_guess = sigma
    if np.size(phi0_guess) == 0:
        phi0 = phi_guess(v0_guess, disp_guess, vmin, dv,n)
    else:
        phi0 = phi0_guess

    builtins.L     = [[] for _ in range(len(box_steps))]
    builtins.gradL = [[] for _ in range(len(box_steps))]
    fmin_it = 0
    for grid_step, n in enumerate(box_steps):
        if grid_step == len(box_steps)-1:
            printing = True
        builtins.grid_step = grid_step
        dv = (vmax-vmin)/n

        Kvals_args = (pvals, rhatvals, vmin, dv, n, N, printing)
        Kvals = Kvals_function_selector(Kvals_args)

        args = (Kvals, N, alpha, dv, n, sigma2)

        phi0r = np.ravel(phi0)  # fmin_cg only takes one-dimensional inputs for the initial guess

        ### This is where the minimization occurs
        print('Started fmin_cg on (%s, %s, %s) grid... ' % (n[0],n[1],n[2]),end='', flush=True)
        builtins.L[grid_step].append(-1*get_neg_L(phi0r,Kvals, N, alpha, dv, n, sigma2))
        builtins.gradL[grid_step].append(np.linalg.norm(get_grad_neg_L(phi0r,Kvals, N, alpha, dv, n, sigma2)))

        mxl, fopt, fcalls, gcalls, flag, phi_all = fmin_cg(get_neg_L, phi0r, fprime=get_grad_neg_L, gtol=1e-6, args=args, retall=True, disp=False, full_output=True, callback=callback_mg)
        fmin_it += np.shape(phi_all)[0] - 1

        print(colored('Finished!','green',attrs=['bold','underline']),end='')
        print(' fopt : %s' % fopt)
        ##########################################

        if grid_step == len(box_steps)-1:
            fmin_cg_output(fopt, fcalls, gcalls, flag, fmin_it)
        else:
            phi0 = zoomed_mxl(mxl.reshape(n))

    builtins.n = n
    builtins.dv = dv

    mxlnew = mxl.reshape(n)

    return mxlnew, fmin_it


#@profile
def Kvals_function_selector(args):
    """This function determines if Kvals can be calculated using a faster numpy pre-allocation method or if it
    requires a block-step method based on the systems available RAM. The block-step method calculates Kvals for
    a block on Nblock stars before making those Kvals sparse, it then proceeds with the next."""

    pvals, rhatvals, vmin, dv, n, N, printing = args

    # We allow only 80% of the available RAM to be used to allow other processes to run simulataneously.
    AvMem = psutil.virtual_memory().available*0.8

    MaxFloats = AvMem / 8  # Number of floats we can handle assuming 8 bytes per float
    Nblock = int(np.floor(MaxFloats / np.prod(n)))  # Largest possible block size
    MemReq = 8*N*np.prod(n)/1e9
    AvMem /= 1e9
    if printing == True:
        print('Allocated RAM: %.2f GB  \nRequired RAM : %.2f GB | Block size = %s' % (AvMem, MemReq, Nblock))

    if Nblock > N:
        if printing == True:
            print('Fast Numpy Kvals run possible, running...\n')
            Kvals = KvalsNumpyMethod(pvals, rhatvals, vmin, dv, n)
            print('Finished kvals.\n')
        elif printing == False:
            Kvals = KvalsNumpyMethod(pvals, rhatvals, vmin, dv, n)
    else:
        if printing == True:
            print('Slightly slower Block Kvals run needed, running...\n')
            Kvals = KvalsBlockMethod(pvals, rhatvals, vmin, dv, n, Nblock)
            print('Finished kvals.\n')
        elif printing == False:
            Kvals = KvalsBlockMethod(pvals, rhatvals, vmin, dv, n, Nblock)

    return Kvals


#
#@profile
def KvalsNumpyMethod(pvals, rhatvals, vmin, dv, n):
    nx, ny, nz = n
    N = len(pvals)

    K0 = np.ravel(calc_K(pvals[0], rhatvals[0], vmin, dv, n))


    Kvals = np.zeros((N, len(K0)))
    Kvals[0] = K0

    for i in range(1, N):
        Kvals[i] = np.ravel(calc_K(pvals[i], rhatvals[i], vmin, dv, n))

    Kvals = scisp.csc_matrix(Kvals)

    return Kvals


#
#@profile
def KvalsBlockMethod(pvals, rhat, vmin, dv, n, nblock):
    n_stars = len(pvals)
    kvals_block = np.zeros((nblock, np.prod(n)))
    kvals_block_shape = kvals_block.shape
    last_block = (n_stars % nblock)

    print("Kvals_block size: ", getsize(kvals_block)/1e9, "GB")
    print("Available RAM: ", psutil.virtual_memory().available/1e9, "GB")

    # Run kvals for for all the blocks and stacking them sparsely after n_stars.
    kvals = None
    for i in range(n_stars - last_block):
        kvals_block[i % nblock] = np.ravel(calc_K(pvals[i], rhat[i], vmin, dv, n))

        if (i + 1) % nblock == 0:
            kvals_coo = scisp.coo_matrix(kvals_block)
            kvals = scisp.vstack((kvals, kvals_coo))
            del kvals_block
            kvals_block = np.zeros(kvals_block_shape)

    kvals_block = np.zeros((last_block, np.prod(n)))
    for i in range(last_block):
        kvals_block[i] = np.ravel(calc_K(pvals[i], rhat[i], vmin, dv, n))

    # Performing the final stack
    kvals_coo = scisp.coo_matrix(kvals_block)
    kvals = scisp.vstack((kvals, kvals_coo))
    kvals = scisp.csc_matrix(kvals)

    return kvals


#@profile
def opt_alpha(stdscr, alpha0, M, N, sample, vmin, dv, n, opt_tol=0.01, mise_tol=0.01, noniso=False):
    """Function that finds the optimal value of alpha for a given sample of stars.
    Given an initial guess of alpha, alpha0, it will draw M samples of size N from the resulting
    distribution f(v) computed using max_L. For each sample we perform the maximisation scheme
    to find the f(v) distribution for range of 10 alpha values between some upper and lower bounds.
    The mean integrated square error (MISE) is then computed for each alpha value. The value with the lowest
    MISE is the alpha_opt value and a new range of 10 alpha values is taken with upper and lower bounds set by
    the alpha values above and below opt_alpha in the previous range. The iterator repeats until the difference
    |alpha_opt_previous-alpha_opt| falls below the mise tolerance at which point the alpha_opt value is set as the new
    alpha_0 initial guess. The process then starts from a new initial range centered on the new alpha0. The iteration
    ends when the difference between |alpha_0_previous-alpha_opt| falls below the optimization tolerance

    alpha0: The initial guess of alpha
    M: The number of samples to compute the MISE for
    N: Number of stars in each sample
    pvals: Array with the tangential velocity vectors for all the  original stars in our sample
    rhatvals: Array with the unit vector for each sample star
    vmin: Vector indicating the anchor of v-space
    dv: The dims of each cell
    n: The dims of our box
    mise_tol: The desired logarithmic tolerance for the minimisation scheme. default value is 0.01.
    opt_tol: The desired logarithmic tolerance for the optimization scheme. default value is 0.01.

    """
    ### Preparing counters
    ti = time.time()

    opt_diff = 10
    opt_it = 0

    mise_diff = 10
    mise_it = 0
    mise_it_tot = 0

    fmin_calls = 0
    fmin_it = 0

    ###### Iterator below ######
    # Setting up the initial range and guess
    logalpha_0 = np.log10(alpha0)
    logalpha_opt = logalpha_0
    logalphavals = np.linspace(logalpha_0 - 1, logalpha_0 + 1, 10)  # The initial set of alpha values

    vv, rrind, cartcoords, coordinds, pvals, rhatvals = read_params(sample, vmin, dv,
                                                                    n)  # Generate the needed params from input

    while opt_diff >= opt_tol:
        if mise_diff < mise_tol or mise_it == 0:
            if 'phi0' in locals():
                phi0, _, its = max_L(alpha0, pvals, rhatvals, vmin, dv, n, phi0_guess=phi0, disp=0, noniso=noniso)
                fmin_calls += 1
                fmin_it += its
            else:
                # print("phi0 not defined, using phi_guess for first phi0")
                phi0, _, its = max_L(alpha0, pvals, rhatvals, vmin, dv, n, disp=0, noniso=noniso)
                fmin_calls += 1
                fmin_it += its
            stdscr.addstr(2, 1, '%s' % "".join(([' '] * 42)))
            stdscr.addstr(5, 1, '%s' % "".join(([' '] * 42)))
            stdscr.addstr(2, 1, str('opt it       : %s' % opt_it), curses.color_pair(0) | curses.A_BOLD)
            stdscr.addstr(5, 1, str('time elapsed : %s hrs' % (np.around((time.time() - ti) / 3600, 3))),
                          curses.color_pair(0) | curses.A_BOLD)
            stdscr.refresh()

        vxvals, vyvals, vzvals, smpcoords, fv0 = sample_pdf(phi0, M, N, rrind, cartcoords, coordinds, vv,
                                                            n)  # Generate M pseudosamples of size N from phi0

        ise = np.zeros((M, len(logalphavals)))  # container for the integrated square errors

        for i in range(M):

            pspvals, psrhatvals = pseudosample_pvals(vxvals, vyvals, vzvals, smpcoords, i)  # Get a pseudosample

            for j in range(len(logalphavals)):
                alpha = 10 ** (logalphavals[j])

                phi, phiall, its = max_L(alpha, pspvals, psrhatvals, vmin, dv, n, phi0_guess=phi0, disp=0,
                                         noniso=noniso)
                fmin_calls += 1
                fmin_it += its

                fv = np.exp(phi)

                ise[i][j] = np.sum((fv - fv0) ** 2)

        mise = np.mean(ise, axis=0)

        # Finds alpha_opt using the mise and sets up new range
        xrange_old = logalphavals
        yrange_old = mise
        logalpha_opt, logalphavals, mise_diff, step = tenstep_mise(logalphavals, mise, logalpha_opt)
        xrange = logalphavals
        ### For the debugger:
        t = (time.time() - ti) / (3600)
        mise_it_tot += 1
        mise_it += 1
        make_string(stdscr, logalpha_0, opt_it, mise_it, mise_it_tot, mise_diff, xrange, xrange_old, yrange_old, t,
                    step, 'tenstep')

        if mise_diff < mise_tol:
            opt_diff = abs(logalpha_opt - logalpha_0)

            logalpha_0 = logalpha_opt  # We set the optimised alpha value to be our new initial guess
            alpha0 = 10 ** (logalpha_0)
            logalphavals = np.linspace(logalpha_0 - 1, logalpha_0 + 1, 10)

            builtins.Nps = M

            opt_it += 1  # One opt iteration completed

    alpha_fin = 10 ** (logalpha_0)

    stdscr.addstr(0, 1, str('The optimal value for alpha is : %.2f' % alpha_fin), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(1, 1, str('To run took : %s opt alpha iterations' % opt_it), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(2, 1, str('            : %s tenstep iterations' % mise_it_tot), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(3, 1, str('            : %s fmin_cg calls' % fmin_calls), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(4, 1, str('            : %s fmin_cg iterations' % fmin_it), curses.color_pair(0) | curses.A_BOLD)
    stdscr.refresh()


    return alpha_fin


#@profile
def opt_alpha_ternary(stdscr, alpha0, M, N, sample, vmin, dv, n, opt_tol=0.01, mise_tol=0.01, noniso=False):
    """Function that finds the optimal value of alpha for a given sample of stars using a ternary search algorithm.
    Given an initial guess of alpha, alpha0, it will draw M samples of size N from the resulting
    distribution f(v) computed using max_L. For each sample we perform the maximisation scheme
    to find the f(v) distribution for the ternary search left, left_third, right_third, and right alpha
    values. The mean integrated square error (MISE) is then computed for each ternary alpha  value to select the
    new ternary range.

    This process repeats until the difference between |left-right| falls below the mise tolerance at which point
    the center of the ternary bounds is the new initial guess alpha_0. The ternary search then starts over until
    it finds |alpha_0_previous - alpha_0| is below the optimization tolerance

    alpha0: The initial guess of alpha
    M: The number of samples to compute the MISE for
    N: Number of stars in each sample
    pvals: Array with the tangential velocity vectors for all the  original stars in our sample
    rhatvals: Array with the unit vector for each sample star
    vmin: Vector indicating the anchor of v-space
    dv: The dims of each cell
    n: The dims of our box
    mise_tol: The desired logarithmic tolerance for the minimisation scheme. default value is 0.01.
    opt_tol: The desired logarithmic tolerance for the optimization scheme. default value is 0.01.

    """
    tree_file = make_treefile()
    with open(tree_file, 'a') as txt:
        txt.write('# Array shape: opt_it x log(alpha) x MISE\n')
    ### Preparing counters
    ti = time.time()

    opt_diff = 10
    opt_it = 0

    mise_diff = 10
    mise_it = 0
    mise_it_tot = 0

    fmin_calls = 0
    fmin_it = 0

    builtins.Nps = M
    ###### Iterator below ######
    # Setting up the initial range and guess

    logalpha_0 = np.log10(alpha0)
    xrange = np.linspace(logalpha_0 - 1, logalpha_0 + 1, 4)
    yrange = np.array([.0, .0, .0, .0])
    ind = np.array([True, True, True, True])

    logalphavals = xrange
    builtins.min_in_range = False
    vv, rrind, cartcoords, coordinds, pvals, rhatvals = read_params(sample, vmin, dv,
                                                                    n)  # Generate the needed params from input

    while opt_diff >= opt_tol:
        if mise_diff < mise_tol or mise_it == 0:
            if 'phi0' in locals():
                phi0, _, its = max_L(alpha0, pvals, rhatvals, vmin, dv, n, phi0_guess=phi0, disp=0, noniso=noniso)
                fmin_calls += 1
                fmin_it += its
            else:
                # print("phi0 not defined, using phi_guess for first phi0")
                phi0, _, its = max_L(alpha0, pvals, rhatvals, vmin, dv, n, disp=0, noniso=noniso)
                fmin_calls += 1
                fmin_it += its
            stdscr.addstr(2, 1, '%s' % "".join(([' '] * 42)))
            stdscr.addstr(5, 1, '%s' % "".join(([' '] * 42)))
            stdscr.addstr(2, 1, str('opt it       : %s' % opt_it), curses.color_pair(0) | curses.A_BOLD)
            stdscr.addstr(5, 1, str('time elapsed : %s hrs' % (np.around((time.time() - ti) / 3600, 3))),
                          curses.color_pair(0) | curses.A_BOLD)
            stdscr.refresh()

        vxvals, vyvals, vzvals, smpcoords, fv0 = sample_pdf(phi0, M, N, rrind, cartcoords, coordinds, vv,
                                                            n)  # Generate M pseudosamples of size N from phi0

        ise = np.zeros((Nps, 4))  # container for the integrated square errors

        for i in range(Nps):

            pspvals, psrhatvals = pseudosample_pvals(vxvals, vyvals, vzvals, smpcoords, i)  # Get a pseudosample

            for j in np.where(ind == True)[0]:
                alpha = 10 ** (logalphavals[j])

                phi, phiall, its = max_L(alpha, pspvals, psrhatvals, vmin, dv, n, phi0_guess=phi0, disp=0,
                                         noniso=noniso)
                fmin_calls += 1
                fmin_it += its

                fv = np.exp(phi)
                ise[i][j] = np.sum((fv - fv0) ** 2)

        mise = np.mean(ise, axis=0)

        np.copyto(yrange, mise, where=ind)

        with open(tree_file, 'a') as txt:
            data = np.ones(4, 3)
            data[:, 0] = np.ones(4, 1) * opt_it
            data[:, 1] = xrange
            data[:, 2] = yrange
            np.savetxt(txt, data)

        xrange_old = xrange
        yrange_old = yrange
        logalphavals, mise_diff, xrange, yrange, ind, step = ternary_mise(xrange, yrange)  # Get the new ternary range

        ### For the debugger:
        t = (time.time() - ti) / (3600)
        mise_it_tot += 1
        mise_it += 1
        make_string(stdscr, logalpha_0, opt_it, mise_it, mise_it_tot, mise_diff, xrange, xrange_old, yrange_old, t,
                    step, 'ternary')

        if mise_diff < mise_tol:
            logalpha_opt = (xrange[0] + xrange[3]) / 2

            opt_diff = abs(logalpha_opt - logalpha_0)

            logalpha_0 = logalpha_opt
            alpha0 = 10 ** (logalpha_0)

            xrange = np.linspace(logalpha_0 - 1, logalpha_0 + 1, 4)
            yrange = np.array([.0, .0, .0, .0])
            ind = np.array([True, True, True, True])

            logalphavals = xrange
            builtins.min_in_range = False
            builtins.Nps = M
            opt_it += 1

    alpha_fin = 10 ** (logalpha_0)

    stdscr.addstr(0, 1, str('The optimal value for alpha is : %.2f' % alpha_fin), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(1, 1, str('To run took : %s opt alpha iterations' % opt_it), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(2, 1, str('            : %s ternary iterations' % mise_it_tot), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(3, 1, str('            : %s fmin_cg calls' % fmin_calls), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(4, 1, str('            : %s fmin_cg iterations' % fmin_it), curses.color_pair(0) | curses.A_BOLD)
    stdscr.refresh()


    return alpha_fin


#@profile
def opt_alpha_gss(stdscr, alpha0, M, N, sample, vmin, dv, n, opt_tol=0.01, mise_tol=0.01, noniso=False):
    """Function that finds the optimal value of alpha for a given sample of stars using a golden-section search
    algorithm. Given an initial guess of alpha, alpha0, it will draw M samples of size N from the resulting
    distribution f(v) computed using max_L. For each sample we perform the maximisation scheme
    to find the f(v) distribution for the golden-section four values. They are:

    x1: Leftmost corner
    x2: 2nd leftmost point
    x3: 2nd rightmost point
    x4: Rightmost point

    They are specified by the intervals x2-x1:x3-x2:x4-x3 having the intervals widths in ratio 2-g:2g-3:2-g, where g
    is the golden ratio. Similarly to a ternary search. if f(x3) < f(x2) the minimum is within x2 to x4 with new range
    x1',x2',x3',x4' = x2, x3, ?, x4. If f(x3) > f(x2) the minimum is within x1 to x3 with new range
    x1',x2',x3',x4' = x1, ?, x2, x3. We calculate a new x3' or x2' point and estimate f() for it. The new range is
    chosen such that the widths follow the previous ratio.

    The mean integrated square error (MISE) is our f() and the process repeats until the new range limit difference
    falls below the mise tolerance, i.e. either |x2 - x4| or |x1 - x3|. At that point the central value (x3 or x2 in
    the previous cases) is set as the new guess for alpha0. The process starts over until it finds that
    |alpha_0_previous - alpha_0| is below the optimization tolerance when the search concludes.

    alpha0: The initial guess of alpha
    M: The number of samples to compute the MISE for
    N: Number of stars in each sample
    pvals: Array with the tangential velocity vectors for all the  original stars in our sample
    rhatvals: Array with the unit vector for each sample star
    vmin: Vector indicating the anchor of v-space
    dv: The dims of each cell
    n: The dims of our box
    mise_tol: The desired logarithmic tolerance for the minimisation scheme. default value is 0.01.
    opt_tol: The desired logarithmic tolerance for the optimization scheme. default value is 0.01.

    """
    tree_file = make_treefile()
    try:
        s = open(tree_file, 'r').read().split('\n')
    except FileNotFoundError:
        with open(tree_file, 'a') as txt:
            txt.write('# Array shape: opt_it x log(alpha) x MISE\n')
    else:
        os.remove(tree_file)
        with open(tree_file, 'a') as txt:
            txt.write('# Array shape: opt_it x log(alpha) x MISE\n')
        stdscr.addstr(0, 51, ("Found an incomplete tree_file, removed it"),
                      curses.color_pair(4) | curses.A_BOLD | curses.A_UNDERLINE)
        stdscr.refresh()
    ### Preparing counters
    ti = time.time()

    opt_diff = 10
    opt_it = 0

    mise_diff = 10
    mise_it = 0
    mise_it_tot = 0

    fmin_calls = 0
    fmin_it = 0

    builtins.Nps = M
    ###### Iterator below ######
    # Setting up the initial range and guess

    logalpha_0 = np.log10(alpha0)

    invphi = (np.sqrt(5) - 1) / 2  # 1 / phi
    invphi2 = (3 - np.sqrt(5)) / 2  # 1 / phi^2

    x1 = logalpha_0 - 1
    x4 = logalpha_0 + 1
    x2 = x1 + invphi2 * (x4 - x1)
    x3 = x1 + invphi * (x4 - x1)
    xrange = np.array([x1, x2, x3, x4])
    yrange = np.array([.0, .0, .0, .0])
    ind = np.array([True, True, True, True])

    logalphavals = xrange
    builtins.min_in_range = False
    vv, rrind, cartcoords, coordinds, pvals, rhatvals = read_params(sample, vmin, dv,
                                                                    n)  # Generate the needed params from input

    while opt_diff >= opt_tol:
        if mise_diff < mise_tol or mise_it == 0:
            if 'phi0' in locals():
                phi0, _, its = max_L(alpha0, pvals, rhatvals, vmin, dv, n, phi0_guess=phi0, disp=0, noniso=noniso)
                fmin_calls += 1
                fmin_it += its
            else:
                # print("phi0 not defined, using phi_guess for first phi0")
                phi0, _, its = max_L(alpha0, pvals, rhatvals, vmin, dv, n, disp=0, noniso=noniso)
                fmin_calls += 1
                fmin_it += its
            stdscr.addstr(2, 1, '%s' % "".join(([' '] * 42)))
            stdscr.addstr(5, 1, '%s' % "".join(([' '] * 42)))
            stdscr.addstr(2, 1, str('opt it       : %s' % opt_it), curses.color_pair(0) | curses.A_BOLD)
            stdscr.addstr(5, 1, str('time elapsed : %s hrs' % (np.around((time.time() - ti) / 3600, 3))),
                          curses.color_pair(0) | curses.A_BOLD)
            stdscr.refresh()
        vxvals, vyvals, vzvals, smpcoords, fv0 = sample_pdf(phi0, Nps, N, rrind, cartcoords, coordinds, vv,
                                                            n)  # Generate M pseudosamples of size N from phi0
        ise = np.zeros((Nps, 4))  # container for the integrated square errors
        for i in range(Nps):

            pspvals, psrhatvals = pseudosample_pvals(vxvals, vyvals, vzvals, smpcoords, i)  # Get a pseudosample

            for j in np.where(ind == True)[0]:
                alpha = 10 ** (logalphavals[j])

                phi, phiall, its = max_L(alpha, pspvals, psrhatvals, vmin, dv, n, phi0_guess=phi0, disp=0,
                                         noniso=noniso)
                fmin_calls += 1
                fmin_it += its

                fv = np.exp(phi)
                ise[i][j] = np.sum((fv - fv0) ** 2)

        mise = np.mean(ise, axis=0)

        try:
            np.copyto(yrange, mise, where=ind)
        except TypeError:
            raise Exception('\nyrange : ' + str(yrange) + ', ' + str(yrange.dtype) +
                            '\nmise   : ' + str(mise) + ', ' + str(mise.dtype) +
                            '\nind    : ' + str(ind) + ', ' + str(ind.dtype))

        with open(tree_file, 'a') as txt:
            data = np.ones((4, 3))
            data[:, 0] = np.ones(4) * opt_it
            data[:, 1] = xrange
            data[:, 2] = yrange
            np.savetxt(txt, data)

        xrange_old = xrange
        yrange_old = yrange
        logalphavals, mise_diff, xrange, yrange, ind, step = gss_mise(xrange, yrange)  # Get the new gss range

        ### For the debugger:
        t = (time.time() - ti) / (3600)
        mise_it_tot += 1
        mise_it += 1
        make_string(stdscr, logalpha_0, opt_it, mise_it, mise_it_tot, mise_diff, xrange, xrange_old, yrange_old, t,
                    step, 'gss')

        if mise_diff < mise_tol:
            logalpha_opt = (xrange[0] + xrange[3]) / 2

            opt_diff = abs(logalpha_opt - logalpha_0)

            logalpha_0 = logalpha_opt
            alpha0 = 10 ** (logalpha_0)

            x1 = logalpha_0 - 1
            x4 = logalpha_0 + 1
            x2 = x1 + invphi2 * (x4 - x1)
            x3 = x1 + invphi * (x4 - x1)
            xrange = np.array([x1, x2, x3, x4])
            yrange = np.array([.0, .0, .0, .0])
            ind = np.array([True, True, True, True])

            logalphavals = xrange
            builtins.min_in_range = False
            builtins.Nps = M
            opt_it += 1

    alpha_fin = 10 ** (logalpha_0)

    with open(tree_file, 'a') as txt:
        txt.write('#FINISHED\n')

    stdscr.erase()
    stdscr.addstr(0, 1, str('The optimal value for alpha is : %.2f' % alpha_fin), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(1, 1, str('To run took : %s opt alpha iterations' % opt_it), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(2, 1, str('            : %s gss iterations' % mise_it_tot), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(3, 1, str('            : %s fmin_cg calls' % fmin_calls), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(4, 1, str('            : %s fmin_cg iterations' % fmin_it), curses.color_pair(0) | curses.A_BOLD)
    stdscr.addstr(6, 1, "Press any key 3 times to log output and close program!",
                  curses.color_pair(2) | curses.A_BOLD | curses.A_BLINK)
    stdscr.addstr(7, 1, "Press any key 3 times to log output and close program!",
                  curses.color_pair(2) | curses.A_BOLD | curses.A_BLINK)
    stdscr.addstr(8, 1, "Press any key 3 times to log output and close program!",
                  curses.color_pair(2) | curses.A_BOLD | curses.A_BLINK)
    stdscr.refresh()
    stdscr.getkey()
    stdscr.getkey()
    stdscr.getkey()

    return alpha_fin


#@profile
def read_params(sample, vmin, dv, n):
    pvals, rhatvals = calc_p_rhat(sample)
    vxmin, vymin, vzmin = vmin
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    # We want to compute the rhat values for every sample
    # So we just draw M sets of size N of the coordinates from the original
    # sample that we will use for each iteration of the sample.

    sample.set_representation_cls(coord.CartesianRepresentation)
    cartcoords = np.zeros((len(sample), 3))
    cartcoords[:, 0] = sample.u
    cartcoords[:, 1] = sample.v
    cartcoords[:, 2] = sample.w
    coordinds = np.linspace(0, len(cartcoords) - 1, len(cartcoords))

    vxmax, vymax, vzmax = vxmin + nx * dvx, vymin + ny * dvy, vzmin + nz * dvz

    vx_bins = np.linspace(vxmin, vxmax, nx+1)
    vy_bins = np.linspace(vymin, vymax, ny+1)
    vz_bins = np.linspace(vzmin, vzmax, nz+1)

    vxc = (vx_bins[1:] + vx_bins[:-1]) / 2
    vyc = (vy_bins[1:] + vy_bins[:-1]) / 2
    vzc = (vz_bins[1:] + vz_bins[:-1]) / 2

    vxx, vyy, vzz = np.meshgrid(vxc, vyc, vzc, indexing='ij')  # The centre values of each cell in our v-space box
    vv = (vxx, vyy, vzz)
    ind = np.indices((nx, ny, nz))

    rind = np.ravel_multi_index(ind, (nx, ny, nz))  # An array containing the 3D coordinates of our box

    rrind = np.ravel(rind)

    return vv, rrind, cartcoords, coordinds, pvals, rhatvals


#@profile
def sample_pdf(phi0, M, N, rrind, cartcoords, coordinds, vv, n):
    nx, ny, nz = n

    fv0 = np.exp(phi0)  # The pdf f(v) given our inital guess alpha0
    fv0s = np.sum(fv0)
    prob = np.ravel(fv0 / fv0s)  # The normalised probability
    np.random.seed(0)
    smp = np.random.choice(rrind, (M, N), p=prob)  # Creation of M samples of size n given f(v)
    smpcoordinds = np.random.choice(coordinds, (M, N)).astype(int)  # We also draw random positions for these stars

    smpcoords = cartcoords[smpcoordinds]

    smpx, smpy, smpz = np.asarray(np.unravel_index(smp, (nx, ny, nz)))

    vxx, vyy, vzz = vv
    vxvals = vxx[smpx, smpy, smpz]
    vyvals = vyy[smpx, smpy, smpz]
    vzvals = vzz[smpx, smpy, smpz]

    return vxvals, vyvals, vzvals, smpcoords, fv0


#@profile
def pseudosample_pvals(vxvals, vyvals, vzvals, smpcoords, i):
    smpvx, smpvy, smpvz = vxvals[i], vyvals[i], vzvals[i]  # For every pseudosample we get velocities

    coordx, coordy, coordz = smpcoords[i, :, 0], smpcoords[i, :, 1], smpcoords[i, :, 2]  # ... as well as coordinates

    psample = coord.Galactic(u=coordx * u.pc, v=coordy * u.pc, w=coordz * u.pc, U=smpvx * (u.km / u.s),
                             V=smpvy * (u.km / u.s), W=smpvz * (u.km / u.s),
                             representation_type=coord.CartesianRepresentation,
                             differential_type=coord.CartesianDifferential)

    psample.set_representation_cls(coord.SphericalRepresentation, coord.SphericalCosLatDifferential)

    pspvals, psrhatvals = calc_p_rhat(psample)

    return pspvals, psrhatvals


#@profile
def tenstep_mise(logalphavals, mise, logalpha_opt_former):
    '''Function that performs the ten step zoom iteration and returns the new alpha0 and alpha range in log.
    Since our initial guess is very broad and covers 10 orders of magnitude, we have to narrow it down.
    This is done by taking the alpha values to the right and left of our optimal value to define the new range.
    If the optimal value is at the edges of our array, then we take the missing upper or lower bound
    to be 2 magnitudes larger than alpha_opt

    We also return the mise_diff to see if we are satisfied with the minimum. If we are at the edge however,
    this is set to 10 simply to ensure the iteration continues.'''

    optind = np.argwhere(mise == np.amin(mise))[0][0]  # The index of the optimimal alpha
    logalpha_opt = logalphavals[optind]
    mise_diff = abs(logalpha_opt_former - logalpha_opt)  # Check the improvement in the minimum from former best guess

    if (mise[:-1] < mise[1:]).all():
        lower = logalpha_opt - 1
        upper = logalphavals[optind + 1]
        mise_diff = 10
        step = 'Down'
    elif (mise[:-1] > mise[1:]).all():
        lower = logalphavals[optind - 1]
        upper = logalpha_opt + 1
        mise_diff = 10
        step = 'Up'
    elif (mise[:optind] > mise[1:optind + 1]).all() and (mise[optind:-1] < mise[optind + 1:]).all():
        lower = logalphavals[optind - 1]
        upper = logalphavals[optind + 1]
        step = 'Zoom'
    else:
        M = Nps
        builtins.Nps = M + 5
        lower = logalphavals[-1]
        upper = logalphavals[0]
        step = 'Again'


    return logalpha_opt, np.linspace(lower, upper, 10), mise_diff, step


#@profile
def ternary_mise(logalphavals, mise):
    '''Function that performs the ternary iteration and returns the new alpha0 and alpha range in log. If the
    lower(upper) edge has the smallest function value, we step 1 in log down(up) from left(right). I.e
    the new range will be, using the old range values, left-1 <-> left_third (right_third <-> right+1)'''

    class Tern:
        #@profile
        def __init__(self, mise, alpha):
            self.x = alpha
            self.y = mise

    left, left_third, right_third, right = [Tern(mise[i], logalphavals[i]) for i in range(4)]

    # This is if either edge has the smallest mise and we need to expand the range to ensure we have the minimum
    # in our range of left <-> right.
    if (left.y < left_third.y < right_third.y < right.y):
        if not (min_in_range):
            xrange = np.linspace(left.x - 1, left_third.x, 4)
            yrange = np.array([0, 0, 0, left_third.y])
            logalphavals = np.array(list(xrange[1:3]) + [0])
            step = 'Down'

            return logalphavals, 10, xrange, yrange, np.array([True, True, True, False]), step
    elif (left.y > left_third.y > right_third.y > right.y):
        if not (min_in_range):
            xrange = np.linspace(right_third.x, right.x + 1, 4)
            yrange = np.array([right_third.x, 0, 0, 0])
            logalphavals = np.array([0] + list(xrange[:1]))
            step = 'Up'

            return logalphavals, 10, xrange, yrange, np.array([False, True, True, True]), step

    # This is if we have the minmum in the range left <-> right. Script only enters here in this case.
    elif ((left_third.y < left.y) and (right_third.y < right.y)):
        builtins.min_in_range = True
        if left_third.y < right_third.y:
            xrange = np.linspace(left.x, right_third.x, 4)
            yrange = np.array([left.y, 0, 0, right_third.y])
            logalphavals = np.array([0] + list(xrange[1:3]) + [0])
            step = 'Left'

            return logalphavals, abs(left.x - right_third.x), xrange, yrange, np.array([False, True, True, False]), step

        else:
            xrange = np.linspace(left_third.x, right.x, 4)
            yrange = np.array([left_third.y, 0, 0, right.y])
            logalphavals = np.array([0] + list(xrange[1:3]) + [0])
            step = 'Right'

            return logalphavals, abs(left_third.x - right.x), xrange, yrange, np.array([False, True, True, False]), step

    # If we go here the y-values imply a non-unimodal function, to which we must increase our
    # number of pseudosamples M and recalculate the MISE of the current range.
    else:
        M = Nps
        builtins.Nps = M + 5
        xrange = np.array([left.x, left_third.x, right_third.x, right.x])
        yrange = np.array([0, 0, 0, 0])
        step = 'Again'
        logalphavals = xrange

        return logalphavals, abs(xrange[0] - xrange[3]), xrange, yrange, np.array([True, True, True, True]), step


#@profile
def gss_mise(xrange, yrange):
    '''Function that performs the Golden section iteration and returns which new alpha needs to be evaluated, as well
    as its location within the array. If x1(x4) has the smallest function value, we step 1 in log down(up) from
    x1(x4). I.e the new range will be, using the old range values, x1-1 <-> x2 (x3 <-> x4+1). Calculates mise_diff.'''

    x1, x2, x3, x4 = xrange
    y1, y2, y3, y4 = yrange
    # This is if either edge has the smallest mise and we need to expand the range to ensure we have the minimum
    # in our range of x1 <-> x4.
    if (y1 < y2 < y3 < y4):
        if not (min_in_range):
            x1 = x1 - 1;
            x4 = x2
            xrange = np.array(
                [x1, (x1 + (3 - np.sqrt(5)) / 2 * (x4 - x1)), (x1 + (np.sqrt(5) - 1) / 2 * (x4 - x1)), x4])
            yrange = np.array([0, 0, 0, y2],dtype=float)
            logalphavals = np.array(list(xrange[1:]) + [0])
            step = 'Down'

            return logalphavals, 10, xrange, yrange, np.array([True, True, True, False]), step
    elif (y1 > y2 > y3 > y4):
        if not (min_in_range):
            x1 = x3;
            x4 = x4 + 1
            xrange = np.array(
                [x1, (x1 + (3 - np.sqrt(5)) / 2 * (x4 - x1)), (x1 + (np.sqrt(5) - 1) / 2 * (x4 - x1)), x4])
            yrange = np.array([y3, 0, 0, 0],dtype=float)
            logalphavals = np.array([0] + list(xrange[1:]))
            step = 'Up'

            return logalphavals, 10, xrange, yrange, np.array([False, True, True, True]), step

    # This is if we have the minmum in the range left <-> right. Script only enters here in this case.
    elif ((y2 < y1) and (y3 < y4)):
        builtins.min_in_range = True
        if y2 < y3:
            xrange = np.array([x1, (x1 + (3 - np.sqrt(5)) / 2 * (x4 - x1)), x2, x3])
            yrange = np.array([y1, 0, y2, y3],dtype=float)
            step = 'Left'
            logalphavals = np.array([0, xrange[1], 0, 0])

            return logalphavals, abs(xrange[0] - xrange[3]), xrange, yrange, np.array([False, True, False, False]), step

        else:
            xrange = np.array([x2, x3, (x1 + (np.sqrt(5) - 1) / 2 * (x4 - x1)), x4])
            yrange = np.array([y2, y3, 0, y4])
            step = 'Right'
            logalphavals = np.array([0, 0, xrange[2], 0],dtype=float)

            return logalphavals, abs(xrange[0] - xrange[3]), xrange, yrange, np.array([False, False, True, False]), step

    # If we go here the y-values imply a non-unimodal function, to which we must increase our
    # number of pseudosamples M and recalculate the MISE of the current range.
    else:
        M = Nps
        builtins.Nps = M + 5
        xrange = np.array([x1, x2, x3, x4])
        yrange = np.array([0, 0, 0, 0],dtype=float)
        step = 'Again'
        logalphavals = xrange

        return logalphavals, abs(xrange[0] - xrange[3]), xrange, yrange, np.array([True, True, True, True]), step
