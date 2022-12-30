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
from scipy.optimize import fmin_cg, minimize
from scipy.interpolate import interpn
from scipy.ndimage import zoom
from datetime import date
from decimal import Decimal
from IPython.display import clear_output
from types import ModuleType, FunctionType
from gc import get_referents


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


def zoomed_mxl(mxl):
    '''Takes a mxl result and upscales it with interpolation'''
    phi0_guess = zoom(mxl, zoom=2, order=3) - np.log(8)

    return phi0_guess


def callback_mg(x):
    i = builtins.grid_step
    builtins.L[i].append(builtins.current_L)
    builtins.gradL[i].append(builtins.current_gradL)


def callback(x):
    builtins.L.append(builtins.current_L)
    builtins.gradL.append(builtins.current_gradL)


def make_polar(sample, pvals, rhat):
    sample_gc = sample.galactocentric
    x = sample_gc.x
    y = sample_gc.y
    z = sample_gc.z

    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z/r)
    phi = np.arctan2(y, x)

    vsun = sample_gc.galcen_v_sun.d_xyz

    R = np.array([[np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)],
                  [np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -np.sin(theta)],
                  [-np.sin(phi), np.cos(phi), np.zeros(len(phi))]]).transpose(2, 0, 1)

    pvals_polar = ( R @ (vsun.value + pvals)[..., np.newaxis]).squeeze()
    rhat_polar = ( R @ rhat[..., np.newaxis]).squeeze()

    return pvals_polar, rhat_polar


def calc_sigma2(pvals, rhat, give_vmean=False, noniso=False):
    """Function that applies equation 12 of DB98 for a set of stars from their tangential velocities and unit vectors.
    Returns the velocity dispersion tensor.

    pvals: array of the pk vectors for all N stars in our sample. Should have dims (N,3)
    rhat: array of N unit vectors, one for each star
    give_vmean: if True, returns the computed mean velocity vector value for the given sample
    noniso: if True, we no longer assume that the sample is isotropic"""

    pmean = pvals.mean(axis=-2)

    # Fast way for outer product eq (4) of DB98
    A = np.identity(3) - rhat[..., np.newaxis] * rhat[..., np.newaxis, :]


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

    pos = np.stack(np.meshgrid(vxc,vyc,vzc, indexing='ij'),axis=3)

    fA = st.multivariate_normal(mean=v0, cov=np.diag(dispA**2))
    fB = st.multivariate_normal(mean=v0, cov=np.diag(dispB**2))

    phi = np.log((fA.pdf(pos) + fB.pdf(pos))/2)
    return phi


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

    exphi = np.exp(phi)
    Kphi = (Kvals @ exphi.T) # Order all Kphi values in 1D arrays and compute the sum of exp(phi)*K(k|l) for each star
    Kphi_sum_tot = np.log(Kphi[Kphi != 0]).sum() # To make sure we don't get infinities and .sum() gives the double sum in the first term

    phi_unr = np.reshape(phi, n)
    phixhi_sum = (sec_der(phi_unr, sigma2, dv) ** 2).sum()

    t1 = Kphi_sum_tot / N
    t2 = exphi.sum()
    t3 = ((alpha * dv[0] * dv[1] * dv[2]) / 2) * phixhi_sum

    L_tilde = t1 - t2 - t3 # eq. 31 in DB98
    neg_L = -1 * L_tilde  # Since we want to maximize L_tilde, we should minimize -L_tilde
    builtins.current_L = L_tilde
    return neg_L


def get_grad_neg_L(phi, *args):
    """In this function we compute the gradient of L. We compute the derivative for each cell and return a
    1D array of length (nx*ny*nz).

    args: see get_L

    """
    Kvals, N, alpha, dv, n, sigma2 = args
    exphi = np.exp(phi)

    Kphi_sum = (Kvals @ exphi.T)
    Kphi_sum[Kphi_sum != 0] = Kphi_sum[Kphi_sum != 0]**(-1)
    K_term = exphi*(Kphi_sum.T @ Kvals) # The final array with the first term for each cell

    phi_unr = np.reshape(phi,n)
    dphixhi = grad_sec_der(phi_unr, sigma2, dv)

    t1 = K_term/N
    t2 = exphi
    t3 = ((alpha * dv[0] * dv[1] * dv[2]) / 2) * dphixhi.ravel()

    grad_L = np.asarray(t1 - t2 - t3).reshape(len(phi), )

    neg_grad_L = -1 * grad_L
    builtins.current_gradL = np.linalg.norm(neg_grad_L)

    return neg_grad_L


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


def max_L(alpha, pvals, rhatvals, vmin, dv, n, phi0_guess=[], v0_guess=[], disp_guess=[], noniso=False, printing=True, polar=False):
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
    if polar:
        sigma2 = disp_guess**2

    Kvals = KvalsSparseMethod(pvals, rhatvals, vmin, dv, n)

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


def multigrid_max_L(alpha, pvals, rhatvals, vmin, dv, n, phi0_guess=[], v0_guess=[], disp_guess=[], noniso=False, printing=False, polar=False):
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
    if polar:
        sigma2 = disp_guess**2

    builtins.L     = [[] for _ in range(len(box_steps))]
    builtins.gradL = [[] for _ in range(len(box_steps))]
    fmin_it = 0

    mxl_dict = {}
    for grid_step, n in enumerate(box_steps):
        if grid_step == len(box_steps)-1:
            printing = True
        builtins.grid_step = grid_step
        dv = (vmax-vmin)/n

        print(f'Starting Kvals after {(time.time() - builtins.ti)/60:.2f}mins...')
        Kvals = KvalsSparseMethod(pvals, rhatvals, vmin, dv, n)

        args = (Kvals, N, alpha, dv, n, sigma2)
        phi0r = np.ravel(phi0)  # fmin_cg only takes one-dimensional inputs for the initial guess

        ### This is where the minimization occurs
        print(f'Started fmin_cg on {(n[0],n[1],n[2])} grid after {(time.time() - builtins.ti)/60:.2f} mins...', end='', flush=True)
        builtins.L[grid_step].append(-1*get_neg_L(phi0r,Kvals, N, alpha, dv, n, sigma2))
        builtins.gradL[grid_step].append(np.linalg.norm(get_grad_neg_L(phi0r,Kvals, N, alpha, dv, n, sigma2)))

#        mxl, fopt, fcalls, gcalls, flag, phi_all = fmin_cg(get_neg_L, phi0r, fprime=get_grad_neg_L, gtol=1e-6, args=args, retall=True, disp=False, full_output=True, callback=callback_mg)
        optr = minimize(fun=get_neg_L, x0=phi0r, method='CG', jac=get_grad_neg_L, args=args, callback=callback_mg, options={'gtol': 1e-6, 'disp': False, 'return_all': True})
        mxl, fopt, fcalls, gcalls, flag = optr.x, optr.fun, optr.nfev, optr.njev, optr.status
        fmin_it += optr.nit

        mxl_dict['x'.join(n.astype(str))] = mxl.reshape(n)

        print(colored('Finished!','green',attrs=['bold','underline']),end='')
        print(' fopt : %s' % fopt)

        if grid_step == len(box_steps)-1:
            fmin_cg_output(fopt, fcalls, gcalls, flag, fmin_it)
        else:
            phi0 = zoomed_mxl(mxl.reshape(n))

    builtins.n = n
    builtins.dv = dv

    mxlnew = mxl.reshape(n)

    return mxl_dict, fmin_it


def KvalsSparseMethod(pvals, rhat, vmin, dv, n):
    n_stars = len(pvals)
    allocated_memory = psutil.virtual_memory().available*0.8/1e9
    required_memory = 8*n_stars*np.ceil(np.linalg.norm(n)).astype(int)/1e9

    if required_memory > allocated_memory:
        raise MemoryError(f'Required memory for {required_memory:.2f} GB exceeds available memory {allocated_memory:.2f} GB')

    n_stars = len(pvals)
    lil = scisp.lil_matrix((n_stars, np.prod(n)))

    for i in range(n_stars):
        coords, vals = calc_sparse_K(pvals[i], rhat[i], vmin, dv, n)
        lil[i, coords] = vals

    return lil.tocsc()


def calc_sparse_K(pk, rhat, vmin, dv, n):
    '''Calculate the values of K simultaneously for all bins for a given star with tangential velocity
        vector pk and a unit vector rhat'''


    vmax = vmin + dv*n

    """We now solve the line equation v = pk + vr*rhat, where v are the intersections of the tangential velocity line
    and the boundaries of each bin"""

    vbins = (np.linspace(vmin[i], vmax[i], n[i] + 1) for i in range(3))

    vr = [(np.array(vb) - pk[i])/rhat[i] for i, vb in enumerate(vbins)]

    """After solving the line equation for each dim we remove the vr values which solve the equation
    for bins outside of our specified box."""

    vrmax = min(map(np.max, vr))
    vrmin = max(map(np.min, vr))

    vr = [arr[(arr <= vrmax) & (arr >= vrmin)] for arr in vr]

    vr = np.concatenate(vr)
    vr.sort()  # We obtain an ordered list with all vr values for intersections between entry and exit points

    if len(vr) == 0:
        return np.array([]), np.array([])

    vr_prime = (vr[:-1] + vr[1:]) / 2
    pks = np.ones((len(vr_prime), len(pk))) * pk[np.newaxis]
    rhats = np.ones((len(vr_prime), len(rhat))) * rhat[np.newaxis]
    vmins = np.ones((len(vr_prime), len(vmin))) * vmin[np.newaxis]
    vr_primestack = np.ones((len(vr_prime), 3))*vr_prime[:, np.newaxis]

    """We now solve the line equation again for values in the middle of each bin with a line segment in it.
    This gives us the coordinates for each relevant bin, given in line_bins.
    Finally we compute the length of each segment and output the multi_index and values."""

    v_prime = pks + vr_primestack * rhats

    line_bins = ((v_prime - vmins) // dv).astype(int)

    line_len = vr[1:] - vr[:-1]  # Gets the length of each line
    non_zero = np.nonzero(line_len)
    line_len = line_len[non_zero]  # Removes duplicate intersections
    line_bins = line_bins[non_zero].T

    vals = (line_len / np.prod(dv))
    return np.ravel_multi_index(line_bins, n), vals


# # # # THE FOLLOW WAS MADE OBSOLITE BY # # # #
# # # pvals = sample.velocity.d_xyz.T.unmasked.value
# # # rhatvals = sample.spherical.unit_vectors()['distance'].xyz.T.unmasked.value

def calc_p_rhat(sample, polar=False):
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
    A = 0*A.unit
    B = 0*B.unit

    if isinstance(sample, Table):
        b = sample['b']
        l = sample['l']
        mul_obs = sample['pml_cosb'].to(u.km / (u.s * u.kpc), equivalencies=u.dimensionless_angles())
        mub_obs = sample['pmb'].to(u.km / (u.s * u.kpc), equivalencies=u.dimensionless_angles())
        p = sample['parallax'].to(u.kpc, equivalencies=u.parallax())

    elif isinstance(sample, (coord.builtin_frames.galactic.Galactic, coord.sky_coordinate.SkyCoord)):
        b = sample.b
        l = sample.l
        mul_obs = sample.pm_l_cosb.to(u.km / (u.s * u.kpc), equivalencies=u.dimensionless_angles())
        mub_obs = sample.pm_b.to(u.km / (u.s * u.kpc), equivalencies=u.dimensionless_angles())
        p = sample.distance.to(u.kpc)

    cosl = np.cos(l)
    cosb = np.cos(b)
    sinl = np.sin(l)
    sinb = np.sin(b)

    if not polar:
        mul = mul_obs - A * np.cos(2 * l) - B
        mub = mub_obs + A * np.sin(2 * l) * cosb * sinb

        pvals = (p * np.vstack((-sinl*mul - cosl*sinb*mub,
                                cosl*mul - sinl*sinb*mub,
                                cosb*mub))).T

        rhatvals = np.array([cosb * cosl, cosb * sinl, sinb]).T

        return pvals, rhatvals

    elif polar:
        mul = mul_obs
        mub = mub_obs
        pvals = (p * np.vstack((-sinl*mul - cosl*sinb*mub,
                                cosl*mul - sinl*sinb*mub,
                                cosb*mub))).T
        rhatvals = np.array([cosb * cosl, cosb * sinl, sinb]).T

        # Get polar coordinates. We start by creating a Galactocentric frame copy of the sample
        sample_galcen = sample.transform_to(frame='galactocentric')
        x = sample_galcen.cartesian.x.value
        y = sample_galcen.cartesian.y.value
        z = sample_galcen.cartesian.z.value

        # Get solar velocity in Galactocentric frame
        v_sun = np.array([sample_galcen.galcen_v_sun.d_x.value,
                          sample_galcen.galcen_v_sun.d_y.value,
                          sample_galcen.galcen_v_sun.d_z.value])*u.km/u.s

        r = np.sqrt(x**2 + y**2 + z**2)
        theta = np.arccos(z/r)
        phi = np.arctan2(y, x)

        R = np.array([[np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)],
                     [np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -np.sin(theta)],
                     [-np.sin(phi), np.cos(phi), np.zeros_like(phi)]]).transpose(2, 0, 1)
        # Get polar pvals
        pvals += v_sun # We add the motion of the sun
        pvals_polar = (R @ pvals[..., np.newaxis]).squeeze()
        rhatvals_polar = (R @ rhatvals[..., np.newaxis]).squeeze()
        return pvals_polar, rhatvals_polar
