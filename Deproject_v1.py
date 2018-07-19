import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import astropy.units as u
import astropy.coordinates as coord
import scipy.stats as st
from astropy.table import Table
from astropy.io import ascii
from scipy import sparse
from scipy.optimize import minimize
from scipy.optimize import fmin_cg
import sparse as sp

def calc_K1(pk,rhat,vmin,dv,n):
    '''Calculate the values of K simultaneously for all bins for a given star with p, rhat'''
    
    vxmin, vymin, vzmin = vmin
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    pkx, pky, pkz = pk
    rhatx, rhaty, rhatz = rhat
    
    vxmax, vymax, vzmax = vxmin+nx*dvx,vymin+ny*dvy,vzmin+nz*dvz
    
    """We now solve the line equation v = pk + vr*rhat, where v are the intersections of the tangential velocity line 
    and the boundaries of each bin"""
    
    vx_bins = np.arange(vxmin, vxmax+dvx, dvx)
    vy_bins = np.arange(vymin, vymax+dvy, dvy)
    vz_bins = np.arange(vzmin, vzmax+dvz, dvz)
    
    vrx = (vx_bins-pkx)/rhatx
    vry = (vy_bins-pky)/rhaty
    vrz = (vz_bins-pkz)/rhatz
    
    """After solving the line equation for each dim we remove the vr values which solve the equation
    for bins outside of our specified box."""
    
    vrmax = min(max(vrx),max(vry),max(vrz))
    vrmin = max(min(vrx),min(vry),min(vrz))
    
    vrx = vrx[(vrx<=vrmax) & (vrx>=vrmin)]
    vry = vry[(vry<=vrmax) & (vry>=vrmin)]
    vrz = vrz[(vrz<=vrmax) & (vrz>=vrmin)]
    vr = np.concatenate((vrx,vry,vrz))
    vr.sort() #We obtain an ordered list with all vr values for intersections between entry and exit points
    
    if len(vr)==0:
        return K
    
    vr_prime =(vr[:-1] + vr[1:]) / 2

    pks = np.zeros((len(vr_prime),len(pk)))
    pks[:] = pk
    rhats = np.zeros((len(vr_prime),len(rhat)))
    rhats[:] = rhat
    vmins = np.zeros((len(vr_prime),len(vmin)))
    vmins[:] = vmin
    vr_primestack = np.ones((len(vr_prime),3))
    vr_primestack *= vr_prime[:,None] 

    """We now solve the line equation again for values in the middle of each bin with a line segment in it.
    This gives us the coordinates for each relevant bin, given in line_bins.
    Finally we compute the length of each segment and add said value to the relevant box in our K-space."""
    v_prime = pks + vr_primestack*rhats
    line_bins = np.floor((v_prime-vmins)/ dv)
      
    line_bins = line_bins.astype(int)
    
    line_len = vr[1:]-vr[:-1] #Gets the length of each line
    non_zero = np.nonzero(line_len)
    line_len = line_len[non_zero] #Removes duplicate intersections
    line_bins = line_bins[non_zero]
    
    Kvals = line_len/(dvx*dvy*dvz)

    K = sp.COO(line_bins.T,Kvals,shape=(nx,ny,nz))
    
    return K

def sec_der1(phi,sigma2,dv):
    
    """Estimates the second deriative for ln(f(v_l)) given a sample of stars (eq. 30 of D98).
    Takes contributions at the phi values of adjacent bins for each bin l.
    
    We create a new, larger box with dimensions n+2 centred on our phi-space box.
    This allows us to disregard any issues at the bins at the boundaries of our phi-space."""
    
    nx, ny, nz = np.shape(phi)
    dvx, dvy, dvz = dv
    dv2 = dv**2
    sigma2x, sigma2y, sigma2z = sigma2
    
    nxx, nyy, nzz = nx+2, ny+2, nz+2 

    phip = np.zeros((nxx,nyy,nzz)) #new larger box
    
#     phip = sp.COO([],shape=(nxx,nyy,nzz))

    phip[1:-1,1:-1,1:-1] = phi #puts the phi-box in the centre of our larger box
    
    kappa = sigma2/dv2
    
    kappa_sum = -2*(sigma2x/dvx**2+sigma2y/dvy**2+sigma2z/dvx**2)
    
    """Here we compute the contributions from all the adjacent bins simultaneously.
    In every dimension we sum the phi values of box l-1 and l+1 and multiply with the relevant factor"""
    
    phi_fac = np.array([phip[0:nxx-2,1:-1,1:-1]+phip[2:nxx,1:-1,1:-1],
                           phip[1:-1,0:nyy-2,1:-1]+phip[1:-1,2:nyy,1:-1],
                           phip[1:-1,1:-1,0:nzz-2]+phip[1:-1,1:-1,2:nzz]])

    phi_arrx = (sigma2[0]/dv2[0])*phi_fac[0]
    phi_arry = (sigma2[1]/dv2[1])*phi_fac[1]
    phi_arrz = (sigma2[2]/dv2[2])*phi_fac[2]
    
    """We sum all contributions from adjacent boxes and finally add the terms for each box l. 
    Yields a box with the same dimensions as phi, containing the second derivative values for each bin."""
    
    phi_arr = phi_arrx+phi_arry+phi_arrz+kappa_sum*phi 
    
    phi_arrs = sp.COO.from_numpy(phi_arr)

    return phi_arrs

def get_negL1(phi,*args):
    
    """The function that we wish to optimise."""
    
    Kvals, N, alpha, dv, n, sigma2 = args
    nx, ny, nz = n
    dvx, dvy, dvz = dv
    
    """We regain the original shape of our phi guess and proceed to compute the various quantities needed from our functions."""
    
    phi_unr = np.reshape(phi,n)
    
    phixhi = sec_der1(phi_unr,sigma2,dv) #last term
    phixhi_sum = sp.COO.sum(phixhi**2)
    exphi = np.exp(phi_unr)
    
    exphi_sp = sp.COO.from_numpy(exphi)
    
    Kphi = exphi_sp*Kvals
    Kphiord = Kphi.reshape((Kphi.shape[0],nx*ny*nz)) #Order all Kphi values in 1D arrays for each star
    Kphi_sum_sp = sp.COO.sum(Kphiord,axis=1) #We compute the sum of exp(phi)*K(k|l) for each star
    Kphi_sum = sp.COO.todense(Kphi_sum_sp)
    
    notzero = Kphi_sum != 0
    Kphi_sum[notzero] = np.log(Kphi_sum[notzero]) #To make sure we don't get infinities
    Kphi_sum_tot = np.sum(Kphi_sum) #Gives the double sum in the first term
    
    L_tilde = Kphi_sum_tot/N - np.sum(exphi)-((alpha*dvx*dvy*dvz)/2)*phixhi_sum #eq. 31 in DB98
    
    negL = -1*L_tilde #Since we want to maximize L_tilde, we should minimize -L_tilde 
    
    return negL

def get_grad_negL(phi,*args):
    
    """In this function we compute the gradient of L. We compute the derivative for each cell and return a 
    1D array of length (nx*ny*nz)."""
    
    Kvals, N, alpha, dv, n, sigma2 = args
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    sigma2x, sigma2y, sigma2z = sigma2
    
    phi_unr = np.reshape(phi,n)
    exphir = np.exp(phi)
    exphi = exphir.reshape(phi_unr.shape)
    
    Kphi = exphi*Kvals
    Kphiord = Kphi.reshape(len(Kphi),nx*ny*nz) #Order all Kphi values in 1D arrays for each star
    Kphi_sum = np.sum(Kphiord,axis=1) #We compute the sum of exp(phi)*K(k|l) for each star
    Kphistack = np.zeros((N,nx*ny*nz))
    
    Kphistack[:] = Kphi_sum[:,None]
    
    #We first create an array of shape (N,nx*ny*nz) and then sum over all of the stars to obtain the first term
    
    K_term0 = Kphiord / Kphistack
    K_term = np.sum(K_term0,axis=0) #The final array with the first term for each cell
    
    kappa_sum = -2*(sigma2x/dvx**2+sigma2y/dvy**2+sigma2z/dvx**2)
    
    dphixhi = sec_der(phi_unr,sigma2,dv)*kappa_sum #last term
    dphixhi_rav = np.ravel(dphixhi) #We ravel to obtain a 1D array of length nx*ny*nz

    grad_L = (K_term/N-exphir-(alpha*dvx*dvy*dvz)*dphixhi_rav)
    
    return -1*grad_L

def max_L(alpha, pvals, rhatvals, vmin, dv, n,v0_guess=[],disp_guess=[], disp=1):
    
    """Function that employs scipy.optimize.fmin_cg to maximise the function get_negL().
    It takes guesses of the distribution (currently only supports Gaussian guesses) and the relevant data from the
    star sample for which the velocity distribution is to be estimated."""
    
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    
    N = len(pvals)
    
    sigma2, vmean = calc_sigma2(pvals,rhatvals,give_vmean=True) 
    
    sigma = np.sqrt(sigma2)
    
    if v0_guess == []:
        v0_guess = vmean
    if disp_guess == []:
        disp_guess = sigma
    
    phi0 = phi_guess(v0_guess,disp_guess,vmin,dv,n) #We obtain phi given our initial guess of the velocity distribution
    
    Klist = []
    
    for i in range(N):
        
        K = calc_K(pvals[i],rhatvals[i],vmin,dv,n)
        
        Klist.append(K)
        
    Kvals = sp.stack(Klist,axis=0)
        
    args = (Kvals, N, alpha, dv, n, sigma2)
    
    phi0r = np.ravel(phi0) #fmin_cg only takes one-dimensional inputs for the initial guess
    
    mxl, phi_all = fmin_cg(get_negL, phi0r, fprime = get_grad_negL, gtol=5e-4, args=args, retall=True,disp=disp)
    
    mxlnew = mxl.reshape(n)

    return mxlnew, phi_all