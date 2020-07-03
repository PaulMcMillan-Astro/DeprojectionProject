import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import astropy.units as u
import astropy.coordinates as coord
import scipy.stats as st
from astropy.table import Table
from astropy.io import ascii
from scipy import sparse as scisp
from scipy.optimize import fmin_cg
# Version adjusted by Paul to fix the edge issues with the gradient

def calc_p_rhat(sample):

    """Function that computes the values of rhat and the tangential velocity for a given sample"""

#Oort constant values from Bovy (2018)
    A = (15.3*(u.km/(u.s*u.kpc))).to(1/u.yr)
    B = (-11.9*(u.km/(u.s*u.kpc))).to(1/u.yr)

    bvals = sample.b.to(u.deg)
    lvals = sample.l.to(u.deg)

    mul_obs = sample.pm_l_cosb.to(1/u.yr,equivalencies = u.dimensionless_angles())
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

    b = np.deg2rad(bvals).value
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

    return pvals, rhatvals

@profile
def max_L(alpha, pvals, rhatvals, vmin, dv, method, n, phi0_guess = [],v0_guess=[],disp_guess=[], disp=1, noniso=False):

    """Function that employs scipy.optimize.fmin_cg to maximise the function get_negL().
    It takes guesses of the distribution (currently only supports Gaussian guesses) and the relevant data from the
    star sample for which the velocity distribution is to be estimated."""

    dvx, dvy, dvz = dv
    nx, ny, nz = n

    N = len(pvals)

    sigma2, vmean = calc_sigma2(pvals,rhatvals,True,noniso)

    sigma = np.sqrt(sigma2)

    if phi0_guess == []:
        phi0 = phi_guess(v0_guess,disp_guess,vmin,dv,n) #We obtain phi given our initial guess of the velocity distribution
        phi0[1:-1,1:-1,1:-1] += np.random.uniform(-1,1,size=(n-2))*5
    else:
        phi0 = phi0_guess
        
    Ks = calc_K(pvals[0],rhatvals[0],vmin,dv,n)
    K0 = np.ravel(Ks)
    
    if method == 1:
        ## SCIPY VSTACK METHOD:
        Kvals = scisp.coo_matrix(K0)
        
        for i in range(1,N): #Loop that yield a sparse array of N K-values
    
            K = np.ravel(calc_K(pvals[i],rhatvals[i],vmin,dv,n))
            K_coo = scisp.coo_matrix(K)
            
            Kvals = scisp.vstack((Kvals,K_coo))

    elif method == 2:
    ## SCIPY PREALLOCATE METHOD
        Kvals = scisp.dok_matrix((N,len(K0)))
        Kvals[0] = scisp.coo_matrix(K0)
        
        for i in range(1,N):
            K = np.ravel(calc_K(pvals[i],rhatvals[i],vmin,dv,n))
            Kvals[i] = scisp.coo_matrix(K)

    
    Kvals = scisp.csc_matrix(Kvals)
    
    """
    INCLUDE THIS FOR NORMAL RUN 
    args = (Kvals, N, alpha, dv, n, sigma2)

    phi0r = np.ravel(phi0) #fmin_cg only takes one-dimensional inputs for the initial guess

    mxl, phi_all = fmin_cg(get_negL, phi0r, fprime = get_grad_negL, gtol=1e-4, args=args, retall=True,disp=disp)

#    mxl, phi_all = fmin_cg(get_negL, phi0r, gtol=5e-4, args=args, retall=True,disp=disp)

    mxlnew = mxl.reshape(n)

    return mxlnew, phi_all
    """
    return

def calc_sigma2(pvals,rhat,give_vmean=False,noniso=False):

    """Function that applies equation 12 of DB98 for a set of stars from their tangential velocities and unit vectors.
    Returns the velocity dispersion tensor.

    pvals: array of the pk vectors for all N stars in our sample. Should have dims (N,3)
    rhat: array of N unit vectors, one for each star
    give_vmean: if True, returns the computed mean velocity vector value for the given sample
    noniso: if True, we no longer assume that the sample is isotropic"""

    pmean = np.mean(pvals, axis=0)

    rhat_outer = rhat[:,:,None]*rhat[:,None,:] #Fast way of getting the outer product for each rhat with itself.

    iden = np.identity(3)

    A0 = np.zeros((len(rhat_outer),3,3))

    A0[:] = iden

    A = A0-rhat_outer

    A_mean = np.mean(A,axis=0)
    A_mean_inv = np.linalg.inv(A_mean)
    v_mean = np.dot(A_mean_inv, pmean)

    pp = pvals - np.dot(A,v_mean) #Computes p' from equation (6) in DB98

    pp2mean = np.mean(np.square(pp),axis=0)
    
    if noniso:
        #In our main method, we rely on built-in tensor algebra functions
        #that perform these computations for us. We set up B by computing the
        #mean of the outer product of the peculiar tangential velocities
        #p' with itself.
        
        B = np.mean(pp[:,:,None]*pp[:,None,:],axis=0)
        
        #Next, we use the tensordot function of numpy to obtain T for each star.
        #We resort to a simple list comprehension operation to obtain these 
        #values
        
        #We could alternatively use np.einsum here
        T_stack = np.asarray([np.tensordot(i,i.T,axes=0) for i in A])
        
        T = np.mean(T_stack,axis=0)
        
        #With the non-singular A tensor at hand we can easily solve for D
        #and obtain the velocity dispersion tensor
        
        D = np.linalg.tensorsolve(T,B,(0,2))
        
        sigma2 = np.diag(D)
        
    else:
        
        B = np.array([[9,-1,-1],[-1,9,-1],[-1,-1,9]])

        sigma2 = (3/14)*np.dot(B,pp2mean) #The velocity dispersion tensor

    if give_vmean == True:

        return sigma2, v_mean

    else:

        return sigma2

def phi_guess(v0,disp0,vmin,dv,n):

    """Provides an initial guess of the phi values in each bin given an assumed distribution f(v).
    For now only allows for a Gaussian type guess given arrays with mean velocities and dispersions for each dimension.

    v0: A vector containing all the mean velocity values for the Gaussian distribution

    disp0: A vector with the velocity dispersions for each Gaussian in x, y, z

    vmin: The anchor point of our velocity space box

    """

    vxmin, vymin, vzmin = vmin
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    v0x, v0y, v0z = v0
    dispx, dispy, dispz = disp0

    vxmax, vymax, vzmax = vxmin+nx*dvx,vymin+ny*dvy,vzmin+nz*dvz

    vx_bins = np.arange(vxmin, vxmax+dvx, dvx)
    vy_bins = np.arange(vymin, vymax+dvy, dvy)
    vz_bins = np.arange(vzmin, vzmax+dvz, dvz)

    vxc = (vx_bins[1:]+vx_bins[:-1])/2
    vyc = (vy_bins[1:]+vy_bins[:-1])/2
    vzc = (vz_bins[1:]+vz_bins[:-1])/2

    """Given the velocities of each bin we compute the 3D Gaussian value."""

    gx = st.norm(v0x,dispx) #Gets the normal distribution for dim x given the mean velocities v0x and dispersion dispx
    gy = st.norm(v0y,dispy)
    gz = st.norm(v0z,dispz)

    gxp = gx.logpdf(vxc) #Computes the log of the probability of our bin values being in the Gaussian
    gyp = gy.logpdf(vyc)
    gzp = gz.logpdf(vzc)

    phi = np.meshgrid(gxp,gyp,gzp,indexing='ij')
    #The meshgrid function couples each phi value for our 3D scenario to the relevant box

    phi_sum = np.sum(phi,axis=0) #We sum the contributions to phi from x,y,z in each box

    return phi_sum
    
def calc_K(pk,rhat,vmin,dv,n):
    '''Calculate the values of K simultaneously for all bins for a given star with tangential velocity
    vector pk and a unit vector rhat'''

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

    K = np.zeros((nx,ny,nz))

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

    K[line_bins[:,0],line_bins[:,1],line_bins[:,2]] = line_len/(dvx*dvy*dvz)

    return K

def get_negL(phi,*args):

    """The function that we wish to optimise. Corresponds to eq. 31 in D98.

    N: Number of stars in our sample

    Kvals: Array of dimensions (N,nx,ny,nz) containing the K-values for each star in our sample.

    alpha: Smoothing parameter that can be found using the function opt_alpha
    """

    Kvals, N, alpha, dv, n, sigma2 = args
    nx, ny, nz = n
    dvx, dvy, dvz = dv

    """We regain the original shape of our phi guess and proceed to compute the
    various quantities needed from our functions."""

    phi_unr = np.reshape(phi,n)

    phixhi = sec_der(phi_unr,sigma2,dv)**2 #last term
    phixhi_sum = np.sum(phixhi)
    exphir = np.exp(phi)

    """We use sparse matrices in our computation of L_tilde because our K- arrays have low
    density. Therefore we use the scipy.sparse package to convert our arrays to sparse arrays.
    Using coo-matrices is faster when building the matrix, but the csc-matrices are faster for
    arithmetic operations"""

    exphi_coo = scisp.coo_matrix(exphir)
    exphi_csc = exphi_coo.tocsc()

    Kphi = Kvals.multiply(exphi_csc) #Order all Kphi values in 1D arrays for each star
    Kphi_sum = Kphi.sum(axis=1) #We compute the sum of exp(phi)*K(k|l) for each star

    Kphi_sumlog = np.log(Kphi_sum[Kphi_sum != 0]) #To make sure we don't get infinities
    Kphi_sum_tot = np.sum(Kphi_sumlog) #Gives the double sum in the first term

    L_tilde = Kphi_sum_tot/N - exphi_csc.sum() - ((alpha*dvx*dvy*dvz)/2)*phixhi_sum #eq. 31 in DB98

    negL = -1*L_tilde #Since we want to maximize L_tilde, we should minimize -L_tilde

    return negL

def sec_der(phi,sigma2,dv):

    """Estimates the second deriative for ln(f(v_l)) given a sample of stars (eq. 30 of D98).
    Takes contributions at the phi values of adjacent bins for each bin l.

    phi: An array of dimensions (nx,ny,nz) that gives the logarithm of the probability f(v) for a given
        velocity to be in the different cells of v-space

    sigma2: The velocity dispersion tensor, should be a vector

    dv: The dimensions of each cell in km/s


    We create a new, larger box with dimensions n+2 centred on our phi-space box.
    This allows us to disregard any issues at the bins at the boundaries of our phi-space."""

    nx, ny, nz = phi.shape
    dvx, dvy, dvz = dv
    dv2 = dv**2
    sigma2x, sigma2y, sigma2z = sigma2

    nxx, nyy, nzz = nx+2, ny+2, nz+2

    phip = np.zeros((nxx,nyy,nzz)) #new larger box

    phip[1:-1,1:-1,1:-1] = phi #puts the phi-box in the centre of our larger box

    """Here we compute the contributions from all the adjacent bins simultaneously.
    In every dimension we sum the phi values of box l-1 and l+1 and multiply with the relevant factor"""
    phi_fac = np.array([phip[0:nxx-2,1:-1,1:-1]+phip[2:nxx,1:-1,1:-1],
                           phip[1:-1,0:nyy-2,1:-1]+phip[1:-1,2:nyy,1:-1],
                           phip[1:-1,1:-1,0:nzz-2]+phip[1:-1,1:-1,2:nzz]])
    
    phi_arrx = (sigma2[0]/dv2[0])*(phi_fac[0]-2*phi)
    phi_arry = (sigma2[1]/dv2[1])*(phi_fac[1]-2*phi)
    phi_arrz = (sigma2[2]/dv2[2])*(phi_fac[2]-2*phi)

    """We sum all contributions from adjacent boxes and finally add the terms for each box l.
    Yields a box with the same dimensions as phi, containing the second derivative values for each bin."""

    phi_arr = phi_arrx+phi_arry+phi_arrz

    phi_arr[0,:,:] = phi_arr[:,0,:] = phi_arr[:,:,0] = 0 #Due to lack of smoothness at edges, we set these values to 0
    phi_arr[-1,:,:] = phi_arr[:,-1,:] = phi_arr[:,:,-1] = 0

    return phi_arr


def get_grad_negL(phi,*args):

    """In this function we compute the gradient of L. We compute the derivative for each cell and return a
    1D array of length (nx*ny*nz).

    args: see get_L

    """

    Kvals, N, alpha, dv, n, sigma2 = args
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    sigma2x, sigma2y, sigma2z = sigma2

    phi_unr = np.reshape(phi,n)

    exphir = np.exp(phi)

    exphi_coo = scisp.coo_matrix(exphir)
    exphi_csc = exphi_coo.tocsc()

    Kphi = exphi_csc.multiply(Kvals)

    Kphi_sum = Kphi_suminv = Kphi.sum(axis=1)

    Kphi_suminv[Kphi_sum.nonzero()] = 1/Kphi_sum[Kphi_sum.nonzero()] #We compute the sum of exp(phi)*K(k|l) for each star

    K_term0 = Kphi.multiply(Kphi_suminv)
    K_term = K_term0.sum(axis=0) #The final array with the first term for each cell

    dphixhi = grad_sec_der(phi_unr,sigma2,dv,n)
    dphixhi_rav = dphixhi.reshape((phi.shape[0],1)).T
    dphixhi_coo = scisp.coo_matrix(dphixhi_rav)
    dphixhi_csc = dphixhi_coo.tocsc()

    grad_L = np.asarray(K_term/N-exphi_csc-((alpha*dvx*dvy*dvz)/2)*dphixhi_csc).reshape(nx*ny*nz,)

    return -1*grad_L


def grad_sec_der(phi,*args):

    sigma2, dv, n = args

    nx, ny, nz = n
    dv2 = dv**2


    """We sum all contributions from adjacent boxes and finally add the terms for each box l.
    Yields a box with the same dimensions as phi, containing the second derivative values for each bin."""

    eta = sec_der(phi,sigma2,dv)

    eta_big = np.zeros((nx+2,ny+2,nz+2))
    eta_big[1:-1,1:-1,1:-1] = eta
    eta_fac_rs = np.array([eta_big[0:-2,1:-1,1:-1]+eta_big[2:,1:-1,1:-1],
                           eta_big[1:-1,0:-2,1:-1]+eta_big[1:-1,2:,1:-1],
                           eta_big[1:-1,1:-1,0:-2]+eta_big[1:-1,1:-1,2:]])


    eta_arrx = (sigma2[0]/dv2[0])*(eta_fac_rs[0]-2*eta)
    eta_arry = (sigma2[1]/dv2[1])*(eta_fac_rs[1]-2*eta)
    eta_arrz = (sigma2[2]/dv2[2])*(eta_fac_rs[2]-2*eta)

    eta_arr = 2*(eta_arrx+eta_arry+eta_arrz)

    return np.ravel(eta_arr)