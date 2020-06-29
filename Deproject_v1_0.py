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

    pvals = s*np.array([-sinl*mul - cosl*sinb*mub,
                     cosl*mul - sinl*sinb*mub,
                     cosb*mub])/u.yr

    rhatvals = np.array([cosb*cosl, cosb*sinl, sinb]).T
    pvals = pvals.to(u.km/u.s).value.T

    return pvals, rhatvals

def model_sample(N):
    """Generates a simple model solar neighbourhood star sample in a Galactic frame of reference assuming a
    velocity distribution that is a sum of three predefined Gaussians.
    
    If one wishes to use other mean velocities, change these in mu0, mu1 and mu2, while the dispersions are changed
    in disp0, disp1 and disp2. The relative weights, w1, w2, w3 determine how many of the stars that belong to
    each Gaussian. As for now, they're just set to be ~1/3.
    
    The stars have distances from the Sun that are in the range [10,100] pc.
    
    Takes the following arugments:
    
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
    
    psx = (np.random.rand(N)*(xmax-xmin)+xmin)*u.pc
    psy = (np.random.rand(N)*(ymax-ymin)+ymin)*u.pc
    psz = (np.random.rand(N)*(zmax-zmin)+zmin)*u.pc
    
    from_g0 = np.where(w_sample < w0) #We get the indices of the stars that belong to the first Gaussian
    from_g1 = np.where((w0 <= w_sample) & (w_sample < (w0+w1)))
    from_g2 = np.where(w_sample > (1-w2))
    
    scale0 = np.random.randn(len(from_g0[0]),3)
    scale1 = np.random.randn(len(from_g1[0]),3)
    scale2 = np.random.randn(len(from_g2[0]),3)
    
    psvels = np.zeros((N,3))
    
    psvels[from_g0] = mu0+scale0*disp0 #We exchange our empty velocity values with the ones obtained from each Gaussian
    psvels[from_g1] = mu1+scale1*disp1
    psvels[from_g2] = mu2+scale2*disp2
    
    psvx, psvy, psvz = psvels.T
    
    #We use Astropy's coord class which makes it easy to keep track of units and conversions
    
    psample = coord.Galactic(u=psx,v=psy,w=psz,U=psvx*(u.km/u.s),
                            V=psvy*(u.km/u.s),W=psvz*(u.km/u.s),representation_type=coord.CartesianRepresentation,differential_type=coord.CartesianDifferential)
    
    psample.set_representation_cls(coord.SphericalRepresentation,coord.SphericalCosLatDifferential)
    
    return psample

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

def calc_sigma2(pvals,rhat,give_vmean=False):

    """Function that applies equation 12 of DB98 for a set of stars from their tangential velocities and unit vectors.
    Returns the velocity dispersion tensor.

    pvals: array of the pk vectors for all N stars in our sample. Should have dims (N,3)
    rhat: array of N unit vectors, one for each star
    give_vmean: if True, returns the computed mean velocity vector value for the given sample"""

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

    B = np.array([[9,-1,-1],[-1,9,-1],[-1,-1,9]])

    sigma2 = (3/14)*np.dot(B,pp2mean) #The velocity dispersion tensor

    if give_vmean == True:

        return sigma2, v_mean

    else:

        return sigma2

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

def max_L(alpha, pvals, rhatvals, vmin, dv, n,phi0_guess = [],v0_guess=[],disp_guess=[], disp=1):

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

    if phi0_guess == []:
        phi0 = phi_guess(v0_guess,disp_guess,vmin,dv,n) #We obtain phi given our initial guess of the velocity distribution
    else:
        phi0 = phi0_guess

    K0 = np.ravel(calc_K(pvals[0],rhatvals[0],vmin,dv,n))

    Kvals = scisp.coo_matrix(K0)

    for i in range(1,N): #Loop that yield a sparse array of N K-values

        K = np.ravel(calc_K(pvals[i],rhatvals[i],vmin,dv,n))
        K_coo = scisp.coo_matrix(K)

        Kvals = scisp.vstack((Kvals,K_coo))

    Kvals_csc = Kvals.tocsc()

    args = (Kvals_csc, N, alpha, dv, n, sigma2)

    phi0r = np.ravel(phi0) #fmin_cg only takes one-dimensional inputs for the initial guess

    mxl, phi_all = fmin_cg(get_negL, phi0r, fprime = get_grad_negL, 
                            gtol=5e-4, args=args, retall=True,disp=disp)

#    mxl, phi_all = fmin_cg(get_negL, phi0r, gtol=5e-4, args=args, retall=True,disp=disp)

    mxlnew = mxl.reshape(n)

    return mxlnew, phi_all

def opt_alpha(alpha0,M,N,sample,vmin,dv,n,tol=0.01):
    
    """Function that finds the optimal value of alpha for a given sample of stars.
    Given an initial guess of alpha, alpha0, it will draw M samples of size N from the resulting
    distribution f(v) computed using max_L. For each sample we perform the maximisation scheme
    to find the f(v) distribution for a set of alpha values. The mean integrated square error (MISE) is then
    computed for each alpha value and the alpha with the lowest MISE, alpha_opt is then our new initial guess.
    
    This process repeats until the difference |alpha0-alpha_opt| falls below a certain threshold.
    
    alpha0: The initial guess of alpha
    M: The number of samples to compute the MISE for
    N: Number of stars in each sample
    pvals: Array with the tangential velocity vectors for all the  original stars in our sample
    rhatvals: Array with the unit vector for each sample star
    vmin: Vector indicating the anchor of v-space
    dv: The dims of each cell
    n: The dims of our box
    tol: The desired logarithmic tolerance for the minimisation scheme. Default value is 0.01.
    
    """
    
    pvals, rhatvals = calc_p_rhat(sample)
    
    #We want to compute the rhat values for every sample
    #So we just draw M sets of size N of the coordinates from the original
    #sample that we will use for each iteration of the sample.
    
    sample.set_representation_cls(coord.CartesianRepresentation)
    cartcoords = np.array([sample.u,sample.v,sample.w]).T
    coordinds = np.linspace(0,len(cartcoords)-1,len(cartcoords))
    
    diff = 1
    s = 0

    logalpha0 = np.log10(alpha0)
    logalphavals = np.linspace(logalpha0-2,logalpha0+2,10) #The initial set of alpha values

    vxmin, vymin, vzmin = vmin
    dvx, dvy, dvz = dv
    nx, ny, nz = n

    vxmax, vymax, vzmax = vxmin+nx*dvx,vymin+ny*dvy,vzmin+nz*dvz

    vx_bins = np.arange(vxmin, vxmax+dvx, dvx)
    vy_bins = np.arange(vymin, vymax+dvy, dvy)
    vz_bins = np.arange(vzmin, vzmax+dvz, dvz)

    vxc = (vx_bins[1:]+vx_bins[:-1])/2
    vyc = (vy_bins[1:]+vy_bins[:-1])/2
    vzc = (vz_bins[1:]+vz_bins[:-1])/2

    vxx, vyy, vzz = np.meshgrid(vxc,vyc,vzc,indexing='ij') #The centre values of each cell in our v-space box

    ind = np.indices((nx,ny,nz))

    rind = np.ravel_multi_index(ind,(nx,ny,nz)) #An array containing the 3D coordinates of our box

    rrind = np.ravel(rind)

    while diff>=tol:

        phi0, phiall = max_L(alpha0, pvals, rhatvals, vmin, dv, n, disp=0)

        fv0 = np.exp(phi0) #The pdf f(v) given our inital guess alpha0
        fv0s = np.sum(fv0)
        prob = np.ravel(fv0/fv0s) #The normalised probability

        smp = np.random.choice(rrind,(M,N),p=prob) #Creation of M samples of size n given f(v)
        smpcoordinds = np.random.choice(coordinds,(M,N)).astype(int) #We also draw random positions for these stars

        smpcoords = cartcoords[smpcoordinds]

        smpx, smpy, smpz = np.asarray(np.unravel_index(smp,(nx,ny,nz)))

        vxvals = vxx[smpx,smpy,smpz].T
        vyvals = vyy[smpx,smpy,smpz].T
        vzvals = vzz[smpx,smpy,smpz].T

        smpvels = np.asarray([vxvals, vyvals, vzvals]).T

        ise = np.zeros((M,len(logalphavals))) #container for the integrated square errors

        for i in range(M):
            
            smpvx, smpvy, smpvz = smpvels[i].T #For every pseudosample we get velocities
            
            coordx, coordy, coordz = smpcoords[i].T#... as well as coordinates
            
            psample = coord.Galactic(u=coordx*u.pc,v=coordy*u.pc,w=coordz*u.pc,U=smpvx*(u.km/u.s),
                            V=smpvy*(u.km/u.s),W=smpvz*(u.km/u.s),representation_type=coord.CartesianRepresentation,
                            differential_type=coord.CartesianDifferential)
    
            psample.set_representation_cls(coord.SphericalRepresentation,coord.SphericalCosLatDifferential)
            
            pspvals, psrhatvals = calc_p_rhat(psample)

            for j in range(len(logalphavals)):

                alpha = 10**(logalphavals[j])

                phi, phiall = max_L(alpha, pspvals, psrhatvals, vmin, dv, n, disp=0)

                fv = np.exp(phi)

                ise[i][j] = np.sum((fv-fv0)**2)

        mise = np.mean(ise,axis=0)

        minmise = np.amin(mise)

        optind = np.argwhere(mise==minmise)[0][0] #The index of the optimimal alpha

        logalpha_opt = logalphavals[optind]

        diff = abs(logalpha0 - logalpha_opt)

        logalpha0 = logalpha_opt #We set the optimised alpha value to be our new initial guess
        
        """Since our initial guess is very broad and covers 10 orders of magnitude, we have to narrow it down.
        This is done by taking the alpha values to the right and left of our optimal value to define the new range.
        If the optimal value is at the edges of our array, then we take the missing upper or lower bound
        to be 2 magnitudes larger than alpha_opt"""

        if optind == 0:
            lower = logalpha0-2
            upper = logalphavals[optind+1]
        elif optind == (len(logalphavals)-1):
            lower = logalphavals[optind-1]
            upper = logalpha0+2
        else:
            lower = logalphavals[optind-1]
            upper = logalphavals[optind+1]

        logalphavals = np.linspace(lower,upper,10)

        s+=1

    alpha_fin = 10**(logalpha0)

    print('It took',s,'iterations')
    print('The optimal value for alpha is',alpha_fin)
    
    return