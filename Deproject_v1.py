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

def model_sample(N,v0,disp0): 

    """Generates a simple model solar neighbourhood star sample in a Galactic frame of reference assuming a 
    Gaussian velocity distribution. The stars have distances from the Sun that are in the range [10,100] pc.
    
    Takes the following arugments:
    
    N: Number of stars in the sample
    v0: 3d array specifying the mean velocities of the Gaussian distributions in x, y, z
    disp: 3d array specifying the velocity dispersions of the distributions in x, y, z"""
    
    xmax, ymax, zmax = np.array([100,100,100]) / np.sqrt(3)
    xmin, ymin, zmin = -xmax,-ymax,-zmax
    
    psx = (np.random.rand(N)*(xmax-xmin)+xmin)*u.pc
    psy = (np.random.rand(N)*(ymax-ymin)+ymin)*u.pc
    psz = (np.random.rand(N)*(zmax-zmin)+zmin)*u.pc
    
    scale = np.random.randn(N,3)
    psvels = v0 + scale*disp0
    
    psvx, psvy, psvz = psvels.T
    
    psample = coord.Galactic(u=psx,v=psy,w=psz,U=psvx*(u.km/u.s),
                             V=psvy*(u.km/u.s),W=psvz*(u.km/u.s),representation_type=coord.CartesianRepresentation,differential_type=coord.CartesianDifferential)
    
    psample.set_representation_cls(coord.SphericalRepresentation)
    
    return psample

def calc_K(pk,rhat,vmin,dv,n):
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
    
    Kvals = line_len/(dvx*dvy*dvz)

    K[line_bins[:,0],line_bins[:,1],line_bins[:,2]] = line_len/(dvx*dvy*dvz)

    return K

def calc_sigma2(pvals,rhat,give_vmean=False):
    
    """Function that applies equation 12 of DB98 for a set of stars from their tangential velocities and unit vectors.
    Returns the velocity dispersion tensor."""
    
    pmean = np.mean(pvals, axis=0)
    
    rhat_outer = rhat[:,:,None]*rhat[:,None,:] #Fast way of getting the outer product for each rhat with itself.

    iden = np.identity(3)
    
    A0 = np.zeros((len(rhat_outer),3,3))
    
    A0[:] = iden 
    
    A = A0-rhat_outer
    
    A_mean = np.mean(A,axis=0)
    A_mean_inv = np.linalg.inv(A_mean)
    v_mean = np.dot(A_mean_inv, pmean)
    
    pp = pvals - np.dot(A,v_mean)
    
    pp2mean = np.mean(np.square(pp),axis=0)
    
    B = np.array([[9,-1,-1],[-1,9,-1],[-1,-1,9]])
    
    sigma2 = (3/14)*np.dot(B,pp2mean)
    
    if give_vmean == True:
    
        return sigma2, v_mean
    
    else:
    
        return sigma2

def sec_der(phi,sigma2,dv):
    
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

    phip[1:-1,1:-1,1:-1] = phi #puts the phi-box in the centre of our larger box
    
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

    return phi_arr

def phi_guess(v0,disp0,vmin,dv,n):
    
    """Provides an initial guess of the phi values in each bin given an assumed distribution f(v). 
    For now only allows for a Gaussian type guess given arrays with mean velocities and dispersions for each dimension."""
    
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
    
    phi = np.meshgrid(gxp,gyp,gzp,indexing='ij') #The meshgrid function couples each phi value for our 3D scenario to the relevant box 
    
    phi_sum = np.sum(phi,axis=0) #We sum the contributions to phi from x,y,z in each box
    
    return phi_sum

def get_negL(phi,*args):
    
    """The function that we wish to optimise."""
    
    Kvals, N, alpha, dv, n, sigma2 = args
    nx, ny, nz = n
    dvx, dvy, dvz = dv
    
    """We regain the original shape of our phi guess and proceed to compute the various quantities needed from our functions."""
    
    phi_unr = np.reshape(phi,n)
    
    phixhi = sec_der(phi_unr,sigma2,dv) #last term
    phixhi_sum = np.sum(phixhi**2)
    exphir = np.exp(phi)
    
    exphi_coo = scisp.coo_matrix(exphir)
    exphi_csc = exphi_coo.tocsc()
    
    Kphi = Kvals.multiply(exphi_csc) #Order all Kphi values in 1D arrays for each star
    Kphi_sum = Kphi.sum(axis=1) #We compute the sum of exp(phi)*K(k|l) for each star
    
    Kphi_sumlog = np.log(Kphi_sum.data) #To make sure we don't get infinities
    Kphi_sum_tot = np.sum(Kphi_sumlog) #Gives the double sum in the first term
    
    L_tilde = Kphi_sum_tot/N - exphi_csc.sum() -((alpha*dvx*dvy*dvz)/2)*phixhi_sum #eq. 31 in DB98
    
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
    
    exphi_coo = scisp.coo_matrix(exphir)
    exphi_csc = exphi_coo.tocsc()
    
    Kphi = exphi_csc.multiply(Kvals)
    
    Kphi_sum = 1/Kphi.sum(axis=1) #We compute the sum of exp(phi)*K(k|l) for each star

    K_term0 = Kphi.multiply(Kphi_sum)
    K_term = K_term0.sum(axis=0) #The final array with the first term for each cell

    kappa_sum = -2*(sigma2x/dvx**2+sigma2y/dvy**2+sigma2z/dvx**2)

    dphixhi = sec_der(phi_unr,sigma2,dv)*kappa_sum #last term
    dphixhi_rav = dphixhi.reshape((phi.shape[0],1)) #We ravel to obtain a 1D array of length nx*ny*nz
    dphixhi_coo = scisp.coo_matrix(dphixhi_rav)
    dphixhi_csc = dphixhi_coo.tocsc().T

    grad_L = np.asarray(K_term/N-exphi_csc-(alpha*dvx*dvy*dvz)*dphixhi_csc).reshape(nx*ny*nz,)
    
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
    
    Kvals = np.zeros((N,nx*ny*nz))

    for i in range(N):

        K = np.ravel(calc_K(pvals[i],rhatvals[i],vmin,dv,n))

        Kvals[i] = K

    Kvals_coo = scisp.coo_matrix(Kvals)
    Kvals_csc = Kvals_coo.tocsc()
        
    args = (Kvals_csc, N, alpha, dv, n, sigma2)
    
    phi0r = np.ravel(phi0) #fmin_cg only takes one-dimensional inputs for the initial guess
    
    mxl, phi_all = fmin_cg(get_negL, phi0r, fprime = get_grad_negL, gtol=5e-4, args=args, retall=True,disp=disp)
    
    mxlnew = mxl.reshape(n)

    return mxlnew, phi_all

def opt_alpha(alpha0,M,N,pvals,rhatvals,vmin,dv,n,tol=0.01):
    
    diff = 1
    s = 0

    logalpha0 = np.log10(alpha0)
    logalphavals = np.linspace(logalpha0-5,logalpha0+5,10) #logarithm

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

    vxx, vyy, vzz = np.meshgrid(vxc,vyc,vzc,indexing='ij')

    ind = np.indices((nx,ny,nz))

    rind = np.ravel_multi_index(ind,(nx,ny,nz))

    rrind = np.ravel(rind)

    while diff>=tol:

        phi0, phiall = max_L(alpha0, pvals, rhatvals, vmin, dv, n, disp=0)

        fv0 = np.exp(phi0)
        fv0s = np.sum(fv0)
        prob = np.ravel(fv0/fv0s)

        smp = np.random.choice(rrind,(M,N),p=prob)

        smpx, smpy, smpz = np.asarray(np.unravel_index(smp,(nx,ny,nz)))

        vxvals = vxx[smpx,smpy,smpz].T
        vyvals = vyy[smpx,smpy,smpz].T
        vzvals = vzz[smpx,smpy,smpz].T

        smpvals = np.asarray([vxvals, vyvals, vzvals]).T

        ise = np.zeros((M,len(logalphavals)))

        for i in range(M):

            for j in range(len(logalphavals)):

                alpha = 10**(logalphavals[j])

                phi, phiall = max_L(alpha, smpvals[i], rhatvals, vmin, dv, n, disp=0)

                fv = np.exp(phi)

                ise[i][j] = np.sum((fv-fv0)**2)

        mise = np.mean(ise,axis=0)

        minmise = np.amin(mise)

        optind = np.argwhere(mise==minmise)[0][0]

        logalpha_opt = logalphavals[optind]

        diff = abs(logalpha0 - logalpha_opt)

        logalpha0 = logalpha_opt

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

    print("The run took", time.time() - start_time, 's')
    print('It took',s,'iterations')
    print('The optimal value for alpha is',alpha_fin)
    
    return