"""Test that for a given model sample checks if the distribution of L is smooth for a number of guesses of the dispersion."""
from Deproject_v0 import *


def test_L(plane, alpha, v0, disp0, pvals, rhatvals, vmin, dv, n):
    
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    disp0x, disp0y, disp0z = disp0
    N = len(pvals)
    
    sigma2 = calc_sigma2(pvals,rhatvals) 
    
    Kvals = np.zeros((N,nx,ny,nz))
    
    for i in range(N):
        K = calc_K(pvals[i],rhatvals[i],vmin,dv,n)
        Kvals[i] += K 
    
    """We generate a set of sigma values for each dimension that will serve as basis for different gueses of phi"""
    
    step = 5
    
    dispx_guess = np.arange(0.6*disp0x,1.4*disp0x)
    dispy_guess = np.arange(0.6*disp0y,1.4*disp0y)
    dispz_guess = np.arange(0.6*disp0z,1.4*disp0z)
    
    if plane == 'xy':
        trueind = np.where(disp0z == dispz_guess)
        disp10, disp20 = disp0x, disp0y
        disp1, disp2 = dispx_guess, dispy_guess
        disp3 = dispz_guess[trueind]
        n1, n2 = nx, ny
        dx, dy = dvx, dvy
        x0, y0 = min(dispx_guess), min(dispy_guess)
        x1, y1 = max(dispx_guess), max(dispy_guess)
        xlab, ylab = '$\sigma_x$', '$\sigma_y$'
    elif plane == 'yz':
        trueind = np.where(disp0x == dispz_guess)
        disp10, disp20 = disp0y, disp0z
        disp1, disp2 = dispy_guess, dispz_guess
        disp3 = dispx_guess[trueind]
        n1, n2 = ny, nz
        dx, dy = dvy, dvz
        x0, y0 = min(dispy_guess), min(dispz_guess)
        x1, y1 = max(dispy_guess), max(dispz_guess)
        xlab, ylab = '$\sigma_y$', '$\sigma_z$'
    elif plane == 'xz':
        trueind = np.where(disp0y == dispz_guess)
        disp10, disp20 = disp0x, disp0z
        disp1, disp2 = dispx_guess, dispz_guess
        disp3 = dispy_guess[trueind]
        n1, n2 = nx, nz
        dx, dy = dvx, dvz
        x0, y0 = min(dispx_guess), min(dispz_guess)
        x1, y1 = max(dispx_guess), max(dispz_guess)
        xlab, ylab = '$\sigma_x$', '$\sigma_z$'
    
    L_vals = np.zeros((len(disp1),len(disp2)))
    
    for i in range(len(disp1)):
        
        for j in range(len(disp2)):
            
            disp_guess = np.array([disp1[i],disp2[j],disp3])
    
            phi0 = phi_guess(v0,disp_guess,vmin,dv,n)

            phi0r = np.ravel(phi0)

            L_vals[i][j] += -get_L(phi0r, Kvals, N, alpha, dv, n, sigma2)
            
    fig, ax = plt.subplots(figsize=(8,6))
    
    dispg1, dispg2 = np.meshgrid(disp1,disp2)
    
    ax.set_title('Contour plot of $\mathcal{L}$ values')
    ax.set_xlim(x0,x1)
    ax.set_ylim(y0,y1)
    ax.set_xlabel(xlab+' [km s$^{-1}$]',size='large')
    ax.set_ylabel(ylab+' [km s$^{-1}$]',size='large')
    
    #contf = ax.contourf(dispg1,dispg2,L_vals.T,cmap = plt.cm.get_cmap('jet'),vmin=-40,vmax=10)
    
    extent = [disp1[0], disp1[-1], disp2[0], disp2[-1]]
    im = ax.imshow(L_vals.T,origin='lower',interpolation='bilinear',cmap = plt.cm.get_cmap('jet'),extent=extent)
    
    cb = plt.colorbar(im,orientation='vertical', extend='max')
    cb.set_label('Value of $\mathcal{L}$',size='large')
    
    ax.scatter(disp10,disp20,c='black')
    
    plt.show()
    
    return

def sanity_check(pvals,rhatvals,phi,vmin,dv,n):
    
    """Test function that for a given pseudosample computes the Gaussian distribution and the corresponding distribution
    found through maximization using the conjugate gradient algorithm."""
    
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    vxmin, vymin, vzmin = vmin
    vxmax, vymax, vzmax = vxmin+nx*dvx,vymin+ny*dvy,vzmin+nz*dvz
    
    vx_bins = np.arange(vxmin, vxmax+dvx, dvx)
    vy_bins = np.arange(vymin, vymax+dvy, dvy)
    vz_bins = np.arange(vzmin, vzmax+dvz, dvz)
    
    vxc = (vx_bins[1:]+vx_bins[:-1])/2
    vyc = (vy_bins[1:]+vy_bins[:-1])/2
    vzc = (vz_bins[1:]+vz_bins[:-1])/2

    sigma2, vmean0 = calc_sigma2(pvals,rhatvals,True)
    
    sigma = np.sqrt(sigma2)
    
    phi0 = phi_guess(vmean0,sigma,vmin,dv,n)
    fv0 = np.exp(phi0)
    
    fv = np.exp(phi)
    
    """Here we make sure to multipy the probability distribution f(v) with the centre values of the velocity bins in
    each dimension. """
    
    vxmean = np.sum(vxc[None,:,None]*fv)/np.sum(fv)
    vymean = np.sum(vyc[:,None,None]*fv)/np.sum(fv) 
    vzmean = np.sum(vzc[None,None,:]*fv)/np.sum(fv)
    vmean = np.array([vxmean,vymean,vzmean])
    
    dispx = np.sum((vxc[None,:,None]**2*fv))/np.sum(fv)-vxmean**2
    dispy = np.sum((vyc[:,None,None]**2*fv))/np.sum(fv)-vymean**2
    dispz = np.sum((vzc[None,None,:]**2*fv))/np.sum(fv)-vzmean**2
    disp = np.sqrt(np.array([dispx,dispy,dispz]))
    
    #vxmean = np.sum(vxcg*fv)/np.sum(fv)
    #vymean = np.sum(vycg*fv)/np.sum(fv)
    #vzmean = np.sum(vzcg*fv)/np.sum(fv)
    #vmean = np.array([vxmean,vymean,vzmean])
    
    #dispx = np.sum((vxcg**2*fv))/np.sum(fv)-vxmean**2
    #dispy = np.sum((vycg**2*fv))/np.sum(fv)-vymean**2
    #dispz = np.sum((vzcg**2*fv))/np.sum(fv)-vzmean**2
    #disp = np.sqrt(np.array([dispx,dispy,dispz]))
    
    return vmean0, sigma, vmean, disp
	
def sanity_check(pvals,rhatvals,phi,vmin,dv,n):
    
    """Test function that for a given pseudosample computes the Gaussian distribution and the corresponding distribution
    found through maximization using the conjugate gradient algorithm."""
    
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    vxmin, vymin, vzmin = vmin
    vxmax, vymax, vzmax = vxmin+nx*dvx,vymin+ny*dvy,vzmin+nz*dvz
    
    vx_bins = np.arange(vxmin, vxmax+dvx, dvx)
    vy_bins = np.arange(vymin, vymax+dvy, dvy)
    vz_bins = np.arange(vzmin, vzmax+dvz, dvz)
    
    vxc = (vx_bins[1:]+vx_bins[:-1])/2
    vyc = (vy_bins[1:]+vy_bins[:-1])/2
    vzc = (vz_bins[1:]+vz_bins[:-1])/2

    sigma2, vmean0 = calc_sigma2(pvals,rhatvals,True)
    
    sigma = np.sqrt(sigma2)
    
    phi0 = phi_guess(vmean0,sigma,vmin,dv,n)
    fv0 = np.exp(phi0)
    
    fv = np.exp(phi)
    
    """Here we make sure to multipy the probability distribution f(v) with the centre values of the velocity bins in
    each dimension. """
    
    vxmean = np.sum(vxc[None,:,None]*fv)/np.sum(fv)
    vymean = np.sum(vyc[:,None,None]*fv)/np.sum(fv) 
    vzmean = np.sum(vzc[None,None,:]*fv)/np.sum(fv)
    vmean = np.array([vxmean,vymean,vzmean])
    
    dispx = np.sum((vxc[None,:,None]**2*fv))/np.sum(fv)-vxmean**2
    dispy = np.sum((vyc[:,None,None]**2*fv))/np.sum(fv)-vymean**2
    dispz = np.sum((vzc[None,None,:]**2*fv))/np.sum(fv)-vzmean**2
    disp = np.sqrt(np.array([dispx,dispy,dispz]))
    
    return vmean0, sigma, vmean, disp

"""Test that for a given model sample checks if the distribution of L is smooth for a number of guesses of the dispersion."""

def test_L(plane, alpha, v0, disp0, pvals, rhatvals, vmin, dv, n):
    
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    disp0x, disp0y, disp0z = disp0
    N = len(pvals)
    
    sigma2 = calc_sigma2(pvals,rhatvals) 
    
    Kvals = np.zeros((N,nx,ny,nz))
    
    for i in range(N):
        K = calc_K(pvals[i],rhatvals[i],vmin,dv,n)
        Kvals[i] += K 
    
    """We generate a set of sigma values for each dimension that will serve as basis for different gueses of phi"""
    
    step = 5
    
    dispx_guess = np.arange(0.6*disp0x,1.4*disp0x)
    dispy_guess = np.arange(0.6*disp0y,1.4*disp0y)
    dispz_guess = np.arange(0.6*disp0z,1.4*disp0z)
        
    #disp_xx, disp_yy, disp_zz = np.meshgrid(dispx_guess,dispy_guess,dispz_guess)
    
    if plane == 'xy':
        trueind = np.where(disp0z == dispz_guess)
        disp10, disp20 = disp0x, disp0y
        disp1, disp2 = dispx_guess, dispy_guess
        disp3 = dispz_guess[trueind]
        n1, n2 = nx, ny
        dx, dy = dvx, dvy
        x0, y0 = min(dispx_guess), min(dispy_guess)
        x1, y1 = max(dispx_guess), max(dispy_guess)
        xlab, ylab = '$\sigma_x$', '$\sigma_y$'
    elif plane == 'yz':
        trueind = np.where(disp0x == dispz_guess)
        disp10, disp20 = disp0y, disp0z
        disp1, disp2 = dispy_guess, dispz_guess
        disp3 = dispx_guess[trueind]
        n1, n2 = ny, nz
        dx, dy = dvy, dvz
        x0, y0 = min(dispy_guess), min(dispz_guess)
        x1, y1 = max(dispy_guess), max(dispz_guess)
        xlab, ylab = '$\sigma_y$', '$\sigma_z$'
    elif plane == 'xz':
        trueind = np.where(disp0y == dispz_guess)
        disp10, disp20 = disp0x, disp0z
        disp1, disp2 = dispx_guess, dispz_guess
        disp3 = dispy_guess[trueind]
        n1, n2 = nx, nz
        dx, dy = dvx, dvz
        x0, y0 = min(dispx_guess), min(dispz_guess)
        x1, y1 = max(dispx_guess), max(dispz_guess)
        xlab, ylab = '$\sigma_x$', '$\sigma_z$'
    
    L_vals = np.zeros((len(disp1),len(disp2)))
    
    for i in range(len(disp1)):
        
        for j in range(len(disp2)):
            
            disp_guess = np.array([disp1[i],disp2[j],disp3])
    
            phi0 = phi_guess(v0,disp_guess,vmin,dv,n)

            phi0r = np.ravel(phi0)

            L_vals[i][j] += -get_L(phi0r, Kvals, N, alpha, dv, n, sigma2)
            
    fig, ax = plt.subplots(figsize=(8,6))
    
    dispg1, dispg2 = np.meshgrid(disp1,disp2)
    
    ax.set_title('Contour plot of $\mathcal{L}$ values')
    ax.set_xlim(x0,x1)
    ax.set_ylim(y0,y1)
    ax.set_xlabel(xlab+' [km s$^{-1}$]',size='large')
    ax.set_ylabel(ylab+' [km s$^{-1}$]',size='large')
    
    #contf = ax.contourf(dispg1,dispg2,L_vals.T,cmap = plt.cm.get_cmap('jet'),vmin=-40,vmax=10)
    
    extent = [disp1[0], disp1[-1], disp2[0], disp2[-1]]
    im = ax.imshow(L_vals.T,origin='lower',interpolation='bilinear',cmap = plt.cm.get_cmap('jet'),extent=extent)
    
    cb = plt.colorbar(im,orientation='vertical', extend='max')
    cb.set_label('Value of $\mathcal{L}$',size='large')
    
    ax.scatter(disp10,disp20,c='black')
    
    plt.show()
    
    return