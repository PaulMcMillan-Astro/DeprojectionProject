from Deproject_v1_0 import *
import matplotlib.pyplot as plt
plt.style.use('classic')

def plot_fv(phi,plane,vmin,dv,n):
    
    """Function that plots the velocity distribution f(v) that corresponds to
    a given array of phi-values."""
    
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    vxmin, vymin, vzmin = vmin
    vxmax, vymax, vzmax = vxmin+nx*dvx,vymin+ny*dvy,vzmin+nz*dvz
    
    if plane == 'xy':
        axsum = 2
        n1, n2 = nx, ny
        dx, dy = dvx, dvy
        x0, y0 = vxmin, vymin
        x1, y1 = vxmax, vymax
        xlab, ylab = '$v_x$', '$v_y$'
    elif plane == 'yz':
        axsum = 0
        n1, n2 = ny, nz
        dx, dy = dvy, dvz
        x0, y0 = vymin, vzmin
        x1, y1 = vymax, vzmax
        xlab, ylab = '$v_y$', '$v_z$'
    elif plane == 'xz':
        axsum = 1
        n1, n2 = nx, nz
        dx, dy = dvx, dvz
        x0, y0 = vxmin, vzmin
        x1, y1 = vxmax, vzmax
        xlab, ylab = '$v_x$', '$v_z$'
        
    fv = np.exp(phi)
    
    fvlog = np.log10(fv+1)
    
    twodfv = np.sum(fv,axis=axsum)
    
    xbins = np.arange(x0,x1+dx,dx)
    ybins = np.arange(y0,x1+dy,dy)
    
    xc = (xbins[1:]+xbins[:-1])/2
    yc = (ybins[1:]+ybins[:-1])/2
    
    fig, ax = plt.subplots(figsize=(8,6))
    
    ax.set_title('Contour plot of velocity distribution')
    ax.set_xlim(-100,100)
    ax.set_ylim(-100,100)
    ax.set_xlabel(xlab+' [km s$^{-1}$]',size='large')
    ax.set_ylabel(ylab+' [km s$^{-1}$]',size='large')
    
    extent = [xbins[0],xbins[-1],ybins[0],ybins[-1]]
#    im = ax.imshow(twodfv.T,origin='lower',interpolation='bilinear',vmin=0,vmax=twodfv.max(),cmap = plt.cm.get_cmap('plasma'),extent=extent)
    im = ax.contourf(twodfv.T,10,origin='lower',cmap = plt.cm.get_cmap('bone_r'),extent=extent,linestyles='solid')
    ax.contour(twodfv.T,10,origin='lower',extent=extent,colors='k',linewidths=0.5)
    
    cb = plt.colorbar(im,orientation='vertical', extend='max')
    plt.setp(cb.ax.get_yticklabels(), visible=False)
    cb.set_label('Number density of stars',size='large')
    
    plt.show()
    
    return 
def plot_L(phi_all,pvals,rhatvals,vmin,dv,n,alpha):
    
    """Function that plots all the values found using L_collector."""
    
    L_vals = L_collector(phi_all,pvals,rhatvals,vmin,dv,n,alpha)
    
    fig, ax = plt.subplots()
    ax.set_title(r'Maximization of $\tilde{\mathcal{L}}$')
    ax.set_xlabel('Iterations')
    ax.set_ylabel(r'Value of $\tilde{\mathcal{L}}$')
    
    x = np.arange(0,len(L_vals),1)
    
    ax.plot(x,L_vals)
    
    plt.show()

    return

def L_collector(phi_all,pvals,rhatvals,vmin,dv,n,alpha):
    
    """Function that takes an array with all values of phi from the maximisation
    scheme and computes the L value for each iteration."""
    
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    vxmin, vymin, vzmin = vmin
    N = len(pvals)
    
    sigma2 = calc_sigma2(pvals,rhatvals) 
    
    L_all = np.zeros(len(phi_all))
    
    K0 = np.ravel(calc_K(pvals[0],rhatvals[0],vmin,dv,n))
    
    Kvals = scisp.coo_matrix(K0)

    for i in range(1,N): #Loop that yield a sparse array of N K-values

        K = np.ravel(calc_K(pvals[i],rhatvals[i],vmin,dv,n))
        K_coo = scisp.coo_matrix(K)

        Kvals = scisp.vstack((Kvals,K_coo))

    Kvals_csc = Kvals.tocsc()
    
    for i in range(len(phi_all)):
        phi = np.ravel(phi_all[i])
        L_all[i] += get_negL(phi, Kvals_csc, N, alpha, dv, n, sigma2)
        
    return -L_all

def plot_K(plane, pvals, rhatvals, vmin, dv, n):
    
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    vxmin, vymin, vzmin = vmin
    vxmax, vymax, vzmax = vxmin+nx*dvx,vymin+ny*dvy,vzmin+nz*dvz
    
    N = len(pvals)

    Kvals = np.zeros((N,nx,ny,nz))
    
    for i in range(N):
        K = calc_K(pvals[i],rhatvals[i],vmin,dv,n)
        Kvals[i] += K 
    
    K_sum = np.sum(Kvals,axis=0)
    
    if plane == 'xy':
        axsum = 2
        n1, n2 = nx, ny
        dx, dy = dvx, dvy
        x0, y0 = vxmin, vymin
        x1, y1 = vxmax, vymax
        xlab, ylab = '$v_x$', '$v_y$'
    elif plane == 'yz':
        axsum = 0
        n1, n2 = ny, nz
        dx, dy = dvy, dvz
        x0, y0 = vymin, vzmin
        x1, y1 = vymax, vzmax
        xlab, ylab = '$v_y$', '$v_z$'
    elif plane == 'xz':
        axsum = 1
        n1, n2 = nx, nz
        dx, dy = dvx, dvz
        x0, y0 = vxmin, vzmin
        x1, y1 = vxmax, vzmax
        xlab, ylab = '$v_x$', '$v_z$'
    
    K_proj = np.sum(K_sum, axis=axsum)
    
    xbins = np.arange(x0,x1+dx,dx)
    ybins = np.arange(y0,x1+dy,dy)
    
    xc = (xbins[1:]+xbins[:-1])/2
    yc = (ybins[1:]+ybins[:-1])/2
    
    fig, ax = plt.subplots(figsize=(8,6))
    
    ax.set_title('Projection of $\sum K(k|l)$ onto'+' ${}$'.format(plane))
    ax.set_xlim(x0,x1)
    ax.set_ylim(y0,y1)
    ax.set_xlabel(xlab+' [km s$^{-1}$]',size='large')
    ax.set_ylabel(ylab+' [km s$^{-1}$]',size='large')
    
    extent = [xbins[0], xbins[-1], ybins[0], ybins[-1]]
    im = ax.imshow(K_proj.T,origin='lower',interpolation='bilinear',vmin=0,vmax=K_proj.max(),cmap = plt.cm.get_cmap('jet'),extent=extent)
    
    cb = plt.colorbar(im,orientation='vertical', extend='max')
    #plt.setp(cb.ax.get_yticklabels(), visible=False)
    cb.set_label('Value of $\sum K(k|l)$',size='large')
    
    plt.show()
    
def plot_grad_L(phi, plane, alpha, pvals, rhatvals, vmin, dv, n):
    
    """Plots the gradient in the middle of a box depending on what plane
    one has provided as input"""
    
    from scipy.optimize import approx_fprime
    
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    vxmin, vymin, vzmin = vmin
    vxmax, vymax, vzmax = vxmin+nx*dvx,vymin+ny*dvy,vzmin+nz*dvz
    
    N = len(pvals)
    
    sigma2 = calc_sigma2(pvals,rhatvals) 
    
    K0 = np.ravel(calc_K(pvals[0],rhatvals[0],vmin,dv,n))
    
    Kvals = scisp.coo_matrix(K0)

    for i in range(1,N): #Loop that yield a sparse array of N K-values

        K = np.ravel(calc_K(pvals[i],rhatvals[i],vmin,dv,n))
        K_coo = scisp.coo_matrix(K)

        Kvals = scisp.vstack((Kvals,K_coo))

    Kvals_csc = Kvals.tocsc()
    
    if plane == 'xy':
        axsum = 2
        n1, n2 = nx, ny
        dx, dy = dvx, dvy
        x0, y0 = vxmin, vymin
        x1, y1 = vxmax, vymax
        xlab, ylab = '$v_x$', '$v_y$'
    elif plane == 'yz':
        axsum = 0
        n1, n2 = ny, nz
        dx, dy = dvy, dvz
        x0, y0 = vymin, vzmin
        x1, y1 = vymax, vzmax
        xlab, ylab = '$v_y$', '$v_z$'
    elif plane == 'xz':
        axsum = 1
        n1, n2 = nx, nz
        dx, dy = dvx, dvz
        x0, y0 = vxmin, vzmin
        x1, y1 = vxmax, vzmax
        xlab, ylab = '$v_x$', '$v_z$'
        
    phir = np.ravel(phi)
        
    gL_vals = -get_grad_negL(phir, Kvals_csc, N, alpha, dv, n, sigma2)
    
    #gL_vals = approx_fprime(phir,get_negL,1e-5,Kvals,N,alpha,dv,n,sigma2)
    
    gL_unr = np.reshape(gL_vals, (nx,ny,nz))
    
    gLproj = np.sum(gL_unr,axis=axsum)
    
    xbins = np.arange(x0,x1+dx,dx)
    ybins = np.arange(y0,x1+dy,dy)
    
    xc = (xbins[1:]+xbins[:-1])/2
    yc = (ybins[1:]+ybins[:-1])/2
    
    fig, ax = plt.subplots(figsize=(8,6))
    
    ax.set_title(r'Contour plot of $\nabla \mathcal{L}$')
    ax.set_xlim(x0,x1)
    ax.set_ylim(y0,y1)
    ax.set_xlabel(xlab+' [km s$^{-1}$]',size='large')
    ax.set_ylabel(ylab+' [km s$^{-1}$]',size='large')
    
    extent = [xc[0], xc[-1], yc[0], yc[-1]]
    im = ax.imshow(gLproj.T,origin='lower',interpolation='bilinear',vmin=gLproj.min(),vmax=gLproj.max(),cmap = plt.cm.get_cmap('jet'),extent=extent)
    
    cb = plt.colorbar(im,orientation='vertical', extend='max')
    #plt.setp(cb.ax.get_yticklabels(), visible=False)
    cb.set_label(r'Value of $\nabla\mathcal{L}$',size='large')
    
    plt.show()
    
    return 
	
def plot_sec_der(phi, plane, pvals, rhatvals, vmin, dv, n):
    
    from scipy.optimize import approx_fprime
    
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    vxmin, vymin, vzmin = vmin
    vxmax, vymax, vzmax = vxmin+nx*dvx,vymin+ny*dvy,vzmin+nz*dvz
    
    N = len(pvals)
    
    sigma2 = calc_sigma2(pvals,rhatvals) 
    
    secd = sec_der(phi,sigma2,dv)
    
    if plane == 'xy':
        axsum = 2
        n1, n2 = nx, ny
        dx, dy = dvx, dvy
        x0, y0 = vxmin, vymin
        x1, y1 = vxmax, vymax
        xlab, ylab = '$v_x$', '$v_y$'
    elif plane == 'yz':
        axsum = 0
        n1, n2 = ny, nz
        dx, dy = dvy, dvz
        x0, y0 = vymin, vzmin
        x1, y1 = vymax, vzmax
        xlab, ylab = '$v_y$', '$v_z$'
    elif plane == 'xz':
        axsum = 1
        n1, n2 = nx, nz
        dx, dy = dvx, dvz
        x0, y0 = vxmin, vzmin
        x1, y1 = vxmax, vzmax
        xlab, ylab = '$v_x$', '$v_z$'
    
    secdproj = np.sum(secd,axis=axsum)
    
    xbins = np.arange(x0,x1+dx,dx)
    ybins = np.arange(y0,x1+dy,dy)
    
    xc = (xbins[1:]+xbins[:-1])/2
    yc = (ybins[1:]+ybins[:-1])/2
    
    fig, ax = plt.subplots(figsize=(8,6))
    
    ax.set_title(r'Contour plot of $\nabla^2 \log f(v)$')
    ax.set_xlim(x0,x1)
    ax.set_ylim(y0,y1)
    ax.set_xlabel(xlab+' [km s$^{-1}$]',size='large')
    ax.set_ylabel(ylab+' [km s$^{-1}$]',size='large')
    
    extent = [xbins[0], xbins[-1], ybins[0], ybins[-1]]
    im = ax.imshow(secdproj.T,origin='lower',interpolation='bilinear',vmin=secdproj.min(),vmax=secdproj.max(),cmap = plt.cm.get_cmap('jet'),extent=extent)
    
    cb = plt.colorbar(im,orientation='vertical', extend='max')
    #plt.setp(cb.ax.get_yticklabels(), visible=False)
    cb.set_label(r'Value of $\nabla^2 \log f(v)$',size='large')
    
    plt.show()
    
    return 