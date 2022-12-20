from Deproject_v1_0 import *
import matplotlib.pyplot as plt

def DM_plt_prefs():
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color'];
    plt.rcParams.update({'font.size'           : 20,
                         'axes.linewidth'      : 1.2,
                         'axes.grid'           : False,
                         'axes.axisbelow'      : False,
                         'figure.figsize'      : [6, 6],
                         'figure.frameon'      : True,
                         'axes.labelsize'      : 20,
                         'xtick.labelsize'     : 16,
                         'xtick.direction'     :'in',
                         'xtick.top'           : True,
                         'xtick.bottom'        : True,
                         'xtick.major.size'    : 7.5,
                         'xtick.minor.size'    : 3.5,
                         'xtick.major.width'   : 1,
                         'xtick.minor.width'   : 1,
                         'xtick.minor.visible' : True,
                         'xtick.major.top'     : True,
                         'xtick.major.bottom'  : True,
                         'xtick.minor.top'     : True,
                         'xtick.minor.bottom'  : True,
                         'ytick.direction'     :'in',
                         'ytick.labelsize'     : 16,
                         'ytick.minor.visible' : True,
                         'ytick.major.right'   : True,
                         'ytick.major.left'    : True,
                         'ytick.minor.right'   : True,
                         'ytick.minor.left'    : True,
                         'ytick.right'         : True,
                         'ytick.left'          : True,
                         'ytick.major.size'    : 7.5,
                         'ytick.minor.size'    : 3.5,
                         'ytick.major.width'   : 1,
                         'ytick.minor.width'   : 1,
                        })
    return

def plot_fv(phi,plane,vmin,dv,n,folder,logging=0, polar=False):
    
    """Function that plots the velocity distribution f(v) that corresponds to
    a given array of phi-values."""
    
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    vxmin, vymin, vzmin = vmin
    vxmax, vymax, vzmax = vxmin+nx*dvx,vymin+ny*dvy,vzmin+nz*dvz
    
    while plane != 'xy' and plane != 'yz' and plane != 'xz':    
        plane = input('Incorrect entry, try again! ')
    if plane == 'xy':
        axsum = 2
        n1, n2 = nx, ny
        dx, dy = dvx, dvy
        x0, y0 = vxmin, vymin
        x1, y1 = vxmax, vymax
        if not polar:
           xlab, ylab = '$v_x$', '$v_y$'
        elif polar:
           xlab, ylab = '$v_r$', '$v_\\theta$'
    elif plane == 'yz':
        axsum = 0
        n1, n2 = ny, nz
        dx, dy = dvy, dvz
        x0, y0 = vymin, vzmin
        x1, y1 = vymax, vzmax
        if not polar:
           xlab, ylab = '$v_y$', '$v_z$'
        elif polar:
           xlab, ylab = '$v_\\theta$', '$v_\phi$'
    elif plane == 'xz':
        axsum = 1
        n1, n2 = nx, nz
        dx, dy = dvx, dvz
        x0, y0 = vxmin, vzmin
        x1, y1 = vxmax, vzmax
        xlab, ylab = '$v_x$', '$v_z$'
        if not polar:
           xlab, ylab = '$v_x$', '$v_z$'
        elif polar:
           xlab, ylab = '$v_r$', '$v_\phi$'
        
    fv = np.exp(phi)
    
    fvlog = np.log10(fv+1)
    
    twodfv = np.sum(fv,axis=axsum)
    
    xbins = np.linspace(x0,x1,n1+1)
    ybins = np.linspace(y0,y1,n2+1)
    
    xc = (xbins[1:]+xbins[:-1])/2
    yc = (ybins[1:]+ybins[:-1])/2
    
    [X,Y] = np.meshgrid(xc,yc);
    fig, ax = plt.subplots(figsize=(8,6))

    
    ax.set_title('Contour plot of velocity distribution')
    ax.set_xlabel(xlab+' [km s$^{-1}$]',size='large')
    ax.set_ylabel(ylab+' [km s$^{-1}$]',size='large')
#    im = ax.imshow(twodfv.T,origin='lower',interpolation='bilinear',vmin=0,vmax=twodfv.max(),cmap = plt.cm.get_cmap('plasma'),extent=extent)
    im = ax.contourf(X,Y,twodfv.T,40,origin='lower',cmap = plt.cm.get_cmap('bone_r'),linestyles='solid')
    ax.contour(X,Y,twodfv.T,10,origin='lower',colors='k',linewidths=0.5)
    
    cb = plt.colorbar(im,orientation='vertical')
    plt.setp(cb.ax.get_yticklabels(), visible=False)
    cb.set_label('Number density of stars',size='large')
    
    if logging:
        plt.savefig('RUNS/' + folder + '/fvplot_' + plane + '.pdf',format='pdf')
    if not builtins.autoplot:
        plt.show()
    
    return 

def plot_L_and_dL(folder, logging=0):
    f,ax = plt.subplots(2,1,figsize=(7,9),sharex=True,frameon=False,gridspec_kw={'height_ratios':[3,1]})
    
    L = builtins.L
    gradL = builtins.gradL    
    did_mg = any(isinstance(entry, list) for entry in L)
    if did_mg:
        mg_string = "Iterations per grid:\n" + " -> ".join([str(len(l)) for l in L])

        L = sum(L,[])
        gradL = sum(gradL,[])
    
    ax[0].plot(range(1,len(gradL)+1),gradL,'k-',linewidth=2)
    ax[0].set_yscale('log')
    ax[0].set_xscale('log')
    ax[0].set_ylabel('$|\\nabla\\widetilde{\mathscr{L}}_\\alpha(\\varphi)|$')

    ax[1].set_ylabel('$\mathscr{L}_\\alpha(\\varphi)$')
    ax[1].plot(range(1,len(L)+1),L,'k-',linewidth=2)
    ax[1].set_xlabel('Iterations')
    ax[1].set_xscale('log')
    
    if did_mg:
        bbox_props = dict(boxstyle="square", fc="white", ec="k", lw=1)
        ax[0].text(0.012,1.02, mg_string, transform=ax[0].transAxes,fontsize=14, bbox=bbox_props)
    
    plt.tight_layout(h_pad=0.02)
    if logging:
        plt.savefig('RUNS/' + folder + '/LanddLplot' + '.pdf',format='pdf')
    if not builtins.autoplot:
        plt.show()

    return

def plot_L(phi_all,pvals,rhatvals,vmin,dv,n,alpha,logging,folder=''):
    
    """ ---This function is now obsolete and replaced by plot_L_and_dL()---
    Function that plots all the values found using L_collector."""
    
    L_vals = L_collector(phi_all,pvals,rhatvals,vmin,dv,n,alpha)
    
    fig, ax = plt.subplots()
    ax.set_title(r'Maximization of $\tilde{\mathcal{L}}$')
    ax.set_xlabel('Iterations')
    ax.set_ylabel(r'Value of $\tilde{\mathcal{L}}$')
    
    x = np.arange(0,len(L_vals),1)
    
    ax.plot(x,L_vals)
    if logging:
        plt.savefig('RUNS/' + folder + '/Lplot' + '.pdf',format='pdf')
    plt.show()

    return

def L_collector(phi_all,pvals,rhatvals,vmin,dv,n,alpha):
    
    """---This function is now obsolete and replaced by plot_L_and_dL()---
    Function that takes an array with all values of phi from the maximisation
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
