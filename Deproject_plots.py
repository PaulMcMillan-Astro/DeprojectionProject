from Deproject_v0 import *
plt.style.use('classic')

def plot_fv(phi,plane,vmin,dv,n):
    
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
    
    twodfv = np.sum(fv,axis=axsum)
    
    xbins = np.arange(x0,x1+dx,dx)
    ybins = np.arange(y0,x1+dy,dy)
    
    xc = (xbins[1:]+xbins[:-1])/2
    yc = (ybins[1:]+ybins[:-1])/2
    
    fig, ax = plt.subplots(figsize=(8,6))
    
    ax.set_title('Contour plot of velocity distribution')
    ax.set_xlim(x0,x1)
    ax.set_ylim(y0,y1)
    ax.set_xlabel(xlab+' [km s$^{-1}$]',size='large')
    ax.set_ylabel(ylab+' [km s$^{-1}$]',size='large')
    
    extent = [xc[0], xc[-1], yc[0], yc[-1]]
    im = ax.imshow(twodfv.T,origin='lower',interpolation='bilinear',vmin=0,vmax=twodfv.max(),cmap = plt.cm.get_cmap('jet'),extent=extent)
    
    cb = plt.colorbar(im,orientation='vertical', extend='max')
    #cb.ax.xaxis.get_ticklabels().set_visible(False)
    cb.set_label('Number density of stars',size='large')
    
    plt.show()
    
    return 

def plot_L(phi_all,pvals,rhatvals,vmin,dv,n,alpha):
    
    L_vals = L_collector(phi_all,pvals,rhatvals,vmin,dv,n,alpha)
    
    fig, ax = plt.subplots()
    ax.set_title('Maximization of $\mathcal{L}$')
    ax.set_xlabel('Iterations')
    ax.set_ylabel('Value of $\mathcal{L}$')
    
    x = np.arange(0,len(L_vals),1)
    
    ax.plot(x,L_vals)
    
    plt.show()

    return

def L_collector(phi_all,pvals,rhatvals,vmin,dv,n,alpha):
    
    dvx, dvy, dvz = dv
    nx, ny, nz = n
    vxmin, vymin, vzmin = vmin
    N = len(pvals)
    
    sigma2 = calc_sigma2(pvals,rhatvals) 
    
    Kvals = np.zeros((N,nx,ny,nz))
    
    L_all = np.zeros(len(phi_all))
    
    for i in range(N):
        K = calc_K(pvals[i],rhatvals[i],vmin,dv,n)
        Kvals[i] += K 
    
    for i in range(len(phi_all)):
        phi = np.ravel(phi_all[i])
        L_all[i] += get_L(phi, Kvals, N, alpha, dv, n, sigma2)
        
    return -L_all