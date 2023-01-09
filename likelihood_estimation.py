import numpy as np
import scipy.stats as st
import psutil
import builtins
import time
from scipy import sparse as scisp
from scipy.ndimage import zoom
from scipy.optimize import fmin_cg, minimize
from termcolor import colored


class Deprojection:
    def __init__(self, alpha, pvals, rhatvals, vmin, dv, n, phi0_guess=[], v0_guess=[], disp_guess=[], noniso=True, polar=False):
        '''
        Constructs the attributes for the Deprojection from the input. Is used to set up a 3D grid in velocity space upon which we estimate the likelihood

        Parameters
        ----------
        alpha : float
            The value of the chosen smoothing parameter from Dehnen 1998 eq. 38, third term.

        pvals : Array of shape  (N, 3)
            The 3-D components of the transverse velocity. Either cartesian or polar.

        rhatvals : Array of shape  (N, 3)
            The unit vectors in the direction of the stars on the celestial sphere in 3-D coordinates. Either cartesian or polar.

        vmin : Array of shape (3,)
            The bottom corners of the grid in each component of velocity space.

        dv : Array of shape (3,)
            The width of each cell if the grid in each component of velocity space.

        n : Array of shape (3,)
            The dimensions of the grid.

        phi0_guess : Array of shape (n), optional
            The probability is f = exp(phi), so this is a guess of the probability distribution in logspace. If phi0_guess is not provided, it will be estimated using phi_guess().

        v0_guess : Array of shape (3,), optional
            Provides an initial guess for the average velocity in each dimension. If not provided it will be estimated by calc_sigma2().

        disp_guess : Array of shape (3,), optional
            Provides an initial guess for the velocity dispersion of the distribution in each velocity component. If not provided it will be stimated by calc_sigma2().

        non_iso : bool, default=True
            Specifies if the determination of vmean and sigma2 from calc_sigma2() should use the non-isotropic method or isotropic method (i.e the Dehnen & Binney 1998) method.

        polar : bool, default=False
            Set to true in .ini file if you require polar coordinates for your velocities. Ignores the calc_sigma2() determinations as it does not work with polar. Instead requires an educated guess of disp_guess.
        '''
        self.alpha = alpha
        self.pvals = pvals
        self.rhatvals = rhatvals
        self.vmin = vmin
        self.dv = dv
        self.n = n
        self.phi0_guess = phi0_guess
        self.v0_guess = v0_guess
        self.disp_guess = disp_guess
        self.non_iso = noniso
        self.polar = polar


    def max_L(self):
        '''
        Function that employs scipy.optimize.minimize with the conjugate gradient method to maximise the function get_neg_L().
        It takes a gaussian guess of the distribution and the relevant data from the star sample for which the velocity distribution is to be estimated.

        Returns
        -------
        mxl : Array of shape (n).
            The final and optimal state of the phi grid that maximises the likelihood. The exponential of this mxl will give the probability distribution.

        fmin_it : int
            The number of iterations required to complete the run.
        '''

        builtins.L     = []
        builtins.gradL = []

        dvx, dvy, dvz = self.dv
        nx, ny, nz = self.n
        N = len(self.pvals)

        self.sigma2, vmean = self.calc_sigma2(give_vmean=True)
        if np.size(self.v0_guess) == 0:
            v0_guess = vmean
        if np.size(self.disp_guess) == 0:
            sigma = np.sqrt(sigma2)
            disp_guess = sigma
        if np.size(self.phi0_guess) == 0:
            phi0 = self.phi_guess(v0_guess, disp_guess)
        else:
            phi0 = self.phi0_guess
        if polar:
            self.sigma2 = self.disp_guess**2

        self.Kvals = self.KvalsSparseMethod(self.dv, self.n)

        phi0r = np.ravel(phi0)  # fmin_cg only takes one-dimensional inputs for the initial guess

        print('Started minimization... ',end='')

        builtins.L.append(-1*self.get_neg_L(phi0r, self.dv, self.n))
        builtins.gradL.append(np.linalg.norm(self.get_grad_neg_L(phi0r, self.dv, self.n)))

        fmin_it = 0
        optr = minimize(fun=self.get_neg_L,
                        x0=phi0r,
                        jac=get_grad_neg_L,
                        args=(self.dv, self.n),
                        callback=self.callback,
                        method='CG',
                        options={'gtol': 1e-6, 'disp': False, 'return_all': True})

        mxl, fopt, fcalls, gcalls, flag = optr.x, optr.fun, optr.nfev, optr.njev, optr.status
        fmin_it += optr.nit
        print(colored('Finished!','green',attrs=['bold','underline']))

        self.fmin_cg_output(fopt, fcalls, gcalls, flag, fmin_it)

        builtins.n = self.n
        builtins.dv = self.dv

        mxl = mxl.reshape(self.n)

        return mxl, fmin_it


    def multigrid_max_L(self):
        '''
        Function that employs scipy.optimize.minimize with the conjugate gradient method to maximise the function get_neg_L().
        It takes a gaussian guess of the distribution and the relevant data from the star sample for which the velocity distribution is to be estimated.

        This function differs from max_L() in that it starts with a low-res box (nx,ny,nz) and itertively increases size to reach the final box size. This will significantly improve runtime and is used in Dehnen 1998.

        Returns
        -------
        mxl_dict : dict
            A dictionary containing the final and optimal state of the phi grid that maximises the likelihood for each grid step.

        fmin_it : int
            The number of iterations required to complete the run.

        '''

        self.N = len(self.pvals)
        self.vmax = self.vmin + self.dv*self.n

        box_steps = self.multigrid_steps()

        self.n = box_steps[0]
        dv = (vmax-self.vmin)/self.n

        self.sigma2, vmean = self.calc_sigma2(give_vmean=True)
        if np.size(self.v0_guess) == 0:
            v0_guess = vmean
        if np.size(self.disp_guess) == 0:
            sigma = np.sqrt(sigma2)
            disp_guess = sigma
        if np.size(self.phi0_guess) == 0:
            phi0 = self.phi_guess(v0_guess, disp_guess)
        else:
            phi0 = self.phi0_guess
        if polar:
            self.sigma2 = self.disp_guess**2

        builtins.L     = [[] for _ in range(len(box_steps))]
        builtins.gradL = [[] for _ in range(len(box_steps))]

        fmin_it = 0
        mxl_dict = {}
        for grid_step, n in enumerate(box_steps):
            builtins.grid_step = grid_step
            dv = (vmax-self.vmin)/n

            print(f'Starting Kvals after {(time.time() - builtins.ti)/60:.2f}mins...')
            self.Kvals = KvalsSparseMethod(dv, n)

            phi0r = np.ravel(phi0)  # fmin_cg only takes one-dimensional inputs for the initial guess

            ### This is where the minimization occurs
            print(f'Started minimization on {(n[0],n[1],n[2])} grid after {(time.time() - builtins.ti)/60:.2f} mins...', end='', flush=True)
            builtins.L[grid_step].append(-1*get_neg_L(phi0r, dv, n))
            builtins.gradL[grid_step].append(np.linalg.norm(get_grad_neg_L(phi0r, dv, n)))

            optr = minimize(fun=get_neg_L,
                            x0=phi0r,
                            method='CG',
                            jac=get_grad_neg_L,
                            args=(dv, n),
                            callback=callback_mg,
                            options={'gtol': 1e-6, 'disp': False, 'return_all': True})

            mxl, fopt, fcalls, gcalls, flag = optr.x, optr.fun, optr.nfev, optr.njev, optr.status
            fmin_it += optr.nit

            mxl_dict['x'.join(n.astype(str))] = mxl.reshape(n)

            print(colored('Finished!','green',attrs=['bold','underline']),end='')
            print(' fopt : %s' % fopt)

            if grid_step == len(box_steps)-1:
                self.fmin_cg_output(fopt, fcalls, gcalls, flag, fmin_it)
            else:
                phi0 = self.zoomed_mxl(mxl.reshape(n))

        builtins.n = n
        builtins.dv = dv

        mxlnew = mxl.reshape(n)

        return mxl_dict, fmin_it


    def multigrid_steps(self, step=5):
        '''This function determines how many multigrids steps can be performed and returns the box size
        in each step. If the starting box dimensions are not divisible by 2**step the final box dimensions
        are changed to accomodate this.

        Parameters
        ----------
        step : int
            The max amount of multigrid steps we are aiming for. We do not divide dimensions below 10.


        Returns
        -------
        array
            An array of shape (steps, ndim).
        '''

        while any(np.round(self.n/(2**step)) < 10):
            step -= 1
        box_steps = (2**np.linspace(0,step,step+1).reshape(step+1,1) * np.round(self.n/(2**step))).astype(int)

        if not all(box_steps[-1] == self.n):
            print('To oct-split box the recommended %s times, dimensions were changed from (%s, %s, %s) to (%s, %s, %s)\n'
                  % (len(box_steps),self.n[0],self.n[1],self.n[2],box_steps[-1,0],box_steps[-1,1],box_steps[-1,2]))

        return box_steps.astype(int)


    def calc_sigma2(self, give_vmean=False):
        """Function that applies equation 12 of DB98 for a set of stars from their tangential velocities and unit vectors. Returns the velocity dispersion tensor.

        Parameters
        ----------
        give_vmean : bool, default=False
            Specifies whether or not the mean velocities should also be output.

        Returns
        -------
        sigma2 : Array of shape (3,)
            The velocity dispersion squared in each component.

        vmean : Array of shape (3,), optional
            The average velocity in each component.
        """

        pmean = self.pvals.mean(axis=-2)

        # Fast way for outer product eq (4) of DB98
        A = np.identity(3) - self.rhat[..., np.newaxis] * self.rhat[..., np.newaxis, :]


        A_mean_inv = np.linalg.inv(A.mean(axis=-3))
        v_mean = np.einsum('...ij,...j->...i', A_mean_inv, pmean)
        del A_mean_inv

        # Computes p' from equation (6) in DB98
        pp = self.pvals - np.einsum('...ikj,...j->...ik', A, v_mean)

        pp2mean = np.mean(np.square(pp), axis=-2)
        if self.noniso:
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


    def phi_guess(self, v0, disp0):
        """Create an initial guess for the phi values in each bin based on the assumed initial guess of average velocity and dispersion.

        The initial guess is the sum of two multivariate Gaussians fA and fB. fA uses the inital guess for v0 and disp0, while fB uses a dispersion twice as large in each dimension. The sum of two Gaussians is critical as testing reveals that the code breaks down for a single Gaussian.

        Parameters
        ----------
        v0 : Array of shape (3,)
            The average velocity in each component.

        disp0 : Array of shape (3,)
            The velocity dispersion in each component.

        Returns
        -------
        phi : Array of shape (n)
            The Gaussian guess phi grid of logarithmic probabilities.
        """

        vxmin, vymin, vzmin = self.vmin
        dvx, dvy, dvz = self.dv
        nx, ny, nz = self.n
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

        # Gets the 3D coordinates of each grid cell in the correct format.
        pos = np.stack(np.meshgrid(vxc,vyc,vzc, indexing='ij'),axis=3)

        # Sets up both of the multivariate Guassians
        fA = st.multivariate_normal(mean=v0, cov=np.diag(dispA**2))
        fB = st.multivariate_normal(mean=v0, cov=np.diag(dispB**2))

        # Get the values of the Guassians at our cell positions and join them.
        phi = np.log((fA.pdf(pos) + fB.pdf(pos))/2)
        return phi


    def KvalsSparseMethod(self, dv, n):
        '''
        Checks if there is RAM avaiable to create an array with n_stars*nx*ny*nz floats. Allows for the use of up to 80% of available RAM. This function is necessary since realistically, someone else might be using/start using the RAM of the server while running, in which case we need to terminate before running into real system memory errors.

        The function subsequently calls calc_sparse_K() for each of the stars and stores it in a scipy.sparse array, of type lil_matrix. This is converted to a sparse csc_matrix for later convenience.

        Parameters
        ----------
        dv : Array of shape (3,)
            The width of the cells. Not self.dv since in the multigrid scheme, coarse grids have different dv.

        n : Array of shape (3,)
            The shape of the grid. Not self.n since in the multigrid scheme, coarse grids have different n.

        Returns
        -------
        scipy.sparse.csc_matrix of shape (n_stars, n)
            The sparse array containing output grid of K for each star
        '''

        n_stars = len(self.pvals)
        allocated_memory = psutil.virtual_memory().available*0.8/1e9
        required_memory = 8*n_stars*np.ceil(np.linalg.norm(n)).astype(int)/1e9

        if required_memory > allocated_memory:
            raise MemoryError(f'Required memory for {required_memory:.2f} GB exceeds available memory {allocated_memory:.2f} GB')

        n_stars = len(pvals)
        lil = scisp.lil_matrix((n_stars, np.prod(n)))

        for i in range(n_stars):
            coords, vals = self.calc_sparse_K(self.pvals[i], self.rhat[i], dv, n)
            lil[i, coords] = vals

        return lil.tocsc()


    def calc_sparse_K(self, pk, rhat, dv, n):
        '''
        The value of K(k|l) on eq. 28 of Dehnen 1998 is defined as the length of the line /v = pk + vr*rhat/ that lies in cell l divided by the cell widths. Since full 3D velocity is a point in our grid, the unknown vr creates a line. This is the line equation above. Here we calculate the grid K of the same size as our grid n. The values in the grid correspond then to the length of the line crossing the cells. Since most of the cells are zero (i.e the line doesn't cross most cells) it is eventually put into a sparse matrix.

        This calculation, by nature, is extremely difficult to vectorise and we do it for a single star (and loop over all stars in KvalsSparseMethod(). We can however, calculate the line length for each cell simultaneously.

        Parameters
        ----------
        pk : Array of shape  (3,)
            The 3-D components of the transverse velocity of a given star. Either cartesian or polar.

        rhatvals : Array of shape  (3,)
            The unit vectors in the direction of the star on the celestial sphere in 3-D coordinates. Either cartesian or polar.

        dv : Array of shape (3,)
            The width of the cells. Not self.dv since in the multigrid scheme, coarse grids have different dv.

        n : Array of shape (3,)
            The shape of the grid. Not self.n since in the multigrid scheme, coarse grids have different n.

        Returns
        -------
        Array of shape (l,)
            Array containing the flattened cell number of the cells that the line crosses through.

        Array of shape (l,)
            Array containing the line lengths through each of the cells that the line crosses through.
        '''

        vmax = self.vmin + dv*n

        # We now solve the line equation v = pk + vr*rhat, where v are the intersections of the tangential velocity line and the boundaries of each bin

        vbins = (np.linspace(self.vmin[i], vmax[i], n[i] + 1) for i in range(3))

        vr = [(np.array(vb) - pk[i])/rhat[i] for i, vb in enumerate(vbins)]

        # After solving the line equation for each dim we remove the vr values which solve the equation for bins outside of our specified box.

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
        vmins = np.ones((len(vr_prime), len(self.vmin))) * self.vmin[np.newaxis]
        vr_primestack = np.ones((len(vr_prime), 3))*vr_prime[:, np.newaxis]

        # We now solve the line equation again for values in the middle of each bin with a line segment in it. This gives us the coordinates for each relevant bin, given in line_bins. Finally we compute the length of each segment and output the multi_index and values.

        v_prime = pks + vr_primestack * rhats

        line_bins = ((v_prime - vmins) // dv).astype(int)

        line_len = vr[1:] - vr[:-1]  # Gets the length of each line
        non_zero = np.nonzero(line_len)
        line_len = line_len[non_zero]  # Removes duplicate intersections
        line_bins = line_bins[non_zero].T

        vals = (line_len / np.prod(dv))
        return np.ravel_multi_index(line_bins, n), vals


    def get_neg_L(self, phi, dv, n):
        """The function that we wish to optimise. Corresponds to eq. 31 in Dehnen 1998. We calculate the three separate terms of the quation, sum them up, and multiply it with -1 in order to use scipy.optimize.minimize() for a maxmimum likelihood.

        Parameters
        ----------
        phi : Array of shape (np.prod(n),)
            Flattened version of the logarithmic probabilities. It is flattened since scipy.maximise requires this.

        dv : Array of shape (3,)
            The width of the cells. Not self.dv since in the multigrid scheme, coarse grids have different dv.

        n : Array of shape (3,)
            The shape of the grid. Not self.n since in the multigrid scheme, coarse grids have different n.

        Returns
        -------
        float
            The negative value of the likelihood in this instance.

        """


        exphi = np.exp(phi)
        Kphi = (self.Kvals @ exphi.T) # Order all Kphi values in 1D arrays and compute the sum of exp(phi)*K(k|l) for each star
        Kphi_sum_tot = np.log(Kphi[Kphi != 0]).sum() # To make sure we don't get infinities and .sum() gives the double sum in the first term

        phi_unr = np.reshape(phi, n)
        phixhi_sum = (self.sec_der(phi_unr, dv)**2).sum()

        t1 = Kphi_sum_tot / self.N
        t2 = exphi.sum()
        t3 = ((self.alpha * dv[0] * dv[1] * dv[2]) / 2) * phixhi_sum

        L_tilde = t1 - t2 - t3 # eq. 31 in DB98
        neg_L = -1 * L_tilde  # Since we want to maximize L_tilde, we should minimize -L_tilde
        builtins.current_L = L_tilde
        return neg_L


    def sec_der(self, phi, dv):
        '''
        Calculates eq. 29 (and in turn eq. 30) of WD98. Each cell /m/ has a contribution
        of -2*phi_m, phi_(m+e_i), and phi_(m-e_i) from itself, its neighbor one step up in the i
        direction and one step down in the i direction respectively.

        The contributions are added onto a zero array of the same shape as phi and only when
        both neighbours are available in the given dimension. e.g. there is no x term when
        treating on [1, :, :] or [:, :, -1]. Or in other words, we do not calculate the value of eq. 30 when the cell /m/ lies along the edges of the grid.

        Parameters
        ----------
        phi : Array of shape (n)
            The logarithmic probabilities.

        dv : Array of shape (3,)
            The width of the cells. Not self.dv since in the multigrid scheme, coarse grids have different dv.

        Returns
        -------
        phi_arr : Array of shape (n)
            The contributions to each cell, equiv. to eq 29 in WD98
        '''

        phi_arr = np.zeros(phi.shape)

        phi_arr[1:-1,:,:] += (self.sigma2[0] / (dv[0]*dv[0])) * (-2*phi[1:-1, :, :] + phi[2:, :, :] + phi[:-2, :, :])
        phi_arr[:,1:-1,:] += (self.sigma2[1] / (dv[1]*dv[1])) * (-2*phi[:, 1:-1, :] + phi[:, 2:, :] + phi[:, :-2, :])
        phi_arr[:,:,1:-1] += (self.sigma2[2] / (dv[2]*dv[2])) * (-2*phi[:, :, 1:-1] + phi[:, :, 2:] + phi[:, :, :-2])

        return phi_arr


    def get_grad_neg_L(self, phi, dv, n):
        """
        The Jacobian, or the derivative/gradient of L. This took us ages to figure out. Requires some new stuff as well as a new function grad_sec_der().

        If you want to ask me (Daniel) questions about this, expect something like this:
        https://youtu.be/2Ll4iSvosDo?t=3

        Parameters
        ----------
        phi : Array of shape (np.prod(n),)
            Flattened version of the logarithmic probabilities. It is flattened since scipy.maximise requires this.

        dv : Array of shape (3,)
            The width of the cells. Not self.dv since in the multigrid scheme, coarse grids have different dv.

        n : Array of shape (3,)
            The shape of the grid. Not self.n since in the multigrid scheme, coarse grids have different n.

        Returns
        -------
        Array
            The negative value of the likelihood in this instance.

        """

        exphi = np.exp(phi)

        Kphi_sum = (self.Kvals @ exphi.T)
        Kphi_sum[Kphi_sum != 0] = Kphi_sum[Kphi_sum != 0]**(-1)
        K_term = exphi*(Kphi_sum.T @ self.Kvals) # The final array with the first term for each cell

        phi_unr = np.reshape(phi, n)
        dphixhi = self.grad_sec_der(phi_unr, dv)

        t1 = K_term/self.N
        t2 = exphi
        t3 = ((self.alpha * dv[0] * dv[1] * dv[2]) / 2) * dphixhi.ravel()

        grad_L = np.asarray(t1 - t2 - t3).reshape(len(phi), )

        neg_grad_L = -1 * grad_L
        builtins.current_gradL = np.linalg.norm(neg_grad_L)

        return neg_grad_L


    def grad_sec_der(self, phi, dv):
        '''
        Here we calculate the equivalent factor to sec_der for the gradients third term, it requires a call to the previous sec_der as well as it factors into the gradient.

        Parameters
        ----------
        phi : Array of shape (n)
            The logarithmic probabilities.

        dv : Array of shape (3,)
            The width of the cells. Not self.dv since in the multigrid scheme, coarse grids have different dv.

        Returns
        -------
        phi_arr : Array of shape (n)
            The contributions to each cell, equiv. to eq 29 in WD98

        '''

        phi_arr = np.zeros(phi.shape)

        phi_arr_loc = self.sec_der(phi, dv) # This gives us the matrix A_m for all m = (i,j,k) cells

        # The x contribution
        phi_arr[:-2,:,:]  += 2 * (phi_arr_loc[1:-1,:,:]) * self.sigma2[0]/(dv[0]*dv[0]) # Adds A_(m-1) contribution
        phi_arr[2:,:,:]   += 2 * (phi_arr_loc[1:-1,:,:]) * self.sigma2[0]/(dv[0]*dv[0]) # Adds A_(m+1) contribution
        phi_arr[1:-1,:,:] -= 4 * (phi_arr_loc[1:-1,:,:]) * self.sigma2[0]/(dv[0]*dv[0]) # Adds A_m contribution

        # The y contribution
        phi_arr[:,:-2,:]  += 2 * (phi_arr_loc[:,1:-1,:]) * self.sigma2[1]/(dv[1]*dv[1]) # Adds A_(m-1) contribution
        phi_arr[:,2:,:]   += 2 * (phi_arr_loc[:,1:-1,:]) * self.sigma2[1]/(dv[1]*dv[1]) # Adds A_(m+1) contribution
        phi_arr[:,1:-1,:] -= 4 * (phi_arr_loc[:,1:-1,:]) * self.sigma2[1]/(dv[1]*dv[1]) # Adds A_m contribution

        # The z contribution
        phi_arr[:,:,:-2]  += 2 * (phi_arr_loc[:,:,1:-1]) * self.sigma2[2]/(dv[2]*dv[2]) # Adds A_(m-1) contribution
        phi_arr[:,:,2:]   += 2 * (phi_arr_loc[:,:,1:-1]) * self.sigma2[2]/(dv[2]*dv[2]) # Adds A_(m+1) contribution
        phi_arr[:,:,1:-1] -= 4 * (phi_arr_loc[:,:,1:-1]) * self.sigma2[2]/(dv[2]*dv[2]) # Adds A_m contribution

        return phi_arr

    def callback(self, x):
        '''
        A callback function that scipy.optimize.minimize() calls for every iteration it does. Its purpose is simply to log the current values of the likelihood and the norm of its derivative so these can be plotted at the end for information about the run.
        '''
        builtins.L.append(builtins.current_L)
        builtins.gradL.append(builtins.current_gradL)


    def fmin_cg_output(self, fopt, fcalls, gcalls, flag, fmin_it):
        '''
        This function collects information about the minimization run and prints a helpful finishing messages.

        Parameters
        ----------
        fopt : float
            The final value of the likelihood.

        fcalls : int
            The number of calls to get_neg_L()

        gcalls : int
            The number of calls to get_grad_neg_L()

        flag : int
            An integer corresponding to different types of termination causes. If not 0, you probably want to fix something.

        fmin_it : int
            The number of iterations required to complete.
        '''
        l1 = ['Optimization terminated sucessfully.\n',
              colored('Warning','red',attrs=['bold'])+': Maximum number of iterations has been exceeded.\n',
              colored('Warning','red',attrs=['bold'])+': Desired error not necessarily achieved due to precision loss.\n',
              colored('Warning','red',attrs=['bold'])+': NaN result encountered.\n']
        l2 = ('         Current function value : %f\n' % fopt)
        l3 = ('         Iterations             : %s\n' % fmin_it)
        l4 = ('         Function evaluations   : %s\n' % fcalls)
        l5 = ('         Gradient evaluations   : %s\n' % gcalls)
        print(l1[flag] + l2 + l3 + l4 + l5)


    def zoomed_mxl(self, mxl):
        '''
        After a successful MLE on current multigrid step, produce phi0_guess of next step by upscaling MLE with interplation. This is done with scipy.ndimage.zoom(). The -log(8) is equivalent to diving by 2 for each dimension in linear space, ensuring the integral of the probability is one.

        Parameters
        ----------
        mxl : array_like
            The input array, our MLE.
        '''

        phi0_guess = zoom(mxl, zoom=2, order=3) - np.log(8)

        return phi0_guess


    def callback_mg(self, x):
        '''
        A callback function that scipy.optimize.minimize() calls for every iteration it does. Its purpose is simply to log the current values of the likelihood and the norm of its derivative so these can be plotted at the end for information about the run.
        '''
        i = builtins.grid_step
        builtins.L[i].append(builtins.current_L)
        builtins.gradL[i].append(builtins.current_gradL)
