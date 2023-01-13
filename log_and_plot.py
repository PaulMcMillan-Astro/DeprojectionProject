from datetime import date
import os
import string
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import builtins

class Logger:
    def __init__(self, mxl):
        self.mxl = mxl
        plt.style.use('pretty_plots.mplstyle')


    def make_folder(self):
        '''Function that creates a directory with name YYYY-MM-DDx where x
        is the first letter that is a non-existing directory

        Returns
        -------
        string
            The name of the output directory
        '''

        alp = string.ascii_lowercase
        date_str = str(date.today())

        list_str = []
        for letter in alp:
            list_str.append(date_str + letter)

        for letter1 in alp:
            for letter2 in alp:
                list_str.append(date_str + letter1 + letter2)

        list_str = np.array(list_str)


        os_dirs = os.popen('ls RUNS | grep "%s"' % date_str).read().split("\n")[:-1]
        os_dirs.sort(key=lambda x: len(x))
        os_dirs = np.array(os_dirs)

        existing_dirs = np.zeros(len(list_str),dtype='<U12')
        existing_dirs[:len(os_dirs)] = os_dirs

        self.folder=list_str[list_str != existing_dirs][0]
        output = subprocess.getoutput('mkdir RUNS/' + self.folder)
        if output:
            raise OSError(output)
        return


    def log_mle_run(self, fmin_it, endtime, n, vmin, dv, polar, use_guess, non_iso, v_guess, disp_guess, alpha, datafile):
        '''
        Function that saves the output maximum likelihood estimate as well as some information about the run.

        Parameters
        ----------
        fmin_it : int
            The number of iterations required to run the minimization.

        endtime : float
            The time it took to run the minimization in numbers.

        n : Array of shape (3,)
            The dimensions of the grid.

        vmin : Array of shape (3,)
            The bottom corners of the grid in each component of velocity space.

        dv : Array of shape (3,)
            The width of each cell if the grid in each component of velocity space.

        polar : bool
            True if the velocities are in Galactocentric spherical coordinates.

        use_guess : bool
            True if the run uses the input guess values of the .ini file and False if they are estimated from calc_sigma2().

        non_iso : bool
            Specifies if the determination of vmean and sigma2 from calc_sigma2() used the non-isotropic method or isotropic method (i.e the Dehnen & Binney 1998) method.

        v_guess : Array of shape (3,)
            The initial guess for the average velocity in each dimension. If use_guess = False it was estimated by calc_sigma2().

        disp_guess : Array of shape (3,)
            The initial guess for the velocity dispersion in each dimension. If use_guess = False it was estimated by calc_sigma2().

        alpha : float
            The value of the chosen smoothing parameter from Dehnen 1998 eq. 38, third term.

        datafile : string
            The name of the datafile used.

        '''

        self.n = n
        self.vmin = vmin
        self.dv = dv
        self.polar = polar
        self.use_guess = use_guess
        self.non_iso = non_iso
        self.v_guess = v_guess
        self.disp_guess = disp_guess
        self.alpha = alpha
        self.datafile = datafile

        # Create a folder for the run and save mxl data
        np.save('RUNS/' + self.folder + '/mxl_data', self.mxl)

        # Save output
        # RUNS folder identifier
        with open('logs/log_dir_identifier.txt', 'a') as logfile:
            if self.folder[10] == 'a' and len(self.folder) == 11:
                mark = '='
                logfile.write('\n' + mark*120 + '\n')
                logfile.write(mark*55 + self.folder[:10] + mark*55 + '\n')
                logfile.write(mark*120 + '\n')

            logfile.write('\nFolder name : ' + self.folder + '\n')
            logfile.write('Datafile    : ' + self.datafile + '\n')
            logfile.write('fmin its    : ' + str(fmin_it) + '\n')
            logfile.write('Time needed : ' + str(endtime/60) + ' hrs\n')

        # Logfile in RUNS folder
        with open('RUNS/' + self.folder + '/log.txt', 'a') as logfile:
            logfile.write('Datafile    : ' + self.datafile + '\n')
            logfile.write('fmin its    : ' + str(fmin_it) + '\n')
            logfile.write('Time needed : ' + str(endtime/60) + ' hrs\n')
            logfile.write('Labels      : n[1x3], vmin[1x3], dv[1x3], polar, use_guess, noniso, mu_guess[1x3], sigma_guess[1x3], alpha\n')
            value_string=str((str(list(self.n)).replace(",",":").replace(":",""),
                              str(list(self.vmin)).replace(",",":").replace(":",""),
                              str(list(self.dv)).replace(",",":").replace(":",""),
                              int(self.polar),
                              int(self.use_guess),
                              int(self.non_iso),
                              str(list(self.v_guess)).replace(",",":").replace(":",""),
                              str(list(self.disp_guess)).replace(",",":").replace(":",""),
                              self.alpha)).replace("'","")[1:-1]

            logfile.write("Values      : " + value_string + '\n')
        return


    def plot_and_save(self, plane='all'):
        """This function plots the MLE in the specified plane along with contour lines which contain 95%, 90%, 80%, 68%, 50%, 33%, 21%, 12%, 6%, 2% of stars within them, going from outside and inwards.

        Parameters
        ----------
        plane : string, default='all'
            Specifies the plane to plot in, can be one of ['UV', 'UW', 'VW', 'rphi', 'rtheta', 'phitheta', 'all']. If 'all', it will loop over the avaialble dimensions and plot all three combinations.
        """

        # If mxl is a multigrid, extract the final grid for plotting.
        if isinstance(self.mxl, dict):
            self.mxl = self.mxl[list(self.mxl.keys())[-1]]

        planes = ['UV', 'UW', 'VW', 'rphi', 'rtheta', 'phitheta', 'all']
        if plane not in planes:
            raise ValueError('Invalid plane type. Expected one of: %s' % planes)
        elif plane in ['UV', 'UW', 'VW'] and self.polar:
            raise ValueError('Specified plane is cartesian but polar=True, this does not compute.')
        elif plane in ['rphi', 'rtheta', 'phitheta'] and not self.polar:
            raise ValueError('Specified plane is spherical but polar=False, this does not compute.')

        dvx, dvy, dvz = self.dv
        nx, ny, nz = self.n
        self.vmax = self.vmin + self.n*self.dv

        xbins = np.linspace(self.vmin[0], self.vmax[0], self.n[0]+1)
        ybins = np.linspace(self.vmin[1], self.vmax[1], self.n[1]+1)
        zbins = np.linspace(self.vmin[2], self.vmax[2], self.n[2]+1)

        xc = (xbins[1:] + xbins[:-1])/2
        yc = (ybins[1:] + ybins[:-1])/2
        zc = (zbins[1:] + zbins[:-1])/2

        self.centers = np.array([xc, yc, zc], dtype=object)
        self.get_axes_grid()

        if not self.polar:
            xlabels = {'UV': '$U$', 'UW': '$U$', 'VW': '$V$'}
            ylabels = {'UV': '$V', 'UW': '$W$', 'VW': '$W$'}

            planes = planes[:3]
        elif self.polar:
            xlabels = {'rphi': '$v_r$', 'rtheta': '$v_r$', 'phitheta': '$v_\phi$'}
            ylabels = {'rphi': '$v_\phi$', 'rtheta': '$v_\\theta$', 'phitheta': '$v_\\theta$'}

            planes = planes[3:6]

        if plane == 'all':
            for plane in planes:
                twod_fv = self.get_fv(plane)

                plt.figure(figsize=(10, 10))
                ax = plt.gca()
                self.plot_fv(twod_fv, ax, plane)

                ax.set_axisbelow(False)
                ax.grid(which='major', axis='both', alpha=0.2, color='gray')
                ax.set_xlabel(f'{xlabels[plane]} [km s$^{{-1}}$]')
                ax.set_ylabel(f'{ylabels[plane]} [km s$^{{-1}}$]')

                plt.savefig('RUNS/' + self.folder + '/fvplot_' + plane + '.pdf',format='pdf')
        else:
            twod_fv = self.get_fv(plane)
            self.plot_fv(twod_fv, plt.gca(), plane)

            ax.set_axisbelow(False)
            ax.grid(which='major', axis='both', alpha=0.2, color='gray')
            ax.set_xlabel(f'{xlabels[plane]} [km s$^{{-1}}$]')
            ax.set_ylabel(f'{ylabels[plane]} [km s$^{{-1}}$]')

        return


    def get_axes_grid(self):
        '''
        Using the calculated centers of all cells from plot_and_save(), this function creates the meshgrids for all combinations of x and y axes that are relevant, depending on whether or not the data is polar.
        '''
        self.grid_x = {}
        self.grid_y = {}
        if self.polar:
            self.grid_x['rtheta'], self.grid_y['rtheta'] = np.meshgrid(self.centers[0], self.centers[1], indexing='ij')
            self.grid_x['rphi'], self.grid_y['rphi'] = np.meshgrid(self.centers[0], self.centers[2], indexing='ij')
            self.grid_x['phitheta'], self.grid_y['phitheta'] = np.meshgrid(self.centers[2], self.centers[1], indexing='ij')
        elif not self.polar:
            self.grid_x['UV'], self.grid_y['UV'] = np.meshgrid(self.centers[0], self.centers[1], indexing='ij')
            self.grid_x['UW'], self.grid_y['UW'] = np.meshgrid(self.centers[0], self.centers[2], indexing='ij')
            self.grid_x['VW'], self.grid_y['VW'] = np.meshgrid(self.centers[1], self.centers[2], indexing='ij')
        return


    def get_fv(self, plane):
        '''
        Since mxl is 3D, we need to perform the discrete integral over the dimension we are not looking in. This means we sum over one dimension to get a 2D projection of the probability.

        Parameters
        ----------
        plane : string
            The plane that is currently being plotted.

        Returns
        -------
        twodfv : Two-dimensional array
            The 2D probability distribution, integrated over the missing dimension.
        '''
        fv = np.exp(self.mxl)

        if plane in ('UV', 'rtheta'):
            twodfv = np.sum(fv, axis=2)
        elif plane in ('VW', 'phitheta'):
            if plane == 'VW':
                twodfv = np.sum(fv, axis=0)
            else:
                twodfv = np.sum(fv, axis=0).T
        elif plane in ('UW', 'rphi'):
            twodfv = np.sum(fv, axis=1)

        return twodfv


    def plot_fv(self, fv, ax, plane):
        '''
        The function that adds the specific contourf plot of the distribution on the current plot.

        Parameters
        ----------
        fv : Two-dimensional array
            The two-dimensional distribution in the current plane.

        ax : matplotlib.axes._subplots.AxesSubplot
            The axis upon which this distribution is being plotted.

        plane :
            The plane that is currently being plotted.

        '''
        fv_sort = np.sort(fv.ravel())
        levels = [fv_sort[np.searchsorted(np.cumsum(fv_sort), fv_sort.sum()*(100 - percent)/100)]
                  for percent in [95, 90, 80, 68, 50, 33, 21, 12, 6, 2]]

        im = ax.contourf(self.grid_x[plane], self.grid_y[plane], fv, 100, cmap='bone_r', linestyles='solid')
        for c in im.collections:
            c.set_rasterized(True)
            c.set_edgecolor('face')

        cs = ax.contour(self.grid_x[plane], self.grid_y[plane], fv, levels=np.unique(levels), colors=['k']*5 + ['w']*5, linewidths=0.85)
        for c in cs.collections:
            c.set_rasterized(True)



        return


    def plot_L_and_dL(self):
        '''
        This function is run at the end of every run to output the performance information on L and dL. It plots L and norm(dL) as a function of the number of iterations, giving insight into the minimization.
        '''
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
        plt.savefig('RUNS/' + self.folder + '/LanddLplot' + '.pdf',format='pdf')

        return


    def get_vars(self, folder):
        '''
        Function that reads the logfile to reconstruct the attributes of the class. Useful for when you're running things in a notebook rather than in the algorithm.

        Parameters
        ----------
        folder : string
            The directory in the RUNS folder that contains the logfile.

        '''

        self.folder = folder
        lines = open(f'RUNS/{self.folder}/log.txt').read()

        self.datafile = lines.split('\n')[0][14:]

        lines = lines.split('\n')[4].split(',')
        self.n = np.array([int(i) for i in lines[0][14:][1:-1].split()])
        self.vmin = np.array([float(i) for i in lines[1][2:-1].split()])
        self.dv = np.array([float(i) for i in lines[2][2:-1].split()])
        self.polar = int(lines[3])
        self.use_guess = int(lines[4])
        self.non_iso = int(lines[5])
        self.v_guess = np.array([float(i) for i in lines[6][2:-1].split()])
        self.disp_guess = np.array([float(i) for i in lines[7][2:-1].split()])
        self.alpha = float(lines[8])
        self.vmax = self.vmin + self.n*self.dv

        return
