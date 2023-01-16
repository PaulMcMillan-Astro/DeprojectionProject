import numpy as np
from astropy.stats import bootstrap
import astropy.coordinates as coord
from astropy.table import QTable
import astropy.units as u

class DataReader:
    def __init__(self, inifile, resample):
        self.file_ = inifile
        self.b = resample
        self.read_params()


    def resample(self, plx_err, pmra_err, pmdec_err, plx_pmra_corr, plx_pmdec_corr, pmra_pmdec_corr):
        '''Use the uncertainties and correlations between plx, pmra, pmdec to resample the values and replace the originals in the Class.


        Parameters
        ----------
        plx_err : astropy.units.Quantity  (N, )
            The parallax uncertainties

        pmra_err : astropy.units.Quantity  (N, )
            The right ascension proper motion uncertainties

        pmdec_err : astropy.units.Quantity  (N, )
            The declination proper motion uncertainties

        plx_pmra_corr : astropy.units.Quantity  (N, )
            The correlation coefficients between pm_RA and plx

        plx_pmdec_corr : astropy.units.Quantity  (N, )
            The correlation coefficients between pm_dec and plx

        pmra_pmdec_corr : astropy.units.Quantity  (N, )
            The correlation coefficients between pm_RA and pm_DEC

        '''
        rng = np.random.default_rng()

        plx_err = plx_err.value
        pmra_err = pmra_err.value
        pmdec_err = pmdec_err.value

        for i in range(len(self.plx)):
            mean = [self.plx[i].value, self.pm_RA[i].value, self.pm_DEC[i].value]
            cov = np.array([
                           [plx_err[i]**2,
                            plx_pmra_corr[i]*plx_err[i]*pmra_err[i],
                            plx_pmdec_corr[i]*plx_err[i]*pmdec_err[i]],
                           [plx_pmra_corr[i]*plx_err[i]*pmra_err[i],
                            pmra_err[i]**2,
                            pmra_pmdec_corr[i]*pmra_err[i]*pmdec_err[i]],
                           [plx_pmdec_corr[i]*plx_err[i]*pmdec_err[i],
                            pmra_pmdec_corr[i]*pmra_err[i]*pmdec_err[i],
                            pmdec_err[i]**2]
                           ])

            self.plx[i], self.pm_RA[i], self.pm_DEC[i] = rng.multivariate_normal(mean, cov)*np.array([self.plx.unit, self.pm_RA.unit, self.pm_DEC.unit])
        return


    def bootstrap_sample(self, pvals, rhatvals):
        '''Performs a bootstrap resample on pvals and rhatvals. This is used for estimating the 1 sigma uncertainties of the warp estimation in Paper 2.

        Parameters
        ----------
        pvals : Array of shape  (N, 3)
            The 3-D components of the transverse velocity. Either cartesian or polar.

        rhatvals : Array of shape  (N, 3)
            The unit vectors in the direction of the stars on the celestial sphere in 3-D coordinates. Either cartesian or polar.

        Returns
        -------
        pvals : Array of shape  (N, 3)
            Bootstrapped version of pvals.

        rhatvals : Array of shape  (N, 3)
            Bootstrapped version of rhatvals.
        '''
        pvals, rhatvals = bootstrap(np.dstack([pvals, rhatvals]), bootnum=1).transpose(3, 0, 1, 2)
        return np.squeeze(pvals), np.squeeze(rhatvals)


    def process_input(self):
        '''Reads the input .ini file. Removes obsolete lines and structures things.

        Returns
        -------
        file : list
            Contains a list where each entry is one of the lines of input from the .ini file.
        '''
        file = open(self.file_, 'r').read()
        file = file.split('\n') # Splits into list with each line separately
        file = [_ for _ in file if not _.startswith('#')] # Removes lines starting with #
        file = [_ for _ in file if not len(_) == 0] # Removes empty lines
        return file


    def read_params(self):
        '''Takes the input .ini file and places all the values in their correct container. Reads the VOTable and extracts the data. If requested, also calls the resample() function.
        '''
        # Read the .ini file
        vars_ = self.process_input()

        # Assign containers
        self.n = np.array(vars_[0].split(',')[:],dtype=int)
        self.vmin = np.array(vars_[1].split(',')[:],dtype=float)
        self.dv = np.array(vars_[2].split(',')[:],dtype=float)
        self.polar = bool(int(vars_[3]))
        self.use_guess = bool(int(vars_[4]))
        self.non_iso = bool(int(vars_[5]))
        self.v_guess = np.array(vars_[6].split(',')[:],dtype=float)
        self.disp_guess = np.array(vars_[7].split(',')[:],dtype=float)
        self.alpha = float(vars_[8])
        self.datafile = vars_[9]

        # Read the VOTable
        if self.datafile.endswith('.vot'):
            data_raw = QTable.read('DATA/' + self.datafile, format='votable')
            self.RA = data_raw['ra']
            self.DEC = data_raw['dec']
            self.plx = data_raw['parallax']
            self.pm_RA = data_raw['pmra']
            self.pm_DEC = data_raw['pmdec']

            if bool(self.b):
                plx_err = data_raw['parallax_error']
                pm_RA_err = data_raw['pmra_error']
                pm_DEC_err = data_raw['pmdec_error']
                plx_pmra_corr = data_raw['parallax_pmra_corr']
                plx_pmdec_corr = data_raw['parallax_pmdec_corr']
                pmra_pmdec_corr = data_raw['pmra_pmdec_corr']

                self.resample(plx_err, pm_RA_err, pm_DEC_err, plx_pmra_corr, plx_pmdec_corr, pmra_pmdec_corr)

        else:
            raise TypeError('Input data is not a .vot file!!')

        print(f'Sample has {len(data_raw)} stars\n')
        return


    def create_sample(self):
        '''Using astropy.coordinates.SkyCoord we set up a SkyCoord in the 'icrs' frame from our data. By setting all the radial velocities to zero, we are also able to convert it to a galactic coordinate frame. From the galactic frame, we can extract 3D velocities (U, V, W) and unit vectors pointing towards the stars. This corresponds to the bold /p/ and /rhat/ from Eq. 5 of Dehnen 1998.

        Returns
        -------
        pvals : Array of shape  (N, 3)
            The 3-D components of the transverse velocity in U, V, W

        rhatvals : Array of shape  (N, 3)
            The unit vectors in the direction of the stars on the celestial sphere in 3-D cartesian coordinates (in the same directions as U, V, W)
        '''
        self.sample = coord.SkyCoord(ra=self.RA,
                                dec=self.DEC,
                                pm_ra_cosdec=self.pm_RA,
                                pm_dec=self.pm_DEC,
                                distance=coord.Distance(parallax=self.plx),
                                radial_velocity=np.zeros(len(self.plx))*u.km/u.s, # This is necessary to extract pvals from astropy
                                frame='icrs').galactic

        pvals = self.sample.velocity.d_xyz.T.unmasked.value
        rhatvals = self.sample.spherical.unit_vectors()['distance'].xyz.T.unmasked.value
        return pvals, rhatvals


    def make_polar(self, pvals, rhatvals):
        '''Using the astropy.coordinates.SkyCoord from the create_sample() function, we convert it to a Galactocentric cartesian frame (rather than heliocentric). We can then extract Galactocentric coordinates and convert them to spherical coordinates. We create the rotation matrix. We then apply it to our old rhatvals and the sum of the old pvals with the motion of the sun added back in.

        Returns
        -------
        pvals : Array of shape  (N, 3)
            The 3-D components of the transverse velocity in vr, vtheta, vphi

        rhatvals : Array of shape  (N, 3)
            The unit vectors in the direction of the stars on the celestial sphere (seen from the sun) converted to 3-D polar coordinates (in the same directions as vr, vtheta, vphi)
        '''
        sample_gc = self.sample.galactocentric
        x = sample_gc.x
        y = sample_gc.y
        z = sample_gc.z

        r = np.sqrt(x**2 + y**2 + z**2)
        theta = np.arccos(z/r)
        phi = np.arctan2(y, x)

        vsun = sample_gc.galcen_v_sun.d_xyz

        R = np.array([[np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)],
                      [np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -np.sin(theta)],
                      [-np.sin(phi), np.cos(phi), np.zeros(len(phi))]]).transpose(2, 0, 1)

        pvals_polar = ( R @ (vsun.value + pvals)[..., np.newaxis]).squeeze()
        rhat_polar = ( R @ rhatvals[..., np.newaxis]).squeeze()

        return pvals_polar, rhat_polar
