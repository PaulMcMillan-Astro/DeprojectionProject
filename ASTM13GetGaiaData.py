#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 14:55:14 2018

@author: David Hobbs, Lund Observatory
"""

import numpy as np
import coordTransform as coordT

import matplotlib
import matplotlib.pyplot as plt
from astroquery.gaia import Gaia
from astropy.io.votable import parse_single_table
from matplotlib.ticker import (MultipleLocator)

deg2rad = np.pi/180.0


source = 'gaia7G_20180705'
source = 'ESA'

if(source == 'ESA'):
    #coord = SkyCoord(ra=280, dec=-60, unit=(u.degree, u.degree), frame='icrs')
    #radius = u.Quantity(1.0, u.arcminute)
    #j = Gaia.cone_search_async(coord, radius)
    #r = j.get_results()
    #r.pprint()

    tables = Gaia.load_tables(only_names=True)
    #for table in (tables):
    #    print (table.get_qualified_name())

    # An asynchronous query (asynchronous rather than synchronous queries should be performed when 
    #retrieving more than 2000 rows) centred on the Pleides (coordinates: 56.75, +24.1167) with a 
    #search radius of 1 degrees and save the results to a file.

    job = Gaia.launch_job_async("SELECT * FROM gaiadr2.gaia_source \
                                WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS', 56.75, 24.1167, 180))=1 \
                                AND abs(parallax_error/parallax)<0.10 \
                                AND parallax>0.0 \
                                AND abs(pmra_error/pmra)<0.10 \
                                AND abs(pmdec_error/pmdec)<0.10 \
                                AND phot_g_mean_mag<7.0 \
                                AND parallax IS NOT NULL \
                                AND pmra IS NOT NULL \
                                AND pmdec IS NOT NULL \
                                AND phot_g_mean_mag IS NOT NULL \
                                AND phot_bp_mean_mag IS NOT NULL \
                                AND phot_rp_mean_mag IS NOT NULL \
                                AND abs(pmdec)>0 \
                                ;",dump_to_file=True)

    print (job)    

    data = job.get_results()
else:
    table = parse_single_table(source+'.vot')
    data = table.array

print('Ready now')

# Make them proper arrays
c = np.asarray(data['bp_rp']);
#G = −2.5 log(flux) + zeropoint
#  = -2.5*log10(550.00816) + 25.691439
G = np.asarray(data['phot_g_mean_mag']);
BP = np.asarray(data['phot_bp_mean_mag']);
RP = np.asarray(data['phot_rp_mean_mag']);
ra = np.asarray(data['ra']);
ra_error = np.asarray(data['ra_error']);
dec = np.asarray(data['dec']);
dec_error = np.asarray(data['dec_error']);
parallax = np.asarray(data['parallax']);
parallax_error = np.asarray(data['parallax_error']);
pmra = np.asarray(data['pmra']);
pmra_error = np.asarray(data['pmra_error']);                  
pmdec = np.asarray(data['pmdec']);
pmdec_error = np.asarray(data['pmdec_error']);

# Statistical indicators to use for selecting sources
#-----------------------------------------------------
# astrometric_n_good_obs_al
# astrometric_gof_al
# astrometric_chi2_al
# astrometric_excess_noise
# astrometric_excess_noise_sig
# astrometric_primary_flag
# astrometric_matched_observations
# visibility_periods_used

l, b, varpi, mulStar, mub, deltalStar, deltab, deltaVarpi, deltaMulStar, deltaMub = \
      coordT.transformIcrsGal(ra, dec, parallax, pmra, pmdec, \
                              ra_error, dec_error, parallax_error, pmra_error, pmdec_error);
                              
print('\nRetrieved %d sources\n' % len(parallax))                           