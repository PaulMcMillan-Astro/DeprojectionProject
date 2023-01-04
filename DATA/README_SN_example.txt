SN_example is a test sample you can run with 20 000 stars within 100 pc. It is created by the following query:

select TOP 20000 source_id, bp_rp, phot_g_mean_mag, phot_bp_rp_excess_factor, ruwe, ra, dec, parallax, pmra, pmdec, parallax_error, pmra_error, pmdec_error, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, visibility_periods_used, astrometric_chi2_al, astrometric_n_good_obs_al, radial_velocity, radial_velocity_error
from gaiadr3.gaia_source
where parallax_over_error > 20
and parallax > 10
and ruwe < 1.15
and phot_g_mean_flux_over_error > 50
and phot_rp_mean_flux_over_error > 20
and phot_bp_mean_flux_over_error > 20
and visibility_periods_used > 8
and astrometric_chi2_al/(astrometric_n_good_obs_al-5) < 1.44*greatest(1,exp(-0.4*(phot_g_mean_mag-19.5)))
and radial_velocity IS NOT NULL
