#========================================================================================
"""
Generates a nyquist-sampled grid of the data, combining the dithers
to a user-defined lon/lat scale

Usage:
		python gen_nyquist.py

"""
#========================================================================================
# Inputs


# Bands to run scaling on (will process all to the same lon/lat grid)
bands = ['ch1-short', 'ch2-short']

working_dir = 'aug/east/'   # Location of dither files
stage_dir = 'stage3_desaturated/'

out_dir = 'd_combined_0.25'

scale = 0.25       # In degrees

defringe = True # If True, will use de-fringed data


#========================================================================================
# Imports


import numpy as np
from astropy.io import fits
import glob
import os
import shutil
import matplotlib.pyplot as plt


#========================================================================================
# Main

sstring = working_dir + stage_dir + 'd*_fringe/'

dithers = sorted(glob.glob(sstring))
len_d = len(dithers)
len_b = len(bands)

if defringe == False:
	for kk in range(len_d):
		dithers[kk] = dithers[kk].replace('_fringe', '')

print('\n{} bands detected'.format(len_b))
print('{} dither positions detected\n'.format(len_d))

if os.path.exists(working_dir + '/' + out_dir):
	shutil.rmtree(working_dir + '/' + out_dir)
os.mkdir(working_dir + '/' + out_dir)


hdul = fits.open(dithers[0] + 'Level3_' + bands[len_b-1] + '_s3d_nav.fits')
spec_data_pre = hdul['SCI'].data
hdul.close()

print('Using\n' + bands[len_b-1] + '\nTo generate the lon/lat files\n')

dzp, dyp, dxp = spec_data_pre.shape

lon = np.zeros((len_d, dyp, dxp))
lat = np.zeros((len_d, dyp, dxp))

lat_min = np.zeros((len_d))
lat_max = np.zeros((len_d))
lon_min = np.zeros((len_d))
lon_max = np.zeros((len_d))

for kk in range(len_d):

	hdul = fits.open(dithers[kk] + 'Level3_' + bands[len_b-1] + '_s3d_nav.fits')

	lon[kk] = hdul['LON'].data
	lat[kk] = hdul['LAT_GRAPHIC'].data

	hdul.close()

	lat_min[kk] = np.amin(lat[kk])
	lat_max[kk] = np.amax(lat[kk])
	lon_min[kk] = np.amin(lon[kk])
	lon_max[kk] = np.amax(lon[kk])

lat_min_total = int(np.amin(lat_min))-1.0
lat_max_total = int(np.amax(lat_max))+1.0
lon_min_total = int(np.amin(lon_min))-1.0
lon_max_total = int(np.amax(lon_max))+1.0

lat_array = np.arange(lat_min_total, lat_max_total + scale, scale)
lon_array = np.arange(lon_max_total, lon_min_total - scale, np.multiply(scale, -1))

lon_grid, lat_grid = np.meshgrid(lon_array, lat_array)

hdu111 = fits.PrimaryHDU(lon_grid)
hdu111.writeto(working_dir + '/' + out_dir + '/longitude_grid.fits', overwrite=True)

hdu222 = fits.PrimaryHDU(lat_grid)
hdu222.writeto(working_dir + '/' + out_dir + '/latitude_grid.fits', overwrite=True)


dy2, dx2 = lon_grid.shape


for jj in range(len_b):

	print(bands[jj])

	hdul = fits.open(dithers[0] + 'Level3_' + bands[jj] + '_s3d_nav.fits')
	spec_data_pre = hdul['SCI'].data
	hdr1 = hdul[0].header
	hdr2 = hdul['SCI'].header

	wave = np.arange(hdr2['NAXIS3'])*hdr2['CDELT3']+hdr2['CRVAL3']

	dz, dy, dx = spec_data_pre.shape
	hdul.close()

	spec_data = np.zeros((len_d, dz, dy, dx))
	spec_err = np.zeros((len_d, dz, dy, dx))

	lon = np.zeros((len_d, dy, dx))
	lat = np.zeros((len_d, dy, dx))

	sa = np.zeros((len_d, dy, dx))
	ea = np.zeros((len_d, dy, dx))
	aa = np.zeros((len_d, dy, dx))

	rv = np.zeros((len_d, dy, dx))
	dopp = np.zeros((len_d, dy, dx))

	for kk in range(len_d):
		hdul = fits.open(dithers[kk] + 'Level3_' + bands[jj] + '_s3d_nav.fits')
		spec_data[kk] = hdul['SCI'].data
		spec_err[kk] = hdul['ERR'].data

		lon[kk] = hdul['LON'].data
		lat[kk] = hdul['LAT_GRAPHIC'].data

		sa[kk] = hdul['INCIDENCE'].data
		ea[kk] = hdul['EMISSION'].data
		aa[kk] = hdul['AZIMUTH'].data

		rv[kk] = hdul['RADIAL_VELOCITY'].data
		dopp[kk] = hdul['DOPPLER'].data

		hdul.close()

	data_nq = np.empty((dz, dy2, dx2))
	data_nq[:] = np.nan

	err_nq = np.empty((dz, dy2, dx2))
	err_nq[:] = np.nan

	sa_nq = np.empty((dy2, dx2))
	sa_nq[:] = np.nan

	ea_nq = np.empty((dy2, dx2))
	ea_nq[:] = np.nan

	aa_nq = np.empty((dy2, dx2))
	aa_nq[:] = np.nan

	rv_nq = np.empty((dy2, dx2))
	rv_nq[:] = np.nan

	dopp_nq = np.empty((dy2, dx2))
	dopp_nq[:] = np.nan

	for j in range(dy2):
		for i in range(dx2):

			arral = np.where(np.logical_and(\
			np.logical_and(lon >= lon_grid[j, i] - 0.75, lon < lon_grid[j, i] + 0.75),\
			np.logical_and(lat >= lat_grid[j, i] - 0.75, lat < lat_grid[j, i] + 0.75)))

			if len(arral[0]) == 0:

				data_nq[:, j, i] = np.nan
				err_nq[:, j, i] = np.nan

				sa_nq[j, i] = np.nan
				ea_nq[j, i] = np.nan
				aa_nq[j, i] = np.nan

				rv_nq[j, i] = np.nan
				dopp_nq[j, i] = np.nan

			else:

				data_nq[:, j, i] = np.nanmean(spec_data[arral[0], :, arral[1], arral[2]], axis=0)
				err_nq[:, j, i] = np.nanmean(spec_err[arral[0], :, arral[1], arral[2]], axis=0)

				sa_nq[j, i] = np.nanmean(sa[arral], axis=0)
				ea_nq[j, i] = np.nanmean(ea[arral], axis=0)
				aa_nq[j, i] = np.nanmean(aa[arral], axis=0)

				rv_nq[j, i] = np.nanmean(rv[arral], axis=0)
				dopp_nq[j, i] = np.nanmean(dopp[arral], axis=0)

	hdr1['PATT_NUM'] = 'ALL'
	hdr2['NAXIS1'] = dx2
	hdr2['NAXIS2'] = dy2

	primary_hdu = fits.PrimaryHDU(header=hdr1)

	sci_hdu = fits.ImageHDU(data_nq, header=hdr2, name='SCI')
	err_hdu = fits.ImageHDU(err_nq, name='ERR')

	lon_hdu = fits.ImageHDU(lon_grid, name='LON_WEST')
	lat_hdu = fits.ImageHDU(lat_grid, name='LAT_PGR')

	sa_hdu = fits.ImageHDU(sa_nq, name='SOLAR_ANGLE')
	ea_hdu = fits.ImageHDU(ea_nq, name='EMISSION_ANGLE')
	aa_hdu = fits.ImageHDU(aa_nq, name='AZIMUTH_ANGLE')

	rv_hdu = fits.ImageHDU(rv_nq, name='RADIAL_VELOCITY')
	dopp_hdu = fits.ImageHDU(dopp_nq, name='DOPPLER')

	hdul_all = fits.HDUList([primary_hdu, sci_hdu, err_hdu, lon_hdu, lat_hdu, sa_hdu, ea_hdu, aa_hdu, rv_hdu, dopp_hdu])
	hdul_all.writeto(working_dir + '/' + out_dir + '/Level3_' + bands[jj] + '_s3d_nav.fits', overwrite=True)


#========================================================================================


print('End of script\n')
