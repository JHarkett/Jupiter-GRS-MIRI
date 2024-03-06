#===============================================================
"""
Combines the input post-retrieval tiles below into a final mosaic


Usage:
		python build_retrieval_mosaic.py



"""
#===============================================================
# Inputs


tiles = ['west', 'centre'] # In the order west centre east
ret_dir = 'results_stage2' # Directory where the results are stored
file = 'residual.fits'

invert = False # Normally == True, if it does not work, then use False

mos_dir = 'mosaics_stage2' # Location to save the final mosaic


#===============================================================
# Imports


import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import shutil
from astropy.io import fits


#===============================================================
# Functions


# Merges the mosaics in a more sophisticated way by analysing the
# overlap region between two tiles
def remove_outliers(bm, sdf, dz):

	for jj in range(dz):
		print(jj)

		arralc = np.where(np.logical_and(np.isnan(bm[0][jj]) == False, np.isnan(bm[1][jj]) == False))

		fraction = np.empty((dy2, dx2))
		fraction[:] = np.nan

		fraction[arralc] = np.divide(bm[0][jj][arralc], bm[1][jj][arralc])
		frac_med = np.nanmedian(fraction)

		print('Overlap 1 fraction = {}'.format(frac_med))

		if invert == False:
			arralc2 = np.where(fraction < frac_med)
			arralc3 = np.where(fraction >= frac_med)
		if invert == True:
			arralc3 = np.where(fraction < frac_med)
			arralc2 = np.where(fraction >= frac_med)

		sdf[jj][arralc2] = bm[1][jj][arralc2]
		sdf[jj][arralc3] = bm[0][jj][arralc3]

		if len_t == 3:
			arralc = np.where(np.logical_and(np.isnan(bm[1][jj]) == False, np.isnan(bm[2][jj]) == False))

			fraction = np.empty((dy2, dx2))
			fraction[:] = np.nan

			fraction[arralc] = np.divide(bm[1][jj][arralc], bm[2][jj][arralc])
			frac_med = np.nanmedian(fraction)

			print('Overlap 2 fraction = {}'.format(frac_med))

			if invert == False:
				arralc2 = np.where(fraction < frac_med)
				arralc3 = np.where(fraction >= frac_med)
			if invert == True:
				arralc3 = np.where(fraction < frac_med)
				arralc2 = np.where(fraction >= frac_med)

			sdf[jj][arralc2] = bm[2][jj][arralc2]
			sdf[jj][arralc3] = bm[1][jj][arralc3]

	return sdf


#===============================================================
# Main


if not os.path.exists(mos_dir):
	os.mkdir(mos_dir)

len_t = len(tiles)

sd1 = [None]*len_t
sd2 = [None]*len_t
sd3 = [None]*len_t
sd4 = [None]*len_t

se1 = [None]*len_t
se2 = [None]*len_t
se3 = [None]*len_t

lat_min = [None]*len_t
lat_max = [None]*len_t
lon_min = [None]*len_t
lon_max = [None]*len_t

lon_data = [None]*len_t
lat_data = [None]*len_t

dx = [None]*len_t
dy = [None]*len_t
dz = [None]*len_t


for kk in range(len_t):

	print(tiles[kk])

	hdul = fits.open(tiles[kk] + '/retrieval_data/' + ret_dir + '/' + file)

	if file == 'residual.fits':
		sd1[kk] = hdul['RESIDUAL'].data
		se1[kk] = hdul['ERR'].data

		sd2[kk] = hdul['SPECTRAL_DATA'].data
		sd3[kk] = hdul['FIT_DATA'].data
		sd4[kk] = hdul['WAVELENGTH'].data

		dz[kk] = len(sd1[kk])

	if file == '000_retrieved.fits':
		sd1[kk] = hdul['VARIABLE'].data
		se1[kk] = hdul['ERR'].data

		dz[kk] = len(sd1[kk])

	if file != 'residual.fits' and file != '000_retrieved.fits':
		sd1[kk] = hdul['VAL1'].data
		se1[kk] = hdul['ERR1'].data
		sd2[kk] = hdul['VAL2'].data
		se2[kk] = hdul['ERR2'].data

		if file == '-1032_retrieved.fits':
			sd3[kk] = hdul['VAL3'].data
			se3[kk] = hdul['ERR3'].data

		dz[kk] = 1

	hdul.close()

	hdul = fits.open(tiles[kk] + '/retrieval_data/' + ret_dir + '/000_retrieved.fits')
	lon_data[kk] = hdul['LON_WEST'].data
	lat_data[kk] = hdul['LAT_PGR'].data
	hdul.close()

	dy[kk], dx[kk] = lon_data[kk].shape

	lon_min[kk] = np.amin(lon_data[kk])
	lon_max[kk] = np.amax(lon_data[kk])
	lat_min[kk] = np.amin(lat_data[kk])
	lat_max[kk] = np.amax(lat_data[kk])

lat_min_total = int(np.amin(lat_min))
lat_max_total = int(np.amax(lat_max))
lon_min_total = int(np.amin(lon_min))
lon_max_total = int(np.amax(lon_max))


# Changed the size of the lon/lat grids so they are the same for july/aug
lat_array = np.arange(-31.0, -11.0, 0.5)
lon_array = np.arange(320.0, 273.0, -0.5)

lon_grid, lat_grid = np.meshgrid(lon_array, lat_array)

dy2, dx2 = lon_grid.shape

print('\n')

bm1 = np.empty((len_t, dz[0], dy2, dx2))
bm2 = np.empty((len_t, dz[0], dy2, dx2))
bm3 = np.empty((len_t, dz[0], dy2, dx2))
bm4 = np.empty((len_t, dz[0], dy2, dx2))

be1 = np.empty((len_t, dz[0], dy2, dx2))
be2 = np.empty((len_t, dz[0], dy2, dx2))
be3 = np.empty((len_t, dz[0], dy2, dx2))

bm1[:] = np.nan
bm2[:] = np.nan
bm3[:] = np.nan
bm4[:] = np.nan

be1[:] = np.nan
be2[:] = np.nan
be3[:] = np.nan

for kk in range(len_t):
	arral = np.where(np.logical_and(lon_grid == lon_max[kk], lat_grid == lat_min[kk]))

	if file == 'residual.fits':
		bm1[kk, :, arral[0][0]:arral[0][0]+dy[kk], arral[1][0]:arral[1][0]+dx[kk]] = sd1[kk]
		be1[kk, :, arral[0][0]:arral[0][0]+dy[kk], arral[1][0]:arral[1][0]+dx[kk]] = se1[kk]

		bm2[kk, :, arral[0][0]:arral[0][0]+dy[kk], arral[1][0]:arral[1][0]+dx[kk]] = sd2[kk]
		bm3[kk, :, arral[0][0]:arral[0][0]+dy[kk], arral[1][0]:arral[1][0]+dx[kk]] = sd3[kk]
		bm4[kk, :, arral[0][0]:arral[0][0]+dy[kk], arral[1][0]:arral[1][0]+dx[kk]] = sd4[kk]

	if file == '000_retrieved.fits':
		bm1[kk, :, arral[0][0]:arral[0][0]+dy[kk], arral[1][0]:arral[1][0]+dx[kk]] = sd1[kk]
		be1[kk, :, arral[0][0]:arral[0][0]+dy[kk], arral[1][0]:arral[1][0]+dx[kk]] = se1[kk]

	if file != 'residual.fits' and file != '000_retrieved.fits':
		bm1[kk, :, arral[0][0]:arral[0][0]+dy[kk], arral[1][0]:arral[1][0]+dx[kk]] = sd1[kk]
		be1[kk, :, arral[0][0]:arral[0][0]+dy[kk], arral[1][0]:arral[1][0]+dx[kk]] = se1[kk]
		bm2[kk, :, arral[0][0]:arral[0][0]+dy[kk], arral[1][0]:arral[1][0]+dx[kk]] = sd2[kk]
		be2[kk, :, arral[0][0]:arral[0][0]+dy[kk], arral[1][0]:arral[1][0]+dx[kk]] = se2[kk]

		if file == '-1032_retrieved.fits':
			bm3[kk, :, arral[0][0]:arral[0][0]+dy[kk], arral[1][0]:arral[1][0]+dx[kk]] = sd3[kk]
			be3[kk, :,  arral[0][0]:arral[0][0]+dy[kk], arral[1][0]:arral[1][0]+dx[kk]] = se3[kk]

sdf1 = np.nanmedian(bm1, axis=0)
sdf2 = np.nanmedian(bm2, axis=0)
sdf3 = np.nanmedian(bm3, axis=0)
sdf4 = np.nanmedian(bm4, axis=0)

sef1 = np.nanmedian(be1, axis=0)
sef2 = np.nanmedian(be2, axis=0)
sef3 = np.nanmedian(be3, axis=0)


if file == 'residual.fits':
	sdf1 = remove_outliers(bm1, sdf1, dz[0])
	sef1 = remove_outliers(be1, sef1, dz[0])

	sdf2 = remove_outliers(bm2, sdf2, dz[0])
	sdf3 = remove_outliers(bm3, sdf3, dz[0])

	primary_hdu = fits.PrimaryHDU()
	residual_hdu = fits.ImageHDU(sdf1, name='RESIDUAL')
	wave_err_hdu = fits.ImageHDU(sef1, name='ERR')
	spec_data_hdu = fits.ImageHDU(sdf2, name='SPECTRAL_DATA')
	fit_data_hdu =  fits.ImageHDU(sdf3, name='FIT_DATA')
	wave_hdu = fits.ImageHDU(sdf4, name='WAVELENGTH')
	lon_hdu = fits.ImageHDU(lon_grid, name='LON_WEST')
	lat_hdu = fits.ImageHDU(lat_grid, name='LAT_PGR')

	hdul = fits.HDUList([primary_hdu, residual_hdu, wave_err_hdu, spec_data_hdu, fit_data_hdu, wave_hdu, lon_hdu, lat_hdu])
	hdul.writeto(mos_dir + '/residual.fits', overwrite=True)


if file == '000_retrieved.fits':
	sdf1 = remove_outliers(bm1, sdf1, dz[0])
	sef1 = remove_outliers(be1, sef1, dz[0])

	primary_hdu = fits.PrimaryHDU()
	data_hdu = fits.ImageHDU(sdf1, name='VARIABLE')
	err_hdu = fits.ImageHDU(sef1, name='ERR')
	lon_hdu = fits.ImageHDU(lon_grid, name='LON_WEST')
	lat_hdu = fits.ImageHDU(lat_grid, name='LAT_PGR')

	hdul = fits.HDUList([primary_hdu, data_hdu, err_hdu, lon_hdu, lat_hdu])
	hdul.writeto(mos_dir + '/000_retrieved.fits', overwrite=True)

if file != 'residual.fits' and file != '000_retrieved.fits':
	sdf1 = remove_outliers(bm1, sdf1, dz[0])
	sef1 = remove_outliers(be1, sef1, dz[0])
	sdf2 = remove_outliers(bm2, sdf2, dz[0])
	sef2 = remove_outliers(be2, sef2, dz[0])

	primary_hdu = fits.PrimaryHDU()
	hdu1 = fits.ImageHDU(sdf1, name='VAL1')
	hdu2 = fits.ImageHDU(sef1, name='ERR1')
	hdu3 = fits.ImageHDU(sdf2, name='VAL2')
	hdu4 = fits.ImageHDU(sef2, name='ERR2')

	lon_hdu = fits.ImageHDU(lon_grid, name='LON_WEST')
	lat_hdu = fits.ImageHDU(lat_grid, name='LAT_PGR')

	if file == '-1032_retrieved.fits':
		sdf3 = remove_outliers(bm3, sdf3, dz[0])
		sef3 = remove_outliers(be3, sef3, dz[0])

		hdu5 = fits.ImageHDU(sdf3, name='VAL3')
		hdu6 = fits.ImageHDU(sef3, name='ERR3')

		hdul = fits.HDUList([primary_hdu, hdu1, hdu2, hdu3, hdu4, hdu5, hdu6, lon_hdu, lat_hdu])

	else:
		hdul = fits.HDUList([primary_hdu, hdu1, hdu2, hdu3, hdu4, lon_hdu, lat_hdu])

	hdul.writeto(mos_dir + '/' + file, overwrite=True)


#===============================================================


print('End of script\n')
