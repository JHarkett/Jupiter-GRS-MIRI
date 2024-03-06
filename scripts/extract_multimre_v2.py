#=====================================================
# Inputs


spec_dir = 'stage2'                # Location of NEMESIS outputs
out_dir = 'results_stage2'         # Desired name of output directory

ref_file = 'jupiter_v2022miri.ref' # Ref file used for retrievals

xsc_file = 'grs_single.xsc'
aero_file = 'aerosol.ref'
prior_file = 'jupiter_v2022miri.ref'
ref_skip = 25

aero_unit = 'opacity/bar'

display_figures = True


#=====================================================
# Description
"""
Extracts and plots variables from an mre file
(Optional) compares to the mirisim input and prior profile
Also generates a new apr file of the results for further retrievals

Usage

	Run extract_chisq_v2.py first

	[Input variables above]

     	python extract_multimre_v2.py


"""
#=====================================================
# Imports


import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits
import sys

sys.path.insert(1, '/data/nemesis/jh852/python/subroutines')
from manifold import *


#=====================================================
# Functions




#=====================================================
# Main


id_data = np.loadtxt(spec_dir + '/id_conserve.txt')
id = id_data[:,0]
lon = id_data[:,1]
lat = id_data[:,2]
id_length = len(id)

ref_pix = int(id[int(id_length/2)])
mre_file = spec_dir + '/core{}/core_1/nemesis.mre'.format(ref_pix)

nwave, nvar, id_name, id_array, v_start, v_length = detect_variables(mre_file)

#plot_wave(nwave, mre_file)

hdul = fits.open(out_dir + '/lon_lat.fits')
lon_grid = hdul['LON_WEST'].data
lat_grid = hdul['LAT_PGR'].data
hdul.close()

lon_arr, lat_arr = gen_lonlat_arr(lon_grid, lat_grid)

dy, dx = lon_grid.shape


wave_res = np.empty((nwave, dy, dx))
wave_res[:] = np.nan

wave_err = np.empty((nwave, dy, dx))
wave_err[:] = np.nan

residual = np.empty((nwave, dy, dx))
residual[:] = np.nan

spec_data = np.empty((nwave, dy, dx))
spec_data[:] = np.nan

fit_data = np.empty((nwave, dy, dx))
fit_data[:] = np.nan

# Extract the residuals
for kk in range(id_length):
	print(kk+1, lon[kk], lat[kk])

	try:
		res_data = np.loadtxt(spec_dir + '/core{}/core_1/nemesis.mre'.format(kk+1), skiprows=5, max_rows=nwave, usecols=(1, 2, 3, 5))

		wave_res[:, np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = res_data[:,0]
		residual[:, np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = res_data[:,3] - res_data[:,1]
		wave_err[:, np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = res_data[:,2]
		spec_data[:, np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = res_data[:,1]
		fit_data[:, np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = res_data[:,3]


	except:
		print('WARNING\nResidual data not found\n')

primary_hdu = fits.PrimaryHDU()
residual_hdu = fits.ImageHDU(residual, name='RESIDUAL')
wave_err_hdu = fits.ImageHDU(wave_err, name='ERR')
spec_data_hdu = fits.ImageHDU(spec_data, name='SPECTRAL_DATA')
fit_data_hdu =	fits.ImageHDU(fit_data, name='FIT_DATA')
wave_hdu = fits.ImageHDU(wave_res, name='WAVELENGTH')
lon_hdu = fits.ImageHDU(lon_grid, name='LON_WEST')
lat_hdu = fits.ImageHDU(lat_grid, name='LAT_PGR')

hdul = fits.HDUList([primary_hdu, residual_hdu, wave_err_hdu, spec_data_hdu, fit_data_hdu, wave_hdu, lon_hdu, lat_hdu])
hdul.writeto(out_dir + '/residual.fits', overwrite=True)


# Extract the variables

alt_data = np.loadtxt(prior_file, skiprows=ref_skip, usecols=(0, 1))
H = alt_data[:,0] # km
P = alt_data[:,1] # atm
P *= 1.01325 # bar

for jj in range(nvar):

	print(id_name[jj])

	grs_var_ret = np.empty((120, dy, dx))
	grs_var_ret[:] = np.nan

	grs_var_ret1 = np.empty((120, dy, dx))
	grs_var_ret1[:] = np.nan

	grs_var_ret2 = np.empty((120, dy, dx))
	grs_var_ret2[:] = np.nan

	grs_ret_err = np.empty((120, dy, dx))
	grs_ret_err[:] = np.nan


	grs_mre_val1 = np.empty((dy, dx))
	grs_mre_val1[:] = np.nan

	grs_mre_val2 = np.empty((dy, dx))
	grs_mre_val2[:] = np.nan

	grs_mre_val3 = np.empty((dy, dx))
	grs_mre_val3[:] = np.nan

	grs_mre_val4 = np.empty((dy, dx))
	grs_mre_val4[:] = np.nan

	grs_mre_val5 = np.empty((dy, dx))
	grs_mre_val5[:] = np.nan

	grs_mre_val6 = np.empty((dy, dx))
	grs_mre_val6[:] = np.nan


	for kk in range(id_length):
		print(kk+1, lon[kk], lat[kk])

		try:

			if id_array[jj][0] == '-1':
				mro_data = np.loadtxt(spec_dir + '/core{}/core_1/nemesis.mre'.format(kk+1), skiprows=v_start[jj], max_rows=v_length[jj], usecols=(4, 5))

				grs_mre_val1[np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = mro_data[:,0][0]
				grs_mre_val2[np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = mro_data[:,0][1]
				grs_mre_val3[np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = mro_data[:,0][2]
				grs_mre_val4[np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = mro_data[:,1][0]
				grs_mre_val5[np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = mro_data[:,1][1]
				grs_mre_val6[np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = mro_data[:,1][2]


			if id_array[jj][0] == '0':
				ret_data = np.loadtxt(spec_dir + '/core{}/core_1/nemesis.mre'.format(kk+1), skiprows=v_start[jj], max_rows=v_length[jj], usecols=(4, 5))

				grs_var_ret[:, np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = ret_data[:,0]
				grs_ret_err[:, np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = ret_data[:,1]

			if id_array[jj][0] == '11':
				mro_data = np.loadtxt(spec_dir + '/core{}/core_1/nemesis.mre'.format(kk+1), skiprows=v_start[jj], max_rows=v_length[jj], usecols=(4, 5))

				grs_mre_val1[np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = mro_data[:,0][0]
				grs_mre_val2[np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = mro_data[:,0][1]
				grs_mre_val3[np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = mro_data[:,1][0]
				grs_mre_val4[np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = mro_data[:,1][1]

			if id_array[jj][0] == '28':
				mro_data = np.loadtxt(spec_dir + '/core{}/core_1/nemesis.mre'.format(kk+1), skiprows=v_start[jj], max_rows=v_length[jj], usecols=(4, 5))

				grs_mre_val1[np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = mro_data[:,0][0]
				grs_mre_val2[np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = mro_data[:,0][1]
				grs_mre_val3[np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = mro_data[:,1][0]
				grs_mre_val4[np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = mro_data[:,1][1]
				print(mro_data[:,0][0])

		except:
			print('WARNING\nNo data found\n')


	if id_array[jj][0] == '0':
		primary_hdu = fits.PrimaryHDU()
		data_hdu = fits.ImageHDU(grs_var_ret, name='VARIABLE')
		err_hdu = fits.ImageHDU(grs_ret_err, name='ERR')
		lon_hdu = fits.ImageHDU(lon_grid, name='LON_WEST')
		lat_hdu = fits.ImageHDU(lat_grid, name='LAT_PGR')

		hdul = fits.HDUList([primary_hdu, data_hdu, err_hdu, lon_hdu, lat_hdu])

		hdul.writeto(out_dir + '/{}_retrieved.fits'.format(id_name[jj]), overwrite=True)

	if id_array[jj][0] == '-1'or id_array[jj][0] == '-2':
		primary_hdu = fits.PrimaryHDU()
		hdu1 = fits.ImageHDU(grs_mre_val3, name='VAL1')
		hdu2 = fits.ImageHDU(grs_mre_val6, name='ERR1')
		hdu3 = fits.ImageHDU(grs_mre_val1, name='VAL2')
		hdu4 = fits.ImageHDU(grs_mre_val4, name='ERR2')
		hdu5 = fits.ImageHDU(grs_mre_val2, name='VAL3')
		hdu6 = fits.ImageHDU(grs_mre_val5, name='ERR3')
		hdul = fits.HDUList([primary_hdu, hdu1, hdu2, hdu3, hdu4, hdu5, hdu6])

		hdul.writeto(out_dir + '/{}_retrieved.fits'.format(id_name[jj]), overwrite=True)

	if id_array[jj][0] == '28' or id_array[jj][0] == '11':
		primary_hdu = fits.PrimaryHDU()
		hdu1 = fits.ImageHDU(grs_mre_val1, name='VAL1')
		hdu2 = fits.ImageHDU(grs_mre_val3, name='ERR1')
		hdu3 = fits.ImageHDU(grs_mre_val2, name='VAL2')
		hdu4 = fits.ImageHDU(grs_mre_val4, name='ERR2')
		hdul = fits.HDUList([primary_hdu, hdu1, hdu2, hdu3, hdu4])

		hdul.writeto(out_dir + '/{}_retrieved.fits'.format(id_name[jj]), overwrite=True)


#=====================================================


print('End of script\n')
