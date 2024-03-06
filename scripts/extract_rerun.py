#=========================================================================
# Description
"""
Detects spaxels that have been re-run and inserts into the final data cube

Usage

		[Input parameters below]

		python extract_rerun.py


"""
#=========================================================================
# Inputs


ret_dir = 'stage2'
results_dir = 'results_stage2'


#=========================================================================
# Imports


import numpy as np
from astropy.io import fits
import sys
import os
import shutil
from subprocess import call

sys.path.insert(1, '/data/nemesis/jh852/python/subroutines')
from manifold import *


#=========================================================================
# Main


print('\n')
nwave, nvar, id_name, id_array, v_start, v_length = detect_variables(ret_dir + '/rerun/core1/core_1/nemesis.mre')

hdul = fits.open(results_dir + '/lon_lat.fits')
lon_grid = hdul['LON_WEST'].data
lat_grid = hdul['LAT_PGR'].data
hdul.close()

lon_arr, lat_arr = gen_lonlat_arr(lon_grid, lat_grid)

id_data = np.loadtxt(ret_dir + '/rerun/id_conserve_aux.txt', skiprows=1)

if len(id_data.shape) != 1:
	id_pre = id_data[:,0]
	id_new = id_data[:,1]
	lon_id = id_data[:,2]
	lat_id = id_data[:,3]

else:
	id_pre = [id_data[0]]
	id_new = [id_data[1]]
	lon_id = [id_data[2]]
	lat_id = [id_data[3]]

len_i = len(id_pre)

print('Extracting spectral fits\n')

hdul = fits.open(results_dir + '/residual.fits')
residual = hdul['RESIDUAL'].data
wave_err = hdul['ERR'].data
spec_data = hdul['SPECTRAL_DATA'].data
fit_data = hdul['FIT_DATA'].data
wave_res = hdul['WAVELENGTH'].data

print('\n')

for jj in range(len_i):
	print('core{}'.format(jj+1))

	try:
		shutil.copy(ret_dir + '/rerun/core{}/core_1/nemesis.mre'.format(int(id_new[jj])), ret_dir + '/core{}/core_1/nemesis.mre'.format(int(id_pre[jj])))
		shutil.copy(ret_dir + '/rerun/core{}/core_1/nemesis.prf'.format(int(id_new[jj])), ret_dir + '/core{}/core_1/nemesis.prf'.format(int(id_pre[jj])))
		shutil.copy(ret_dir + '/rerun/core{}/core_1/aerosol.prf'.format(int(id_new[jj])), ret_dir + '/core{}/core_1/aerosol.prf'.format(int(id_pre[jj])))

		res_data = np.loadtxt(ret_dir + '/rerun/core{}/core_1/nemesis.mre'.format(jj+1), skiprows=5, max_rows=nwave, usecols=(1, 2, 3, 5))

		x_pos = np.where(lon_arr == lon_id[jj])[0][0]
		y_pos = np.where(lat_arr == lat_id[jj])[0][0]

		wave_res[:, y_pos, x_pos] = res_data[:,0]
		residual[:, y_pos, x_pos] = res_data[:,3] - res_data[:,1]
		wave_err[:, y_pos, x_pos] = res_data[:,2]
		spec_data[:, y_pos, x_pos] = res_data[:,1]
		fit_data[:, y_pos, x_pos] = res_data[:,3]

	except:
		print('No data available\n')

hdul['RESIDUAL'].data = residual
hdul['ERR'].data = wave_err
hdul['SPECTRAL_DATA'].data = spec_data
hdul['FIT_DATA'].data = fit_data
hdul['WAVELENGTH'].data = wave_res

hdul.writeto(results_dir + '/residual.fits', overwrite=True)

print('Extracting species data\n')

for kk in range(nvar):
	print(id_name[kk])

	hdul = fits.open(results_dir + '/' + id_name[kk] + '_retrieved.fits')

	if id_array[kk][0] == '-1':
		val1 = hdul['VAL1'].data
		err1 = hdul['ERR1'].data
		val2 = hdul['VAL2'].data
		err2 = hdul['ERR2'].data
		val3 = hdul['VAL3'].data
		err3 = hdul['ERR3'].data

	if id_array[kk][0] == '0':
		val1 = hdul['VARIABLE'].data
		err1 = hdul['ERR'].data

	if id_array[kk][2] == '20':
		val1 = hdul['VAL1'].data
		err1 = hdul['ERR1'].data
		val2 = hdul['VAL2'].data
		err2 = hdul['ERR2'].data

	for jj in range(len_i):
		#print('({aa}, {bb})'.format(aa=lon_id[jj], bb=lat_id[jj]))
		print('core{}'.format(jj+1))

		try:

			mro_data = np.loadtxt(ret_dir + '/rerun/core{}/core_1/nemesis.mre'.format(jj+1), skiprows=v_start[kk], max_rows=v_length[kk], usecols=(4, 5))

			x_pos = np.where(lon_arr == lon_id[jj])[0][0]
			y_pos = np.where(lat_arr == lat_id[jj])[0][0]

			if id_array[kk][0] == '-1':
				val1[y_pos, x_pos] = mro_data[:,0][2]
				err1[y_pos, x_pos] = mro_data[:,1][2]
				val2[y_pos, x_pos] = mro_data[:,0][0]
				err2[y_pos, x_pos] = mro_data[:,1][0]
				val3[y_pos, x_pos] = mro_data[:,0][1]
				err3[y_pos, x_pos] = mro_data[:,1][1]

			if id_array[kk][0] == '0':
				val1[:, y_pos, x_pos] = mro_data[:,0]
				err1[:, y_pos, x_pos] = mro_data[:,1]

			if id_array[kk][2] == '20':
				val1[y_pos, x_pos] = mro_data[:,0][0]
				err1[y_pos, x_pos] = mro_data[:,1][0]
				val2[y_pos, x_pos] = mro_data[:,0][1]
				err2[y_pos, x_pos] = mro_data[:,1][1]

		except:
			print('No data available\n')

	if id_array[kk][0] == '-1':
		hdul['VAL1'].data = val1
		hdul['ERR1'].data = err1
		hdul['VAL2'].data = val2
		hdul['ERR2'].data = err2
		hdul['VAL3'].data = val3
		hdul['ERR3'].data = err3

	if id_array[kk][0] == '0':
		hdul['VARIABLE'].data = val1
		hdul['ERR'].data = err1

	if id_array[kk][2] == '20':
		hdul['VAL1'].data = val1
		hdul['ERR1'].data = err1
		hdul['VAL2'].data = val2
		hdul['ERR2'].data = err2

	hdul.writeto(results_dir + '/' + id_name[kk] + '_retrieved.fits', overwrite=True)


#=========================================================================

print('End of script\n')
