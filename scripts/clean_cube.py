#===============================================================
# Description
"""
Takes the average of a cube of pixels up to 2 pixel lengths,
returns a (hopefully) cleaner cube

Usage:

		python clean_cube.py [file_name]

	where [file_name] is the name of the file to be cleaned


"""
#===============================================================
# Imports


import numpy as np
from astropy.io import fits
import sys


#===============================================================
# Main


file = sys.argv[1]

if not '000' in file:
	print('\nmre file detected\n')

	hdul = fits.open(file)
	val1 = hdul['VAL1'].data
	err1 = hdul['ERR1'].data
	val2 = hdul['VAL2'].data
	err2 = hdul['ERR2'].data

	dy, dx = val1[0].shape

	print('\nShape of cubes = {aa} X {bb} pixels\n'.format(aa=dx, bb=dy))

	val1_clean = np.empty((dy, dx))
	err1_clean = np.empty((dy, dx))
	val2_clean = np.empty((dy, dx))
	err2_clean = np.empty((dy, dx))

	val1_clean[:] = np.nan
	err1_clean[:] = np.nan
	val2_clean[:] = np.nan
	err2_clean[:] = np.nan

	if '-' in file:
		print('\nAerosol layer detected\n')

		val3 = hdul['VAL3'].data
		err3 = hdul['ERR3'].data

		val3_clean = np.empty((dy, dx))
		err3_clean = np.empty((dy, dx))

		val3_clean[:] = np.nan
		err3_clean[:] = np.nan

	for j in range(2, dy-2):
		for i in range(2, dx-2):
			val1_clean[j, i] = np.nanmedian(val1[0, j-1:j+2, i-1:i+2])
			err1_clean[j, i] = np.nanmedian(err1[0, j-1:j+2, i-1:i+2])
			val2_clean[j, i] = np.nanmedian(val2[0, j-1:j+2, i-1:i+2])
			err2_clean[j, i] = np.nanmedian(err2[0, j-1:j+2, i-1:i+2])

			if '-' in file:
				val3_clean[j, i] = np.nanmedian(val3[0, j-1:j+2, i-1:i+2])
				err3_clean[j, i] = np.nanmedian(err3[0, j-1:j+2, i-1:i+2])

	val1_clean[np.where(np.isnan(val1[0]) == True)] = np.nan
	err1_clean[np.where(np.isnan(err1[0]) == True)] = np.nan
	val2_clean[np.where(np.isnan(val2[0]) == True)] = np.nan 
	err2_clean[np.where(np.isnan(err2[0]) == True)] = np.nan

	if '-' in file:
		val3_clean[np.where(np.isnan(val3[0]) == True)] = np.nan 
		err3_clean[np.where(np.isnan(err3[0]) == True)] = np.nan

	hdul['VAL1'].data = val1_clean
	hdul['ERR1'].data = err1_clean
	hdul['VAL2'].data = val2_clean
	hdul['ERR2'].data = err2_clean

	if '-' in file:
		hdul['VAL3'].data = val3_clean
		hdul['ERR3'].data = err3_clean

	hdul.writeto(file.replace('.fits', '_cleaned.fits'), overwrite=True)


else:
	print('\nUsing normal method\n')

	hdul = fits.open(file)
	var = hdul['VARIABLE'].data
	err = hdul['ERR'].data

	dz, dy, dx = var.shape
	print('\nShape of cube = {aa} X {bb} X {cc} pixels\n'.format(aa=dx, bb=dy, cc=dz))

	var_clean = np.empty((dz, dy, dx))
	err_clean = np.empty((dz, dy, dx))

	var_clean[:] = np.nan
	err_clean[:] = np.nan

	for k in range(dz):
		print(k)
		for j in range(2, dy-2):
			for i in range(2, dx-2):
				var_clean[k, j , i] = np.nanmedian(var[k, j-1:j+2, i-1:i+2])
				err_clean[k, j , i] = np.nanmedian(err[k, j-1:j+2, i-1:i+2])


	var_clean[np.where(np.isnan(var) == True)] = np.nan
	err_clean[np.where(np.isnan(err) == True)] = np.nan

	hdul['VARIABLE'].data = var_clean
	hdul['ERR'].data = err_clean

	hdul.writeto(file.replace('.fits', '_cleaned.fits'), overwrite=True)



#===============================================================


print('End of script')
