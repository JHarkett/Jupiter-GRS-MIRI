#=============================================================================
# Description
"""
Locates the wavelength positions where all spaxels are np.nan (complete
saturation)

Usage
		[Input parameters below]

		python analyse_nans.py


"""
#=============================================================================
# Inputs


tile_dirs = ['july/centre', 'july/east', 'aug/centre', 'aug/west']
dir_in = 'd_combined_0.5'
band = 'ch1-medium'


#=============================================================================
# Imports


import numpy as np
from astropy.io import fits
import sys

sys.path.insert(1, '/data/nemesis/jh852/python/subroutines')
from manifold import *


#=============================================================================
# Main


len_t = len(tile_dirs)

print('\nAnalysing {} tile data-sets\n'.format(len_t))

total = []

for kk in range(len_t):
	print(tile_dirs[kk])

	hdul = fits.open(tile_dirs[kk] + '/d_combined_0.5/Level3_' + band + '_s3d_nav.fits')
	spec_cube = hdul['SCI'].data
	hdul.close()

	if kk == 0:
		wave = gen_wave_arr(tile_dirs[kk] + '/d_combined_0.5/Level3_' + band + '_s3d_nav.fits')

	spec_arr_pre = np.nanmean(spec_cube, axis=2)
	spec_arr = np.nanmean(spec_arr_pre, axis=1)

	arral = np.where(np.isnan(spec_arr) == True)
	print(arral[0])

	for jj in range(len(arral[0])):
		if not arral[0][jj] in total:
			total.append(arral[0][jj])

	print('\n')

print(sorted(total))
print('\n')



#=============================================================================

print('End of script\n')
