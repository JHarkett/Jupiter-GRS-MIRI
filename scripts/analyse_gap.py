#=======================================================================================
# Description
"""
Loads a cube (input below), quantifies the size of the gap (in terms of physical distance
and scaling factor). Saves the results as a cube.

Usage
		[Input parameters below]

		python analyse_gap.py


"""
#=======================================================================================
# Inputs


epoch = 'july'
tile = 'centre'

bands1 = ['ch2-long', 'ch3-short'] # 1nd comparison band
bands2 = ['ch3-short', 'ch3-medium'] # 2rd comparison band

wave_bound1 = [11.55, 13.34] # 1th wavelength
wave_bound2 = [11.70, 13.47] # 2st wavelength


#=======================================================================================
# Inputs


import numpy as np
from astropy.io import fits
import sys
import matplotlib.pyplot as plt

sys.path.insert(1, '/data/nemesis/jh852/python/subroutines')
from manifold import *


#=======================================================================================
# Functions


# Opens defined file and returns important parameters
def open_hdul(file):
	hdul = fits.open(file)
	spec_cube = hdul['SCI'].data
	
	hdul.close()

	return spec_cube


#=======================================================================================
# Main


len_b = len(bands1)
print('\n{} gaps detected\n'.format(len_b))

primary_hdu = fits.PrimaryHDU()

for kk in range(len_b):
	print('Analysing\n' + bands1[kk] + ' vs ' + bands2[kk] + '\n')

	# Loading the data
	spec_cube1 = open_hdul(epoch + '/' + tile + '/d_combined_0.5/Level3_' + bands1[kk] + '_s3d_nav.fits')
	spec_cube2 = open_hdul(epoch + '/' + tile + '/d_combined_0.5/Level3_' + bands2[kk] + '_s3d_nav.fits')

	wave1 = gen_wave_arr(epoch + '/' + tile + '/d_combined_0.5/Level3_' + bands1[kk] + '_s3d_nav.fits')
	wave2 = gen_wave_arr(epoch + '/' + tile + '/d_combined_0.5/Level3_' + bands2[kk] + '_s3d_nav.fits')

	wave_idx_low1= find_nearest(wave1, wave_bound1[kk])
	wave_idx_high1 = find_nearest(wave1, wave_bound2[kk])
	wave_idx_low2 = find_nearest(wave2, wave_bound1[kk])
	wave_idx_high2 = find_nearest(wave2, wave_bound2[kk])


	# Constraining to the key wavelengths

	spec_cube1 = spec_cube1[wave_idx_low1:wave_idx_high1]
	spec_cube2 = spec_cube1[wave_idx_low2:wave_idx_high2]

	spec_frame1 = np.nanmean(spec_cube1, axis=0)
	spec_frame2 = np.nanmean(spec_cube2, axis=0)

	subtractor = np.subtract(spec_frame2, spec_frame1)
	divider = np.divide(spec_frame2, spec_frame1)





#=======================================================================================

print('End of script\n')
