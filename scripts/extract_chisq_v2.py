#=====================================================
# Inputs


spec_dir = 'stage2'                # Location of NEMESIS outputs
out_dir = 'results_stage2'         # Desired name of output directory

ref_file = 'jupiter_v2022miri.ref' # Ref file used for retrievals

show = True                        # Show results?
chisq_toll = 1000.0                # Max allowed chisq


#=====================================================
# Description
"""
Extracts chisq for the retrieval and saves as a fits file input to 
extract_multimre_v2.py

Usage

		[Alter inputs above]

		python extract_chisq_v2.py

"""
#=====================================================
# Imports


import numpy as np
import os
import glob
from astropy.io import fits
import matplotlib.pyplot as plt
import shutil
import sys
sys.path.insert(1, '/data/nemesis/jh852/python/subroutines')
from manifold import *


#=====================================================
# Main


if not os.path.exists(out_dir):
	os.mkdir(out_dir)

id_data = np.loadtxt(spec_dir + '/id_conserve.txt')
id = id_data[:,0]
lon = id_data[:,1]
lat = id_data[:,2]
id_length = len(id)

min_lon = np.nanmin(lon) - 0.5
max_lon = np.nanmax(lon) + 0.5
min_lat = np.nanmin(lat) - 0.5
max_lat = np.nanmax(lat) + 0.5

print('\nLon range = {aa} - {bb} degrees'.format(aa=max_lon, bb=min_lon))
print('Lat range = {aa} - {bb} degrees\n'.format(aa=min_lat, bb=max_lat))

lon_arr = np.arange(max_lon, min_lon-0.1, -0.5)
lat_arr = np.arange(min_lat, max_lat+0.1, 0.5)

lon_grid, lat_grid = np.meshgrid(lon_arr, lat_arr)

primary_hdu = fits.PrimaryHDU()
lon_hdu = fits.ImageHDU(lon_grid, name='LON_WEST')
lat_hdu = fits.ImageHDU(lat_grid, name='LAT_PGR')

hdul = fits.HDUList([primary_hdu, lon_hdu, lat_hdu])
hdul.writeto(out_dir + '/lon_lat.fits', overwrite=True)


string = 'chisq/ny is equal to :'

print('Extracting data\n')

dy, dx = lon_grid.shape

chi_array = np.empty((dy, dx))
chi_array[:] = np.nan

print('Building cube\n')

for kk in range(id_length):
	print(kk+1, lon[kk], lat[kk])

	try:
		spec = np.loadtxt(spec_dir + '/core{aa}/core_1/nemesis.mre'.format(aa=kk+1), skiprows=5, max_rows=1, usecols=(5))
		if np.isnan(spec) == True:
			raise TypeError('NaN values not permitted')

		with open(spec_dir + '/core{aa}/core_1/log_{aa}'.format(aa=kk+1)) as f:
			data = f.readlines()
		f.close()

		index = [idx for idx, s in enumerate(data) if 'chisq/ny is equal to :' in s][0]

		chi_array[np.where(lat_arr == lat[kk])[0][0], np.where(lon_arr == lon[kk])[0][0]] = float(data[index].replace('chisq/ny is equal to :', '').replace(' ', ''))

	except:
		print('WARNING\nChisq not found')


hdu_sav = fits.PrimaryHDU(chi_array)
hdu_sav.writeto(out_dir + '/chisq.fits',overwrite=True)

xtick_array, x_label_array, ytick_array, y_label_array = define_labels(lon_arr, lat_arr)

plt.imshow(chi_array,cmap='viridis',origin='lower')
plt.colorbar()

plt.xticks(xtick_array, x_label_array)
plt.yticks(ytick_array, y_label_array)

plt.savefig(out_dir + '/chisq.png')
if show == True:
	plt.show()
plt.clf()


#=====================================================

print('End of script\n')
