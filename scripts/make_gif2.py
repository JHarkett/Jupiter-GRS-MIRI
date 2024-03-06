#==================================================================
"""
Generates a moving gif of the July and August temperature data

Usage

		Input parameters below

		python make_gif.py


"""
#==================================================================   
# Inputs


assign_axes = False
lon_range = [319,275]
lat_range = [-30,-12]

# If it is False: fill out below
# Range to extend the displayed longitude around the GRS
lon_bounder = 20.0


#==================================================================   
# Imports


import numpy as np
from astropy.io import fits
import glob
import shutil
from PIL import Image
import matplotlib.pyplot as plt
import os
import imageio
from matplotlib.patches import Ellipse
import sys




#==================================================================   
# Main




out_name = '000_animation_stage2.gif'

hdul = fits.open('july/mosaics_stage2/000_retrieved_cleaned.fits')
spec_data = hdul['VARIABLE'].data
lon = hdul['LON_WEST'].data
lat = hdul['LAT_PGR'].data
hdul.close()
dz, dy, dx = spec_data.shape

hdul2 = fits.open('aug/mosaics_stage2/000_retrieved_cleaned.fits')
spec_data2 = hdul2['VARIABLE'].data
lon2 = hdul2['LON_WEST'].data
lat2 = hdul2['LAT_PGR'].data
hdul2.close()
dz2, dy2, dx2 = spec_data2.shape


pressure = np.loadtxt('nemesis.prf',skiprows=24,usecols=(1))
pressure *= 1013.25

spec_data = spec_data[::-1]
spec_data2 = spec_data2[::-1]
pressure = pressure[::-1]

idx_top = find_nearest(pressure,0.01)
idx_bottom = find_nearest(pressure,1000)


spec_data_trim = spec_data[idx_top:idx_bottom+1]
spec_data_trim2 = spec_data2[idx_top:idx_bottom+1]
pressure_trim = pressure[idx_top:idx_bottom+1]

idx_slow = find_nearest(pressure_trim, 100)

ny, nx = lat.shape

lon_arr = [None]*dx
lat_arr = [None]*dy

for i in range(dx):
	lon_arr[i] = lon[0][i]
for j in range(dy):
	lat_arr[j] = lat[j][0]

lon_arr2 = [None]*dx2
lat_arr2 = [None]*dy2

for i in range(dx2):
	lon_arr2[i] = lon2[0][i]
for j in range(dy2):
	lat_arr2[j] = lat2[j][0]

x_label_array = np.arange((round(np.amax(lon_arr)/5)*5)+5,(round(np.amin(lon_arr)/5)*5)-5,-5)
lengthx = len(x_label_array)
xtick_array = [None]*lengthx
for i in range(lengthx):
	xtick_array[i] = find_nearest(lon_arr,x_label_array[i])

y_label_array = np.arange((round(np.amin(lat_arr)/5)*5)-5,(round(np.amax(lat_arr)/5)*5)+5,5)
lengthy = len(y_label_array)
ytick_array = [None]*lengthy
for i in range(lengthy):
	ytick_array[i] = find_nearest(lat_arr,y_label_array[i])


x_label_array2 = np.arange((round(np.amax(lon_arr2)/5)*5)+5, (round(np.amin(lon_arr2)/5)*5)-5,-5)
lengthx2 = len(x_label_array2)
xtick_array2 = [None]*lengthx2
for i in range(lengthx2):
	xtick_array2[i] = find_nearest(lon_arr2,x_label_array2[i])

y_label_array2 = np.arange((round(np.amin(lat_arr2)/5)*5)-5,(round(np.amax(lat_arr2)/5)*5)+5,5)
lengthy2 = len(y_label_array2)
ytick_array2 = [None]*lengthy2
for i in range(lengthy2):
	ytick_array2[i] = find_nearest(lat_arr2,y_label_array2[i])


temp = [None]*dz
temp2 = [None]*dz
spec_crop = [None]*dz

for kk in range(120):
	temp[kk] = np.nanmean(spec_data[kk])
	temp2[kk] = np.nanmean(spec_data2[kk])

temp_trim = temp[idx_top:idx_bottom+1]
temp_trim2 = temp[idx_top:idx_bottom+1]

filenames = []

for k in range(len(pressure_trim)):
	fig = plt.subplots(figsize=(13, 10))

	plt.subplot2grid((2, 3), (0, 0), colspan=1, rowspan=1)
	plt.plot(temp, pressure,'k-',linewidth=0.8)
	plt.xlabel('Mean Temperature (K)',fontsize=8)
	plt.ylabel('Pressure (mbar)',fontsize=8)
	plt.grid()
	plt.hlines(pressure_trim[k],xmin=np.amin(temp),xmax=np.amax(temp),color='red')
	plt.ylim(max(pressure),min(pressure))
	plt.yscale('log')

	min_arr = [np.nanmin(spec_data_trim[k]), np.nanmin(spec_data_trim2[k])]
	max_arr = [np.nanmax(spec_data_trim[k]), np.nanmax(spec_data_trim2[k])]

	plt.subplot2grid((2, 3), (0, 1), colspan=2, rowspan=1)

	plt.imshow(spec_data_trim[k],cmap='inferno',origin='lower', vmin=min(min_arr), vmax=max(max_arr))
	plt.gca().add_patch(Ellipse((find_nearest(lon_arr, 295.639), find_nearest(lat_arr, -20.7)), width=(11.09979*2), height=(8.880474*2), fill=None, linestyle='--', edgecolor='red', linewidth=2.0))

	plt.xlabel('West Longitude ($^\circ$)',fontsize=8)
	plt.ylabel('Planetographic Latitude ($^\circ$)',fontsize=8)
	plt.xticks(xtick_array,x_label_array)
	plt.yticks(ytick_array,y_label_array)
	plt.xlim(find_nearest(lon_arr, 295.639 + lon_bounder),find_nearest(lon_arr, 295.639 - lon_bounder))
	plt.ylim(find_nearest(lat_arr,lat_range[0]),find_nearest(lat_arr,lat_range[1]))
	plt.title('July\nMean Temperature: {bb} K'.format(bb=round(temp_trim[k],2)))

	plt.subplot2grid((2, 3), (1, 1), colspan=2, rowspan=1)

	plt.imshow(spec_data_trim2[k],cmap='inferno',origin='lower', vmin=min(min_arr), vmax=max(max_arr))
	plt.gca().add_patch(Ellipse((find_nearest(lon_arr, 301.0), find_nearest(lat_arr, -20.7)), width=50, height=45, fill=None, linestyle='--', edgecolor='red', linewidth=2.0))

	plt.xlabel('West Longitude ($^\circ$)',fontsize=8)
	plt.ylabel('Planetographic Latitude ($^\circ$)', fontsize=8)
	plt.xticks(xtick_array2,x_label_array2)
	plt.yticks(ytick_array2,y_label_array2)
	plt.xlim(find_nearest(lon_arr2, 301.0 + lon_bounder),find_nearest(lon_arr2, 301.0 - lon_bounder))
	plt.ylim(find_nearest(lat_arr2, lat_range[0]), find_nearest(lat_arr2, lat_range[1]))
	plt.title('August\nMean Temperature: {bb} K'.format(bb=round(temp_trim2[k],2)))

	plt.subplots_adjust(wspace=0.3)
	plt.suptitle('Pressure: {aa} mbar'.format(aa=round(pressure_trim[k],2)), weight='bold')

	filename = 'figures/{}.png'.format(k)
	filenames.append(filename)
	plt.savefig(filename)
	plt.close()


with imageio.get_writer('figures/' + out_name, mode='I') as writer:
	for filename in filenames:
		image = imageio.imread(filename)
		name1 = filename.replace('.png','')
		number = int(name1.replace('figures/',''))
		writer.append_data(image)
		if number >= idx_slow:
			writer.append_data(image)

for filename in set(filenames):
	os.remove(filename)


print('End of script\n')

