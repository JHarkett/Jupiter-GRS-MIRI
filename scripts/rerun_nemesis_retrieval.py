#======================================================================
# Description
"""
Detects retrievals that have not run properly and re-runs them

Usage
		[Input parameters below]

		python rerun_nemesis_retrieval.py


"""
#======================================================================
# Inputs


ret_dir = 'stage2'
results_dir = 'results_stage2'

ref_file = 'prior.ref'
xsc_file = 'grs_single.xsc'
aero_file = 'aerosol.ref'

walltime = '20:00:00'
vmem = '16G'


#======================================================================
# Imports


import numpy as np
from astropy.io import fits
import sys
import os
import shutil
from subprocess import call

sys.path.insert(1, '/data/nemesis/jh852/python/subroutines')
from manifold import *


#======================================================================
# Main


print('\nDetecting missed spaxels\n')

hdul = fits.open(results_dir + '/000_retrieved.fits')
temp_cube = hdul['VARIABLE'].data
lon_grid = hdul['LON_WEST'].data
lat_grid = hdul['LAT_PGR'].data
hdul.close()

temp = temp_cube[0]

lon_arr, lat_arr = gen_lonlat_arr(lon_grid, lat_grid)
xtick_array, x_label_array, ytick_array, y_label_array = define_labels(lon_arr, lat_arr)

id_data = np.loadtxt(ret_dir + '/id_conserve.txt', skiprows=1)
id = id_data[:,0]
lon_list = id_data[:,1]
lat_list = id_data[:,2]


arral = np.where(np.isnan(temp) == True)

len_n = len(arral[0])

x_pos = []
y_pos = []
core_pre = []
lon_pos = []
lat_pos = []

for kk in range(len_n):
	arral2 = np.where(np.logical_and(lon_list == lon_grid[arral[0][kk], arral[1][kk]], lat_list == lat_grid[arral[0][kk], arral[1][kk]]))

	if len(arral2[0]) == 1:
		x_pos.append(arral[1][kk])
		y_pos.append(arral[0][kk])

		core_pre.append(int(id[arral2[0][0]]))

		lon_pos.append(lon_grid[arral[0][kk], arral[1][kk]])
		lat_pos.append(lat_grid[arral[0][kk], arral[1][kk]])

if len(x_pos) == 0:
	print('No missed spaxels detected\n')

else:
	print('Missed spaxels detected\n')

	# The example plot
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 7))

	ax.imshow(temp, cmap='inferno', origin='lower')
	ax.scatter(x_pos, y_pos, s=10, c='green')

	ax.set_xlabel('West Longitude ($\degree$)')
	ax.set_ylabel('PGR Latitude ($\degree$)')
	ax.set_xticks(xtick_array, x_label_array)
	ax.set_yticks(ytick_array, y_label_array)

	plt.show()


	print('Generating files\n')

	if os.path.exists(ret_dir + '/rerun'):
		shutil.rmtree(ret_dir + '/rerun')
	os.mkdir(ret_dir + '/rerun')

	len_p = len(x_pos)

	f = open(ret_dir + '/rerun/id_conserve_aux.txt', 'a')

	f.write("#ID_pre	ID_new	lon	lat\n")

	for kk in range(len_p):
		print('{aa}/{bb}'.format(aa=kk+1, bb=len_p))

		os.mkdir(ret_dir + '/rerun/core{}'.format(kk+1))

		shutil.copy(ret_dir + '/core{}/aerosol.ref'.format(core_pre[kk]), ret_dir + '/rerun/core{}/aerosol.ref'.format(kk+1))
		shutil.copy(ret_dir + '/core{}/fmerror_data.txt'.format(core_pre[kk]), ret_dir + '/rerun/core{}/fmerror_data.txt'.format(kk+1))
		shutil.copy(ret_dir + '/core{}/grs_single.xsc'.format(core_pre[kk]), ret_dir + '/rerun/core{}/grs_single.xsc'.format(kk+1))
		shutil.copy(ret_dir + '/core{}/'.format(core_pre[kk]) + ref_file, ret_dir + '/rerun/core{}/'.format(kk+1) + ref_file)
		shutil.copy(ret_dir + '/core{}/'.format(core_pre[kk]) + ref_file.replace('.ref', '.prf'), ret_dir + '/rerun/core{}/'.format(kk+1) + ref_file.replace('.ref', '.prf'))
		shutil.copy(ret_dir + '/core{}/jwst_mrs.spx'.format(core_pre[kk]), ret_dir + '/rerun/core{}/jwst_mrs.spx'.format(kk+1))
		shutil.copy(ret_dir + '/core{}/preNemesis.py'.format(core_pre[kk]), ret_dir + '/rerun/core{}/preNemesis.py'.format(kk+1))
		shutil.copy(ret_dir + '/core{}/tempapr.dat'.format(core_pre[kk]), ret_dir + '/rerun/core{}/tempapr.dat'.format(kk+1))

		os.chdir(ret_dir + '/rerun/core{}'.format(kk+1))
		call('python preNemesis.py',shell=True)
		os.chdir('../../..')

		f.write("{aa}	{bb}	{cc}	{dd}\n".format(aa=core_pre[kk], bb=kk+1, cc=lon_pos[kk], dd=lat_pos[kk]))

	f.close()

	print('Generating qsub file\nwalltime = {pp}\nvmem = {pn}\n'.format(pp=walltime, pn=vmem))

	f = open(ret_dir + '/rerun/submitjob','a')
	f.write('#!/bin/bash\n')
	f.write('#\n')
	f.write('#SBATCH --job-name=nemesis\n')
	f.write('#SBATCH --account=nemesis\n')
	f.write('#SBATCH --time=' + walltime + '\n')
	f.write('#SBATCH --mem=' + vmem + '\n')
	f.write('#SBATCH --export=NONE\n')
	f.write('#SBATCH --cpus-per-task=1\n')
	f.write('#SBATCH --array=1-{}\n'.format(len_p))
	f.write('export PATH=~/bin/ifort:$PATH\n')
	f.write('cd $SLURM_SUBMIT_DIR\n')
	f.write('inputdir=core${SLURM_ARRAY_TASK_ID}/core_1\n')
	f.write('cd $inputdir\n')
	f.write('Nemesis < nemesis.nam > log_${SLURM_ARRAY_TASK_ID}\n')
	f.close()


#======================================================================

print('End of script\n')
