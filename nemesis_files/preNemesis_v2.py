#================================================================
"""
Generates the required files for a Nemesis retrieval/forward model

Usage

		Set input parameters below

		python preNemesis.py

Updates

Version 2
	- Updated inputs for scattering runs
	- Removed second retrieval feature


"""
#================================================================
# Main inputs (will likely need setting every time)


planet = 5

ref_file = 'prior.ref' # Not needed if second_retrieval == True
prf_file = 'prior.prf'
ref_skip = 25 # For the ref file

xscfile='grs_double.xsc'
aerfile='../../aerosol_double.ref' # Not needed if second_retrieval == True
spxfile = 'jwst_mrs.spx'

phase_files = ['../../PHASE1.DAT', '../../PHASE2.DAT'] # fmt = PHASE[n].dat
hgphase_files = ['../../hgphase1.dat', '../../hgphase2.dat'] # fmt = hgphase[n].dat

fcloud_file = '../../fcloud.ref'

var=['-1 0 32', '-2 0 32', '11 0 20', '28 0 20', '01 0 20']

apr_override = False
var_override = ['']
vals_override = ['']

temp_override = False # True if retrieving temperature
temp_apr = 'tempapr.dat'

synth = 0       # Set to 1 to compute a forward model
nit=15          # Number of Nemesis iterations

walltime = '22:00:00'
vmem = '15gb'

gen_submitjob = False
run_retrieval = False


#================================================================
# Other useful inputs


obs = 'strat'
fmerror = 'file'
fmerror_data = 'fmerror_data.txt'

ciafile='dimers_ff_fb_bb_dnu2.0.tab'
dnu_cia=2.0 # Wavenumber step of CIA file
numpara=24  # Number of para-fractions in CIA file (12 or 24).

solfile = 'sun_spec_echo.dat'

ktab_path = '/data/nemesis/specdata/ktables/miri/mrs2024_combi/'
ktab_fmt = '[species].combi.kta'
species = ['nh3_all', 'ph3', 'h2', 'ash3', 'geh4', 'co', 'h2o', 'c2h6_all', 'ch3d', 'ch4', 'ch4_13C']


ispace = 1    # ISPACE=1 (wavelength), 0 (wavenumber)
iscat = 1     # ISCAT =0 (thermal) =1 (scattering), Also sets lower boundary in .set file to be Lambert
ilbl = 0    # Line by line calculations

sun = 1       # Sunlight on(1) or off(0)

inormal = 0   # Equilibrium ortho/para (0) or normal (3:1, 1)
iray = 1      # Rayleigh optical depths off (0), gas giant (1),CO2(2),N2-O2 (>2)
ih2o = 0      # Additional H2O continuum on (1) or off (0)
ich4 = 0      # Additional CH4 continuum on (1) or off (0)
io3 = 0       # Additional O3 continuum on (1) or off (0)
inh3 = 0
iptf = 0      # Only set to 1 if using high-T CH4 partition function
imie = 0      # Set to 0 to use hgphase.dat, 1 to compute from Mie theory.
iuv = 0

layht = -80.0000                     # Bottom layer for *set file, always >= base layer in PRF file.
laytyp = 1


#================================================================
# Imports


import numpy as np
import shutil
import glob
import os
import time
import sys

sys.path.insert(1, '/data/nemesis/jh852/python/subroutines')
from manifold import *


#================================================================
# Main


print('Creating working directory\n')

retval = os.getcwd().replace('/lustre/alice3','')
if synth == 1:
	nit = 0

if os.path.exists('core_1'):
	shutil.rmtree('core_1')

os.mkdir('core_1')


#---------- PHASE[n].DAT file ----------


for kk in range(len(phase_files)):
	shutil.copy(phase_files[kk], 'core_1/' + phase_files[kk].replace('..', '').replace('/', ''))


#---------- hgphase[n].dat file ----------


for kk in range(len(phase_files)):
	shutil.copy(hgphase_files[kk], 'core_1/' + hgphase_files[kk].replace('..', '').replace('/', ''))


#---------- nemesis.sol file ----------


f = open('core_1/nemesis.sol','a')
f.write(solfile)
f.close()


#---------- Aerosol reference file ----------


shutil.copy(aerfile,'core_1/aerosol.ref')


#---------- Log file ----------


dt = time.asctime(time.localtime(time.time()))

f = open('core_1/logfile','a')
f.write(obs + '       1\n' + retval + '/' + spxfile + '\n' + dt)
f.close()


#---------- CIA file ----------


f = open('core_1/nemesis.cia','a')
f.write(ciafile + '\n      {aa}\n      {bb}'.format(aa=dnu_cia, bb=numpara))
f.close()


#---------- kls file ----------


length_s = len(species)

f = open('core_1/nemesis.kls','a')

for kk in range(length_s):
	f.write(ktab_path + ktab_fmt.replace('[species]',species[kk]) + '\n')
f.close()


#---------- ref file ----------


shutil.copy(ref_file,'core_1/nemesis.ref')

ref_data = np.loadtxt('core_1/nemesis.ref',skiprows=ref_skip,usecols=(0,2))
height = ref_data[:,0]
temp = ref_data[:,1]
len_layers = len(temp)
base_temp = temp[find_nearest(height,layht)]

print('Reference file contains: {} pressures\n'.format(len_layers))


#---------- prf file ----------


shutil.copy(prf_file, 'core_1/nemesis.prf')


#---------- xsc file ----------


shutil.copy(xscfile,'core_1/nemesis.xsc')


#---------- fcloud file ----------


shutil.copy(fcloud_file, 'core_1/fcloud.ref')


#---------- abort file ----------


shutil.copy('/data/nemesis/jh852/python/nemesis/template_core/nemesis.abo','core_1/nemesis.abo')


#---------- flag file ----------


f = open('core_1/nemesis.fla','a')
f.write('       {aa}	! Inormal (0=eqm, 1=normal)\n       {bb}	! Iray (0=off, 1=on)\n       {cc}	! IH2O (1 = turn extra continuum on)\n       {dd}	! ICH4 (1 = turn extra continuum on)\n       {ee}	! IO3 (1 = turn extra continuum on)\n       {ff}	! INH3 (1 = turn extra continuum on)\n       {gg}	! Iptf (0=default, 1=CH4 High-T)\n       {hh}	! IMie\n       {ii}	! UV Cross-sections\n       0	! INLTE (0=LTE)'.format(aa=inormal, bb=iray, cc=ih2o, dd=ich4, ee=io3, ff=inh3, gg=iptf, hh=imie, ii=iuv))
f.close()


#---------- name file ----------


print("Name of the retrieval is currently by default set to be 'nemesis'\n")

shutil.copy('/data/nemesis/jh852/python/nemesis/template_core/nemesis.nam','core_1/nemesis.nam')


#---------- setup file ----------


print('Generating set file')

if planet == 5:
	dist=5.2
if planet == 6:
	dist=9.546
if planet == 11:
	dist=9.546 # Titan
if planet == 7:
	dist=19.2
if planet == 8:
	dist=30.0
if planet == 13:
	dist=0.045

print('Distance from sun = {}'.format(dist))

f = open('core_1/nemesis.set','a')
f.write('*********************************************************\n Number of zenith angles :  5\n  0.165278957666387	  0.327539761183898\n  0.477924949810444	  0.292042683679684\n  0.738773865105505	  0.224889342063117\n  0.919533908166459	  0.133305990851069\n   1.00000000000000	  2.222222222222220E-002\n')
f.write(' Number of fourier components :  1\n Number of azimuth angles for fourier analysis : 100\n')
f.write(' Sunlight on(1) or off(0) :  {aaa}\n Distance from Sun (AU) :   {bbb}\n Lower boundary cond. Thermal(0) Lambert(1) :  {ccc}\n Ground albedo :   0.000\n Surface temperature :  {ddd}\n'.format(aaa=sun, bbb=dist, ccc=iscat, ddd=base_temp))
f.write('*********************************************************\n Alt. at base of bot.layer (not limb) :   {aaaa:.4f}\n Number of atm layers : 39\n Layer type :  {cccc}\n Layer integration :  1\n*********************************************************'.format(aaaa=layht, cccc=laytyp))
f.close()


#---------- parah2 ref file ----------


shutil.copy('/data/nemesis/jh852/python/nemesis/template_core/parah2.ref','core_1/parah2.ref')


#---------- spx file ----------


shutil.copy(spxfile,'core_1/nemesis.spx')
with open('core_1/nemesis.spx') as f:
	data_s = f.readlines()

fwhm,lat,lon,ngeom = data_s[0].split()
ngeom = int(ngeom)

tally = 1

nconv = [None]*ngeom
wave_spx = [None]*ngeom

for kk in range(ngeom):

	nconv[kk] = int(data_s[tally].split()[0])

	wave_spx[kk] = np.loadtxt('core_1/nemesis.spx',skiprows=tally+3,usecols=(0), max_rows=nconv[kk])

	tally = tally + nconv[kk] + 3

wave_comp = np.concatenate((wave_spx))
length_wc = len(wave_comp)


#---------- fm error file ----------


if fmerror == 'file':
	fm_in = np.loadtxt(fmerror_data)
	wave_fm = fm_in[:,0]
	fmerror_vals = fm_in[:,1]

	if not len(wave_fm) == len(wave_comp):
		print('ERROR\nfmerror file wavelengths are different to spx file wavelengths\n')
		wave_fm = wave_comp
		fmerror_vals = [0.000000]*length_wc

else:
	wave_fm = wave_comp
	fmerror_vals = [fmerror]*length_wc

cccp = [None]*2
cccp[0] = wave_fm
cccp[1] = fmerror_vals

pccc = np.transpose(cccp)

f = open('core_1/fmerror.dat','a')
f.write('{}\n'.format(length_wc))

for jk in range(length_wc):
	f.write('{qq:6e}	{ww:.6e}\n'.format(qq=pccc[jk][0], ww=pccc[jk][1]))

f.close()


#---------- apr file ----------


print('Generating apr file')

length_v = len(var)
id_array = [None]*length_v

print('{} variables being retrieved\n'.format(length_v))

print('Using default values\n')

if apr_override == True:
	print('Overriding the default values for the following:')

	len_ov = len(var_override)

	for ae in range(len_ov):
		print(var_override[ae])

	print('\n')

g = open('core_1/nemesis.apr','a')
g.write('*******Apriori File*******\n           {}\n'.format(length_v))

id_array2 = [['2', '2', '2']]


for kk in range(length_v):
	print(var[kk])
	id_array[kk] = var[kk].split()

	with open('/data/nemesis/jh852/python/nemesis/subroutines/apr_data.txt') as f:
		apr_data = f.readlines()
	length_aprdat = len(apr_data)
	f.close()

	for jji in range(length_aprdat):
		if apr_data[jji][0:1] == id_array[kk][0] or apr_data[jji][0:2] == id_array[kk][0]:
			if apr_data[jji][2] == id_array[kk][1] or apr_data[jji][3] == id_array[kk][1]:
				species, iso, name, name2, p_knee, err_p_knee, p_zero, err_p_zero, const_vmr, err_const_vmr, fsh, err_fsh, err_pro0, err_pro3, pro_fac3, op, err_op, width, dwidth = apr_data[jji].split()

	g.write(var[kk].replace('01', '1') + ' - ' + name2.replace('_',' ') + '\n')

	n = -1

	for i in range(len(id_array2)):
		if id_array2[i] == id_array[kk]:
			n = v_start2[i]

	if id_array[kk] in id_array2:
		n = v_start2[id_array2.index(id_array[kk])]

	if n != -1:
		print('Retrieving results of previous retrieval')

		pre_data = np.loadtxt(prior_dir + '/core_1/nemesis.mre', skiprows=n, max_rows=120, usecols=(4, 5))

		pressure = np.loadtxt('core_1/nemesis.ref',skiprows=ref_skip,usecols=(1))
		prof = pre_data[:,0]
		delta = pre_data[:,1]

		g.write(name + 'apr.dat\n')

		a = [None]*3
		a[0] = pressure
		a[1] = prof
		a[2] = delta
		b = np.transpose(a)

		np.savetxt('core_1/' + name + 'apr.dat', b, fmt='%.5e', delimiter='	', comments='', header='120    1.50000')


	if temp_override == False and n == -1:
		if id_array[kk][2] == '0':
			g.write(name + 'apr.dat\n')

			with open('core_1/nemesis.ref') as h:
				ref_data = h.readlines()
			length_h = len(ref_data)
			h.close()

			for iij in range(length_h):
				if ref_data[iij][2:4] == id_array[kk][0] and ref_data[iij][8:9] == id_array[kk][1]:
					colu = iij

			if species == '0':
				colu = 2

			pre_data = np.loadtxt('core_1/nemesis.ref',skiprows=ref_skip,usecols=(1,colu))
			pressure = pre_data[:,0]
			prof = pre_data[:,1]

			length_red = len(pressure)

			if species == '0':
				delta = [4.00000]*length_red

			else:
				delta = prof*float(err_pro0)

			a = [None]*3
			a[0] = pressure
			a[1] = prof
			a[2] = delta
			b = np.transpose(a)

			np.savetxt('core_1/' + name + 'apr.dat',b,fmt='%.5e', delimiter='	', comments='', header='120	1.50000')


	if temp_override == True and id_array[kk][2] == '0':
		g.write(name + 'apr.dat\n')
		shutil.copy(temp_apr,'core_1/tempapr.dat')

	if apr_override == True and var[kk] in var_override:
		idx_ov = var_override.index(var[kk])

		g.write(vals_override[idx_ov])

	else:
		if id_array[kk][2] == '3':
			g.write(pro_fac3 + ' ' + err_pro3 + '\n')

		if id_array[kk][2] == '4':
			g.write(p_knee + ' ' + err_p_knee + '\n' + const_vmr + ' ' + err_const_vmr + '\n' + fsh + ' ' + err_fsh + '\n')

		if id_array[kk][2] == '20':
			g.write(p_knee + ' ' + p_zero + '\n' + const_vmr + ' ' + err_const_vmr + '\n' + fsh + ' ' + err_fsh + '\n')

		if id_array[kk][2] == '32':
			g.write(p_knee + ' ' + err_p_knee + '\n' + op + ' ' + err_op + '\n' + fsh + ' ' + err_fsh + '\n')

		if id_array[kk][2] == '47':
			g.write(op + ' ' + err_op + '\n' + p_knee + ' ' + err_p_knee + '\n' + width + ' ' + dwidth + '\n')

		if id_array[kk][2] == '48':
			g.write(p_knee + ' ' + err_p_knee + '\n' + p_zero + ' ' + err_p_zero + '\n' + op + ' ' + err_op + '\n' + fsh + ' ' + err_fsh + '\n')

g.close()


print('\nSuccessfully generated apr file\n')


#---------- input file ----------


subs = 0

f = open('core_1/nemesis.inp','a')
f.write('       {aap}       {abp}       {acp}              ! ispace,iscat\n0.0               ! Wavenumber offset to add to measured spectrum\nfmerror.dat\n      {adp}  ! Number of iterations\n'.format(aap=ispace, abp=iscat, acp=ilbl, adp=nit))
f.write('0.1               ! Minimum % change in phi before terminating\n1 1   ! Total spectra to fit and starting ID\n       {}  ! Dont substitute previous retrievals\n       0   ! Output format'.format(subs))
f.close()


#---------- qsub file ----------


if gen_submitjob == True:
	if os.path.exists('submitjob'):
		os.remove('submitjob')

	f = open('submitjob','a')
	f.write('#!/bin/bash\n')
	f.write('#\n')
	f.write('#SBATCH --job-name=nemesis\n')
	f.write('#SBATCH --account=nemesis\n')
	f.write('#SBATCH --time=' + walltime + '\n')
	f.write('#SBATCH --mem=' + vmem + '\n')
	f.write('#SBATCH --export=NONE\n')
	f.write('#SBATCH --cpus-per-task=1\n')
	f.write('#SBATCH --array=1-1\n')
	f.write('export PATH=~/bin/ifort:$PATH\n')
	f.write('cd $SLURM_SUBMIT_DIR\n')
	f.write('inputdir=core_${SLURM_ARRAY_TASK_ID}\n')
	f.write('cd $inputdir\n')
	f.write('Nemesis < nemesis.nam > log_${SLURM_ARRAY_TASK_ID}\n')
	f.close()


#---------- running retrieval (optional) ----------


if run_retrieval == True:
	print('Starting retrieval')

	os.system('sbatch submitjob')


#================================================================


print('End of script\n')
