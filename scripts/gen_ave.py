import numpy as np
from astropy.io import fits



hdul = fits.open('results_stage1/-1032_retrieved.fits')
alt = hdul['VAL1'].data
op = hdul['VAL2'].data
fsh = hdul['VAL3'].data
hdul.close()


print('\nMedian altitude = {}'.format(np.nanmedian(alt)))
print('Median fsh = {}\n'.format(np.nanmedian(fsh)))


print('End of script\n')
