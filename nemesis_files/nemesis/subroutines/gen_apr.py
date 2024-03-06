import numpy as np
import sys


var1 = int(sys.argv[1])
var2 = int(sys.argv[2])
var3 = int(sys.argv[3])
retval = sys.argv[4]


if not os.path.exists(retval + '/core_1/nemesis.apr'):
	f = open(retval + '/core_1/nemesis.apr','a')
	f.write('*******Apriori File*******\n\n')
else:
	f = open(retval + '/core_1/nemesis.apr','a')	

if var1 == 0:
	name = 'temp'
	if var3 == 0:
		




print('End of script\n')
