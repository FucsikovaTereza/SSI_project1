import numpy as np
import numpy.matlib
np.seterr(divide='ignore', invalid='ignore') # ignore errors

import functions as f

###   run simulation  ###

for i in range(1):
    X, V, i_T1, i_T2 = f.run_simulation(f.xA)
    


