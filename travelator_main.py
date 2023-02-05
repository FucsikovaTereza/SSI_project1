import numpy as np
import numpy.matlib
np.seterr(divide='ignore', invalid='ignore') # ignore errors

from setting_variables import N, P
from functions import run_simulation
from plot_functions import hist


###   run simulation  ###

i_T1 = []
i_T2 = []
for i in range(5): 
    X, V = run_simulation(i_T1, i_T2, print_stats = False, plot_env = False)

hist(i_T1, i_T2)
    


