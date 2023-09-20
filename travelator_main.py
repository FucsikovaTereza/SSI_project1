import numpy as np
import numpy.matlib
from scipy import stats
import csv
import datetime
from itertools import zip_longest
from collections import defaultdict
np.seterr(divide='ignore', invalid='ignore') # ignore errors

from setting_variables import N, P
from functions import run_simulation, load_data
from plot_functions import hist


###   run simulation  ###

for i in range(1000): 
    i_last = []
    i_T1 = []
    i_T2 = []
    X, V = run_simulation(i_last, i_T1, i_T2, print_stats = False, plot_env = False, same_speed = False)
    
    # Open the CSV file in write mode
    csv_file_path = f'F80S20_{i+1}.csv'
    with open(csv_file_path, mode='a', newline='') as file:
        writer = csv.writer(file)
    
        # Use zip_longest to iterate over the arrays and pad with None for differing lengths
        for data_row in zip_longest(i_T1, i_T2, i_last, fillvalue=None):
            writer.writerow(data_row)

                
i_T1_arr, i_T2_arr, i_last_arr = load_data('F20S80')     
hist(i_T1_arr, i_T2_arr, file_path = None)