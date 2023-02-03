###   setting variables   ###

import numpy as np
from math import floor

# dimensions of the corridor and travelators
corridor_x = 60
corridor_y = 70
travelators_x = 10
travelators_y = 65 
    
    
# pedestrians setting
N = 20                                                  # number of pedestrians
xini = np.zeros([N,2])                                 # initial position
xini[:,0] = np.random.randint(corridor_x, size = N)
xini[:,1] = np.random.randint(corridor_y/2, size = N)
vini = 0 * xini 
    
# plot setting                                       # initial speed
P = [corridor_x, corridor_x/3, corridor_x*2/3, corridor_y, travelators_x, travelators_x/2, travelators_y, corridor_x/2]   # possible values for plot
Points = np.array([[(0, 0), (0, P[3])], 
                 [(P[0], P[0]), (0, P[3])], 
                 [(0, P[1]-P[5]), (P[3], P[3])], 
                 [(P[1]+P[5], P[2]-P[5]), (P[3], P[3])],
                 [(P[2]+P[5], P[0]), (P[3], P[3])],
                 [(P[1]-P[5], P[1]-P[5]), (P[3], P[3]+P[6])],
                 [(P[1]+P[5], P[1]+P[5]), (P[3], P[3]+P[6])],
                 [(P[2]-P[5], P[2]-P[5]), (P[3], P[3]+P[6])],
                 [(P[2]+P[5], P[2]+P[5]), (P[3], P[3]+P[6])]
                 ])
m, n = zip(*Points)


#SFM variables
    
tmax = 70
h = 0.1                        # diskretisation step
nmax = floor(tmax/h)
v0 = 2 * np.ones((N, 1))       # optimal velocity [m/s]
xi = 1.1                       # reach of pedestrians [m] 
xw = 0.9
tau = 0.1                      # reaction time [s]
U0 = 55                        # force of interaction between agents [N]
U0w = 50
lamb = 0
epsilon = travelators_x/2      # distance from the exit [m]
