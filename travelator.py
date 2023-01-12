###   importing libraries   ###
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm
from math import floor
import numpy.matlib

###   setting variables   ###

# dimensions of the corridor and travelators
corridor_x = 60
corridor_y = 70
travelators_x = 10
travelators_y = 65 


# pedestrians setting
N = 5                                                  # number of pedestrians
xini = np.zeros([N,2])                                 # initial position
xini[:,0] = np.random.randint(corridor_x, size = N)
xini[:,1] = np.random.randint(corridor_y/2, size = N)
vini = 0 * xini                                        # initial speed
P = [corridor_x, corridor_x/3, corridor_x*2/3, corridor_y, travelators_x, travelators_x/2, travelators_y]   # possible values for plot
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



###   plot   ###

def plot_environment(P, m, n, X):
    figure, axes = plt.subplots()
    if corridor_x <= 2*travelators_x:
        raise ValueError('Invalid values for x axis')
        
    plt.xlim([0 - 0.5*P[0], 1.5*P[0]])
    plt.ylim([0 - 0.5*P[4], P[6] + P[3]])
    
    for i in range(Points.shape[0]):
        plt.plot(m[i], n[i], '-k')

    axes.set_aspect('equal', adjustable='box')
    
    
    for j in range(len(X)):
        axes.add_artist(plt.Circle((X[j, :]), 1.5, color='r'))
        #axes.add_artist(plt.Circle((xA[j, :]), 0.8, color='b'))
    
    
    plt.axis('off')
    #plt.savefig("425.pdf", format="pdf")
    plt.show() 
    
plot_environment(P, m, n, xini)


    
###   move   ###

# navigation to travelators
xA = np.zeros([N,2])                                  # positions of atractors
r = np.zeros([N,2])                                   # directions to atractors
for i in range(0, N):
    v = [(corridor_x/3), corridor_y] - (xini[i,:])
    u = [(corridor_x*2/3), corridor_y] - (xini[i,:])
    
    if norm(v) < norm(u):
        xA[i,:] = [(corridor_x/3), corridor_y]
        r[i,:] = v
    else:
        xA[i,:] = [(corridor_x*2/3), corridor_y]
        r[i,:] = u
        
nr = norm(r, axis = 1)


# forces
def F(x, v):
    global tau, v0, xA, U0, xi, lamb
    
    if type(x) != numpy.ndarray:
        x = np.array(x)
    if type(v) != numpy.ndarray:
        v = np.array(v)
    
    # Fmot - attraction to an attractor
    nr2 = np.matlib.repmat(nr, 2, 1)
    s_alphaA = r / nr2.transpose()
    v1 = np.matlib.repmat(v0, 1, 2)
    f = ((s_alphaA * v1) - v ) / tau
    Fmot = f
    
    # Fint - interactions between agents
    f2 = np.zeros(np.shape(f))
    nv = norm(v, axis = 1)
    for i in range(N):
        xBA = np.matlib.repmat(x[i, :], N, 1) - x
        nxBA = norm(xBA, axis = 1)
        nxBA1 = np.zeros((N, 2))
        for j in range(N):
            nxBA1[j, :] = nxBA[j]
        s = xBA / nxBA1
        cosPhi = - np.sum(np.matlib.repmat(v[i, :], N, 1) * xBA, axis = 1) / (nv[i] * nxBA)
        f2now = U0/xi * np.exp(-nxBA1/xi) * s 
        e = (1 - lamb) * (1 + cosPhi)/2 + lamb
        e1 = np.zeros((N, 2))
        for j in range(N):
            e1[j, :] = e[j]
        f2now = e1 * f2now
        f2now[np.isnan(f2now)] = 0
        f2[i, :] = sum(f2now, 1)
    
    # Fenv - interaction with the environment
    
    # finding closest point from the walls
    wA = np.zeros((N,2))
    w = np.zeros((N,2))
    for i in range(N):
        A = [0, x[i,1]]
        B = [P[0], x[i,1]]
        if 0 <= x[i,0] <= Points[2,0,1] or Points[3,0,0] <= x[i,0] <= Points[3,0,1] or Points[4,0,0] <= x[i,0] <= Points[4,0,1]:
            C = [x[i,0], P[3]]
            XC = C - x[i,:]
        else: 
            XC = np.nan     
        XA = A - x[i,:]
        XB = B - x[i,:] 
        short = min(norm(XA), norm(XB), norm(XC))
        if short == norm(XA):
            wA[i,:] = A
            w[i,:] = XA
        elif short == norm(XB):
            wA[i,:] = B
            w[i,:] = XB
        else:
            wA[i,:] = C
            w[i,:] = XC
        
    f3 = np.zeros(np.shape(f))
    nw = norm(w, axis = 1)
    for i in range(N):
        xWall = np.matlib.repmat(wA[i, :], N, 1) - x
        nxWall = norm(xWall, axis = 1)
        nxWall1 = np.zeros((N, 2))
        for j in range(N):
            nxWall1[j, :] = nxWall[j]
        sw = xWall / nxWall1
        cosPhiw = - np.sum(np.matlib.repmat(w[i, :], N, 1) * xWall, axis = 1) / (nw[i] * nxWall)
        f3now = U0w/xw * np.exp(-nxWall1/xw) * sw 
        ew = (1 - lamb) * (1 + cosPhiw)/2 + lamb
        ew1 = np.zeros((N, 2))
        for j in range(N):
            ew1[j, :] = ew[j]
        f3now = ew1 * f3now
        f3now[np.isnan(f3now)] = 0
        f3[i, :] = sum(f3now, 1)
    
    f = f + f2 #+ f3
    Fint = f2
    Fenv = f3
    
    return f, Fmot, Fint, Fenv


###   initialisation   ###

tmax = 70
h = 0.1                        # diskretisation step
nmax = floor(tmax/h)
v0 = 2 * np.ones((N, 1))       # optimal velocity [m/s]
X = np.zeros((N, 2, nmax+1))
V = np.zeros((N, 2, nmax+1))
X[:, :, 0] = xini              # Nth pedestrian, 0 - x position and 1 - y position, iteration
V[:, :, 0] = vini
xi = 1.1                       # reach of pedestrians [m] 
xw = 0.9
tau = 0.1                      # reaction time [s]
U0 = 55                        # force of interaction between agents [N]
U0w = 50
lamb = 0
epsilon = travelators_x/2      # distance from the exit [m]
Yi = np.zeros((N))



###   simulation   ###

end_y = corridor_y + travelators_y
active_cor = [k for k in range(N)]
active_T = [k for k in range(N)]
i = 0
xini_T = np.zeros((N, 2))
v_T = 3 * np.ones((N, 1))

while len(active_T) > 0:
    
    
    if len(active_cor) != 0: 
        # calling forces
        Ftot, Fmot, Fint, Fenv = F(X[:, :, i], V[:, :, i])
        
        # Euler
        X[active_cor, :, i+1] = X[active_cor, :, i] + h * V[active_cor, :, i]
        V[active_cor, :, i+1] = V[active_cor, :, i] + h * Ftot[active_cor]
    
        
        for j in active_cor:
            print(j, norm(X[j, :, i] - xA[j]))
            if norm(X[j, :, i] - xA[j]) < epsilon:
                xini_T[j,:] = X[j, :, i]
                active_cor.remove(j)
                print("na travelatoru:", j, X[j, :, i],)
    
    if xini_T != np.array([]):

        for j in active_T:
            if X[j, 1, i] != np.nan:
                print(j, "y:",  X[j,1,i])
                if xini_T[j, 0] > corridor_x/2:
                    X[j, 0, i+1] = X[j, 0, i]
                    X[j, 1, i+1] = X[j, 1, i] + h * v_T[j]
                    V[j, 1, i+1] = V[j, 1, i]
                    
                elif xini_T[j,0] > 1 and xini_T[j,0] < corridor_x/2:
                    X[j, 0, i+1] = X[j, 0, i]
                    X[j, 1, i+1] = X[j, 1, i] + h * (v_T[j] + v0[j])
                    V[j, 1, i+1] = V[j, 1, i]
                        
                if X[j, 1, i] >= end_y:
                    active_T.remove(j)
                    print("odstranen:", j, X[j, :, i])
                    Yi[j] = i
                    X[j, :, i+1:nmax+1] = np.nan
                   
    # plot
    plot_environment(P, m, n, X[active_T, :, i])
                
    i += 1

#X = np.save('X.npy', X)

###   Statistics for specific model   ###

#tensor of positions from the first to the last step on the travelators
#X_T = X[:,:,int(min(i_T[:,0])):int(max(i_T[:,1]))]
#xA_i = xA[:,0]
#xA_i[xA_i == 20] = 1
#xA_i[xA_i == 40] = 2
#x_TA = np.hstack((X_T, xA_i))

#X = np.load('X.npy', X)

#plt.plot(X[0,1,0:int(i_T[0,0])], 'b--', label = "agent 1 koridor")
#plt.plot(X[1,1,0:int(i_T[1,0])], 'r--', label = "agent 2 koridor")
#plt.plot(X[2,1,0:int(i_T[2,0])], 'g--', label = "agent 3 koridor")
#plt.plot([int(i_T[0,0]), int(i_T[0,1])], [int(X[0,1,int(i_T[0,0])]), int(X[0,1,int(i_T[0,1])])], 'b-', label = "agent 1 travelátor 1")
#plt.plot([int(i_T[1,0]), int(i_T[1,1])], [int(X[1,1,int(i_T[1,0])]), int(X[1,1,int(i_T[1,1])])], 'r-', label = "agent 2 travelátor 2")
#plt.plot([int(i_T[2,0]), int(i_T[2,1])], [int(X[2,1,int(i_T[2,0])]), int(X[2,1,int(i_T[2,1])])], 'g-', label = "agent 3 travelátor 2")
#plt.xlabel("počet iterací")
#plt.ylabel("y-ová osa pro agenty")
#plt.legend()
#plt.savefig("rychlejsi.pdf", format="pdf")
#plt.show()


