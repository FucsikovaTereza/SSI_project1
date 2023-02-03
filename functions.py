import numpy as np
from numpy.linalg import norm
from plot_functions import plot_environment
from setting_variables import N, P, Points, xini, vini, tmax, h, nmax, v0, tau, lamb, epsilon, U0, xi, U0w, xw

###   choosing atractor   ###

# navigation to travelators 1
def closer_atractor_navigation():
    global N, P, xini
    xA = np.zeros([N,2])                                  # positions of atractors
    r = np.zeros([N,2])                                   # directions to atractors
    for i in range(0, N):
        v = [P[1], P[3]] - (xini[i,:])
        u = [P[2], P[3]] - (xini[i,:])
        
        if norm(v) < norm(u):
            xA[i,:] = [P[1], P[3]]
            r[i,:] = v
        else:
            xA[i,:] = [P[2], P[3]]
            r[i,:] = u
    return r, xA
   
# navigation to travelators 2
def random_atractor_navigation():
    global N, P, xini
    xA = np.zeros([N,2])                                  # positions of atractors
    r = np.zeros([N,2])                                   # directions to atractors
    for i in range(0, N):
        v = [P[1], P[3]] - (xini[i,:])
        u = [P[2], P[3]] - (xini[i,:])
        
        if norm(v) < norm(u):
            xA[i,:] = [P[1], P[3]]
            r[i,:] = v
        else:
            xA[i,:] = [P[2], P[3]]
            r[i,:] = u
    return r, xA


# choose the way to find an atractor; 0 is for choosing closer atractor, 1 for random selection
w = 0
if w == 0:
    r, xA = closer_atractor_navigation()
else:
    r, xA = random_atractor_navigation()

###   forces   ###

# Fmot - attraction to an attractor
def F_mot(v):
    global r, v0, tau
    nr = norm(r, axis = 1)
    nr2 = np.matlib.repmat(nr, 2, 1)
    s_alphaA = r / nr2.transpose()
    v1 = np.matlib.repmat(v0, 1, 2)
    f = ((s_alphaA * v1) - v ) / tau
    return f

# Fint - interactions between agents
def F_int(x, v, f, U, z):
    global N, lamb
    f = np.zeros(np.shape(f))
    nv = norm(v, axis = 1)
    
    for i in range(N):
        xBA = np.matlib.repmat(x[i, :], N, 1) - x
        nxBA = norm(xBA, axis = 1)
        nxBA1 = np.zeros((N, 2))
        
        for j in range(N):
            nxBA1[j, :] = nxBA[j]
        s = xBA / nxBA1
        cosPhi = - np.sum(np.matlib.repmat(v[i, :], N, 1) * xBA, axis = 1) / (nv[i] * nxBA)
        fnow = U0/xi * np.exp(-nxBA1/xi) * s 
        e = (1 - lamb) * (1 + cosPhi)/2 + lamb
        e1 = np.zeros((N, 2))
        
        for j in range(N):
            e1[j, :] = e[j]
        fnow = e1 * fnow
        fnow[np.isnan(fnow)] = 0
        f[i, :] = sum(fnow, 1)
        return f


# finding closest point from the walls
def find_closest_wallpoint(x):
    global N, Points, P
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
    return w, wA

def F_all(x, v):  
    global U0, xi, U0w, xw
    
    if type(x) != np.ndarray:
        x = np.array(x)
    if type(v) != np.ndarray:
        v = np.array(v)
    
    Fmot = F_mot(v)
    Fint = F_int(x, v, Fmot, U0, xi)
    
    # Fenv - interaction with the environment
    w, wA = find_closest_wallpoint(x)
    Fenv = F_int(wA, w, Fmot, U0w, xw)
    
    f = Fmot + Fint - Fenv
    return f, Fmot, Fint, Fenv


###   simulation   ###

def run_simulation(xA):
    global N, P, nmax, xini, vini, h, epsilon
    X = np.zeros((N, 2, nmax+1))
    V = np.zeros((N, 2, nmax+1))
    X[:, :, 0] = xini              # Nth pedestrian, 0 - x position and 1 - y position, iteration
    V[:, :, 0] = vini
    end_y = P[3] + P[6]
    active_cor = [k for k in range(N)]
    active_T = [k for k in range(N)]
    i = 0
    xini_T = np.zeros((N, 2))
    v_T = 3 * np.ones((N, 1))
    
    i_T1 = []
    i_T2 = []
    
    for k in range(10):
        while len(active_T) > 0:
            
            if len(active_cor) != 0: 
                # calling forces
                Ftot, Fmot, Fint, Fenv = F_all(X[:, :, i], V[:, :, i])
                    
                # Euler
                X[active_cor, :, i+1] = X[active_cor, :, i] + h * V[active_cor, :, i]
                V[active_cor, :, i+1] = V[active_cor, :, i] + h * Ftot[active_cor]
                
                    
                for j in active_cor:
                    # print(j, norm(X[j, :, i] - xA[j]))
                    if norm(X[j, :, i] - xA[j]) < epsilon:
                        xini_T[j,:] = X[j, :, i]
                        active_cor.remove(j)
                        print("na travelatoru:", j, X[j, :, i],)
                
            if len(active_cor) < N:
                    # finished?
                for j in active_T:
                    if X[j, 1, i] >= end_y:
                        active_T.remove(j)
                        print("odstranen:", j, X[j, :, i])
                        if xini_T[j, 0] > P[7]:
                            i_T1.append(i)
                        elif xini_T[j,0] > 1 and xini_T[j,0] < P[7]:
                            i_T2.append(i)
            
                for j in active_T:
                    # print(j, "y:",  X[j,1,i])
                    if xini_T[j, 0] > P[7]:  # if in a slow travelator
                        X[j, 0, i+1] = X[j, 0, i]
                        X[j, 1, i+1] = X[j, 1, i] + h * v_T[j]
                        V[j, 1, i+1] = V[j, 1, i]
                            
                    elif xini_T[j,0] > 1 and xini_T[j,0] < P[7]:  # fast travelator
                        X[j, 0, i+1] = X[j, 0, i]
                        X[j, 1, i+1] = X[j, 1, i] + h * (v_T[j] + v0[j])
                        V[j, 1, i+1] = V[j, 1, i]
                              
            # plot
            plot_environment(X[active_T, :, i])
                            
            i += 1
    return X, V, i_T1, i_T2
   