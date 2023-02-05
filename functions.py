import numpy as np
from numpy.linalg import norm
import random
from setting_variables import *
from plot_functions import plot_environment



###   choosing atractor   ###

# navigation to travelators 1
def closer_atractor_navigation(xini):
    global N, P
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
def random_atractor_navigation(xini):
    global N, P
    xA = np.zeros([N,2])                                  # positions of atractors
    r = np.zeros([N,2])                                   # directions to atractors
    for i in range(0, N):
        k = random.randint(0,1)
        if k == 0:
            r[i,:] = [P[1], P[3]] - (xini[i,:])
            xA[i,:] = [P[1], P[3]]
        else:
            r[i,:] = [P[2], P[3]] - (xini[i,:])
            xA[i,:] = [P[2], P[3]]
    return r, xA

# navigation to travelators 3
def navigate_concrete_amount_of_peds(xini, faster_peds):
    global N, P
    xA = np.zeros([N,2])                                  # positions of atractors
    r = np.zeros([N,2])                                   # directions to atractors
    for i in range(faster_peds):
        xA[i,:] = [P[1], P[3]]
        r[i,:] = [P[1], P[3]] - (xini[i,:])
    for i in range(faster_peds, N):
        xA[i,:] = [P[2], P[3]]
        r[i,:] = [P[2], P[3]] - (xini[i,:])
    return r, xA


###   forces   ###

# Fmot - attraction to an attractor
def F_mot(v, r):
    global v0, tau
    nr = norm(r, axis = 1)
    nr2 = np.matlib.repmat(nr, 2, 1)
    s_alphaA = r / nr2.transpose()
    v1 = np.matlib.repmat(v0, 1, 2)
    f = ((s_alphaA * v1) - v ) / tau
    return f

# Fint - interactions between agents
def F_int(x, v, f):
    global N, lamb, U0, xi
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
    return f2

# Fenv - interaction with the environment
def F_env(wA, x, w, f):
    global N, lamb, U0w, xw
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
    return f3


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


def F_all(x, v, xini, vini, r):  
    global U0, xi, U0w, xw   
    if type(x) != np.ndarray:
        x = np.array(x)
    if type(v) != np.ndarray:
        v = np.array(v)
    
    Fmot = F_mot(v, r)
    Fint = F_int(x, v, Fmot)
    w, wA = find_closest_wallpoint(x)
    Fenv = F_env(wA, x, w, Fmot)
    
    f = Fmot + Fint - Fenv
    return f, Fmot, Fint, Fenv


###   simulation   ###

def run_simulation(i_T1, i_T2, print_stats = True, plot_env = True):
    global N, P, nmax, h, epsilon, choice, faster_peds
    
    xini = np.zeros([N,2])                                 # initial position
    xini[:,0] = np.random.randint(P[0], size = N)
    xini[:,1] = np.random.randint(P[8], size = N)
    vini = 0 * xini
    
    if choice == 0:
        r, xA = closer_atractor_navigation(xini)
    elif choice == 1:
        r, xA = random_atractor_navigation(xini)
    elif choice == 2:
        r, xA = navigate_concrete_amount_of_peds(xini, faster_peds)
    
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
    
    for k in range(10):
        while len(active_T) > 0:    
            if len(active_cor) != 0: 
                # calling forces
                Ftot, Fmot, Fint, Fenv = F_all(X[:, :, i], V[:, :, i], xini, vini, r)
                # Euler
                X[active_cor, :, i+1] = X[active_cor, :, i] + h * V[active_cor, :, i]
                V[active_cor, :, i+1] = V[active_cor, :, i] + h * Ftot[active_cor]
                    
                for j in active_cor:
                    # print(j, norm(X[j, :, i] - xA[j]))
                    if norm(X[j, :, i] - xA[j]) < epsilon:
                        xini_T[j,:] = X[j, :, i]
                        active_cor.remove(j)
                        if print_stats == True:
                            print("na travelatoru:", j, X[j, :, i],)
            if len(active_cor) < N:
                    # finished?
                for j in active_T:
                    if X[j, 1, i] >= end_y:
                        active_T.remove(j)
                        if print_stats == True:
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
            if plot_env == True:
                plot_environment(X[active_T, :, i])
                            
            i += 1
    return X, V