###   plot   ###
import matplotlib.pyplot as plt
from numpy import arange
from setting_variables import P, m, n
    
def plot_environment(X, file_path = False):
    global P, m, n
    figure, axes = plt.subplots()
    if P[0] <= 2*P[4]:
        raise ValueError('Invalid values for x axis')       
    plt.xlim([0 - 0.5*P[0], 1.5*P[0]])
    plt.ylim([0 - 0.5*P[4], P[6] + P[3]])
    
    for i in range(9):
        plt.plot(m[i], n[i], '-k')
    axes.set_aspect('equal', adjustable='box') 
    for j in range(len(X)):
        axes.add_artist(plt.Circle((X[j, :]), 1.5, color='r'))
        #axes.add_artist(plt.Circle((xA[j, :]), 0.8, color='b'))
    plt.axis('off')
    
    if file_path == True:
        plt.savefig(file_path, format="pdf")
    plt.show() 
        
###   Statistics for specific model   ###

#plot graph for 3 pedestrians

def plot_3ped(X, i_T1, i_T2, file_path = False):
    i_T = i_T1.append(i_T2)
    
    plt.plot(X[0,1,0:int(i_T[0,0])], 'b--', label = "agent 1")
    plt.plot(X[1,1,0:int(i_T[1,0])], 'r--', label = "agent 2")
    plt.plot(X[2,1,0:int(i_T[2,0])], 'g--', label = "agent 3")
    
    plt.plot([int(i_T[0,0]), int(i_T[0,1])], [int(X[0,1,int(i_T[0,0])]), int(X[0,1,int(i_T[0,1])])], 'b-', label = "agent 1 travelátor 1")
    plt.plot([int(i_T[1,0]), int(i_T[1,1])], [int(X[1,1,int(i_T[1,0])]), int(X[1,1,int(i_T[1,1])])], 'r-', label = "agent 2 travelátor 2")
    plt.plot([int(i_T[2,0]), int(i_T[2,1])], [int(X[2,1,int(i_T[2,0])]), int(X[2,1,int(i_T[2,1])])], 'g-', label = "agent 3 travelátor 2")
    
    plt.xlabel("počet iterací")
    plt.ylabel("y-ová osa pro agenty")
    plt.legend()
    
    if file_path == True:
        plt.savefig(file_path, format="pdf")   
    plt.show()


def hist(i_T1, i_T2, file_path = False):
    plt.hist(i_T2, alpha=0.7, label='travelátor 1', rwidth=0.8, color = 'green')
    plt.hist(i_T1, alpha=0.8, label='travelátor 2', rwidth=0.8, color = 'orange')
    
    plt.xlabel("čas v systému")
    plt.ylabel("počet agentů")
    #plt.yticks(arange(451, step=50))
    plt.legend()
    
    if file_path == True:
        plt.savefig(file_path, format="pdf")
    plt.show()