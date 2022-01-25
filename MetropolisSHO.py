from typing import no_type_check_decorator
import numpy as np
import math
import matplotlib.pyplot as plt
global Ncf
Ncf = 20000
global N_cor 
N_cor = 20
global N
N = 20
global a
a = .5
global binsize
binsize = 1
HighPrecision = False

def findSigma(G,sigma): #find uncertainty estimate(input is an array of all configuration)
    global Ncf
    global binsize
    newG = bin(G,binsize)
    for n in range(0,N):
        sigmaSquared = 0
        for i in range(0,int(Ncf/binsize)):
            sigmaSquared += newG[i][n]**2
        sigmaSquared = sigmaSquared/(Ncf/binsize)
        myVar = 0
        for i in range(0,int(Ncf/binsize)):  #find expectation value of Gamma(GammaBar)
            myVar += newG[i][n]
        myVar = myVar/(Ncf/binsize)
        myVar = myVar**2
        sigmaSquared = sigmaSquared - myVar
        sigmaSquared = sigmaSquared/(Ncf/binsize)
        sigma[n] = math.sqrt(sigmaSquared)

def S(i,x):#finds the action
    global a
    global N
    if(HighPrecision == True):
        return x[i]*(x[i]-x[(i+1)%N]-x[(i-1)%N])/(a) + x[i]*(-x[(i+2)%N]+4*x[(i+1)%N]-3*x[i]+4*x[(i-1)%N]-x[(i-2)%N])/(12*a) + a*x[i]**2/2
    else:
        return x[i]*(x[i]-x[(i+1)%N]-x[(i-1)%N])/(a) + a*x[i]**2/2

def update(x):#updates the configuration
    eps = 1.4
    for i in range(0,N):
        old_x = x[i]
        old_Si = S(i,x)
        x[i] += np.random.uniform(-eps, eps)
        dS = S(i,x) - old_Si
        if dS>0 and math.exp(-dS)<np.random.uniform(0,1):
            x[i] = old_x

def compute_G(x,n):
    global N
    g = 0
    for i in range(0,N):
        g = g + x[i]*x[(i+n)%N]
    return g/N

def Thermalize(x):#initializes the lattice
    for i in range(0, 5*N_cor):
        update(x)
 
def Run(x,G,y):#runs the MC simulation
    for alpha in range(0,Ncf):
        for j in range(0,N_cor):
            update(x)
        for n in range(0,N):
            y[alpha] = x
            G[alpha][n] = compute_G(x,n)
    
def MCaverage(x,G,avg_G):
    global N
    global binsize
    newG = bin(G,binsize)
    for n in range(0,N):
        avg_G[n] = 0
        for alpha in range(0,int(Ncf/binsize)):
            avg_G[n] += newG[alpha][n]
        avg_G[n] = avg_G[n]/(Ncf/binsize)
        #print("G(%d) = %g" % (n,avg_G[n]))

def DeltaE(avg_G,sigma):#throws error when i+1 is negative
    global N
    global Ncf
    global a
    DelE = [0 for n in range(0,N)]
    for n in range(0,N):
        DelE[n] = math.log(abs(avg_G[n]/avg_G[(n+1)%N]))/a
    fig, ax = plt.subplots()
    ax.errorbar([x/2 for x in range(N)], DelE, yerr=sigma, fmt='go')
    #ax.set_xlim([0, 2])
    ax.set_xlabel('t')
    ax.set_ylabel('Energy')
    ax.set_title('Quantum SHO Simulation of Delta E = (E_0 - E_1), N_cf = %d' % (Ncf))
    plt.show()

def boostrap(G):#creates bootstrap copies of ensembles to estimate error
    Ncf = len(G)
    G_bootstrap = []
    for i in range(0,Ncf):
        alpha = int(np.random.uniform(0,Ncf))
        G_bootstrap.append(G[alpha])
    return G_bootstrap

def bin(G,binsize):#bins G to save space
    G_binned = []
    for i in range(0,len(G),binsize):
        G_avg = [0 for k in range(0,N)]
        for j in range(0,binsize):
            for n in range(0,N):
                G_avg[n] = G_avg[n] + G[i+j][n]/binsize
        G_binned.append(G_avg)
    return G_binned

def avg_y(y):#averages y to attempt to find the approximate wavefunction
    avg = [0 for i in range(0,N)]
    for n in range(0,N):
        for alpha in range(0,Ncf):
            avg[n] += y[alpha][n]
        avg[n] = avg[n]/Ncf
    return avg

def main():
    y = [[0 for n in range(0,N)] for alpha in range(0,Ncf)]
    x = [0 for i in range(0,N)]
    G = [[0 for n in range(0,N)] for alpha in range(0,Ncf)]
    avg_G = [0 for n in range(0,N)]
    sigma = [0 for n in range(0,N)]
    Thermalize(x)
    Run(x,G,y)
    if(1==1):#finds Delta_E
        MCaverage(x,G, avg_G)
        findSigma(G,sigma)
        DeltaE(avg_G,sigma)
    else: #shows average wavefunction found
        fig, ax = plt.subplots()
        ax.errorbar([x/2 for x in range(N)], avg_y(y), fmt='go')
        ax.set_xlabel('t')
        ax.set_ylabel('Energy')
        ax.set_title('Quantum SHO Simulation of Delta E = (E_0 - E_1), N_cf = %d' % (Ncf))
        plt.show()

main()