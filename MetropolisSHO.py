from typing import no_type_check_decorator
import numpy as np
import math
import matplotlib.pyplot as plt
global Ncf
Ncf = 10000
global N_cor 
N_cor = 20
global N
N = 20
global a
a = .5

def findSigma(G,sigma): #find sigma^2(input is an array of all configuration)
    global Ncf
    for n in range(0,N):
        sigmaSquared = 0
        for i in range(Ncf):
            sigmaSquared += G[i][n]**2
        sigmaSquared = sigmaSquared/Ncf
        myVar = 0
        for i in range(Ncf):  #find expectation value of Gamma(GammaBar)
            myVar += G[i][n]
        myVar = myVar/Ncf
        myVar = myVar**2
        sigmaSquared = sigmaSquared - myVar
        sigmaSquared = sigmaSquared/Ncf
        sigma[n] = math.sqrt(sigmaSquared)
    print(sigma)

def S(i,x):
    global a
    global N
    return a*x[i]**2/2 + x[i]*(x[i]-x[(i+1)%N]-x[(i-1)%N])/a

def update(x):
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

def Thermalize(x):
    for i in range(0, 5*N_cor):
        update(x)
 
def Run(x,G,y):
    for alpha in range(0,Ncf):
        for j in range(0,N_cor):
            update(x)
        for n in range(0,N):
            y[alpha] = x
            G[alpha][n] = compute_G(x,n)
    
def MCaverage(x,G,avg_G):
    global N
    for n in range(0,N):
        avg_G[n] = 0
        for alpha in range(0,Ncf):
            avg_G[n] += G[alpha][n]
        avg_G[n] = avg_G[n]/Ncf
        #print("G(%d) = %g" % (n,avg_G[n]))

def DeltaE(avg_G,sigma):#throws error when i+1 is negative
    global N
    global Ncf
    global a
    DelE = [0 for n in range(0,N)]
    for n in range(0,N):
        DelE[n] = math.log(avg_G[n]/avg_G[(n+1)%N])/a
    fig, ax = plt.subplots()
    ax.errorbar([x/2 for x in range(N)], DelE, yerr=sigma, fmt='go')
    ax.set_xlim([0, 2])
    ax.set_xlabel('t')
    ax.set_ylabel('Energy')
    ax.set_title('Quantum SHO Simulation of Delta E = (E_0 - E_1), N_cf = %d' % (Ncf))
    plt.show()

def main():
    y = [[0 for n in range(0,N)] for alpha in range(0,Ncf)]
    x = [0 for i in range(0,N)]
    G = [[0 for n in range(0,N)] for alpha in range(0,Ncf)]
    avg_G = [0 for n in range(0,N)]
    sigma = [0 for n in range(0,N)]
    Thermalize(x)
    Run(x,G,y)
    MCaverage(x,G, avg_G)
    findSigma(G,sigma)
    DeltaE(avg_G,sigma)

main()