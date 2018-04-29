from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import math
import os
os.chdir('../result/exhaustion/model_4_4_2')
from data import lenx;
from data import leny;
from data import states as evenoddstates;
os.chdir('../../worm')
from observers_355_0 import m0
from observers_355_0 import g0
from observers_355_0 import p0 as single_fermion_propagator

C = (0,0,0,0)
states = {}
statesnum = 0
for key in evenoddstates:
    if ((1<<key[0])&C[0]) != 0:
        continue
    if ((1<<key[1])&C[1]) != 0:
        continue
    if ((1<<key[2])&C[2]) != 0:
        continue
    if ((1<<key[3])&C[3]) != 0:
        continue
    if(states.has_key((key[4],key[5],key[6]))):
        states[(key[4],key[5],key[6])]+= evenoddstates[key]
    else :
        states[(key[4],key[5],key[6])] = evenoddstates[key]
    statesnum += evenoddstates[key]
print statesnum 
area = lenx*leny;
maxc = 2*area;
def prob(c,n1,n2,f1,f2):
    return (0.5**(c/2))*(f1**n1)*(f2**n2)
def check(F1,F2,Z):
    Z0 = F1 - F1
    for key in states:
        i = key[0]
        j = key[1]
        k = key[2]
        Z0 += prob(i,j,k,F1,F2)*states[i,j,k]/Z
    return Z0
def navr(F1,F2,Z):
    nAvr = [F1-F1,F1-F1,F1-F1]
    for key in states:
        i = key[0]
        j = key[1]
        k = key[2]
        nAvr[0] += states[i,j,k]*prob(i,j,k,F1,F2)*(area - j - k)
        nAvr[1] += states[i,j,k]*prob(i,j,k,F1,F2)*j
        nAvr[2] += states[i,j,k]*prob(i,j,k,F1,F2)*k
    return (nAvr[0]/Z,nAvr[1]/Z,nAvr[2]/Z)
def nvar(F1,F2,Z,nAvr):
    nVar = [F1-F1,F1-F1,F1-F1]
    for key in states:
        i = key[0]
        j = key[1]
        k = key[2]
        nVar[0] += states[i,j,k]*prob(i,j,k,F1,F2)*(area - j - k - nAvr[0])**2
        nVar[1] += states[i,j,k]*prob(i,j,k,F1,F2)*(j-nAvr[1])**2
        nVar[2] += states[i,j,k]*prob(i,j,k,F1,F2)*(k-nAvr[2])**2
    return (nVar[0]/Z,nVar[1]/Z,nVar[2]/Z)
def ncovar(F1,F2,Z,nAvr):
    ncoVar = [F1-F1,F1-F1,F1-F1]
    for key in states:
        i = key[0]
        j = key[1]
        k = key[2]
        ncoVar[0] += states[i,j,k]*prob(i,j,k,F1,F2)*(k-nAvr[2])*(j-nAvr[1])
        ncoVar[1] += states[i,j,k]*prob(i,j,k,F1,F2)*(k-nAvr[2])*(area - j - k - nAvr[0])
        ncoVar[2] += states[i,j,k]*prob(i,j,k,F1,F2)*(j-nAvr[1])*(area - j - k - nAvr[0])
    return (ncoVar[0]/Z,ncoVar[1]/Z,ncoVar[2]/Z)
def nijsumvar(F1,F2,Z,nAvr):
    nijSumSquareSum = [F1-F1,F1-F1,F1-F1]
    for key in states:
        i = key[0]
        j = key[1]
        k = key[2]
        nijSumSquareSum[0] += states[i,j,k]*prob(i,j,k,F1,F2)*(k+j)**2
        nijSumSquareSum[1] += states[i,j,k]*prob(i,j,k,F1,F2)*(area - j)**2
        nijSumSquareSum[2] += states[i,j,k]*prob(i,j,k,F1,F2)*(area - k)**2
    return (nijSumSquareSum[0]/Z - (nAvr[1]+nAvr[2])**2,nijSumSquareSum[1]/Z - (nAvr[2]+nAvr[0])**2,nijSumSquareSum[2]/Z - (nAvr[0]+nAvr[1])**2)
def chi(f1,f2,nAvr,nVar):
    return -(f2*nAvr[1]+2*f1*f1*nAvr[0])/f1/area
def Cchi(f1,f2,nAvr,nVar,nijSumVar):
    return -((4*f1**4-2*f1**2*f2)*nVar[0] \
            +(f2**2-2*f1**2*f2)*nVar[1]   \
            +(2*f1**2*f2)*nijSumVar[2]       \
            -(4*f1**4-2*f1**2*f2)*nAvr[0] \
            -(f2*f2)*nAvr[1]              )/f1/f1/area;
def Z_cff(F1,F2):
    Z = F1 - F1
    for key in states:
        i = key[0]
        j = key[1]
        k = key[2]
        Z = Z + prob(i,j,k,F1,F2)*states[i,j,k]
    return Z
def plotXYZ(X,Y,Z,namex,namey,namez):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    ax.set_xlabel(namex)
    ax.set_ylabel(namey)
    ax.set_zlabel(namez)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.savefig(namex+namey+namez+'+0.3'+'.png')
    plt.show()
def plotlnZ():
    msamples = np.arange(-0.5,1.0,0.05)
    gsamples = np.arange(0.0,1.0,0.1)
    M, G = np.meshgrid(msamples, gsamples)
    e1 = 2+M
    e0 = e1*e1+G*G
    F0 = e0/e0
    F1 = e1/e0
    F2 = 1.0/e0
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(M, G, Z_cff(F1,F2), cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    name = 'lnZ'
    ax.set_zlabel(name)
    ax.set_xlabel('g')
    ax.set_ylabel('m bare')
    fig.colorbar(surf, shrink=0.5, aspect=5, )
    #plt.savefig('Analyze-m bare-g-'+name+'.png')
    plt.show()
def print2D(F1,F2,M,G,nAvr,nVar,Chi,CChi):
    for i in np.arange(len(M[0])):
        print F1[0][i]," ",F2[0][i]," ",M[0][i]," ",Chi[0][i]," ",CChi[0][i]
    for i in np.arange(len(M[0])):
        print M[0][i],nAvr[1][0][i],nAvr[2][0][i],nVar[1][0][i],nVar[2][0][i]
def plot2D(M,G,nAvr,nVar,Chi,CChi):
    fig = plt.figure()
    plt.plot(M[0],CChi[0])
    plt.show()
def plot3D(M,G,nAvr,nVar,Chi,CChi):
    plotXYZ(M,G,nAvr[0]+nAvr[1]+nAvr[2],'m_bare','g','nAvr')
    plotXYZ(M,G,nAvr[0],'m_bare','g','nAvr0')
    plotXYZ(M,G,nAvr[1],'m_bare','g','nAvr1')
    plotXYZ(M,G,nAvr[2],'m_bare','g','nAvr2')
    plotXYZ(M,G,nVar[0],'m_bare','g','nVar0')
    plotXYZ(M,G,nVar[1],'m_bare','g','nVar1')
    plotXYZ(M,G,nVar[2],'m_bare','g','nVar2')
    plotXYZ(M,G,Chi,'m_bare','g','chi')
    plotXYZ(M,G,check(F1,F2,Z),'m_bare','g','check')
    plotXYZ(M,G,CChi,'m_bare','g','Cchi')
if __name__ == '__main__':
    
    #msamples = np.arange(-0.3,0.3,0.02)
    #gsamples = np.arange(0.0,1.0,0.05)
    msamples = np.array([m0])
    gsamples = np.array([g0])
    
    M, G = np.meshgrid(msamples, gsamples)
    e1 = 2+M
    e0 = e1*e1+G*G
    F0 = e0/e0
    F1 = e1/e0
    F2 = 1.0/e0
    print M,G
    print F0,F1,F2
    Z = Z_cff(F1,F2)
    nAvr = navr(F1,F2,Z)
    nVar = nvar(F1,F2,Z,nAvr)
    nCoVar = ncovar(F1,F2,Z,nAvr)
    nijSumVar = nijsumvar(F1,F2,Z,nAvr)
    Chi = chi(F1,F2,nAvr,nVar)
    print '<xi(z,t)xi(z,t)>',single_fermion_propagator[0][0]
    print 'Chi',Chi
    CChi=Cchi(F1,F2,nAvr,nVar,nijSumVar)
    print2D(F1,F2,M,G,nAvr,nVar,Chi,CChi)
    plot2D(M,G,nAvr,nVar,Chi,CChi)
    #plot3D(M,G,nAvr,nVar,Chi,CChi)