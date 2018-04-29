from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import math
import os
os.chdir('../../result/worm')
from worm_propagator_300 import propagators
L = 128
count=200000
def print_ratio(resultid,task,name,p):
    for i in np.arange(len(propagators[(resultid,task,name,0)][5])-1):
        flag = False
        for q in p:
            if np.abs(propagators[(resultid,task,name,q)][5][i+1]) <= 1e-20 or np.abs(propagators[(resultid,task,name,q)][5][i]) <= 1e-20:
                flag = True
        if(flag):
            print 'zero'
            continue
        ratelist = []
        for q in p:
            ratelist.append(-np.log(propagators[(resultid,task,name,q)][5][i+1]/propagators[(resultid,task,name,q)][5][i]))
        print ratelist
print_ratio(3,18,"single",[0,2,4])
#print -np.log(propagators[(3,0,"gamma5")][5][10]/propagators[(3,0,"gamma5")][5][9])
#print len(G)
m = 1.20
def func(i):
    return np.exp(-m*i) + np.exp(-m*(L-i))
def plot_log(resultidgroup,task,name,p):
    fig = plt.figure()
    for resultid in resultidgroup:
        for t in task:
            for n in name:
                for q in p:
                    data = propagators[resultid,t,n,q][5]
                    plt.scatter(np.arange(len(data)),np.log(data)/np.log(10), s=5, zorder=3)
    plt.show()
#plot_log([0,1,2,3],[10],["single"],[0])
def fit_m(L,l,data):
    data_avr = np.zeros(l)
    data_abs = np.zeros(l)
    data_var = np.zeros(l)
    data_err = np.zeros(l)
    base_var = np.zeros(l)
    for i in np.arange(l):
        for j in np.arange(L):
            data_avr[i] += data[j][i]
        data_avr[i] /= L
    data_abs = np.abs(data_avr)
    for i in np.arange(l):
        for j in np.arange(L):
            data_var[i] += (data[j][i]-data_avr[i])**2
    minval = data_abs[0]
    for i in np.arange(l):
        if(np.abs(data_abs[i])<1e-20):
            continue
        if(minval>data_abs[i])
            minval = data_abs[i]
    for i in np.arange(l):
        if(np.abs(data_var[i])<1e-40):
            np.abs += minval**2
    # grad down ?
    # func = A1e-m1x(1 + A2e-dmx)
    data_err = np.sqrt(data_var)
    print data_var[i]
def fit_m_array(resultidgroup,task,name,p):
    L = len(resultidgroup)
    m = []
    if L == 0:
        return m
    for t in task:
        for n in name:
            for q in p:
                data = []
                for resultid in resultidgroup:
                    data.append(propagators[resultid,t,n,q][5])
                m.append(fit_m(L,len(data[0]),data))
    return m
fit_m_array([0,1,2,3],[10],["single"],[0])
    #fig = plt.figure()
#plt.scatter(np.arange(len(propagators[3,10,"gamma5",0][5])),np.log(propagators[3,10,"gamma5",0][5])/np.log(10), c='r', s=5, zorder=3)
#plt.scatter(np.arange(len(propagators[3,40,"single"][5])),np.log(propagators[4,40,"single"][5])/np.log(10), c='b', s=5, zorder=3)
#plt.scatter(np.arange(len(propagators[3,40,"single"][5])),np.log(propagators[5,40,"single"][5])/np.log(10), c='b', s=5, zorder=3)
#plt.plot(np.arange(4),(propagator[1][5][0]*func(np.arange(4))-propagator[1][5][0:4]),'b-')
#plt.show()
