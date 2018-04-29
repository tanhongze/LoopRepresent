from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import math
import os
os.chdir('../../result/worm')
from mc_327 import observers as observers
plt.style.use('seaborn-whitegrid')

def jackknife(data):
    num = len(data)
    fsum = np.zeros(num)
    bsum = np.zeros(num)
    jsum = np.zeros(num)
    for i in np.arange(1,num):
        fsum[      i] = fsum[  i-1]+data[  i-1]
        bsum[num-1-i] = bsum[num-i]+data[num-i]
    for i in np.arange(num):
        jsum[i] = (fsum[i]+bsum[i])/(num-1)
    avr = np.mean(data)
    err = np.sqrt((num-1)*np.var(jsum-avr))
    return (avr,err,jsum)
Observers = {}
for keys in observers:
    if Observers.has_key(keys[0:5]):
        Observers[keys[0:5]][1].append(observers[keys][1])
    else :
        Observers[keys[0:5]] = (observers[keys][0],[observers[keys][1]])
graph = {0:([],[],[]),1:([],[],[]),2:([],[],[]),3:([],[],[])}
graph2= []
for i in np.arange(10):
    keys = []
    keys.append((0,i,32,32,1))
    keys.append((0,i,32,32,2))
    keys.append((0,i,32,32,4))
    keys.append((0,i,32,32,8))
    if not Observers.has_key(keys[0]):
        continue
    graph2.append(([],[],[]))
    for j in np.arange(4):
        if not Observers.has_key(keys[j]):
            continue
        (avr,err,jsum) = jackknife(Observers[keys[j]][1])
        graph2[-1][0].append(2**j)
        graph2[-1][1].append(avr)
        graph2[-1][2].append(err)
        if(avr<-1.1):
            continue
        graph[j][0].append(Observers[keys[j]][0])
        graph[j][1].append(avr)
        graph[j][2].append(err)


Avr1=[12,11,7,7,6,5]
Var1=[0.5,0.4,0.3,1,0.3,0.5]
Avr2=[10,8,5,4,3,3]
Var2=[0.4,0.3,0.4,0.6,0.3,0.5]
Time=np.arange(1,7,1)
plt.errorbar(Time,Avr1,yerr=Var1,fmt='.k')
plt.errorbar(Time,Avr2,yerr=Var2,capsize=2)
plt.xlabel('moth')
plt.ylabel('ame/cm') 
plt.show()
        
lines = []
for keys in np.arange(4):
    if not graph.has_key(keys):
        continue
    print graph[keys][0]
    print graph[keys][1]
    print graph[keys][2]
    l = plt.errorbar(graph[keys][0],graph[keys][1],yerr=graph[keys][2],capsize=2)
    lines.append(l)
plt.xlabel(r'$g$')
plt.ylabel(r'$am_c$')
plt.legend(lines,[r'N=1',r'N=2',r'N=4',r'N=8'],loc = 'best')
plt.show()

lines = []
for keys in np.arange(1,4):
    l = plt.errorbar(graph2[keys][0],graph2[keys][1],yerr=graph2[keys][2],capsize=2)
    lines.append(l)
plt.xlabel(r'$N_f$')
plt.ylabel(r'$am_c$')
plt.legend(lines,[r'g='+str(0.01),r'N=',r'N=4',r'N=8'],loc = 'best')
plt.show()