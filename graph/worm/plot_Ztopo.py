from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import math
import os
os.chdir('../../result/worm')
from Ztopo_348 import observers as observers
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
def cared_key(g,N,m,X):
    l = []
    for i in np.arange(m):
        l.append((0,g,X,X,N,i))
    return l
record_num = 50
cared_keys = {
0:cared_key(0,1,record_num, 8),
1:cared_key(0,1,record_num,16),
2:cared_key(0,1,record_num,32),
3:cared_key(0,2,record_num, 8),
4:cared_key(0,2,record_num,16),
5:cared_key(0,2,record_num,32)
}
graph = {0:([],[],[])
,1:([],[],[])
,2:([],[],[])
,3:([],[],[])
,4:([],[],[])
,5:([],[],[])}
for i in np.arange(6):
    for keys in cared_keys[i]:
        if observers.has_key(keys):
            graph[i][0].append(observers[keys][1])
            graph[i][1].append(observers[keys][2])
            graph[i][2].append(observers[keys][3])
        else :
            continue
        
lines = []
for keys in np.arange(6):
    l = plt.errorbar(graph[keys][0],graph[keys][1],yerr=graph[keys][2],capsize=2)
    lines.append(l)
plt.xlabel(r'$am$')
plt.ylabel(r'$Z^{0100..00}$')
plt.legend(lines,[r'N=1,Z=T=8',r'N=1,Z=T=16',r'N=1,Z=T=32',r'N=2,Z=T=8',r'N=2,Z=T=16',r'N=2,Z=T=32'],loc = 'best')
plt.show()
