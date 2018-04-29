version = 105
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import math
import os
os.chdir('../../result/worm')
from mc_301 import observers as observers
def var(a,b,c):
    avr = (a+b+c)/3
    val = np.sqrt((a-avr)**2+(b-avr)**2+(c-avr)**2)/2
    #print a,b,c,val
    return val
    
    
if __name__ == '__main__':
    fig = plt.figure()
    plt.yticks(np.linspace(-2.0,0.1,22))
    yval = [[],[],[],[]]
    xval = [[],[],[],[]]
    yerr = [[],[],[],[]]
    for keys in observers:
        print observers[keys]
        xval[keys[0]].append(observers[keys][0])
        yval[keys[0]].append(observers[keys][1])
        yerr[keys[0]].append(observers[keys][2])
    print xval[0]
    print yval[0]
    color = ['r','g','b','y']
    plt.scatter (xval[0], yval[0], c=color[0], s=5, zorder=5,label = 'N = '+str(1))
    plt.errorbar(xval[0], yval[0], yerr=yerr[0], zorder=0,fmt="none",marker="none")
    plt.scatter (xval[1], yval[1], c=color[1], s=5, zorder=6,label = 'N = '+str(2))
    plt.errorbar(xval[1], yval[1], yerr=yerr[1], zorder=1,fmt="none",marker="none")
    plt.scatter (xval[2], yval[2], c=color[2], s=5, zorder=7,label = 'N = '+str(4))
    plt.errorbar(xval[2], yval[2], yerr=yerr[2], zorder=2,fmt="none",marker="none")
    plt.scatter (xval[3], yval[3], c=color[3], s=5, zorder=8,label = 'N = '+str(8))
    plt.errorbar(xval[3], yval[3], yerr=yerr[3], zorder=3,fmt="none",marker="none")
    plt.legend()
    plt.show()