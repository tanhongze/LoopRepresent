version = 105
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import math
import os
os.chdir('../../result/worm')
#from mc_meff_105 import observers
from mc_meff_105 import observers as observers1
from mc_meff_114 import observers as observers2
from mc_meff_117 import observers as observers3
from mc_meff_120 import observers as observers4
def var(a,b,c):
    avr = (a+b+c)/3
    val = np.sqrt((a-avr)**2+(b-avr)**2+(c-avr)**2)/2
    #print a,b,c,val
    return val
def plot105():
    mc = -0.575
    xval = []
    yvalb = []
    yerrb = []
    yvalf = []
    yerrf = []
    for i in range(90):
        xval.append(observers[3,i][3])
        yvalb.append((observers[3,i][6]+observers[4,i][6]+observers[5,i][6])/3)
        yerrb.append(var(observers[3,i][6],observers[4,i][6],observers[5,i][6]))
        yvalf.append((observers[3,i][8]+observers[4,i][8]+observers[5,i][8])/3)
        yerrf.append(var(observers[3,i][8],observers[4,i][8],observers[5,i][8]))
    fig = plt.figure()
    plt.yticks(np.linspace(-0.1,2.1,21))
    plt.scatter(xval, yvalb, c='r', s=5, zorder=3)
    plt.errorbar(xval, yvalb, yerr=yerrb, zorder=0)#, fmt="none",marker="none")
    plt.scatter(xval, yvalf, c='b', s=5, zorder=4)
    plt.errorbar(xval, yvalf, yerr=yerrf, zorder=1)#, fmt="none",marker="none")
    plt.show()
def plot_b(observers,color):
    xval = []
    yval = []
    for key in observers:
        xval.append(observers[key][3])
        yval.append(observers[key][6])
    plt.scatter(xval, yval, c=color, s=5, zorder=3)
def plot_f(observers,color):
    xval = []
    yval = []
    for key in observers:
        xval.append(observers[key][3])
        yval.append(observers[key][8])
    plt.scatter(xval, yval, c=color, s=5, zorder=3)
    
    
if __name__ == '__main__':
    fig = plt.figure()
    plt.yticks(np.linspace(-0.1,1.4,16))
    plot_f(observers4,'b')
    plot_b(observers4,'r')
    plt.show()