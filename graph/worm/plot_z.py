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
from Ztopo_meff_149_0 import observers as observers1
from Ztopo_meff_149_1 import observers as observers2
from Ztopo_meff_149_2 import observers as observers3
def var(a,b,c):
    avr = (a+b+c)/3
    val = np.sqrt((a-avr)**2+(b-avr)**2+(c-avr)**2)/2
    #print a,b,c,val
    return val
def xaxis(observers):
    xval = []
    for key in observers:
        xval.append(observers[key][0])
    plt.plot([np.min(xval),np.max(xval)],[0.0,0.0])
    return xval
def plot_z(xval,observers,color):
    yval = []
    for key in observers:
        yval.append(observers[key][2])
    plt.plot(xval, yval,color)
    
    
if __name__ == '__main__':
    mc = 3.6
    fig = plt.figure()
    xval = xaxis(observers1)
    plt.yticks(np.linspace(-0.3,1.2,16))
    plot_z(xval,observers1,'b-')
    plot_z(xval,observers2,'r-')
    plot_z(xval,observers3,'g-')
    plt.show()