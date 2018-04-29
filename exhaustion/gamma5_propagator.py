from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import math
import os
os.chdir('../result/exhaustion/model_3_3_1')
from check_data import states as check_states
from propagator_data import loops as loops_topo
from propagator_data import worms
from propagator_data import occupy0
os.chdir('../../worm')
import observers_354_0 as observers
def get_ones(len):
    ones = [0,1]
    for i in range(2,len):
        ones.append(ones[i/2]+i%2)
    return ones
ones = get_ones(512)
def check_same():
    for keys in loops_topo:
        check_states[(keys[0],keys[1],keys[2],ones[keys[3]])] -= loops_topo[keys]
    for keys in check_states:
        if check_states[keys] != 0:
            print 'some problems with',keys
def get_loops():
    loops = {}
    for loop in loops_topo:
        key = (loop[2],loop[3])
        if loops.has_key(key):
            loops[key] += loops_topo[loop]
        else :
            loops[key] = loops_topo[loop]
    return loops
loops = get_loops()
def Z(m,g):
    f1 = (2+m)/(((2+m)**2) + g**2)
    f2 = (1.0)/(((2+m)**2) + g**2)
    Z = 0
    for loop1 in loops:
        for loop2 in loops:
            n1 = ones[loop1[1]^loop2[1]]
            n2 = ones[loop1[1]&loop2[1]]
            corner = loop1[0]+loop2[0]
            Z += (0.707**corner)*(f1**n1)*(f2**n2)*loops[loop1]*loops[loop2]
    return Z
def print_occupy(ctx,X,Y):
    for i in range(X):
        line = []
        for j in range(Y):
            if(ctx&(1<<(i*Y+j)))!=0:
                line.append(1)
            else:
                line.append(0)
        print line
def gamma5_propagator(m,g):
    x = []
    y = []
    for worm in worms:
        x.append(worm[0])
        y.append(worm[1])
    X = np.max(x)+1
    Y = np.max(y)+1
    states = {}
    for loop in loops:
        for worm in worms:
            if (loop[1]&worm[3]) != 0:
                continue
            key = (worm[0],worm[1],worm[2],worm[3]|loop[1],loop[0]+worm[4]+worm[5],worm[4]-worm[5])
            if(states.has_key(key)):
                states[key]+= loops[loop]*worms[worm]
            else :
                states[key] = loops[loop]*worms[worm]
    state_count = np.zeros((X,Y))
    for key in states:
        state_count[key[0]][key[1]] += states[key]
    xg,yg= np.meshgrid(np.arange(X),np.arange(Y))
    count = {}
    for state1 in states:
        for state2 in states:
            if state1[0]!=state2[0]:
                continue
            if state1[1]!=state2[1]:
                continue
            if (state1[2]^2)==state2[2]:
                continue
            if (((state1[2]+state1[5])&3)^2) == ((state2[2]+state2[5])&3):
                continue
            headd1 = state1[2]
            headd2 = state2[2]
            taild1 = ((state1[2]+state1[5])&3)^2
            taild2 = ((state2[2]+state2[5])&3)^2
            rot1 = headd1-headd2
            if rot1 == 3:
                rot1 =-1
            if rot1 ==-3:
                rot1 = 1
            rot2 = taild2-taild1
            if rot2 == 3:
                rot2 =-1
            if rot2 ==-3:
                rot2 = 1
            phase = (state1[5]-state2[5]+rot1+rot2)&7
            if phase == 0:
                sign = 1
            elif phase == 4:
                sign = -1
            else :
                print 'error'
            key = (state1[0],state1[1],ones[state1[3]^state2[3]],ones[state1[3]&state2[3]]
                ,state1[4]+state2[4]+(headd1!=headd2)+(taild1!=taild2)
                ,sign)
            if(count.has_key(key)):
                count[key]+= states[state1]*states[state2]
            else :
                count[key] = states[state1]*states[state2]
    for loop1 in loops:
        for loop2 in loops:
            if (loop1[1]&occupy0) != 0:
                continue
            if (loop2[1]&occupy0) != 0:
                continue
            key = (0,0,ones[loop1[1]^loop2[1]],ones[(loop1[1]&loop2[1])|occupy0],loop1[0]+loop2[0],2)
            if(count.has_key(key)):
                count[key]+= loops[loop1]*loops[loop2]
            else :
                count[key] = loops[loop1]*loops[loop2]
    count_count = np.zeros((2,X,Y))
    for key in count:
        count_count[(key[5]+1)/2][key[0]][key[1]] += count[key]
    print count_count
    Z0 = Z(m,g)
    f1 = (2+m)/(((2+m)**2) + g**2)
    f2 = (1.0)/(((2+m)**2) + g**2)
    propagator = np.zeros((X,Y))
    for key in count:
        propagator[key[0]][key[1]] += (f1**key[2])*(f2**key[3])*((0.5**0.5)**key[4])*key[5]*count[key]
    return (xg,yg,propagator/Z0)
def plot_xy(x,y,z):
    surf = ax.plot_surface(x,y,z)
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('psi')
if __name__ == '__main__':
    check_same()
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    (xe,ye,ze) = gamma5_propagator(observers.m0,observers.g0)
    (xw1,yw1,zw1) = (np.array(observers.x0),np.array(observers.y0),np.array(observers.q0))
    xw2 = []
    yw2 = []
    zw2 = []
    for i in np.arange(3):
        for j in np.arange(3):
            xw2.append(xw1[i][j])
            yw2.append(yw1[i][j])
            zw2.append(zw1[i][j])
    surf = ax.plot_surface(xe,ye,ze,cmap='rainbow')
    scat = ax.scatter(xw1,yw2,zw2,c='b',alpha=0.4,s=10)
    #plot_xy(xw1,yw1,zw1)
    print 'worm:'
    print zw1
    print 'exhaustion:'
    print ze
    print 'differ:'
    #print ze-zw1
    plt.show()
    