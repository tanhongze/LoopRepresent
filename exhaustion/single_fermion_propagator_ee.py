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
import observers_216_0 as observers
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
        if loop[0] != 0 or loop[1] !=0:
            continue
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
def get_worms1():
    worms1 = {}
    for key in worms:
        if ((key[4]-key[5])&7) == 2:
            continue
        if ((key[4]-key[5])&7) == 6:
            continue
        corner = key[4]+key[5] + ((key[4]+key[5])&1)
        if ((key[4]-key[5])&7) in [7,0,1]:
            sign = 1
        else :
            sign = -1
        key1 = (key[0],key[1],key[3],corner,sign)
        if worms1.has_key(key1):
            worms1[key1] += worms[key]
        else :
            worms1[key1] = worms[key]
    return worms1
def single_fermion_count():
    worms1 = get_worms1()
    # (x,y,n,corner,sign)
    states = {}
    for worm in worms1:
        for loop in loops:
            occu = loop[1] & worm[2]
            if occu != 0:
                continue
            key = (worm[0],worm[1],loop[1]|worm[2],loop[0]+worm[3],worm[4])
            if states.has_key(key) :
                states[key] += loops[loop] * worms1[worm]
            else :
                states[key] = loops[loop] * worms1[worm]
    for loop in loops:
        occu = loop[1] & occupy0
        if occu != 0:
            continue
        key = (0,0,loop[1]|occupy0,loop[0],2)
        if states.has_key(key) :
            states[key] += loops[loop]
        else :
            states[key] = loops[loop]
    state_count = 0
    for key in states:
        state_count += states[key]
    print 'state count',state_count
    pr = {}
    for state in states:
        for loop in loops:
            key = (state[0],state[1],ones[state[2]^loop[1]],ones[state[2]&loop[1]],loop[0]+state[3],state[4])
            if pr.has_key(key) :
                pr[key]+= loops[loop] * states[state]
            else :
                pr[key] = loops[loop] * states[state]
    #print pr
    return pr
def single_fermion_propagator(m,g):
    pr = single_fermion_count()
    x = []
    y = []
    Z0 = Z(0.0,1.1)
    print Z0
    for key in pr:
        x.append(key[0])
        y.append(key[1])
    X = np.max(x)+1
    Y = np.max(y)+1
    propagator = np.zeros((X,Y))
    x = np.arange(X)
    y = np.arange(Y)
    (xg,yg) = np.meshgrid(x,y)
    f1 = (2+m)/(((2+m)**2) + g**2)
    f2 = (1.0)/(((2+m)**2) + g**2)
    for key in pr:
        propagator[key[0]][key[1]] += (f1**key[2]) *(f2**key[3]) *((0.5**0.5)**key[4]) *key[5] *pr[key]
    return (xg,yg,propagator/Z0)
def plot_xy(x,y,z):
    surf = ax.plot_surface(x,y,z)
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('psi')
if __name__ == '__main__':
    check_same()
    loop_count_full = 0
    for loop in loops_topo:
        loop_count_full+=loops_topo[loop]
    print loop_count_full
    loop_count_cast = 0
    for loop in loops:
        loop_count_cast += loops[loop]
    print loop_count_cast
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    (xe,ye,ze) = single_fermion_propagator(0.0,1.1)
    (xw1,yw1,zw1) = (np.array(observers.x0),np.array(observers.y0),np.array(observers.p0))
    surf = ax.plot_surface(xe,ye,ze,cmap='rainbow')
    plot_xy(xw1,yw1,zw1)
    print 'worm:'
    print zw1
    print 'exhaustion:'
    print ze
    print 'differ:'
    print ze-zw1
    plt.show()
    