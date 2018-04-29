from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import math
import os
os.chdir('../result/exhaustion/model_3_3_2')
from data import lenx;
from data import leny;
from data import states as evenoddstates;
os.chdir('../model_3_3_1')
from check_data import states as singlestates;
C = (0,0,0,0)
states = {}
checker= {}
statesnum = 0
for key in evenoddstates:
    if key[0] != C[0]:
        continue
    if key[1] != C[1]:
        continue
    if key[2] != C[2]:
        continue
    if key[3] != C[3]:
        continue
    if(states.has_key((key[4],key[5]+2*key[6]))):
        states[(key[4],key[5]+2*key[6])] += evenoddstates[key]
    else:
        states[(key[4],key[5]+2*key[6])] = evenoddstates[key]
    statesnum += evenoddstates[key]
for key1 in singlestates:
    for key2 in singlestates:
        if key1[0] != 0 :
            continue
        if key1[1] != 0 :
            continue
        if key2[0] != 0 :
            continue
        if key2[1] != 0 :
            continue
        if(checker.has_key((key1[2]+key2[2],key1[3]+key2[3]))):
            checker[(key1[2]+key2[2],key1[3]+key2[3])] += singlestates[key1]*singlestates[key2]
        else:
            checker[(key1[2]+key2[2],key1[3]+key2[3])]  = singlestates[key1]*singlestates[key2]
#for key1 in singlestates: 
#    if key1[0] != 0 :
#        continue
#    if key1[1] != 0 :
#        continue
#    checker[(key1[2]*2,key1[3]*2)] -= singlestates[key1]**2
for key in states:
    print key,states[key],checker[key]
print "---------"
for key in checker:
    print key,checker[key]