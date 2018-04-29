from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import math
import os

os.chdir('../../result/worm')
from worm_propagator_375 import propagators as propagators30
from worm_propagator_376 import propagators as propagators25
from worm_propagator_377 import propagators as propagators10
propagators = propagators30
from worm_propagator_375 import propagators as propagators0
from worm_propagator_376 import propagators as propagators1
from worm_propagator_377 import propagators as propagators2
from worm_propagator_378 import propagators as propagators3
from worm_propagator_380 import propagators as propagators4
from worm_propagator_385 import propagators as propagators5
plt.style.use('seaborn-whitegrid')
def fit_propagator(input_data,input_err2,cutdown,count,T):
    r = 0.0
    epsilon = 1e-21
    t = np.arange(cutdown,T-cutdown)
    data = np.log(np.abs(input_data[cutdown:T-cutdown])+epsilon)
    A = np.mean(data)
    data/= A
    m = 0.1
    err2 = input_err2[cutdown:T-cutdown]*1.0
    err2/= np.min(err2)
    weight = 1/(1+((data-np.max(data))**2))
    def plot_cmp(curve,data):
        plt.plot(t,curve,'r-')
        plt.plot(t,data,'g+')
        plt.show()
    def func(m,cutdown,T,t):
        val = np.log(np.exp(-m*t)+np.exp(-m*(T-t)))
        A_eff = np.mean(val)
        val = (A_eff**-1)*val+epsilon
        return val
    print data
    print t
    for p in np.arange(4):
        def err(l):
            return np.sum(weight*(data-func(m+l,cutdown,T,t))**2)
        step = 1
        err1 = err(step)
        err0 = err(0)
        print err1,err0,'start'
        while step>1e-10:
            if(err1>=err0):
                step *= 0.5
            else :
                m+= step
                m = m *(1+0.001*(0.5-np.random.rand()))
                err0 = err1
            err1 = err(step)
            print 'err',err1
            print 'm',m+step
            #plot_cmp(func(m+step,cutdown,T,t),data)
            #os.system('pause')
    plot_cmp(func(m,cutdown,T,t),data)
    return (A,m)
def func(A,m,T,t):
    val = np.exp(-m*t)+np.exp(-m*(T-t))
    avr = np.mean(val)
    return A*(avr**-1)*val
def plot_cmp(A,m,data,cutdown,T):
    print A
    print m
    tsamples = np.arange(cutdown,T-cutdown)
    plt.plot(tsamples,func(A,m,T,tsamples),'r-')
    plt.plot(tsamples,data[cutdown:T-cutdown],'g+')
    plt.show()
    
def jackknife_scalar(data):
    num = len(data)
    fsum = np.zeros(shape = (num))
    bsum = np.zeros(shape = (num))
    jsum = np.zeros(shape = (num))
    for i in np.arange(1,num):
        fsum[      i] = (fsum[  i-1]+data[  i-1])
        bsum[num-1-i] = (bsum[num-i]+data[num-i])
    for i in np.arange(num):
        jsum[i] = (fsum[i]+bsum[i])/(num-1)
    avr = (bsum[0] + data[0])/num
    var = 0.0
    for i in np.arange(num):
        var+=(num-1)*((jsum[i]-avr)**2)*1.0/num
    return (avr,var,jsum)
def jackknife(data):
    num = len(data)
    fsum = np.zeros(shape = (num,len(data[0])))
    bsum = np.zeros(shape = (num,len(data[0])))
    jsum = np.zeros(shape = (num,len(data[0])))
    for i in np.arange(1,num):
        fsum[      i] = (fsum[  i-1]+data[  i-1])
        bsum[num-1-i] = (bsum[num-i]+data[num-i])
    for i in np.arange(num):
        jsum[i] = (fsum[i]+bsum[i])/(num-1)
    avr = (bsum[0] + data[0])/num
    var = np.zeros(len(data[0]))
    for i in np.arange(num):
        var+=(num-1)*((jsum[i]-avr)**2)*1.0/num
    return (avr,var,jsum)
def distances(T):
    d = np.arange(T)
    for i in np.arange(T):
        if(T-d[i]<d[i]):
            d[i] = T-d[i]
    return d
def check_data(propagators):
    data = {}
    index= {}
    count = 1
    cutdown = 0
    T = 0
    plt.yscale('log')
    pass_key = 0
    for keys in propagators:
        T = keys[0]
        if(keys!=(T,T,1,keys[3],0,0,keys[6],"gamma5",0)):
        
            continue
        if not data.has_key(keys[3]):
            data[keys[3]]=[]
            index[keys[3]] = propagators[keys][3]
        data[keys[3]].append(propagators[keys][7])
        count = propagators[keys][2]
    sum_num = 0
    sum = np.zeros(T)
    lines = []
    legend= []
    for keys in data:
        num = len(data[keys])
        (data_avr,data_var,jsum) = jackknife(data[keys])
        sum+= data_avr
        sum_num+=1
        l = plt.errorbar(np.arange(cutdown,T-cutdown),np.abs(data_avr),yerr=np.sqrt(data_var),capsize=2)
        lines.append(l)
        legend.append(r'$m_{r} = '+str(0.01*int(100*np.log(index[keys])))+'$')
        
        pass_avr = data_avr
        pass_var = data_var
        pass_count = count
    
    #print lines,legend
    plt.xlabel('t')
    plt.ylabel('C(t)')
    plt.legend(lines,legend,loc = 'best')
    plt.show()
    return (data_avr,data_var,count,T)
def check_cputime(propagators_array):
    time_avr = []
    time_err = []
    time_index = []
    T = 0
    plt.yscale('linear')
    for propagators in propagators_array:
        reweight = 1.0
        data = {}
        index= {}
        for keys in propagators:
            T = keys[0]
            if keys!=(T,T,1,keys[3],0,0,keys[6],"gamma5",0):
                continue
            if(keys[6]<5):
                continue
            if not data.has_key(keys[3]):
                data[keys[3]] = []
                index[keys[3]] = propagators[keys][3]
            data[keys[3]].append(propagators[keys][0]+0.01*propagators[keys][1])
        for keys in data:
            (data_avr,data_var,jsum) = jackknife_scalar(data[keys])
            time_index.append(0.01*int(100*np.log(index[keys])))
            time_avr.append(data_avr)
            time_err.append(np.sqrt(data_var))
    lines = []
    print time_index
    for k in np.arange(len(time_index)):
        l = plt.errorbar(time_index,time_avr,yerr=time_err,capsize=2)
        lines.append(l)
    plt.legend(lines,['Z=T=64,'])
    plt.xlabel(r'$m_r$')
    plt.ylabel(r'$CPU time/s$')
    plt.show()
    return time_index
(data_avr,data_var,count,T) = check_data(propagators5)
check_cputime([propagators5])
cutdown = 2
rates = []
for i in np.arange(len(data_avr)-3):
    try:
       rates.append(np.abs(0.5*np.log(np.abs(data_avr[i+1]*data_avr[i+3]/data_avr[i+2]/data_avr[i]))))
    except:
        break
print rates
print np.mean(rates)
if False:
    (A1,m1) = fit_propagator(data_avr,data_var,cutdown,count,T)
    print '----'
    plt.yscale('log')
    print A1,m1
    plt.plot(np.arange(T),data_avr/np.mean(data_avr),'+')
    plt.plot(np.arange(T),func(np.mean(data_avr),m1,T,np.arange(T)),'r-')
    #plt.plot(np.arange(cutdown,T-cutdown),func(A2,m2,T,np.arange(cutdown,T-cutdown)),'r-')
    #plt.plot(np.arange(cutdown,T-cutdown),func(A3,m3,T,np.arange(cutdown,T-cutdown)),'r-')
    plt.errorbar(np.arange(T),data_avr/np.mean(data_avr),yerr=np.sqrt(data_var)/np.mean(data_avr),capsize=2)
    #plt.plot(np.arange(cutdown,T-cutdown),jsum[3][cutdown:T-cutdown],'g+')
    #plt.plot(np.arange(cutdown,T-cutdown),jsum[5][cutdown:T-cutdown],'g+')
    plt.show()
    