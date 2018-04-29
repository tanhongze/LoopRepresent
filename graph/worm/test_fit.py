from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import math
import os

T = 100
count = 100000
def func(A,m,T,t):
    val = 0.0
    assert len(A) == len(m)
    for i in np.arange(len(A)):
        val += A[i]*np.exp(-m[i]*t)+A[i]*np.exp(-m[i]*(T-t))
    return val
def func_gradA(A,m,T,t,i):
    return np.exp(-m[i]*t)+np.exp(-m[i]*(T-t))
def func_gradm(A,m,T,t,i):
    return A[i]*(-t*np.exp(-m[i]*t)-(T-t)*np.exp(-m[i]*(T-t)))

arg1A = [1.0,10.0]
arg1m = [1.0,1.2]
def plot_cmp(A,m,data):
    tsamples = np.arange(T)
    plt.plot(tsamples,func(arg1A,arg1m,T,tsamples),'b-')
    plt.plot(tsamples,func(A,m,T,tsamples),'r-')
    plt.plot(tsamples,data,'g+')
    plt.show()
def err_hat(count,r,t):
    return (1.0/count)*np.exp(-r*np.abs(t-T*0.5))
def fit_propagator(input_data,input_err2,cutdown,count,T):
    r = 0.0
    t = np.arange(cutdown,T-cutdown)
    data = np.abs(input_data[cutdown:T-cutdown])
    rate = np.max(data)
    data/= np.max(data)   
    err2 = input_err2[cutdown:T-cutdown] + err_hat(count,r,t)
    err2/= np.min(err2)
    print data
    print err2
    weight = data/err2
    A_retval = []
    m_retval = []
    k = 0
    K = 1
    def plot_cmp(curve,data,cutdown,T):
        tsamples = np.arange(cutdown,T-cutdown)
        plt.plot(tsamples,curve,'r-')
        plt.plot(tsamples,data,'g+')
        plt.show()
    def func(A,m,cutdown,T,t):
        val = 0.0
        assert len(A) == len(m)
        for i in np.arange(len(A)):
            val += A[i]*np.exp(-m[i]*(t-cutdown))+A[i]*np.exp(-m[i]*(T-cutdown-t))
        return val
    def func_gradA(A,m,T,t,i):
        return 0.0
    def func_gradm(A,m,T,t,i):
        return A[i]*(-(t-cutdown)*np.exp(-m[i]*(t-cutdown))-(T-cutdown-t)*np.exp(-m[i]*(T-cutdown-t)))
    for k in np.arange(1):
        Atemp = []
        mtemp = []
        for x in A_retval:
            Atemp.append(x)
        for x in m_retval:
            mtemp.append(x)
        if len(m_retval)>0:
            Atemp.append(data[0]-func(Atemp,mtemp,T,cutdown))
        else:
            Atemp.append(data[0])
        if len(m_retval)>0:
            mtemp.append(np.min(m_retval)*0.5)
        else:
            mtemp.append(0.1)
        A = np.array(Atemp)
        m = np.array(mtemp)
        #plot_cmp(func(A,m,cutdown,T,t),data,cutdown,T)
        heat = 100
        for p in np.arange(4):
            gradA = np.zeros(k+1)
            gradm = np.zeros(k+1)
            for q in np.arange(k+1):
                gradA[q] = np.sum((-2.0/err2)*func_gradA(A,m,T,t,q)*(data-func(A,m,cutdown,T,t)))
                gradm[q] = np.sum((-2.0/err2)*func_gradm(A,m,T,t,q)*(data-func(A,m,cutdown,T,t)))
            normgrad = np.sqrt(np.sum(gradA**2)+np.sum(gradm**2))
            
            if(normgrad>1e-12):
                gradA /= normgrad
                gradm /= normgrad
            def derr(l,gradA,gram):
                val = 0.0
                assert len(A) == len(m)
                for i in np.arange(len(A)):
                    val += np.sum((-2.0/err2)*gradA[q]*func_gradA(A+l*gradA,m+l*gradm,T,t,q)*(data-func(A+l*gradA,m+l*gradm,T,t))
                +(-2.0/err2)*gradm[q]*func_gradm(A+l*gradA,m+l*gradm,T,t,q)*(data-func(A+l*gradA,m+l*gradm,T,t)))
                return val
            def err(l,gradA,gram):
                val = func(A+l*gradA,m+l*gradm,cutdown,T,t)
                if(False):
                    tsamples = np.arange(cutdown,T-cutdown)
                    plt.plot(tsamples,func(A,m,cutdown,T,tsamples),'b-')
                    plt.plot(tsamples,data,'r+')
                    plt.show()
                return np.sum(weight*(data-func(A+l*gradA,m+l*gradm,cutdown,T,t))**2)
            step = -2
            l = 0
            err1 = err(l+step,gradA,gradm)
            err0 = err(l,gradA,gradm)
            print err1,err0,'start'
            print gradA,gradm
            while step<-1e-10:
                if(err1>=err0):
                    print err1,'>=',err0,'denied'
                    step *= 0.5
                    if False:
                        print "-- update --"
                        print A,A+l*gradA
                        print m,m+l*gradm
                        print gradA
                        print gradm
                    A += l*gradA
                    m += l*gradm
                    A = A *(1+(err0*0.1+0.001)*(0.5-np.random.rand(len(A))))
                    A = A *(1+(err0*0.1+0.001)*(0.5-np.random.rand(len(m))))
                    l = 0
                    gradA = np.zeros(k+1)
                    gradm = np.zeros(k+1)
                    for q in np.arange(k+1):
                        gradA[q] = np.sum((-2.0/err2)*func_gradA(A,m,T,t,q)*(data-func(A,m,cutdown,T,t)))
                        gradm[q] = np.sum((-2.0/err2)*func_gradm(A,m,T,t,q)*(data-func(A,m,cutdown,T,t)))
                    
                    normgrad = np.sqrt(np.sum(gradA**2)+np.sum(gradm**2))            
                    if(normgrad>1e-12):
                        gradA /= normgrad
                        gradm /= normgrad
                    #plot_cmp(func(A,m,cutdown,T,t),data,cutdown,T)
                else :
                    l += step
                    err0 = err1
                    print err1,'<=',err0,'accept'
                err1 = err(l+step,gradA,gradm)
                #os.system('pause')
            A = A + l*gradA
            m = m + l*gradm
            A = A *(1+(err0*0.1+0.001)*(0.5-np.random.rand(len(A))))
            A = A *(1+(err0*0.1+0.001)*(0.5-np.random.rand(len(m))))
            #plot_cmp(A,m,input_data)
        #submit
        for i in np.arange(len(A)-1):
            A_retval[i]=A[i]
            m_retval[i]=m[i]
        A_retval.append(A[len(A)-1])
        m_retval.append(m[len(m)-1])
    A_retval = np.array(A_retval)*rate
    m_retval = np.array(m_retval)
    #plot_cmp(func(A_retval,m_retval,cutdown,T,t),data,cutdown,T)
    return (A_retval,m_retval)
def fit_cosh(input_data,input_err2,count,r,T,K,func,func_gradA,func_gradm):
    t = np.arange(2,T)
    data = input_data[2:] 
    err2 = input_err2[2:] + err_hat(count,r,t)
    err2/= np.min(err2)
    print err2
    A_retval = []
    m_retval = []
    k = 0
    for k in np.arange(K):
        Atemp = []
        mtemp = []
        for x in A_retval:
            Atemp.append(x)
        for x in m_retval:
            mtemp.append(x)
        if len(m_retval)>0:
            Atemp.append(data[0]-func(Atemp,mtemp,T,2))
        else:
            Atemp.append(data[0])
        if len(m_retval)>0:
            mtemp.append(np.min(m_retval)*0.5)
        else:
            mtemp.append(0.1)
        A = np.array(Atemp)
        m = np.array(mtemp)
        #plot_cmp(A,m,input_data)
        print '----- k = ',k,' start -----'
        print 'A',A
        print 'm',m
        heat = 100
        for p in np.arange(2000):
            gradA = np.zeros(k+1)
            gradm = np.zeros(k+1)
            for q in np.arange(k+1):
                gradA[q] = np.sum((-2.0/err2)*func_gradA(A,m,T,t,q)*(data-func(A,m,T,t)))
                gradm[q] = np.sum((-2.0/err2)*func_gradm(A,m,T,t,q)*(data-func(A,m,T,t)))
            normgrad = np.sqrt(np.sum(gradA**2)+np.sum(gradm**2))
            
            if(normgrad>1):
                gradA /= normgrad
                gradm /= normgrad
            def derr(l):
                val = 0.0
                assert len(A) == len(m)
                for i in np.arange(len(A)):
                    val += np.sum((-2.0/err2)*gradA[q]*func_gradA(A+l*gradA,m+l*gradm,T,t,q)*(data-func(A+l*gradA,m+l*gradm,T,t))
                +(-2.0/err2)*gradm[q]*func_gradm(A+l*gradA,m+l*gradm,T,t,q)*(data-func(A+l*gradA,m+l*gradm,T,t)))
                return val
            def err(l):
                return np.sum(((data-func(A+l*gradA,m+l*gradm,T,t))**2))
            step = -2
            l = 0
            err1 = err(l+step)
            err0 = err(l)
            while step<-0.0001:
                if(err1+0.000001>=err0):
                    step *= 0.5
                else :
                    l += step
                    err0 = err1
                err1 = err(l+step)
            print 'step by',l
            print 'err',err0
            A = A + l*gradA
            A = A *(1+(err0*0.1+0.001)*(0.5-np.random.rand(len(A))))
            m = m + l*gradm
            A = A *(1+(err0*0.1+0.001)*(0.5-np.random.rand(len(m))))
            print 'update argument',(A,m)
            #plot_cmp(A,m,input_data)
        #submit
        for i in np.arange(len(A)-1):
            A_retval[i]=A[i]
            m_retval[i]=m[i]
        A_retval.append(A[len(A)-1])
        m_retval.append(m[len(m)-1])
        print 'A',A_retval
        print 'm',m_retval
        print '----- k = ',k,' end -----'
    plot_cmp(A_retval,m_retval,input_data)
    return (A_retval,m_retval)
if __name__ == '__main__':
    arg2A = [1.0]
    arg2m = [1.0]
    r = 1.0
    tsamples = np.arange(T)
    psamples = func(arg1A,arg1m,T,tsamples)
    randval = np.random.rand(T)-0.5
    #print randval
    p_exp = []
    for i in np.arange(10):
        p_exp.append(np.array(psamples*count+np.sqrt(psamples*count)*0*(np.random.rand(T)-0.5),dtype = 'int'))
    p_sum = np.zeros(T)
    p_avr = np.zeros(T)
    p_var = np.zeros(T)
    for i in np.arange(10):
        p_sum += p_exp[i]
    p_avr = p_sum / 10.0
    for i in np.arange(10):
        p_var += (p_exp[i]-p_avr)**2
    p_var = p_var / 9.0
    p_avr/= count
    p_var/= count
    #print p_avr
    #print p_var
    (A_fit,m_fit) = fit_cosh(p_avr,p_var,count,0.0,T,2,func,func_gradA,func_gradm)
    print 'fit_result',A_fit,m_fit
    #plt.plot(tsamples,np.log(func(arg1A,arg1m,T,tsamples)))
    #plt.plot(tsamples,np.log(func(arg2A,arg2m,T,tsamples)))
    plt.plot(tsamples,func(arg1A,arg1m,T,tsamples),'b-')
    plt.plot(tsamples,func(A_fit,m_fit,T,tsamples),'r+')
    plt.show()