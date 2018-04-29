#ifndef UTILS_STATISTIC_H
#define UTILS_STATISTIC_H
#include<cmath>
template<typename T>
void ndarray_zero(T &t)
{
    t = 0;
}

constexpr unsigned int log2_upper(unsigned int x,int bound)
{
    return (x>(1<<(bound-1)))?bound:log2_upper(x,bound-1);
}
constexpr unsigned int log2_upper(unsigned int x)
{
    return
        (x>(1<<30))?31:
        (x==0)?0:log2_upper(x,30);
}
static_assert(log2_upper(1)==0,"WA");
static_assert(log2_upper(2)==1,"WA");
static_assert(log2_upper(3)==2,"WA");
template<typename ... T>
void ndarray_zero(auto t,int n,T ... other)
{
    for(int i=0;i<n;i++)
    {
        ndarray_zero(t[i],other...);
    }
}
void tensorprint(auto t)
{
    using std::cout;
    cout<<t;
}
template<typename ... T>
void tensorprint(auto t,int length,T ... other)
{
    using std::cout;
    cout<<"{";tensorprint(t[0],other...);
    for(int i=1;i<length;i++)
    {
        cout<<",";
        tensorprint(t[i],other...);
    }
    cout<<"}\n";
}
double exp(double a,unsigned int n)
{
    if(n==1)return a;
    else if(n==0)return 1.0;
    double temp = exp(a,n>>1);
    if(n&1)
        return a*temp*temp;
    else
        return temp*temp;
}
double exp(double a,int n)
{
    if(n<0)
        return exp(1.0/a,(unsigned int)-n);
    else
        return exp(a,(unsigned int)n);
}
double Sigma(double u[],int n)
{
    double s = 0.0;
    for(int i=0;i<n;i++)
    {
        s +=u[i];
    }
    return s;

}
double Var(double u[],double v,int n)
{
    double var = 0.0;
    for(int i=0;i<n;i++)
    {
        var+=(u[i]-v)*(u[i]-v);
    }
    return var;
}
double Var(double u[],double v[],int n)
{
    double var = 0.0;
    for(int i=0;i<n;i++)
    {
        var+=(u[i]-v[i])*(u[i]-v[i]);
    }
    return var;
}
template<int x>
inline double Cosh(double m,int i)
{
    return exp(-m*i)+exp(-m*(x-i));
}
template<int X>
inline double normalizedCosh(double m,int i)
{
    if(m*X<0.01)return (exp(-m*i)+exp(-m*(X-i)))*(m+0.5*m*m)/(2*m*X+m*m*X*X);
    else return (exp(-m*i)+exp(-m*(X-i)))*(1-exp(-m))/(1-exp(-m*X));
}
template<int x,int y,int z,bool ta,bool tb>
void matrixMult(double c[x][z],double a[ta?y:x][ta?x:y],double b[tb?z:y][tb?y:z])
{
    for(int i=0;i<x;i++)
    {
        for(int j=0;j<z;j++)
        {
            c[i][j] = 0.0;
            for(int k=0;k<y;k++)
            {
                     if(!ta&&!tb)c[i][j]+=a[i][k]*b[k][j];
                else if( ta&&!tb)c[i][j]+=a[k][i]*b[k][j];
                else if( ta&& tb)c[i][j]+=a[k][i]*b[j][k];
                else if(!ta&& tb)c[i][j]+=a[i][k]*b[j][k];
            }
        }
    }
}
template<typename T1>
T1 max(T1 t1)
{
    return t1;
}
template<typename T1>
T1 min(T1 t1)
{
    return t1;
}
template<typename T1,typename ... T>
T1 max(T1 t1,T ... other)
{
    auto t2 = max(other...);
    return t1<t2?t2:t1;
}
template<typename T1,typename ... T>
T1 min(T1 t1,T ... other)
{
    auto t2 = min(other...);
    return t1<t2?t1:t2;
}
void fittingShape(int L,int cutdown,double l,double r,double &p,double &error,auto f,auto g)
{
    double G[L];
    double F[L];
    for(int i=cutdown;i<L;i++)
        F[i-cutdown] = f(i);
    double sumF = Sigma(F,L-cutdown);
    for(int i=0;i<L-cutdown;i++)
        F[i] = F[i]/sumF;
    double sumG = 0.0;
    double phi = (sqrt(5.)-1.)/2;
    double m1 = l + (r-l)*phi*phi;
    double m2 = l + (r-l)*phi;
    double m;
    double vl;
    double vm1;
    double vm2;
    double vr;
    double vm;
    #define PRODUCE(p) \
    for(int i=cutdown;i<L;i++)\
        G[i-cutdown] = g(p,i); \
    sumG = Sigma(G,L-cutdown); \
    for(int i=0;i<L-cutdown;i++)\
        G[i] = G[i]/sumG;
    PRODUCE(l)
    vl = Var(F,G,L-cutdown);
    PRODUCE(m1)
    vm1 = Var(F,G,L-cutdown);
    PRODUCE(m2)
    vm2 = Var(F,G,L-cutdown);
    PRODUCE(r)
    vr = Var(F,G,L-cutdown);
    int step = 0;
    while(r-l>1e-6)
    {
        step++;
        if(step>100000)
        {
            std::cout<<"Weird hound"<<std::endl;
            break;
        }
        if(vl<=vm1&&vl<=vm2||vr<=vm1&&vr<=vm2)
        {
            if(vl<vr)
            {
                r=m2;
                vr=vm2;
                m2=m1;
                vm2=vm1;
                m1 = l + (m1-l)*phi;
                PRODUCE(m1)
                vm1 = Var(F,G,L-cutdown);
            }
            else
            {
                l  = m1;
                vl = vm1;
                m1 = m2;
                vm1 = vm2;
                m2 = r;
                vm2 = vr;
                r+=1;
                PRODUCE(r)
                vr = Var(F,G,L-cutdown);
                if(r>10)
                {
                    break;
                }
            }
        }
        else
        {
            if(vm1<vm2)
            {
                m = l + (m1-l)*phi;
                r  = m2;
                vr =vm2;
                m2 = m1;
                vm2=vm1;
                m1 = m;
                PRODUCE(m1)
                vm1= Var(F,G,L-cutdown);
            }
            else
            {
                m = r + (m2-r)*phi;
                l  = m1;
                vl =vm1;
                m1 = m2;
                vm1=vm2;
                m2 = m ;
                PRODUCE(m2);
                vm2 = Var(F,G,L-cutdown);
            }
        }
    }
    p = (r+l)/2;
    PRODUCE(p);
    if(0)
    {
        fprintf(stderr,"--get fitting--\n");
        for(int i=0;i<L-cutdown;i++)
        {
            fprintf(stderr,"%d:%f %f\n",i,F[i],G[i]);
        }
        system("pause");
    }
    error = Var(F,G,L-cutdown);
    return ;
}
#endif // UTILS_STATISTIC_H
