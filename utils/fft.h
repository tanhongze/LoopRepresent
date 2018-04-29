#ifndef FFT_H
#define FFT_H
#define r(ar,ai,br,bi) (ar*br-ai*bi)
#define i(ar,ai,br,bi) (ar*bi+ai*br)
void fft(double xr[],double xi[],double ar[],double ai[],int logn,bool inv)
{
    int n = 1<<logn;
    double pi = 4*atan(1);
    int rev[n];
    double ur;
    double ui;
    rev[0] = 0;
    rev[1] = 1<<(logn-1);
    for(int k=2;k<=logn;k++)
    {
        int upper = 1<<k;
        int mask = (1<<(k-1))-1;
        for(int j=(1<<(k-1));j<upper;j++)
        {
            rev[j] = rev[j&mask]|(1<<(logn-k));
        }
    }
    for(int k=0;k<n;k++)
    {
        ar[rev[k]] = xr[k];
        ai[rev[k]] = xi[k];
    }
    for(int s = 1;s<=logn;s++)
    {
        int m = 1<<s;
        int mid = m>>1;
        ur =  cos((2*pi)/m);
        if(inv)ui =  sin((2*pi)/m);
        else   ui = -sin((2*pi)/m);
        for(int k = 0;k<n;k+=m)
        {
            double wr = 1.0;
            double wi = 0.0;
            for(int j = 0;j<mid;j++)
            {
                double zr = r(ar[mid+k+j],ai[mid+k+j],wr,wi);
                double zi = i(ar[mid+k+j],ai[mid+k+j],wr,wi);
                ar[mid+k+j] = ar[k+j]-zr;
                ai[mid+k+j] = ai[k+j]-zi;
                ar[    k+j] = ar[k+j]+zr;
                ai[    k+j] = ai[k+j]+zi;
                double rr = r(ur,ui,wr,wi);
                double ri = i(ur,ui,wr,wi);
                wr = rr;
                wi = ri;
            }
        }
    }
    return ;
}
void fft_real_even(double xr[],double ar[],int logn,bool inv)
{
    int n = 1<<logn;
    double pi = 4*atan(1);
    int rev[n];
    double ai[n];
    double ur;
    double ui;
    rev[0] = 0;
    rev[1] = 1<<(logn-1);
    for(int k=2;k<=logn;k++)
    {
        int upper = 1<<k;
        int mask = (1<<(k-1))-1;
        for(int j=(1<<(k-1));j<upper;j++)
        {
            rev[j] = rev[j&mask]|(1<<(logn-k));
        }
    }
    for(int k=0;k<n;k++)
    {
        ar[rev[k]] = xr[k];
        ai[rev[k]] = 0;
    }
    for(int s = 1;s<=logn;s++)
    {
        int m = 1<<s;
        int mid = m>>1;
        ur =  cos((2*pi)/m);
        if(inv)ui =  sin((2*pi)/m);
        else   ui = -sin((2*pi)/m);
        for(int k = 0;k<n;k+=m)
        {
            double wr = 1.0;
            double wi = 0.0;
            for(int j = 0;j<mid;j++)
            {
                double zr = r(ar[mid+k+j],ai[mid+k+j],wr,wi);
                double zi = i(ar[mid+k+j],ai[mid+k+j],wr,wi);
                ar[mid+k+j] = ar[k+j]-zr;
                ai[mid+k+j] = ai[k+j]-zi;
                ar[    k+j] = ar[k+j]+zr;
                ai[    k+j] = ai[k+j]+zi;
                double rr = r(ur,ui,wr,wi);
                double ri = i(ur,ui,wr,wi);
                wr = rr;
                wi = ri;
            }
        }
    }
    return ;
}
#undef r
#undef i
#endif //FFT_H
