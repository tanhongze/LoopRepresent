#include<cstdio>
number c[10];
unsigned long long bond_count(unsigned long long i,unsigned long long j)
{
    typedef unsigned long long number;
    if(j==0)return 1;
    number a = i*(i-1);
    for(number k=1;k<j;k++)
    {
        a /= 2*k;
        a *= (i-2*k)*(i-2*k-1);
    }
    a /= 2*j;
    return a;
}
int main()
{
    fprintf(stderr,"start\n");
    for(number i=0;i<30;i++)
    {
        fprintf(stderr,"%lld:",i);
        for(number j=0;j<=i/2;j++)
        {
            number x = count(i,j);
            fprintf(stderr,"%llu ",x);
        }
        fprintf(stderr,"\n");
    }
    return 0;
}
