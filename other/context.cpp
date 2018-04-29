#include "../common.h"
static const int X = 3;
static const int Y = 3;
static const int Xh= 1;
static const int Xl= 2;
movement<X,4,8,2,unsigned short> movex;
movement<Y,4,1,4,unsigned short> movey;
int loadContext(const int loop[][Y],const int & U,const int & u,const int & v)
    {
        int t[3];
        t[1] = v;
        if(v==0)t[0]=Y-1;
        else t[0] = v-1;
        if(v==Y-1)t[2] = 0;
        else t[2] = v+1;
        int result = 0;
        int shift = u-1;
        if(U==Xh-1&&u==Xl&&Xl==0)
            for(int j=0;j<3;j++)
                result  |=(((loop[U-1][t[j]]>>31)&1)<<(  3*j))
                        | (((loop[U  ][t[j]]    )&1)<<(1+3*j))
                        | (( loop[0  ][t[j]]     &1)<<(2+3*j));
        else if(U==Xh-1&&u==Xl&&Xl!=0)
            for(int j=0;j<3;j++)
                result  |=(((loop[U  ][t[j]]>>((Xl+31)&31))&3)<<(  3*j))
                        | (( loop[0  ][t[j]]     &1)<<(2+3*j));
        else if(U==0&&u==0)
            for(int j=0;j<3;j++)
                result  |=(((loop[Xh-1][t[j]]>>Xl)&1)<<(  3*j))
                        | (( loop[U  ][t[j]]     &3)<<(1+3*j));
        else if(u==31)
            for(int j=0;j<3;j++)
                result  |=(((loop[U  ][t[j]]>>30)&3)<<(  3*j))
                        | (( loop[U+1][t[j]]     &1)<<(2+3*j));
        else if(u==0)
            for(int j=0;j<3;j++)
                result  |=(((loop[U-1][t[j]]>>31)&1)<<(  3*j))
                        | (( loop[U  ][t[j]]     &3)<<(1+3*j));
        else
            for(int j=0;j<3;j++)
                result  |=((loop[U][t[j]]>>shift)&7)<<(  3*j);
        return result;
    }
int main()
{
    int loop[2][(X+31)>>5][Y] = {};
    unsigned char link[2][X][Y] = {};
    int nn[4][X][Y] = {};
    int n[X][Y] = {};
    int m[X][Y] = {};
    loop[0][0][0] = 3;
    loop[1][0][0] = 1;
    loop[1][0][1] = 1;
    loop2link<X,Y,1>(loop,link);
    printContexts<true>(1);
    int o[] = LOOP_OCCUPATION_ARRAY;
    printBits<1,4,true>(1);
    printBits<1,4,true>(o[1]);
    printBits<1,4,true>(o[4]);
    //printBits<1,4,true>(o[8]);
    //printBits<1,4,true>(o[32]);
    printBits<1,4,true>(o[64]);
    printBits<1,4,true>(o[256]);
    unsigned char nchange[16][16]; // nchange[old][new];
    for(int i=0;i<16;i++)
    {
        for(int j=0;j<16;j++)
        {
            int inc = j&(~i);
            int dec = i&(~j);
            int info = 0;
            int base = 1;
            int id = 0;
            for(int k=0;k<4;k++,base*=3)
            {
                if(inc&(1<<k))info |= 2<<(2*k);
                else if(dec&(1<<k))info |= 0;
                else info |= 1<<(2*k);
                if(inc&(1<<k))id += 2*base;
                else if(dec&(1<<k)) id += 0;
                else id += base;
            }
            nchange[i][j] = info;
        }
    }
    int U=0,u=0,v=0;
    int c1 = loadContext(loop[0],U,u,v);
    U=0,u=0,v=0;
    int c2 = loadContext(loop[1],U,u,v);
    printBits<1,4,true>(o[c1]);
    printBits<1,4,true>(o[c2]);
    printBits<1,8,true>(nchange[o[c1]][o[c2]]);
    printf("occupation:\n");
    for(int i=0;i<X;i++)
    {
        int buffer[2][2*Y];
        for(int j=0;j<Y;j++)
        {
            int U = 0;
            int u = i;
            int v = j;
            int ctx = loadContext(loop[0],U,u,v);
            buffer[0][j*2  ] = o[ctx]&1;
            buffer[0][j*2+1] = o[ctx]&4;
            buffer[1][j*2  ] = o[ctx]&2;
            buffer[1][j*2+1] = o[ctx]&8;
        }
        for(int k=0;k<2;k++)
        {
            for(int j=0;j<Y;j++)
                fprintf(stderr,"%d%d",buffer[k][j*2],buffer[k][j*2+1]);
            printf("\n");
        }
    }
    for(int i=0;i<X;i++)
    {
        int buffer[2][2*Y];
        for(int j=0;j<Y;j++)
        {
            int U = 0;
            int u = i;
            int v = j;
            int ctx = loadContext(loop[1],U,u,v);
            buffer[0][j*2  ] = o[ctx]&1;
            buffer[0][j*2+1] = o[ctx]&4;
            buffer[1][j*2  ] = o[ctx]&2;
            buffer[1][j*2+1] = o[ctx]&8;
        }
        for(int k=0;k<2;k++)
        {
            for(int j=0;j<Y;j++)
                fprintf(stderr,"%d%d",buffer[k][j*2],buffer[k][j*2+1]);
            printf("\n");
        }
    }

    printf("occupation compute\n");
    for(int i=0;i<X;i++)
    {
        for(int j=0;j<Y;j++)
        {
            int U = 0;
            int u = i;
            int v = j;
            int ctx0 = loadContext(loop[0],U,u,v);
            int ctx1 = loadContext(loop[1],U,u,v);
            int ctx;
            ctx = ctx0;
            if(o[ctx]&1)nn[0][movex(i,1)][movey(j,2)]++;
            if(o[ctx]&2)nn[1][movex(i,1)][j]++;
            if(o[ctx]&4)nn[2][i][movey(j,2)]++;
            if(o[ctx]&8)nn[3][i][j]++;
            ctx = ctx1;
            if(o[ctx]&1)nn[0][movex(i,1)][movey(j,2)]+=2;
            if(o[ctx]&2)nn[1][movex(i,1)][j]+=2;
            if(o[ctx]&4)nn[2][i][movey(j,2)]+=2;
            if(o[ctx]&8)nn[3][i][j]+=2;
        }
    }
    for(int i=0;i<X;i++)
    {
        for(int j=0;j<Y;j++)
        {
            m[i][j] += (link[0][i][j]!=0);
            m[i][j] += (link[1][i][j]!=0);
        }
    }
    printf("n using id\n");
    for(int k=0;k<4;k++)
    {
        for(int i=0;i<X;i++)
        {
            for(int j=0;j<Y;j++)
            {
                printf("%d ",nn[0][i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }
    printf("n using link\n");
    for(int i=0;i<X;i++)
    {
        for(int j=0;j<Y;j++)
        {
            printf("%d ",m[i][j]);
        }
        printf("\n");
    }
    return 0;
}
