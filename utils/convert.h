#ifndef UTILS_CONVERT_H
#define UTILS_CONVERT_H
template<int X,int Y,int N>
void loop2link(int loop[][(X+31)>>5][Y],unsigned char link[][X][Y])
{
    static int spin2link[16] = {0,6,3,5,12,10,15,9,9,15,10,12,5,3,6,0};
    static const int Xh = (X+31)>>5;
    for(int f=0;f<2*N;f++)
    {
        for(int i=0;i<X;i++)
        {
            int u;
            int v;
            for(int j=0;j<Y;j++)
            {
                unsigned char info = 0;
                u = i;v = j;
                if((loop[f][u>>5][v]>>(u&31))&1)info|=1;
                u = i+1;v = j;
                if(u>=X)u-=X;
                if((loop[f][u>>5][v]>>(u&31))&1)info|=4;
                u = i+1;v = j+1;
                if(u>=X)u-=X;
                if(v>=Y)v-=Y;
                if((loop[f][u>>5][v]>>(u&31))&1)info|=8;
                u = i;v = j+1;
                if(v>=Y)v-=Y;
                if((loop[f][u>>5][v]>>(u&31))&1)info|=2;
                link[f][i][j] = spin2link[info];
            }
        }
    }
}
template<int X,int Y,int N>
void loop2link(unsigned long long bc,int loop[][(X+31)>>5][Y],unsigned char link[][X][Y])
{
    static int spin2link[16] = {0,6,3,5,12,10,15,9,9,15,10,12,5,3,6,0};
    static const int Xh = (X+31)>>5;
    static const int XH = (X+32)>>5;
    int _loop[2*N][XH][Y+1] = {};
    for(int f=0;f<2*N;f++)
    {
        for(int i=0;i<Xh;i++)
        {
            for(int j=0;j<Y;j++)
            {
                _loop[f][i][j] = loop[f][i][j];
            }
        }
        for(int j=0;j<Y;j++)
        {
            if(((bc&(1<<(f<<1)))!=0)^(loop[f][0][j]&1))
                _loop[f][X>>5][j]|=1<<(X&31);
        }
        for(int i=0;i<XH;i++)
        {
            int n = X+1 - i*32;
            unsigned int valid = (n>=32)?0xffffffff:((1<<n)-1);
            if((bc&(2<<(f<<1)))!=0)_loop[f][i][Y]=(~loop[f][i][Y])&valid;
            else _loop[f][i][Y]=loop[f][i][Y]&valid;
        }
    }
    for(int f=0;f<2*N;f++)
    {
        for(int i=0;i<X;i++)
        {
            int u;
            int v;
            for(int j=0;j<Y;j++)
            {
                unsigned char info = 0;
                u = i;v = j;
                if((loop[f][u>>5][v]>>(u&31))&1)info|=1;
                u = i+1;v = j;
                if((loop[f][u>>5][v]>>(u&31))&1)info|=4;
                u = i+1;v = j+1;
                if((loop[f][u>>5][v]>>(u&31))&1)info|=8;
                u = i;v = j+1;
                if((loop[f][u>>5][v]>>(u&31))&1)info|=2;
                link[f][i][j] = spin2link[info];
            }
        }
    }
}
#endif // UTILS_CONVERT_H
