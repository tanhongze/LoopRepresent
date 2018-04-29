#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cmath>
#include<omp.h>
#include<cassert>
#include<list>
#include "../common.h"
    template<int sizeX,int sizeY>
    class propagatorStateNumber
    {
      public:
        static const int X = (sizeX<sizeY)?sizeY:sizeX;
        static const int Y = (sizeX<sizeY)?sizeX:sizeY;
        static const int maxcorner = 2*X*Y;
        static const int maxpoints =   X*Y;
        static const int mask = (1<<Y)-1;
        static const int numindex = (1<<(Y-1));
        static const int numlink = (1<<Y);
        int search_ylink[X+1];
        int search_kinds[X+1];
        //static const bool disableFolding = false;
        int ones[1<<9];
        int evenodd[2][1<<7];
        int evenoddindex[1<<8];
        unsigned int shift[1<<(Y-1)];
        unsigned int singlelayer[2][2][maxpoints+1][1<<maxpoints];
        unsigned int graphnumber[2][2];
        unsigned int occupy[2][2][1<<(Y-1)][1<<(Y-1)];
        unsigned int kinds[2][numindex][numindex];
        unsigned int oddx[2][2][1<<(Y-1)][1<<(Y-1)];
        movement<X,4,8,2,unsigned char> movex;
        movement<Y,4,1,4,unsigned char> movey;
        std::list<unsigned int> paths[X][Y];
        std::list<unsigned int> ns[X][Y];
        #define CORNER(i,j) (ones[i^j])
        void initSingleLayer()
        {
            for(int ex=0;ex<2;ex++)
                for(int ey=0;ey<2;ey++)
                    for(int corner = 0;corner<=maxpoints;corner++)
                        for(int n1 = 0;n1<(1<<maxpoints);n1++)
                                singlelayer[ex][ey][corner][n1] = 0;
            for(int i=0;i<(1<<(Y-1));i++)
            {
                getsinglegraphnumber(0,i,i,X);
                getsinglegraphnumber(1,i,i,X);
            }
        }
        void initBitCount()
        {
            ones[0] = 0;
            ones[1] = 1;
            for(int i=2;i<(1<<9);i++)
            {
                ones[i] = ones[i>>1]+(i%2);
            }
            int evenoddnum[2] = {};
            for(int i=0;i<(1<<8);i++)
            {
                evenodd[ones[i]%2][evenoddnum[ones[i]%2]++] = i;
            }
            for(int i=0;i<(1<<8);i++)
            {
                evenoddindex[evenodd[i%2][i>>1]] = i>>1;
            }
        }
        void initBitShift()
        {
            for(int i=0;i<(1<<(Y-1));i++)
            {
                unsigned differ = evenodd[0][i];
                shift[i] = 0;
                int haslink = 0;
                for(int j=0;j<Y;j++)
                {
                    if(differ&(1<<j))haslink ^= 1;
                    if(haslink)shift[i] |= 1<<j;
                }
            }
        }
        void initBitOccupy()
        {
            for(int i=0;i<numindex;i++)
            {
                for(int j=0;j<numindex;j++)
                {
                    for(int e=0;e<2;e++)
                    {
                        unsigned int up   = evenodd[e][i];
                        unsigned int down = evenodd[e][j];
                        unsigned int differ = up^down;
                        unsigned int index = evenoddindex[differ];
                        unsigned int shifted[2] = {shift[index],shift[index]^mask};
                        unsigned int transvers[2] = {};
                        unsigned int occupying[2] = {};
                        transvers[0] = (shifted[0]|(shifted[0]<<1)|(shifted[0]>>(Y-1)))&mask;
                        transvers[1] = (shifted[1]|(shifted[1]<<1)|(shifted[1]>>(Y-1)))&mask;
                        occupying[0] = (up&down)|transvers[0];
                        occupying[1] = (up&down)|transvers[1];
                        if(up&down)
                        {
                            if(up&down&transvers[0])
                            {
                                if(up&down&transvers[1])
                                {
                                    kinds[e][i][j] = 0;
                                    occupy[0][e][i][j] = 0;
                                    occupy[1][e][i][j] = 0;
                                    oddx[0][e][i][j] = 0;
                                    oddx[1][e][i][j] = 0;
                                }
                                else
                                {
                                    kinds[e][i][j] = 1;
                                    occupy[0][e][i][j] = occupying[1];
                                    occupy[1][e][i][j] = occupying[1];
                                    oddx[0][e][i][j] = shifted[1]&1;
                                    oddx[1][e][i][j] = shifted[1]&1;
                                }
                            }
                            else
                            {
                                kinds[e][i][j] = 1;
                                occupy[0][e][i][j] = occupying[0];
                                occupy[1][e][i][j] = occupying[0];
                                oddx[0][e][i][j] = shifted[0]&1;
                                oddx[1][e][i][j] = shifted[0]&1;
                            }
                        }
                        else
                        {
                            kinds[e][i][j] = 2;
                            occupy[0][e][i][j] = occupying[0];
                            occupy[1][e][i][j] = occupying[1];
                            oddx[0][e][i][j] = shifted[0]&1;
                            oddx[1][e][i][j] = shifted[1]&1;
                        }
                    }
                }
            }
        }
        propagatorStateNumber()
        {
            initBitCount();
            initBitShift();
            initBitOccupy();
            initSingleLayer();
        }
        void getsinglegraphnumber(int ey,int i,int j,int depth)
        {
            if(depth<=1)
            {
                if(kinds[ey][i][j]==0)return;
                int eo;
                int corner;
                int n1;
                search_ylink[depth] = j;
                int num = 0;
                for(int s=0;s<kinds[ey][i][j];s++)
                {
                    eo = 0;
                    corner = 0;
                    n1 = 0;
                    for(int k=X;k>1;k--)
                    {
                        eo ^= oddx[search_kinds[k]][ey][search_ylink[k]][search_ylink[k-1]];
                        corner += CORNER(evenodd[ey][search_ylink[k]],evenodd[ey][search_ylink[k-1]]);
                        n1 = (n1<<Y)|occupy[search_kinds[k]][ey][search_ylink[k]][search_ylink[k-1]];
                    }
                    eo ^= oddx[s][ey][j][i];
                    corner += CORNER(evenodd[ey][j],evenodd[ey][i]);
                    n1 = (n1<<Y)|occupy[s][ey][j][i];
                    singlelayer[eo][ey][corner][n1] ++;
                }
                return;
            }
            else
            {
                search_ylink[depth] = j;
                for(int k=0;k<(1<<(Y-1));k++)
                {
                    if(kinds[ey][j][k]==0)continue;
                    for(int s = 0;s<kinds[ey][j][k];s++)
                    {
                        search_kinds[depth] = s;
                        getsinglegraphnumber(ey,i,k,depth-1);
                    }
                }
                return;
            }
        }
        void print_path(unsigned int path,int n)
        {
            unsigned char link[1][X][Y] = {};
            int _d = 10;
            int x0 = 0;
            int y0 = 0;
            int num = ones[n]-1;
            for(int i=0;i<num;i++)
            {
                fprintf(stderr,"(%d,%d)-%d>",x0,y0,((path>>((num-i-1)*2))&3));
                link[0][x0][y0] = (1<<((path>>((num-i-1)*2))&3))|(1<<_d);
                _d = ((path>>(2*(num-i-1)))&3)^2;
                x0 = movex(x0,(path>>(2*(num-i-1)))&3);
                y0 = movey(y0,(path>>(2*(num-i-1)))&3);
            }
            link[0][x0][y0] = 1<<_d;
            fprintf(stderr,"(%d,%d)\n",x0,y0);
            printLink(link,0,x,y,0,0);

        }
        #define OCCUPY_BIT(x,y) (1<<(x*Y+y))
        void worm_path(unsigned char u,unsigned char v,unsigned char x,unsigned char y,unsigned char l,unsigned char r,unsigned char d,int n,unsigned int path)
        {
            if(u==x&&v==y)
            {
                if(1)
                {
                    printf("reach (%d,%d)\n",x,y);
                    print_path(path,n);
                    system("pause");
                }
                paths[x][y].push_back(path);
                ns[x][y].push_back(n);
                return;
            }
            else
            {
                if(1)
                {
                    printBits<3,3,true>(n);
                    printf("(%d,%d)--->\n",u,v);
                }
                for(int _d=0;_d<4;_d++)
                {
                    if(_d==2)continue;
                    unsigned char nextx = movex(u,(d+_d)&3);
                    unsigned char nexty = movey(v,(d+_d)&3);
                    if(n&OCCUPY_BIT(nextx,nexty))continue;
                    worm_path(nextx,nexty,x,y,l+(_d==1),r+(_d==3),(d+_d)&3,n|OCCUPY_BIT(nextx,nexty),(path<<2)|((d+_d)&3));
                }
            }
        }
        #define OCCUPY_BIT(x,y) (1<<(x*Y+y))
        void worm_goto(unsigned char u,unsigned char v,unsigned char x,unsigned char y,unsigned char l,unsigned char r,unsigned char d,int n,int count[1<<(X*Y)][X*Y][X*Y])
        {
            if(u==x&&v==y)
            {
                if(0)
                {
                    printf("reach (%d,%d)\n",x,y);
                    system("pause");
                }
                count[n][l][r] ++;
                return;
            }
            else
            {
                if(1)
                {
                    printBits<3,3,true>(n);
                    printf("(%d,%d)--->\n",u,v);
                }
                for(int _d=0;_d<4;_d++)
                {
                    if(_d==2)continue;
                    unsigned char nextx = movex(u,(d+_d)&3);
                    unsigned char nexty = movey(v,(d+_d)&3);
                    if(n&OCCUPY_BIT(nextx,nexty))continue;
                    worm_goto(nextx,nexty,x,y,l+(_d==1),r+(_d==3),(d+_d)&3,n|OCCUPY_BIT(nextx,nexty),count);
                }
            }
        }
        void savePropagatorDataTo(char* filename)
        {
            FILE* fp = fopen(filename,"wt");
            fprintf(fp,"loops = {");
            assert(singlelayer[0][0][0][0]==1);
            int first = 1;
            for(int ex=0;ex<2;ex++)
            {
                for(int ey=0;ey<2;ey++)
                {
                    for(int corner = 0;corner<=maxpoints;corner++)
                    {
                        for(int n1 = 0;n1<(1<<maxpoints);n1++)
                        {
                            if(singlelayer[ex][ey][corner][n1]==0)continue;
                            if(!first)
                            {
                                fprintf(fp,"\n,");
                            }
                            first = 0;
                            fprintf(fp,"(%d,%d,%d,%d):%d",ex,ey,corner,n1,singlelayer[ex][ey][corner][n1]);
                        }
                    }
                }
            }
            fprintf(fp,"}\n");
            first = 1;
            fprintf(fp,"occupy0 = %d\n",OCCUPY_BIT(0,0));
            fprintf(fp,"worms = {");
            for(int i=0;i<X;i++)
            {
                for(int j=0;j<Y;j++)
                {
                    if(i==0&&j==0)continue;
                    for(int d=0;d<4;d++)
                    {
                        int count[1<<(X*Y)][X*Y][X*Y]={};
                        int u = movex(0,d);
                        int v = movey(0,d);
                        worm_goto(u,v,i,j,0,0,d,OCCUPY_BIT(0,0)|OCCUPY_BIT(u,v),count);
                        for(int n=0;n<(1<<(X*Y));n++)
                        {
                            for(int r=0;r<X*Y;r++)
                            {
                                for(int l=0;l<X*Y;l++)
                                {
                                    if(count[n][l][r]==0)continue;
                                    if(!first)fprintf(fp,"\n,");
                                    first = 0;
                                    fprintf(fp,"(%d,%d,%d,%d,%d,%d):%d",i,j,d,n,l,r,count[n][l][r]);
                                }
                            }
                        }
                    }
                }
            }
            fprintf(fp,"}\n");
            fclose(fp);
        }
        #undef OCCUPY_BIT
    };
int main()
{
    char filename[256];
    propagatorStateNumber<3,3> counter;
    sprintf(filename,"../result/exhaustion/model_%d_%d_1/propagator_data.py",counter.X,counter.Y);
    //counter.savePropagatorDataTo(filename);
    #define OCCUPY_BIT(x,y) (1<<(x*3+y))
    for(int d=0;d<4;d++)
    {
        fprintf(stderr,"direct:%d\n",d);
        counter.worm_path(counter.movex(0,d),counter.movey(0,d),0,1,0,0,d,OCCUPY_BIT(0,0)|OCCUPY_BIT(counter.movex(0,d),counter.movey(0,d)),d);
    }
    #undef OCCUPY_BIT
    return 0;

}
