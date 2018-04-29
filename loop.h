#ifndef LOOP_H
#define LOOP_H
#include "common.h"
#define STAT_ACC_RATE
namespace assertion{
    template<typename Tx,typename Tl,typename Tr>
    void inRange(Tx x,Tl l,Tr r)
    {
        assert(l<=x);
        assert(x<=r);
    }
};
template<int X,int Y,int N>
class fermionloop
{
  public:
    double g;
    double m;
    int debug;
    std::mt19937 randomizer;
  private:
    typedef unsigned short loctype;
    double p[2*N+1];
    double f[2*N+1];
    double coeff[1+2*N][9];
    unsigned int coeffi[9][2*N+1][2*N+1][2*N+1][2*N+1][81];
    static_assert(X>1,"X<=1 GN-model is meaningless");
    static_assert(Y>1,"Y<=1 GN-model is meaningless");
    static_assert(N>0,"N<=0 GN-model is meaningless");
    static const int Xh = (X+31)>>5;
    static const int Yh = (Y+31)>>5;
    static const int Xl = (X-1)&31;
    static const int Yl = (Y-1)&31;
    static const int V = X*Y;
  // rules
    static movement<X,4,8,2,loctype> movex;
    static movement<Y,4,1,4,loctype> movey;
    static unsigned int forbidden[16];
    static unsigned char ones[512];
    static unsigned char corners[512];
    static unsigned char occupation[512];
    static unsigned char nchange[16][16];
    static unsigned char nchangeid[16][16];
  //runtime
    int ntrace[2*N+1];
    int ctrace;
    int dn[2*N+1];
    int dcorner;
  public:
  //lattice
    int loop[N*2][Xh][Y];
    char occupying[X][Y];
  //statistic
    int count;
    unsigned long long acceptCount;
    unsigned long long randomCount;
    unsigned long long nSum[2*N+1];
    unsigned long long ninjSum[2*N+1][2*N+1];
    double ninjAvr[2*N+1][2*N+1];
    //double ninjVar[2*N+1][2*N+1];
    double nAvr[2*N+1];
    double nVar[2*N+1];
    double ninjCoVar[2*N+1][2*N+1];
    double Chi;
    double Cx;
  private:
    unsigned int fp2uint(double t)
    {
        if(t<=0.0)return 0;
        if(t>=1.0)return 1u<<31;
        return (1u<<31)*t;
    }
  //coefficients compute
    void refresh_coeff()
    {
        double t = 2.0+m;
        p[0] = 1.0;
        p[1] = 2+m;
        for(int i=2;i<=2*N;i++)
        {
            p[i] = (2+m)*p[i-1]+(i-1)*g*g*p[i-2];
        }
        f[0] = 1.0;
        for(int i=1;i<=2*N;i++)
        {
            f[i] = p[2*N-i]/p[2*N];
        }
        for(int i=-4;i<=4;i++)
        {
            for(int j=0;j<=2*N;j++)
                coeff[j][4+i] = exp(f[j],i);
        }
        for(int c=-4;c<=4;c++)
          for(int i=0;i<=2*N;i++)
              for(int j=0;j<=2*N;j++)
                  for(int k=0;k<=2*N;k++)
                      for(int l=0;l<=2*N;l++)
                        for(int id=0;id<81;id++)
                        {
                            int _id = id;
                            int di = _id%3-1;
                            _id/=3;
                            int dj = _id%3-1;
                            _id/=3;
                            int dk = _id%3-1;
                            _id/=3;
                            int dl = _id%3-1;
                            if(i+di<0||i+di>2*N){coeffi[4+c][i][j][k][l][id] = 0;continue;}
                            if(j+dj<0||j+dj>2*N){coeffi[4+c][i][j][k][l][id] = 0;continue;}
                            if(k+dk<0||k+dk>2*N){coeffi[4+c][i][j][k][l][id] = 0;continue;}
                            if(l+dl<0||l+dl>2*N){coeffi[4+c][i][j][k][l][id] = 0;continue;}
                            int occupy[2*N+1] = {};
                            occupy[i]--;
                            occupy[i+di]++;
                            occupy[j]--;
                            occupy[j+dj]++;
                            occupy[k]--;
                            occupy[k+dk]++;
                            occupy[l]--;
                            occupy[l+dl]++;
                            bool first = true;
                            double r = 1.0;
                            for(int p = 0;p<=2*N;p++)
                            {
                                if(occupy[p]==0)continue;
                                if(first)
                                {
                                    r = coeff[p][4+occupy[p]];
                                    first = false;
                                }
                                else
                                {
                                    r *= coeff[p][4+occupy[p]];
                                }
                            }
                            coeffi[4+c][i][j][k][l][id] = fp2uint(r*exp(0.5,c/2));
                        }
    }
  public:
    fermionloop(double _g,double _m,unsigned int seed):g(_g),m(_m),debug(0),randomizer(seed)
    {
        zeroinitlattice();
        refresh_coeff();
        reset_statistic();
        #pragma omp critical
        {
            fprintf(stderr,"initialize with m:%.2f g:%.2f\n",m,g);
        }
    }
    void reset_coeff(double _g,double _m)
    {
        m = _m;
        g = _g;
        refresh_coeff();
        reset_statistic();
    }
  private:
    inline void updateContext(const int loop[Xh][Y] ,unsigned int &cur,const int & U,const int & u,const int & v) const
    {
        if(U==0&&u==0)
        {
            cur = (cur>>3)|(((loop[Xh-1][v]>>Xl)&1)<<6)|((loop[U][v]&3)<<7);
        }
        else if(U==Xh-1&&u==Xl&&Xl==0)
        {
            cur = (cur>>3)|(((loop[U-1][v]>>31)&1)<<6)|(((loop[U][v])&1)<<7)|(((loop[0][v])&1)<<8);
        }
        else if(U==Xh-1&&u==Xl&&Xl!=0)
        {
            cur = (cur>>3)|(((loop[U][v]>>((Xl+31)&31))&3)<<6)|(((loop[0][v])&1)<<8);
        }
        else if(u==0)
        {
            cur = (cur>>3)|(((loop[U-1][v]>>31)&1)<<6)|((loop[U][v]&3)<<7);
        }
        else if(u==31)
        {
            cur = (cur>>3)|(((loop[U][v]>>30)&3)<<6)|(((loop[U+1][v])&1)<<8);
        }
        else
        {
            cur = (cur>>3)|((loop[U][v]>>(u-1))&7)<<6;
        }
        return ;
    }
    inline int loadContext(const int loop[Xh][Y],const int & U,const int & u,const int & v) const
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
    inline unsigned int iProbability(const int cur,const int x,const int y,const int u,const int v)
    {
        int next = cur^16;
        dcorner  = corners[next]-corners[cur];
        return coeffi[4+dcorner]
        [occupying[u][v]][occupying[x][v]]
        [occupying[u][y]][occupying[x][y]]
        [nchangeid[occupation[cur]][occupation[next]]];
    }
  public:
    void zeroinitlattice()
    {
        ctrace = 0;
        ntrace[0] = V;
        for(int f=1;f<=N*2;f++)
        {
            ntrace[f] = 0;
        }
        for(int f=0;f<N*2;f++)
        {
            for(int i=0;i<Xh;i++)
            {
                for(int j=0;j<Y;j++)
                {
                    loop[f][i][j] = 0;
                }
            }
        }
        for(int i=0;i<X;i++)
        {
            for(int j=0;j<Y;j++)
            {
                occupying[i][j] = 0;
            }
        }
    }
    void heatbath(int r[Xh][Y])
    {
        for(int I=0;I<Xh;I++)
              for(int i=0;i<32;i++)
              {
                int x = I*32+i;
                if(x>=X)break;
                #ifdef DEBUG_PRINT
                unsigned char link[2*N][X][Y];
                loop2link<X,Y,N>(loop,link);
                fprintf(stderr,"#\n");
                for(int f=0;f<2*N;f++)
                {
                    fprintf(stderr,"----------Lattice start----------\n");
                    printLinkwithLoop(link,loop,f,0,0,0,0);
                    fprintf(stderr,"----------Lattice  end ----------\n");
                }
                fprintf(stderr,"#\n");
                #endif // DEBUG_PRINT
                unsigned int rcontext = loadContext(r,I,i,Y-1);
                int u = movex(x,1);
                int v = Y-2;
                #ifdef DEBUG_PROB
                static bool showed[512] = {};
                #endif // DEBUG_PROB
                // j = y-1;
                {
                    if(!((forbidden[rcontext>>5]>>(rcontext&31))&1))
                    {
                        unsigned int iprob = iProbability(rcontext,x,Y-1,u,v);
                        unsigned int irand = randomizer()&0x7fffffff;
                        #ifdef DEBUG_PROB
                        if(iprob==0)printContexts<true>(rcontext,rcontext^16);
                        if(!showed[rcontext])
                        {
                            showed[rcontext] = true;
                            printContext(rcontext,rcontext^16);
                            fprintf(stderr,"%d %d\n%d %d\n"
                                ,occupying[movex(x,1)][movey(Y-1,2)],occupying[movex(x,1)][Y-1]
                                ,occupying[x][movey(Y-1,2)],occupying[x][Y-1]);
                            fprintf(stderr,"Prob: %.5f Rand:%.5f",iprob*1.0/(1u<<31),irand*1.0/(1u<<31));
                            system("pause");
                        }
                        #endif // DEBUG_PROB
                        if(iprob>irand)
                        {
                            #ifdef DEBUG_PRINT
                            fprintf(stderr,">>>>>>>>>>Accept Load Step<<<<<<<<<<\n");
                            printContext(rcontext,rcontext^16);
                            fprintf(stderr,"n:");
                            for(int i=0;i<=2*N;i++)
                            {
                                fprintf(stderr,"%3d ",ntrace[i]);
                            }
                            fprintf(stderr,"\n");
                            fprintf(stderr,"%d %d\n%d %d\n"
                                    ,occupying[u][v],occupying[u][Y-1]
                                    ,occupying[x][v],occupying[x][Y-1]);
                            #endif // DEBUG_PRINT
                            ntrace[occupying[u][v]] --;ntrace[occupying[u][Y-1]] --;
                            ntrace[occupying[x][v]] --;ntrace[occupying[x][Y-1]] --;
                            #ifdef DEBUG_PRINT
                            fprintf(stderr,"n:");
                            for(int i=0;i<=2*N;i++)
                            {
                                fprintf(stderr,"%3d ",ntrace[i]);
                            }
                            fprintf(stderr,"\n");
                            #endif // DEBUG_PRINT
                            int changing = nchange[occupation[rcontext]][occupation[rcontext^16]];
                            occupying[u][v  ] += ((changing   )&3)-1;occupying[u][Y-1] += ((changing>>4)&3)-1;
                            occupying[x][v  ] += ((changing>>2)&3)-1;occupying[x][Y-1] += ((changing>>6)&3)-1;
                            ntrace[occupying[u][v  ]] ++;ntrace[occupying[u][Y-1]] ++;
                            ntrace[occupying[x][v  ]] ++;ntrace[occupying[x][Y-1]] ++;
                            #ifdef DEBUG_PRINT
                            printContext(rcontext,rcontext^16);
                            fprintf(stderr,"n:");
                            for(int i=0;i<=2*N;i++)
                            {
                                fprintf(stderr,"%3d ",ntrace[i]);
                            }
                            fprintf(stderr,"\n");
                            fprintf(stderr,"%d %d\n%d %d\n"
                                    ,occupying[u][v],occupying[u][Y-1]
                                    ,occupying[x][v],occupying[x][Y-1]);
                            for(int x=0;x<X;x++)
                            {
                                for(int y=0;y<Y;y++)
                                {
                                    fprintf(stderr,"%2d ",occupying[x][y]);
                                }
                                fprintf(stderr,"\n");
                            }
                            fprintf(stderr,">>>>>>>>>>>>End<<<<<<<<<<<<\n");
                            system("pause");
                            #endif // DEBUG_PRINT
                            #ifdef RUNTIME_CHECK
                            int total = 0;
                            for(int i=0;i<=2*N;i++)
                            {
                                assertion::inRange(ntrace[i],0,V);
                                total += ntrace[i];
                            }
                            assert(total == V);
                            #endif
                            ctrace+=dcorner;
                            r[I][Y-1]^=1<<i;
                            rcontext ^= 16;
                            #ifdef STAT_ACC_RATE
                            acceptCount += 1;
                            #endif
                        }
                        #ifdef STAT_ACC_RATE
                        randomCount += 1;
                        #endif
                    }
                }
                v = Y - 1;
                for(int j=1;j<Y;j++)
                {
                    #ifdef DEBUG_PRINT
                    unsigned char link[2*N][X][Y];
                    loop2link<X,Y,N>(loop,link);
                    fprintf(stderr,"#\n");
                    for(int f=0;f<2*N;f++)
                    {
                        fprintf(stderr,"---------Lattice start---------\n");
                        printLinkwithLoop(link,loop,f,0,0,0,0);
                        fprintf(stderr,"---------Lattice  end ---------\n");
                    }
                    fprintf(stderr,"#\n");
                    #endif // DEBUG_PRINT
                    updateContext(r,rcontext,I,i,j);
                    if(!((forbidden[rcontext>>5]>>(rcontext&31))&1))
                    {
                        unsigned int iprob = iProbability(rcontext,x,j-1,u,v);
                        unsigned int irand = (randomizer()&0x7fffffff);
                        #ifdef DEBUG_PROB
                        if(iprob==0)printContexts<true>(rcontext,rcontext^16);
                        if(!showed[rcontext])
                        {
                            showed[rcontext] = true;
                            printContext(rcontext,rcontext^16);
                            fprintf(stderr,"%d %d\n%d %d\n"
                                ,occupying[movex(x,1)][movey(j-1,2)],occupying[movex(x,1)][j-1]
                                ,occupying[x][movey(j-1,2)],occupying[x][j-1]);
                            fprintf(stderr,"Prob: %.5f Rand:%.5f",iprob*1.0/(1u<<31),irand*1.0/(1u<<31));
                            system("pause");
                        }
                        #endif // DEBUG_PROB
                        if(iprob>irand)
                        {
                            #ifdef DEBUG_PRINT
                            fprintf(stderr,">>>>>>>>>>Accept Update Step<<<<<<<<<<\n");
                            printContext(rcontext,rcontext^16);
                            fprintf(stderr,"n:");
                            for(int i=0;i<=2*N;i++)
                            {
                                fprintf(stderr,"%3d ",ntrace[i]);
                            }
                            fprintf(stderr,"\n");
                            fprintf(stderr,"%d %d\n%d %d\n"
                                    ,occupying[u][v],occupying[u][j-1]
                                    ,occupying[x][v],occupying[x][j-1]);
                            #endif // DEBUG_PRINT
                            ntrace[occupying[u][v]] --;ntrace[occupying[u][j-1]] --;
                            ntrace[occupying[x][v]] --;ntrace[occupying[x][j-1]] --;
                            #ifdef DEBUG_PRINT
                            fprintf(stderr,"n:");
                            for(int i=0;i<=2*N;i++)
                            {
                                fprintf(stderr,"%3d ",ntrace[i]);
                            }
                            fprintf(stderr,"\n");
                            #endif // DEBUG_PRINT
                            int changing = nchange[occupation[rcontext]][occupation[rcontext^16]];
                            occupying[u][v  ] += ((changing   )&3)-1;occupying[u][j-1] += ((changing>>4)&3)-1;
                            occupying[x][v  ] += ((changing>>2)&3)-1;occupying[x][j-1] += ((changing>>6)&3)-1;
                            ntrace[occupying[u][v  ]] ++;ntrace[occupying[u][j-1]] ++;
                            ntrace[occupying[x][v  ]] ++;ntrace[occupying[x][j-1]] ++;
                            #ifdef DEBUG_PRINT
                            fprintf(stderr,"changing:\n%d %d\n%d %d\n"
                                    ,((changing   )&3)-1,((changing>>4)&3)-1
                                    ,((changing>>2)&3)-1,((changing>>6)&3)-1);
                            fprintf(stderr,"n:");
                            for(int i=0;i<=2*N;i++)
                            {
                                fprintf(stderr,"%3d ",ntrace[i]);
                            }
                            fprintf(stderr,"\n");
                            fprintf(stderr,"%d %d\n%d %d\n"
                                    ,occupying[u][v],occupying[u][j-1]
                                    ,occupying[x][v],occupying[x][j-1]);
                            fprintf(stderr,">>>>>>>>>>>>End<<<<<<<<<<<<\n");
                            system("pause");
                            #endif // DEBUG_PRINT
                            ctrace+=dcorner;
                            r[I][j-1]^=1<<i;
                            rcontext ^= 16;
                            #ifdef RUNTIME_CHECK
                            int total = 0;
                            for(int i=0;i<=2*N;i++)
                            {
                                assertion::inRange(ntrace[i],0,V);
                                total += ntrace[i];
                            }
                            assert(total == V);
                            #endif
                            #ifdef STAT_ACC_RATE
                            acceptCount += 1;
                            #endif
                        }
                        #ifdef STAT_ACC_RATE
                        randomCount += 1;
                        #endif
                    }
                    v = j - 1;
                }
                #ifdef DEBUG_PRINT
                for(int x=0;x<X;x++)
                {
                    for(int y=0;y<Y;y++)
                    {
                        fprintf(stderr,"%2d ",occupying[x][y]);
                    }
                    fprintf(stderr,"\n");
                }
                system("pause");
                #endif // DEBUG_PRINT
              }
    }
    void skip(int runs)
    {
        for(int i=0;i<runs;i++)
        {
            for(int j=0;j<2*N;j++)
            {
                heatbath(loop[j]);
            }
        }
    }
    void sample(int runs,int drop)
    {
        for(int i=0;i<runs;i++)
        {
            skip(drop);
            record();
        }
    }
    void record()
    {
        count++;
        unsigned long long n[2*N+1];
        for(int k=0;k<=2*N;k++)
        {
            n[k] = ntrace[k];
            nSum[k] += n[k];
        }
        for(int j=0;j<=2*N;j++)
        {
            for(int k=0;k<=2*N;k++)
            {
                ninjSum[j][k] += n[j]*n[k];
            }
        }
    }
    //statistic
    void reset_statistic()
    {
        count = 0;
        acceptCount = 0;
        randomCount = 0;
        for(int k=0;k<=2*N;k++)
        {
            nSum[k] = 0;
        }
        for(int k=0;k<=2*N;k++)
        {
            nAvr[k] = 0.0;
            nVar[k] = 0.0;
            Chi = 0.0;
            Cx = 0.0;
        }
        for(int j=0;j<=2*N;j++)
        {
            for(int k=0;k<=2*N;k++)
            {
                ninjSum[j][k] =0;
                ninjAvr[j][k] = 0.0;
                ninjCoVar[j][k] = 0.0;
            }
        }
    }
    void getMoreStatistic()
    {
        for(int k=0;k<=2*N;k++)
            nAvr[k]=nSum[k]*1.0/count;
        for(int j=0;j<=2*N;j++)
        {
            for(int k=0;k<=2*N;k++)
            {
                ninjAvr[j][k] = ninjSum[j][k]*1.0/count;
                ninjCoVar[j][k] = ninjAvr[j][k] - nAvr[j]*nAvr[k];
            }
            nVar[j] = ninjAvr[j][j] - nAvr[j]*nAvr[j]; //  = ninjCoVar[j][j]
        }
        Chi = ChiValue();
        Cx  = CchiValue();
    }

    double ChiValue()
    {
        double chi = 0.0;
        for(int i=0;i<2*N;i++)
        {
            chi += nAvr[i]*(2*N-i)*exp(p[1],2*N-i-1)/p[2*N-i];
        }
        return -chi/V;
    }
    double CchiValue()
    {
        double cchi_part1 = 0.0;
        double cchi_part2 = 0.0;
        double cchi_part3 = 0.0;
        double cchi = 0.0;
        for(int i=0;i<2*N-1;i++)
        {
            cchi_part1 += nAvr[i]*(2*N-i)*(2*N-i-1)*exp(p[1],2*N-i-2)/p[2*N-i];
        }
        for(int i=0;i<2*N;i++)
        {
            cchi_part2 += nAvr[i]*(2*N-i)*(2*N-i)*exp(p[1],2*(2*N-i-1))/exp(p[2*N-i],2);
        }
        for(int i=0;i<2*N;i++)
        {
            double cchi_part3_part = 0.0;
            for(int j=0;j<2*N;j++)
            {
                cchi_part3_part += ninjCoVar[i][j]*(2*N-i)*(2*N-j)*exp(p[1],2*N-i-1+2*N-j-1)/p[2*N-i]/p[2*N-j];
            }
            cchi_part3 += cchi_part3_part;
        }
        cchi = cchi_part1-cchi_part2+cchi_part3;
        return -cchi/V;
    }
};
template<int X,int Y,int N>
movement<X,4,8,2,typename fermionloop<X,Y,N>::loctype> fermionloop<X,Y,N>::movex;
template<int X,int Y,int N>
movement<Y,4,1,4,typename fermionloop<X,Y,N>::loctype> fermionloop<X,Y,N>::movey;
template<int X,int Y,int N>
unsigned int fermionloop<X,Y,N>::forbidden[16] = LOOP_FORBIDDEN_ARRAY;
template<int X,int Y,int N>
unsigned char fermionloop<X,Y,N>::ones[512] = ONES_ARRAY_8;
template<int X,int Y,int N>
unsigned char fermionloop<X,Y,N>::corners[512] = LOOP_CORNER_ARRAY;
template<int X,int Y,int N>
unsigned char fermionloop<X,Y,N>::occupation[512] = LOOP_OCCUPATION_ARRAY;
template<int X,int Y,int N>
unsigned char fermionloop<X,Y,N>::nchange[16][16] = LOOP_NCHANGE_ARRAY;
template<int X,int Y,int N>
unsigned char fermionloop<X,Y,N>::nchangeid[16][16] = LOOP_NCHANGE_ID_ARRAY;
#undef STAT_ACC_RATE
#endif
