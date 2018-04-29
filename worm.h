#ifndef FERMION_WORM_H
#define FERMION_WORM_H
#include "common.h"
template<int X,int Y,int N>
class fermionworm
{
public:
    typedef unsigned short loctype;
    double m;
    double g;
    int debug;
    bool debug_quasi_pdf  = false;
    std::mt19937 randomizer;
    std::mt19937_64 randomizer64;
    #define RANDOM (randomizer()&0x7fffffff)
    #define RANDOM64 (randomizer64()&0x3fffffffffffffffull)
    double get_factor(int n,bool normalized)
    {
        if(normalized)return f[n];
        else return p[n];
    }
private:
    double p[2*N+1];
    double f[2*N+1];
    static_assert(X>1,"X<=1 GN-model is meaningless");
    static_assert(Y>1,"Y<=1 GN-model is meaningless");
    static_assert(N>0,"N<=0 GN-model is meaningless");
    static_assert(X<=30000,"X>30000 GN-model is too large");
    static_assert(Y<=30000,"Y>30000 GN-model is too large");
    static_assert(N<=200,"N>200 GN-model is too large");
    static const int debug_reallocate = 0;
    static const int debug_move = 0;
    static const int debug_reallocate_gamma5 = 0;
    static const int debug_move_gamma5 = 0;
    static const int Xh = (X+31)>>5;
    static const int Yh = (Y+31)>>5;
    static constexpr int log2X = log2_upper(X);
    static constexpr int log2Y = log2_upper(Y);
    static constexpr int upperX= 1<<log2X;
    static constexpr int upperY= 1<<log2Y;
    static const int Xl =  X&31;
    static const int Yl =  Y&31;
    static const int V = X*Y;
    static const int width = X>Y?Y:X;
    static const int length= X>Y?X:Y;
    static int iscorner;
    static int needbreak;
    static int isopen;
    static unsigned char ones[512];
    static unsigned char link2direct[16];
    static unsigned char breakphase[4][16][16];
    static unsigned char turnphase[4][16];
    static unsigned char breakopt[2][16];
    static unsigned char gamma5fuse[4][4];
    static unsigned char normalfuse[4][4];
    static unsigned char vproducts[4][4];
    static movement<X,4,8,2,loctype> movex;
    static movement<Y,4,1,4,loctype> movey;
    static distance<X,loctype> distancex;
    static distance<Y,loctype> distancey;
    static relative<X,loctype> relativex;
    static relative<Y,loctype> relativey;
    static unsigned char arbiter[4][16];
    static unsigned char statecorner[8];
    static char statephase[8];
    static const int maxdn = 2;
    static const int maxreweight = 4;
    static const int maxdcorner = 8;
    static const int maxpropagators = 8;
    unsigned long long int phy_prob[2*N+1][2*maxdn+1][2*N+1][2*maxdn+1];
    unsigned long long int geo_prob[maxpropagators][5][2*maxdcorner+1];
    int log_phy_prob[2*N+1][2*maxdn+1][2*N+1][2*maxdn+1];
    int log_geo_prob[maxpropagators][5][2*maxdcorner+1];
    double reweight[maxpropagators];
public:
    double quasipdf_weight;
    //lattice
    unsigned char link[2*N][X][Y];
    unsigned char occupying[X][Y];
    unsigned long long windings;
    //statistics
    long long count;
    unsigned long long acceptCount;
    unsigned long long randomCount;
    long long Ztopo;
    int n[2*N+1];
    unsigned long long  ncount[2*N+1];
    long long  pcount[maxpropagators][8][X][Y];
    unsigned long long runned[maxpropagators];
    double propagator[maxpropagators][upperX][upperY];
    double G[maxpropagators][upperX][upperY];
    double meff[maxpropagators][upperY];
    double fit_error[maxpropagators][upperY];
private:
    //coefficients compute
    #define LOOP_ASSERT(check,f) if(!(check)){assert_fail(f,headx,heady,tailx,taily);assert(check);}
    void assert_fail(int f,int headx,int heady,int tailx,int taily)
    {
        if(debug_move)
        {
            printLink(link,f,headx,heady,tailx,taily);
            for(int i=0;i<X;i++)
            {
                for(int j=0;j<Y;j++)
                {
                    fprintf(stderr,"%d ",occupying[i][j]);
                }
                fprintf(stderr,"\n");
            }
        }
    }
    inline void branch(const char name[])
    {
        //#define DEBUG_BRANCH
        #ifdef DEBUG_BRANCH
        fprintf(stderr,"Entered %s\n",name);
        #endif
    }
    inline void accept()
    {
        #define STAT_ACC_RATE
        //#define PRINT_ACC
        #ifdef STAT_ACC_RATE
        acceptCount++;
        randomCount++;
        #endif
        #ifdef PRINT_ACC
        fprintf(stderr,"Accepted\n");
        #endif // PRINT_ACC
    }
    inline bool decision(unsigned char id,char dcorner,char ddistance,unsigned char n1,char dn1,unsigned char n2,char dn2)
    {
        long long logp = log_phy_prob[n1][maxdn+dn1][n2][maxdn+dn2]+log_geo_prob[id][2+ddistance][maxdcorner+dcorner];
        long long prob = phy_prob[n1][maxdn+dn1][n2][maxdn+dn2]*    geo_prob[id][2+ddistance][maxdcorner+dcorner];
        assert(n1+dn1<=2*N);
        assert(n2+dn2<=2*N);
        assert(n1+dn1>=0);
        assert(n2+dn2>=0);
        //#define DEBUG_PROB
        #ifdef DEBUG_PROB

        #endif // DEBUG_PROB
        return (logp>0)||(prob>RANDOM64);
    }
    void print_prob(unsigned char id,char dcorner,char ddistance,unsigned char n1,char dn1,unsigned char n2,char dn2)
    {
        long long logp = log_phy_prob[n1][maxdn+dn1][n2][maxdn+dn2]+log_geo_prob[id][2+ddistance][maxdcorner+dcorner];
        long long prob = phy_prob[n1][maxdn+dn1][n2][maxdn+dn2]*    geo_prob[id][2+ddistance][maxdcorner+dcorner];
        fprintf(stderr,"id:%d\n",id);
        fprintf(stderr,"d-corner  :%d\n",dcorner  );
        fprintf(stderr,"d-distance:%d\n",ddistance);
        fprintf(stderr,"n1:%d->%d\n",n1,n1+dn1);
        fprintf(stderr,"n2:%d->%d\n",n2,n2+dn2);
        fprintf(stderr,"phy-prob:%f",phy_prob[n1][maxdn+dn1][n2][maxdn+dn2]*1.0/(1ull<<31));
        fprintf(stderr,"geo-prob:%f",geo_prob[id][2+ddistance][maxdcorner+dcorner]*1.0/(1ull<<31));
        fprintf(stderr,"prob:%f(log:%f)\n",prob*1.0/(1ull<<62),logp*1.0/(1<<16));
        if(0)system("pause");
    }
    inline void denied()
    {
        #ifdef STAT_ACC_RATE
        randomCount++;
        #undef STAT_ACC_RATE
        #endif
        #ifdef PRINT_ACC
        fprintf(stderr,"Denied\n");
        #undef PRINT_ACC
        #endif // PRINT_ACC
    }
    unsigned int fp2uint(double t)
    {
        if(t<=0.0)return 0;
        if(t>=1.0)return 1u<<31;
        return (unsigned int)(((1u<<31)*t)+0.5);
    }
    unsigned long long fp2ull(double t)
    {
        if(t<=0.0)return 0;
        if(t>=(1llu<<31))return 1llu<<62;
        return (unsigned long long)(((1llu<<31)*t)+0.5);
    }
    int logfp2int(double t)
    {
        return (int)((1<<16)*log(t));
    }
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
    void refresh_reweight(int i)
    {
        double r2i = sqrt(0.5);
        for(int dcorner = -maxdcorner;dcorner<=maxdcorner;dcorner++)
        {
            for(int j=-2;j<=2;j++)
            {
                geo_prob[i][2+j][maxdcorner+dcorner] = fp2ull(exp(reweight[i],j)*exp(r2i,dcorner));
                log_geo_prob[i][2+j][maxdcorner+dcorner] = logfp2int(exp(reweight[i],j)*exp(r2i,dcorner));
            }
        }
    }
    void refresh_coeff()
    {
        double t = 2.0+m;
        p[0] = 1.0;
        p[1] = 2+m;
        for(int i=2;i<=2*N;i++)
        {
            p[i] = p[i-1]*(2+m)+(i-1)*g*g*p[i-2];
        }
        f[0] = 1.0;
        for(int i=1;i<=2*N;i++)
        {
            f[i] = p[2*N-i]/p[2*N];
        }
        for(int dni=-maxdn;dni<=maxdn;dni++)
        {
            for(int dnj=-maxdn;dnj<=maxdn;dnj++)
            {
                for(int ni=0;ni<=2*N;ni++)
                {
                    for(int nj=0;nj<=2*N;nj++)
                    {
                        if(dni+ni<0||dni+ni>2*N
                            ||dnj+nj<0||dnj+nj>2*N)
                        {
                            phy_prob[ni][maxdn+dni][nj][maxdn+dnj] = 0;
                            log_phy_prob[ni][maxdn+dni][nj][maxdn+dnj] = -(1<<24);
                        }
                        else
                        {
                            phy_prob[ni][maxdn+dni][nj][maxdn+dnj] = fp2ull(f[ni+dni]*f[nj+dnj]/f[ni]/f[nj]);
                            log_phy_prob[ni][maxdn+dni][nj][maxdn+dnj] = logfp2int(f[ni+dni]*f[nj+dnj]/f[ni]/f[nj]);
                        }
                    }
                }
            }
        }
        for(int i=0;i<maxpropagators;i++)
        {
            refresh_reweight(i);
        }
    }
    int countWinding(unsigned char link[X][Y])
    {
        int winding = 0;
        unsigned char linkState = 0;
        for(int i=0;i<X;i++)
        {
            linkState ^= link[i][0];
        }
        if(linkState&4)winding ^= 2;
        linkState = 0;
        for(int i=0;i<Y;i++)
        {
            linkState ^= link[0][i];
        }
        if(linkState&2)winding ^= 1;
        return winding;
    }
public:
    fermionworm(double _g,double _m,unsigned int seed):g(_g),m(_m),randomizer(seed)
    {
        for(int i=0;i<maxpropagators;i++)reweight[i] = 1.0;
        refresh_coeff();
        zeroinitlattice();
        reset_statistic();
    }
    void reset_coeff(double _g,double _m)
    {
        m = _m;
        g = _g;
        refresh_coeff();
        reset_statistic();
    }
    double get_reweight(int id)
    {
        return reweight[id];
    }
    void set_reweight(int id,double w)
    {
        reweight[id] = w;
        refresh_reweight(id);
        reset_statistic();
    }
    void zeroinitlattice()
    {
        for(int f=0;f<N*2;f++)
        {
            for(int i=0;i<X;i++)
            {
                for(int j=0;j<Y;j++)
                {
                    link[f][i][j] = 0;
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
        start_check_occupy();
    }
    void initfromloop(int loop[2*N][Xh][Y])
    {
        loop2link<X,Y,N>(loop,link);
        windings  = 0;
        for(int i=0;i<X;i++)
        {
            for(int j=0;j<Y;j++)
            {
                occupying[i][j] = 0;
            }
        }
        for(int f=0;f<2*N;f++)
        {
            for(int i=0;i<X;i++)
            {
                for(int j=0;j<Y;j++)
                {
                    if(link[f][i][j])
                    {
                        occupying[i][j] += 1;
                    }
                }
            }
        }
        start_check_occupy();
    }
    void initfromloop(unsigned long long windings,int loop[2*N][Xh][Y])
    {
        loop2link<X,Y,N>(windings,loop,link);
        this->windings  = windings;
        for(int i=0;i<X;i++)
        {
            for(int j=0;j<Y;j++)
            {
                occupying[i][j] = 0;
            }
        }
        for(int f=0;f<2*N;f++)
        {
            for(int i=0;i<X;i++)
            {
                for(int j=0;j<Y;j++)
                {
                    if(link[f][i][j])
                    {
                        occupying[i][j] += 1;
                    }
                }
            }
        }
        start_check_occupy();
    }
private:
    template<bool enable_phase = true>
    unsigned char refresh_state(const unsigned char f,const loctype tailx,const loctype taily)
    {
        if(!enable_phase)
        {
            return 0;
        }
        else
        {
            if(ones[link[f][tailx][taily]]!=1)
            {
                return 0;
            }
            unsigned char d = link2direct[link[f][tailx][taily]];
            loctype tempx = movex(tailx,d);
            loctype tempy = movey(taily,d);
            unsigned char state = 0;
            unsigned char cur = link[f][tempx][tempy];
            while(ones[cur]==2)
            {
                cur ^= 1<<(d^2);
                if(ones[cur]!=1)
                {
                    fprintf(stderr,"encounter weird state");
                    printLink(link,f,tempx,tempy,tailx,taily);
                    assert(ones[cur]==1);
                }
                unsigned char _d = link2direct[cur];
                if(((d+1)&3)==_d)state ++;
                else if(((d-1)&3)==_d) state --;
                d = _d;
                tempx = movex(tempx,_d);
                tempy = movey(tempy,_d);
                cur = link[f][tempx][tempy];
            }
            return state&7;
        }
    }
    template<bool enable_phase = true>
    inline unsigned char psibargamma5psi_state(const unsigned char f1,const unsigned char f2,const loctype tailx,const loctype taily)
    {
        if(!enable_phase)
        {
            return 0;
        }
        else
        {
            return (refresh_state(f1,tailx,taily) - refresh_state(f2,tailx,taily)
                + turnphase[2^link2direct[link[f1][tailx][taily]]][link[f2][tailx][taily]])&7;
        }
    }
    bool print_corner_kinds = false;
    int quasipdf_corner(const unsigned char f1,const unsigned char f2,
                       const loctype tailx,const loctype taily,
                       const loctype headx,const loctype heady,
                       const loctype bodyx,const loctype bodyy)
    {
        int corner = 0;
        loctype nextx = headx;
        loctype nexty = heady;
        #define at(x,y,z) ((!(z&1)||(x==tailx&&y==taily))&&((x==headx&&y==heady)||!(z&2))&&((x==bodyx&&y==bodyy)||!(z&4)))
        for(int i=0;i<X;i++)
        {
            for(int j=0;j<Y;j++)
            {
                if(iscorner&(1<<link[f1][i][j]))corner ++;
                if(iscorner&(1<<link[f2][i][j]))corner ++;
                if(!at(i,j,1)&&!at(i,j,2)&&!at(i,j,4)&&(isopen&(1<<link[f1][i][j])))
                {
                    nextx = i;
                    nexty = j;
                }
                for(int f3=0;f3<2*N;f3++)
                {
                    if(f1==f3)continue;
                    if(f2==f3)continue;
                    if(iscorner&(1<<link[f3][i][j]))corner ++;
                }
            }
        }
        if((isopen&(1<<link[f2][headx][heady]))&&(isopen&(1<<link[f1][headx][heady]))
        && (isopen&(1<<link[f2][tailx][taily]))&&(isopen&(1<<link[f1][tailx][taily]))
        &&!(((isopen&(1<<link[f2][bodyx][bodyy]))||(isopen&(1<<link[f1][bodyx][bodyy])))^at(bodyx,bodyy,1)))
        {
            return corner
                +gamma5fuse[link2direct[link[f1][headx][heady]]][link2direct[link[f2][headx][heady]]]
                +gamma5fuse[link2direct[link[f1][tailx][taily]]][link2direct[link[f2][tailx][taily]]];
        }
        if(at(nextx,nexty,2)&&at(bodyx,bodyy,1))
        {
            corner += normalfuse[link2direct[link[f2][headx][heady]]][link2direct[link[f2][tailx][taily]]];
            if(print_corner_kinds)fprintf(stderr,"end points contribute %d corners, mode 1.\n",
                normalfuse[link2direct[link[f2][headx][heady]]][link2direct[link[f2][tailx][taily]]]);
        }
        else if(at(nextx,nexty,2))
        {
            corner +=
                + gamma5fuse[link2direct[link[f1][bodyx][bodyy]]][link2direct[link[f2][headx][heady]]]
                + gamma5fuse[link2direct[link[f1][tailx][taily]]][link2direct[link[f2][tailx][taily]]];
            if(print_corner_kinds)fprintf(stderr,"end points contribute %d&%d corners, mode 2.\n",
                + gamma5fuse[link2direct[link[f1][bodyx][bodyy]]][link2direct[link[f2][headx][heady]]],
                + gamma5fuse[link2direct[link[f1][tailx][taily]]][link2direct[link[f2][tailx][taily]]]);
            if(print_corner_kinds)fprintf(stderr,"where body %d connects head %d &tail %d connects tail %d.\n",
                link[f1][bodyx][bodyy],link[f2][headx][heady],
                link[f1][tailx][taily],link[f2][tailx][taily]);
        }
        else if(at(bodyx,bodyy,1))
        {
            corner +=
                + gamma5fuse[link2direct[link[f2][tailx][taily]]][link2direct[link[f1][nextx][nexty]]]
                + gamma5fuse[link2direct[link[f1][headx][heady]]][link2direct[link[f2][headx][heady]]];
            if(print_corner_kinds)fprintf(stderr,"end points contribute %d&%d corners, mode 3.\n",
                + gamma5fuse[link2direct[link[f1][bodyx][bodyy]]][link2direct[link[f2][headx][heady]]],
                + gamma5fuse[link2direct[link[f1][tailx][taily]]][link2direct[link[f2][tailx][taily]]]);
        }
        else
        {
            corner +=
                + normalfuse[link2direct[link[f1][bodyx][bodyy]]][link2direct[link[f1][nextx][nexty]]]
                + gamma5fuse[link2direct[link[f1][headx][heady]]][link2direct[link[f2][headx][heady]]]
                + gamma5fuse[link2direct[link[f1][tailx][taily]]][link2direct[link[f2][tailx][taily]]];
            if(print_corner_kinds)fprintf(stderr,"end points contribute %d&%d&%d corners, mode 4,next=(%d,%d).\n",
                + normalfuse[link2direct[link[f1][bodyx][bodyy]]][link2direct[link[f1][nextx][nexty]]],
                + gamma5fuse[link2direct[link[f1][headx][heady]]][link2direct[link[f2][headx][heady]]],
                + gamma5fuse[link2direct[link[f1][tailx][taily]]][link2direct[link[f2][tailx][taily]]],
                nextx,nexty);
        }
        #undef at
        return corner;
    }
    void quasipdf_place(unsigned char state,
                       const unsigned char f1,const unsigned char f2,
                       const loctype tailx,const loctype taily,
                       const loctype headx,const loctype heady,
                       const loctype bodyx,const loctype bodyy,
                       const unsigned char d0,unsigned long long quasipdf[2][length][length][length])
    {
        if(debug_quasi_pdf)
            fprintf(stderr,">> Entered quasi-PDF place (%d,%d),%d\n",bodyx,bodyy,d0);
        if(debug_quasi_pdf)
        {
            for(int i=0;i<X;i++)
            {
                for(int j=0;j<Y;j++)
                {
                    fprintf(stderr,"%d ",occupying[i][j]);
                }
                fprintf(stderr,"\n");
            }
            fprintf(stderr,"----- background f=%d-----\n",f2);
            printLink(link,f2,headx,heady,tailx,taily);
            fprintf(stderr,"----- background f=%d-----\n",f2);
            fprintf(stderr,"-----   front f=%d  ------\n",f1);
            printLink(link,f1,headx,heady,tailx,taily);
            fprintf(stderr,"-----   front f=%d  ------\n",f1);
            if(0)system("pause");
        }
        bool debug_prob = false;
        loctype nextx = movex(bodyx,d0);
        loctype nexty = movey(bodyy,d0);
        unsigned char body = link[f1][bodyx][bodyy];
        unsigned char next = link[f1][nextx][nexty];
        const unsigned char taild1 = link2direct[link[f1][tailx][taily]];
        const unsigned char taild2 = link2direct[link[f2][tailx][taily]];
        const unsigned char headd1 = link2direct[link[f1][headx][heady]];
        const unsigned char headd2 = link2direct[link[f2][headx][heady]];
        auto eqmask = [headx,tailx,heady,taily](int x,int y)->unsigned char
            {if(x==headx&&y==heady)return 16;else if(x==tailx&&y==taily)return 32;else return 0;};
        auto direct = [headx,tailx,heady,taily,taild2,headd2](int x,int y,int ctx)->unsigned char
            {if(x==headx&&y==heady)return headd2^2;
            else if(x==tailx&&y==taily)return taild2^2;
            else return (isopen&(1<<ctx))?link2direct[ctx]:0;};
        unsigned char _body = (body^(1<< d0   ))|eqmask(bodyx,bodyy);
        unsigned char _next = (next^(1<<(d0^2)))|eqmask(nextx,nexty);
        char dnbody = (_body!=0)-(body!=0);
        char dnnext = (_next!=0)-(next!=0);
        unsigned char _bodyd = direct(bodyx,bodyy,_body);
        unsigned char _nextd = direct(nextx,nexty,_next);
        const bool nobody = ((bodyx==tailx)&&(bodyy==taily));
        const bool check_occupy = false;
        if(check_occupy)start_check_occupy();
        if(check_occupy)end_check_occupy();
        bool quit;
        bool nonext = ((nextx==headx)&&(nexty==heady));
        bool _nonext = ((nextx==headx)&&(nexty==heady));
        int corner_cnt[20] = {};
        if(true)
        {
            corner_cnt[0] = +((!nobody&&!nonext)?normalfuse[link2direct[_next&0xf]][link2direct[_body&0xf]]:0);
            corner_cnt[1] = +(( nobody&&!nonext)?gamma5fuse[taild2][link2direct[_next&0xf]]:0);
            corner_cnt[2] = +((!nobody&& nonext)?gamma5fuse[headd2][link2direct[_body&0xf]]:0);
            corner_cnt[3] = +(( nobody&& nonext)?normalfuse[headd2][taild2]:0);
            corner_cnt[4] = -( nobody?gamma5fuse[taild2][taild1]:0);
            corner_cnt[5] = -( nonext?gamma5fuse[headd2][headd1]:0);
            corner_cnt[6] = -((iscorner&(1<<(body&0xf)))!=0);
            corner_cnt[7] = -((iscorner&(1<<(next&0xf)))!=0);
        }
        char dcorner =
                +((!nobody&&!nonext)?normalfuse[link2direct[_next&0xf]][link2direct[_body&0xf]]:0)
                +(( nobody&&!nonext)?gamma5fuse[taild2][link2direct[_next&0xf]]:0)
                +((!nobody&& nonext)?gamma5fuse[headd2][link2direct[_body&0xf]]:0)
                +(( nobody&& nonext)?normalfuse[headd2][taild2]:0)
                -( nobody?gamma5fuse[taild2][taild1]:0)
                -( nonext?gamma5fuse[headd2][headd1]:0)
                -((iscorner&(1<<body))!=0)
                -((iscorner&(1<<next))!=0);
        const bool debug_corner = false||debug_quasi_pdf;
        const bool check_corner = true ||debug_quasi_pdf;
        if(debug_quasi_pdf)
        {
            fprintf(stderr," next:%d", next);
            fprintf(stderr,"_next:%d",_next);
        }
        loctype tempx;
        loctype tempy;
        unsigned char temp;
        unsigned char _temp;
        unsigned char tempn;
        char dntemp;
        loctype passx;
        loctype passy;
        unsigned char pass;
        unsigned char _pass;
        unsigned char passn;
        char dnpass;
        unsigned char d;
        unsigned char _d;
        unsigned char _state;
        char ddistance =
                -distancex(bodyx,nextx)
                -distancey(bodyy,nexty);
        unsigned char bodyn = occupying[bodyx][bodyy];
        auto quasipdf_state = [this,f1,f2,tailx,headx,bodyx,taily,heady,bodyy,_body,direct](int x,int y,int ctx)->unsigned char
            {return (
                +refresh_state<true>(f1,tailx,taily)
                +vproducts[2^direct(bodyx,bodyy,_body)][direct(x,y,ctx)]
                +refresh_state<true>(f1,x,y)
                +((isopen&(1<<link[f1][headx][heady]))?vproducts[link2direct[link[f1][headx][heady]]][link2direct[link[f2][headx][heady]]]:0)
                +refresh_state<true>(f2,headx,heady)
                +((isopen&(1<<link[f1][tailx][taily]))?vproducts[link2direct[link[f2][tailx][taily]]][link2direct[link[f1][tailx][taily]]]:0))&7;};
        assert(nonext == ((nextx==headx)&&(nexty==heady)));
        int __local_corner = __corner;
        unsigned char nextn = occupying[nextx][nexty];
        if(debug_prob)print_prob(0,dcorner,ddistance,bodyn,dnbody,nextn,dnnext);
        if(decision(0,dcorner,ddistance,bodyn,dnbody,nextn,dnnext))
        {
            if(check_corner)
            {
                __local_corner = quasipdf_corner(f1,f2,tailx,taily,headx,heady,bodyx,bodyy);
                if(__corner!=__local_corner)
                {
                    fprintf(stderr,"expect %d corners but get %d of them\n",__local_corner,__corner);
                    fprintf(stderr,"-----  back f=%d  ------\n",f2);
                    printLink(link,f2,headx,heady,tailx,taily,bodyx,bodyy,nextx,nexty);
                    fprintf(stderr,"-----  back f=%d  ------\n",f2);
                    fprintf(stderr,"----- front f=%d  ------\n",f1);
                    printLink(link,f1,headx,heady,tailx,taily,bodyx,bodyy,nextx,nexty);
                    fprintf(stderr,"----- front f=%d  ------\n",f1);
                    assert(false);
                }
            }
            accept();
            change_corner(dcorner);
            print_occupy();
            link[f1][bodyx][bodyy] = _body&0xf;
            link[f1][nextx][nexty] = _next&0xf;
            change_occupy(occupying[nextx][nexty],dnnext);
            change_occupy(occupying[bodyx][bodyy],dnbody);
            print_occupy();
            occupying[bodyx][bodyy] += dnbody;
            occupying[nextx][nexty] += dnnext;
            state = quasipdf_state(nextx,nexty,_next);
            next = _next;
            body = _body;
            d = d0;
            if(check_corner)
            {
                __local_corner = quasipdf_corner(f1,f2,tailx,taily,headx,heady,bodyx,bodyy);
                if(__corner!=__local_corner)
                {
                    if(true)
                    {
                        for(int i=0;i<8;i++)
                            fprintf(stderr,"factor %d:%d\n",i,corner_cnt[i]);
                    }
                    fprintf(stderr,"computed: nobody = %d,nonext = %d\n",nobody,nonext);
                    fprintf(stderr,"points:tail = (%d,%d),head = (%d,%d),body = (%d,%d),next = (%d,%d)",
                            tailx,taily,headx,heady,bodyx,bodyy,nextx,nexty);
                    fprintf(stderr,"Error encountered at runs %d in placing quasi-PDF.\n",count);
                    fprintf(stderr,"originally have %d corners \n",__corner-dcorner);
                    fprintf(stderr,"expect %d corners but get %d of them\n",__local_corner,__corner);
                    print_corner_kinds = true;
                    fprintf(stderr,"where containing :\n");
                    quasipdf_corner(f1,f2,tailx,taily,headx,heady,bodyx,bodyy);
                    fprintf(stderr,"-----  back f=%d  ------\n",f2);
                    printLink(link,f2,headx,heady,tailx,taily,bodyx,bodyy,nextx,nexty);
                    fprintf(stderr,"-----  back f=%d  ------\n",f2);
                    fprintf(stderr,"----- front f=%d  ------\n",f1);
                    printLink(link,f1,headx,heady,tailx,taily,bodyx,bodyy,nextx,nexty);
                    fprintf(stderr,"----- front f=%d  ------\n",f1);
                    assert(false);
                }
            }
        }
        else
        {
            if(debug_quasi_pdf)fprintf(stderr,">> Denied quasi-PDF place (%d,%d),%d\n",bodyx,bodyy,d0);
            return;
        }
        while(nextx!=bodyx||nexty!=bodyy)
        {
            d = d^(randomizer()&2);
            tempx = movex(nextx,d);
            tempy = movey(nexty,d);
            if(debug_quasi_pdf)
            {
                print_occupy();
                fprintf(stderr,"----- front f=%d  ------\n",f1);
                printLink(link,f1,headx,heady,tailx,taily,bodyx,bodyy,nextx,nexty);
                fprintf(stderr,"----- front f=%d  ------\n",f1);
                fprintf(stderr,"try to go direct %d\n",d);
                if(0)system("pause");
            }

            if((state&3)==0)
                quasipdf[state>>2][(5&(1<<d0))?relativey(bodyy,nexty):relativex(bodyx,nextx)]
                    [(5&(1<<d0))?relativey(heady,taily):relativex(headx,tailx)]
                    [(5&(1<<d0))?relativex(headx,tailx):relativey(heady,taily)]++;
            if((!nobody)&&(tempx==tailx)&&(tempy==taily))
            {
                if(false)
                {
                    fprintf(stderr,"----- front f=%d  ------\n",f1);
                    printLink(link,f1,headx,heady,tailx,taily,bodyx,bodyy,nextx,nexty);
                    fprintf(stderr,"----- front f=%d  ------\n",f1);
                    fprintf(stderr,"T:(%d,%d)",tailx,taily);
                    fprintf(stderr,"N:(%d,%d)",nextx,nexty);
                    fprintf(stderr,"H:(%d,%d)",headx,heady);
                    fprintf(stderr,"t:(%d,%d)",tempx,tempy);
                    fprintf(stderr,"Forbidden move.\n");
                    if(1)system("pause");
                }
                continue;
            }
            temp = link[f1][tempx][tempy]|eqmask(tempx,tempy);
            _next=(next^(1<< d   ));
            _temp=(temp^(1<<(d^2)));

            if(((isopen&(1<<(temp&0xf)))!=0)&&((temp&(1<<d0))==0)&&((temp&(1<<(2^d0)))==0)&&(headx==tempx)&&(heady==tempy))
            {
                if(false)
                {
                    fprintf(stderr,"----- front f=%d  ------\n",f1);
                    printLink(link,f1,headx,heady,tailx,taily,bodyx,bodyy,nextx,nexty);
                    fprintf(stderr,"----- front f=%d  ------\n",f1);
                    fprintf(stderr,"T:(%d,%d)",tailx,taily);
                    fprintf(stderr,"N:(%d,%d)",nextx,nexty);
                    fprintf(stderr,"H:(%d,%d)",headx,heady);
                    fprintf(stderr,"t:(%d,%d)",tempx,tempy);
                    fprintf(stderr,"Forbidden move.\n");
                    if(1)system("pause");
                }
                continue;
            }
            if((needbreak&(1<<_temp)))
            {
                if(false&&(_temp&(1<<d)))
                {
                    //branch("break relink");
                    if(false)
                    {
                        fprintf(stderr,"break relink\n");
                    }
                    passx = movex(tempx,d);
                    passy = movey(tempy,d);
                    if((!nobody)&&(passx==tailx)&&(passy==taily))
                    {
                        continue;
                    }
                    quit    = (passx==bodyx)&&(passy==bodyy);
                    _nonext = (passx==headx)&&(passy==heady);
                    pass = link[f1][passx][passy];
                    passn= occupying[passx][passy];
                    _temp =  temp^(1<<(d^2))^(1<< d   );
                    _pass = (pass^(1<<(d^2)))|eqmask(passx,passy);
                    if(X<=2||Y<=2)
                    {
                        if(passx==nextx&&passy==nexty)
                        {
                            _next= _pass = pass^(1<<d)^(1<<(_d^2));
                        }
                    }
                    int dcorner;
                    if(debug_corner)
                    {
                        fprintf(stderr,"factor  1:%d\n",
                        + ((!nobody&&!_nonext&&!quit)?normalfuse[link2direct[_pass&0xf]][link2direct[_body&0xf]]:0));
                        fprintf(stderr,"factor  2:%d\n",
                        + (( nobody&&!_nonext&&!quit)?gamma5fuse[taild2][link2direct[_pass&0xf]]:0));
                        fprintf(stderr,"factor  3:%d\n",
                        + ((!nobody&& _nonext&&!quit)?gamma5fuse[headd2][link2direct[_body&0xf]]:0));
                        fprintf(stderr,"factor  4:%d\n",
                        + (( nobody&& _nonext&&!quit)?normalfuse[headd2][taild2]:0));
                        fprintf(stderr,"factor  5:%d\n",
                        + (( nobody&& quit)?gamma5fuse[taild2][link2direct[_pass&0xf]]:0));
                        fprintf(stderr,"factor  6:%d\n",
                        + (( nonext&& quit)?gamma5fuse[headd2][headd1]:0));
                        fprintf(stderr,"factor  7:%d\n",
                        + ((iscorner&(1<<(_pass&0xf)))!=0));
                        fprintf(stderr,"factor  8:%d\n",
                        + ((iscorner&(1<<(_next&0xf)))!=0));
                        fprintf(stderr,"factor  9:%d\n",
                        + ((iscorner&(1<<(_temp&0xf)))!=0));
                        fprintf(stderr,"factor 10:%d\n",
                        - ((!nobody&&!nonext&& quit)?normalfuse[link2direct[ next&0xf]][link2direct[_body&0xf]]:0));
                        fprintf(stderr,"factor 11:%d\n",
                        - (( nobody&&!nonext&& quit)?gamma5fuse[taild2][link2direct[ next&0xf]]:0));
                        fprintf(stderr,"factor 12:%d\n",
                        - ((!nobody&& nonext&& quit)?gamma5fuse[headd2][link2direct[_body&0xf]]:0));
                        fprintf(stderr,"factor 13:%d\n",
                        - (( nobody&& nonext&& quit)?normalfuse[headd2][taild2]:0));
                        fprintf(stderr,"factor 14:%d\n",
                        - ((!nobody&&!nonext&&!quit)?normalfuse[link2direct[ next&0xf]][link2direct[_body&0xf]]:0));
                        fprintf(stderr,"factor 15:%d\n",
                        - (( nobody&&!nonext&&!quit)?gamma5fuse[taild2][link2direct[ next&0xf]]:0));
                        fprintf(stderr,"factor 16:%d\n",
                        - ((!nobody&& nonext&&!quit)?gamma5fuse[headd2][link2direct[_body&0xf]]:0));
                        fprintf(stderr,"factor 17:%d\n",
                        - (( nobody&& nonext&&!quit)?normalfuse[headd2][taild2]:0));
                        fprintf(stderr,"factor 18:%d\n",
                        - ((iscorner&(1<<( next&0xf)))!=0));
                        fprintf(stderr,"factor 19:%d\n",
                        - ((iscorner&(1<<( pass&0xf)))!=0));
                        fprintf(stderr,"factor 20:%d\n",
                        - ((iscorner&(1<<( temp&0xf)))!=0));
                    }
                    dcorner =
                        + ((!nobody&&!_nonext&&!quit)?normalfuse[link2direct[_pass&0xf]][link2direct[_body&0xf]]:0)
                        + (( nobody&&!_nonext&&!quit)?gamma5fuse[taild2][link2direct[_pass&0xf]]:0)
                        + ((!nobody&& _nonext&&!quit)?gamma5fuse[headd2][link2direct[_body&0xf]]:0)
                        + (( nobody&& _nonext&&!quit)?normalfuse[headd2][taild2]:0)
                        + (( nobody&&           quit)?gamma5fuse[taild2][link2direct[_pass&0xf]]:0)
                        + (( nonext&&           quit)?gamma5fuse[headd2][headd1]:0)
                        + (( nonext&&          !quit)?gamma5fuse[headd2][link2direct[_next&0xf]]:0)
                        + ((iscorner&(1<<(_pass&0xf)))!=0)
                        + ((iscorner&(1<<(_next&0xf)))!=0)
                        + ((iscorner&(1<<(_temp&0xf)))!=0)
                        - ((!nobody&&!nonext&& quit)?normalfuse[link2direct[ next&0xf]][link2direct[_body&0xf]]:0)
                        - (( nobody&&!nonext&& quit)?gamma5fuse[taild2][link2direct[ next&0xf]]:0)
                        - ((!nobody&& nonext&& quit)?gamma5fuse[headd2][link2direct[_body&0xf]]:0)
                        - (( nobody&& nonext&& quit)?normalfuse[headd2][taild2]:0)
                        - ((!nobody&&!nonext&&!quit)?normalfuse[link2direct[ next&0xf]][link2direct[_body&0xf]]:0)
                        - (( nobody&&!nonext&&!quit)?gamma5fuse[taild2][link2direct[ next&0xf]]:0)
                        - ((!nobody&& nonext&&!quit)?gamma5fuse[headd2][link2direct[_body&0xf]]:0)
                        - (( nobody&& nonext&&!quit)?normalfuse[headd2][taild2]:0)
                        - ((_nonext&&!quit)?gamma5fuse[headd2][link2direct[ pass&0xf]]:0)
                        - ((iscorner&(1<<( next&0xf)))!=0)
                        - ((iscorner&(1<<( pass&0xf)))!=0)
                        - ((iscorner&(1<<( temp&0xf)))!=0);
                    char dnpass =
                        + (_pass!=0)
                        - ( pass!=0);
                    char dnnext =
                        + (_next!=0)
                        - ( next!=0);
                    char ddistance =
                        +distancex(bodyx,nextx)
                        +distancey(bodyy,nexty)
                        -distancex(bodyx,passx)
                        -distancey(bodyy,passy);
                    if(debug_prob)print_prob(0,dcorner,ddistance,passn,dnpass,nextn,dnnext);
                    if(decision(0,dcorner,ddistance,passn,dnpass,nextn,dnnext))
                    {
                        accept();
                        change_corner(dcorner);
                        print_occupy();
                        change_occupy(occupying[passx][passy],dnpass);
                        change_occupy(occupying[nextx][nexty],dnnext);
                        print_occupy();
                        link[f1][nextx][nexty] = _next&0xf;
                        link[f1][tempx][tempy] = _temp&0xf;
                        link[f1][passx][passy] = _pass&0xf;
                        nextx = passx;
                        nexty = passy;
                        nonext = _nonext;
                        occupying[passx][passy] = passn+dnpass;
                        next  = _pass;
                        nextn = passn+dnpass;
                        state = quasipdf_state(passx,passy,_pass);
                        if(check_corner)
                        {
                            int local_corner__ = quasipdf_corner(f1,f2,tailx,taily,headx,heady,bodyx,bodyy);
                            if(debug_corner)fprintf(stderr,"local corner:%d+%d should = %d\n",__local_corner,dcorner,local_corner__);
                            if(__local_corner+dcorner!=local_corner__)
                            {
                                fprintf(stderr,"Error encountered at runs %d in break relink quasi-PDF.\n",count);
                                fprintf(stderr,"-----  back f=%d  ------\n",f2);
                                printLink(link,f2,headx,heady,tailx,taily,bodyx,bodyy,nextx,nexty);
                                fprintf(stderr,"-----  back f=%d  ------\n",f2);
                                fprintf(stderr,"----- front f=%d  ------\n",f1);
                                printLink(link,f1,headx,heady,tailx,taily,bodyx,bodyy,nextx,nexty);
                                fprintf(stderr,"----- front f=%d  ------\n",f1);
                                assert(false);
                            }
                            __local_corner = local_corner__;
                        }
                    }
                    else
                    {
                        denied();
                    }
                }
                else denied();
            }
            else
            {
                branch("retrieve or occupy or remove");
                _nonext = (tempx==headx)&&(tempy==heady);
                quit    = (tempx==bodyx)&&(tempy==bodyy);
                char dntemp=
                    + (_temp!=0)
                    - ( temp!=0);
                char dnnext =
                    + (_next!=0)
                    - ( next!=0);
                /*
                char dcorner =
                    +((isopen&(1<<(_temp&0xf)))?normalfuse[direct(bodyx,bodyy,_body)][link2direct[_temp&0xf]]:((iscorner&(1<<(_temp&0xf)))!=0))
                    +((isopen&(1<<(_next&0xf)))?gamma5fuse[direct(nextx,nexty,_next)][link2direct[_next&0xf]]:((iscorner&(1<<(_next&0xf)))!=0))
                    -((isopen&(1<<( next&0xf)))?normalfuse[direct(bodyx,bodyy,_body)][link2direct[ next&0xf]]:((iscorner&(1<<( next&0xf)))!=0))
                    -((isopen&(1<<( temp&0xf)))?gamma5fuse[direct(tempx,tempy, temp)][link2direct[ temp&0xf]]:((iscorner&(1<<( temp&0xf)))!=0));
                */
                char dcorner;
                if(debug_corner)
                {
                    fprintf(stderr,"factor 1:%d(%d,%d)\n",
                            (( nobody&& quit)?gamma5fuse[taild2][link2direct[_temp&0xf]]:0),
                            taild2,link2direct[_temp&0xf]);
                    fprintf(stderr,"factor 2:%d(%d,%d)\n",
                            (( nonext&& quit)?gamma5fuse[headd2][link2direct[_next&0xf]]:0),
                            headd2,link2direct[_next&0xf]);

                    fprintf(stderr,"factor 3:%d(%d,%d)\n",
                            ((!nobody&&!_nonext&&!quit)?normalfuse[link2direct[_temp&0xf]][link2direct[_body&0xf]]:0),
                            link2direct[_temp&0xf],link2direct[_body&0xf]);
                    fprintf(stderr,"factor 4:%d(%d,%d)\n",
                            (( nobody&&!_nonext&&!quit)?gamma5fuse[taild2][link2direct[_temp&0xf]]:0),
                            taild2,link2direct[_temp&0xf]);
                    fprintf(stderr,"factor 5:%d(%d,%d)\n",
                            ((!nobody&& _nonext&&!quit)?gamma5fuse[headd2][link2direct[_body&0xf]]:0),
                            headd2,link2direct[_body&0xf]);
                    fprintf(stderr,"factor 6:%d(%d,%d)\n",
                            (( nobody&& _nonext&&!quit)?normalfuse[headd2][taild2]:0),
                            headd2,taild2);
                    fprintf(stderr,"factor 7:%d(%d,%d),%s\n",
                            - ((!nobody&&!nonext&&!quit)?normalfuse[link2direct[ next&0xf]][link2direct[_body&0xf]]:0),
                            link2direct[ next&0xf],link2direct[_body&0xf],
                            (!nobody&&!nonext&&!quit)?"active":"inactive");
                    fprintf(stderr,"factor 8:%d(%d,%d),%s\n",
                            - (( nobody&&!nonext&&!quit)?gamma5fuse[taild2][link2direct[ next&0xf]]:0),
                            taild2,link2direct[ next&0xf],
                            ( nobody&&!nonext&&!quit)?"active":"inactive");
                    fprintf(stderr,"factor 9:%d(%d,%d),%s\n",
                            - ((!nobody&& nonext&&!quit)?gamma5fuse[headd2][link2direct[_body&0xf]]:0),
                            headd2,link2direct[_body&0xf],
                            (!nobody&& nonext&&!quit)?"active":"inactive");
                    fprintf(stderr,"factor 10:%d(%d,%d),%s\n",
                            - (( nobody&& nonext&&!quit)?normalfuse[headd2][taild2]:0),
                            headd2,taild2,
                            ( nobody&& nonext&&!quit)?"active":"inactive");

                    fprintf(stderr,"factor 17:%d(%d)\n",
                            - ((iscorner&(1<<( next&0xf)))!=0), next&0xf);
                    fprintf(stderr,"factor 18:%d(%d)\n",
                            - ((iscorner&(1<<( temp&0xf)))!=0), temp&0xf);
                }
                dcorner =
                    + (( nobody&& quit)?gamma5fuse[taild2][link2direct[_temp&0xf]]:0)
                    + (( nonext&& quit)?gamma5fuse[headd2][link2direct[_next&0xf]]:0)
                    + ((!nobody&&!_nonext&&!quit)?normalfuse[link2direct[_temp&0xf]][link2direct[_body&0xf]]:0)
                    + (( nobody&&!_nonext&&!quit)?gamma5fuse[taild2][link2direct[_temp&0xf]]:0)
                    + ((!nobody&& _nonext&&!quit)?gamma5fuse[headd2][link2direct[_body&0xf]]:0)
                    + (( nobody&& _nonext&&!quit)?normalfuse[headd2][taild2]:0)
                    + (( nonext&&!quit)?gamma5fuse[headd2][link2direct[_next&0xf]]:0)
                    + ((iscorner&(1<<(_temp&0xf)))!=0)
                    + ((iscorner&(1<<(_next&0xf)))!=0)
                    - ((!nobody&&!nonext&& quit)?normalfuse[link2direct[ next&0xf]][link2direct[ temp&0xf]]:0)
                    - (( nobody&&!nonext&& quit)?gamma5fuse[taild2][link2direct[ next&0xf]]:0)
                    - ((!nobody&& nonext&& quit)?gamma5fuse[headd2][link2direct[ temp&0xf]]:0)
                    - (( nobody&& nonext&& quit)?normalfuse[headd2][taild2]:0)
                    - ((!nobody&&!nonext&&!quit)?normalfuse[link2direct[ next&0xf]][link2direct[_body&0xf]]:0)
                    - (( nobody&&!nonext&&!quit)?gamma5fuse[taild2][link2direct[ next&0xf]]:0)
                    - ((!nobody&& nonext&&!quit)?gamma5fuse[headd2][link2direct[_body&0xf]]:0)
                    - (( nobody&& nonext&&!quit)?normalfuse[headd2][taild2]:0)
                    - ((_nonext&&!quit)?gamma5fuse[headd2][link2direct[ temp&0xf]]:0)
                    - ((iscorner&(1<<( next&0xf)))!=0)
                    - ((iscorner&(1<<( temp&0xf)))!=0);
                char ddistance =
                    +distancex(bodyx,nextx)
                    +distancey(bodyy,nexty)
                    -distancex(bodyx,tempx)
                    -distancey(bodyy,tempy);
                unsigned char nextn = occupying[nextx][nexty];
                unsigned char tempn = occupying[tempx][tempy];
                if(false)
                {
                    fprintf(stderr,"context change:\n");
                    fprintf(stderr,"next:");
                    printBits<1,6,false>( next);
                    fprintf(stderr,"->");
                    printBits<1,6, true>(_next);
                    fprintf(stderr,"temp:");
                    printBits<1,6,false>( temp);
                    fprintf(stderr,"->");
                    printBits<1,6, true>(_temp);
                }
                if(debug_prob)print_prob(0,dcorner,ddistance,nextn,dnnext,tempn,dntemp);
                if(decision(0,dcorner,ddistance,nextn,dnnext,tempn,dntemp))
                {
                    accept();
                    change_corner(dcorner);
                    print_occupy();
                    change_occupy(occupying[tempx][tempy],dntemp);
                    change_occupy(occupying[nextx][nexty],dnnext);
                    print_occupy();
                    occupying[tempx][tempy] = tempn+dntemp;
                    occupying[nextx][nexty] = nextn+dnnext;
                    if(false)
                    {
                        fprintf(stderr,"update next(%d,%d):\n",nextx,nexty);
                        printBits<1,4,true>(_next);
                    }
                    link[f1][nextx][nexty] = _next&0xf;
                    link[f1][tempx][tempy] = _temp&0xf;
                    state = quasipdf_state(tempx,tempy,_temp);
                    nextx = tempx;
                    nexty = tempy;
                    next = _temp;
                    nonext = _nonext;
                    runned[0]++;
                    if(check_corner)
                    {
                        int local_corner__ = quasipdf_corner(f1,f2,tailx,taily,headx,heady,bodyx,bodyy);
                        if(debug_corner)fprintf(stderr,"local corner:%d+%d should = %d\n",__local_corner,dcorner,local_corner__);
                        if(__local_corner+dcorner!=local_corner__)
                        {
                            fprintf(stderr,"Error encountered at runs %d in other branch.\n",count);
                            fprintf(stderr,"-----  back f=%d  ------\n",f2);
                            printLink(link,f2,headx,heady,tailx,taily);
                            fprintf(stderr,"-----  back f=%d  ------\n",f2);
                            fprintf(stderr,"----- front f=%d  ------\n",f1);
                            printLink(link,f1,nextx,nexty,bodyx,bodyy);
                            fprintf(stderr,"----- front f=%d  ------\n",f1);
                            assert(false);
                        }
                        __local_corner = local_corner__;
                    }
                }
                else
                {
                    denied();
                }
            }
        }

        if(false)
        {
            bool err = false;
            for(int i=0;i<X;i++)
            {
                for(int j=0;j<Y;j++)
                {
                    if(i==headx&&j==heady)
                        if(!(isopen&(1<<link[f1][i][j]))||!(isopen&(1<<link[f2][i][j])))
                            err = true;
                    if(i==tailx&&j==taily)
                        if(!(isopen&(1<<link[f1][i][j]))||!(isopen&(1<<link[f2][i][j])))
                            err = true;
                }
            }
            printLink(link,f1,headx,heady,tailx,taily);
            if(err)
            {
                printLink(link,f1,headx,heady,tailx,taily);
                assert(false);
            }
        }
        if(check_occupy)
        {
            end_check_occupy();
        }
    }
    void quasipdf_scan(unsigned char state,
                       const unsigned char f1,const unsigned char f2,
                       const loctype tailx,const loctype taily,
                       const loctype headx,const loctype heady,
                       unsigned long long quasipdf[2][length][length][length])
    {
        if(debug_quasi_pdf)fprintf(stderr,">> Entered quasi-PDF scan.\n");
        loctype bodyx = tailx;
        loctype bodyy = taily;
        if(debug_quasi_pdf)fprintf(stderr,">> Entered quasi-PDF scan for f1.\n");
        unsigned char body = link[f1][bodyx][bodyy];
        unsigned char d = link2direct[body];
        unsigned char _d;
        int steps = 0;
        int skips = 0;
        quasipdf_place(state,f1,f2,tailx,taily,headx,heady,bodyx,bodyy,d,quasipdf);
        body = link[f1][bodyx][bodyy];
        d = link2direct[body];
        bodyx = movex(bodyx,d);
        bodyy = movey(bodyy,d);
        body = link[f1][bodyx][bodyy];
        loctype tempx = (iscorner&(1<<body))?bodyx:tailx;
        loctype tempy = (iscorner&(1<<body))?bodyy:taily;
        unsigned char tempd = d;
        bool debug_recover = false;
        if(debug_recover)
        {
            fprintf(stderr,"entered block from(%d,%d).\n",tailx,taily);
        }
        while(bodyx!=headx||bodyy!=heady)
        {
            if(debug_recover)
            {
                fprintf(stderr,"(%d,%d),from last goto %d\n",bodyx,bodyy,d);
                fprintf(stderr,"(check point)(%d,%d),from last goto %d\n",tempx,tempy,tempd);
            }
            if(skips>0)
            {
                body = link[f1][bodyx][bodyy]^(1<<(d^2));
                d = link2direct[body];
                bodyx = movex(bodyx,d);
                bodyy = movey(bodyy,d);
                steps++;
                skips--;
                if(debug_recover&&skips == 0)
                {
                    fprintf(stderr,"recover.\n");
                    system("pause");
                }
            }
            else
            {
                body = link[f1][bodyx][bodyy]^(1<<(d^2));
                if(ones[body]==1&&(link[f1][bodyx][bodyy]))
                {
                    _d = link2direct[body];
                    quasipdf_place(state,f1,f2,tailx,taily,headx,heady,bodyx,bodyy,_d,quasipdf);
                    steps++;
                    body = link[f1][bodyx][bodyy]^(1<<(d^2));
                }
                if(ones[body]!=1||(!link[f1][bodyx][bodyy]))
                {
                    body = link[f1][tailx][taily];
                    d = link2direct[body];
                    bodyx = movex(tailx,d);
                    bodyy = movey(taily,d);
                    skips = steps;
                    steps = 0;
                }
                else
                {
                    d = link2direct[body];
                    bodyx = movex(bodyx,d);
                    bodyy = movey(bodyy,d);
                }
            }
        }
    }
    //#define ACTIVE
    inline void begin_reallocate(const char name[])
    {
        #ifdef ACTIVE
        fprintf(stderr,"------reallocate %s begin------\n",name);
        for(int f = 0;f<2*N;f++)
        {
            printLink(link,f,0,0,0,0);
            fprintf(stderr,"------ f = %d ------\n",f);
        }
        for(int i=0;i<X;i++)
        {
            for(int j=0;j<Y;j++)
            {
                fprintf(stderr,"%d ",occupying[i][j]);
            }
            fprintf(stderr,"\n");
        }
        #undef ACTIVE
        #endif // ACTIVE
    }
    //#define ACTIVE
    inline void end_reallocate(const unsigned char state,const loctype tailx,const loctype taily,const loctype headx,const loctype heady)
    {
        #ifdef ACTIVE
        fprintf(stderr,"state:%d\n",state);
        fprintf(stderr,"------reallocate end ------\n");
        #endif // ACTIVE
    }
    template<typename ... T>
    inline void end_reallocate(const unsigned char state,const loctype tailx,const loctype taily,const loctype headx,const loctype heady,const unsigned char f,T...otherf)
    {
        #ifdef ACTIVE
        printLink(link,f,headx,heady,tailx,taily);
        fprintf(stderr,"------ f = %d ------\n",f);
        end_reallocate(state,tailx,taily,headx,heady,otherf...);
        #undef ACTIVE
        #endif // ACTIVE
    }

    //#define ACTIVE
    inline void end_step(const unsigned char state,const unsigned int runned,const loctype tailx,const loctype taily,const loctype headx,const loctype heady)
    {
        #ifdef ACTIVE
        for(int i=0;i<X;i++)
        {
            for(int j=0;j<Y;j++)
            {
                fprintf(stderr,"%d ",occupying[i][j]);
            }
            fprintf(stderr,"\n");
        }
        fprintf(stderr,"state:%d\n",state);
        fprintf(stderr,"-----step %d end -----\n",runned);
        if(1)system("pause");
        #endif // ACTIVE
    }
    template<typename ... T>
    inline void end_step(const unsigned char state,const unsigned long long runned[maxpropagators],const loctype tailx,const loctype taily,const loctype headx,const loctype heady,const unsigned char f,T...otherf)
    {
        #ifdef ACTIVE
        printLink(link,f,headx,heady,tailx,taily);
        fprintf(stderr,"------ f = %d ------\n",f);
        end_step(state,runned,tailx,taily,headx,heady,otherf...);
        #undef ACTIVE
        #endif // ACTIVE
    }
    #define ACTIVE
    void check_is_loop()
    {
        #ifdef ACTIVE
        for(int i=0;i<X;i++)
        {
            for(int j=0;j<Y;j++)
            {
                int occu = 0;
                for(int f=0;f<2*N;f++)
                {
                    unsigned char ctx = link[f][i][j];
                    assert(ones[ctx]==0||ones[ctx]==2);
                    if(ctx)occu++;
                }
                assert(occu==occupying[i][j]);
            }
        }
        #undef ACTIVE
        #endif // ACTIVE
    }
    int count_corner(unsigned char f)
    {
        int corner = 0;
        for(int i=0;i<X;i++)
        {
            for(int j=0;j<Y;j++)
            {
                if(iscorner&(1<<link[f][i][j]))corner++;
            }
        }
        return corner;
    }
    int count_corner_head_tail(unsigned char f,unsigned char fuse[4][4])
    {
        int corner = 0;
        short head = -1;
        for(int i=0;i<X;i++)
        {
            for(int j=0;j<Y;j++)
            {
                unsigned char ctx = link[f][i][j];
                if(iscorner&(1<<ctx))corner++;
                if(isopen&(1<<ctx))
                {
                    if(head==-1)head = ctx;
                    else if(head>0)corner += fuse[link2direct[ctx]][link2direct[head]];
                    else assert(false);
                }
            }
        }
        return corner;
    }
    int count_corner_fuse_local(unsigned char f1,unsigned char f2,unsigned char fuse[4][4])
    {
        int corner = 0;
        for(int i=0;i<X;i++)
        {
            for(int j=0;j<Y;j++)
            {
                unsigned char ctx1 = link[f1][i][j];
                unsigned char ctx2 = link[f2][i][j];
                if(iscorner&(1<<ctx1))corner++;
                if(iscorner&(1<<ctx2))corner++;
                if(isopen&(1<<ctx1))
                {
                    assert(isopen&(1<<ctx2));
                    corner += fuse[link2direct[ctx1]][link2direct[ctx2]];
                }
            }
        }
        return corner;
    }
    #define ACTIVE
    int __corner;
    inline void start_check_corner()
    {
        #ifdef ACTIVE
        __corner =0;
        for(int f=0;f<2*N;f++)
        {
            __corner += count_corner(f);
        }
        #endif // ACTIVE
    }
    inline void change_corner(int dc)
    {
        #ifdef ACTIVE
        //#define PRINT_DCORNER
        #ifdef PRINT_DCORNER
        fprintf(stderr,"corner:%d->%d(%d)\n",__corner,__corner+dc,dc);
        #endif // PRINT
        __corner += dc;
        #endif // ACTIVE
    }
    inline void end_check_corner()
    {
        #ifdef ACTIVE
        int corner__ = 0;
        for(int f=0;f<2*N;f++)
        {
            corner__ += count_corner(f);
        }
        if(__corner!=corner__)
        {
            for(int f=0;f<2*N;f++)
            {
                printLink(link,f,0,0,0,0);
            }
            fprintf(stderr,"error encountered at runs %d.\n",count);
            assert(__corner==corner__);
        }
        #endif // ACTIVE
    }
    inline void check_corner_head_tail(int f0,loctype headx,loctype heady,loctype tailx,loctype taily)
    {
        #define ACTIVE_INNER
        #ifdef ACTIVE
        #ifdef ACTVIE_INNER
        int corner__ = count_corner_head_tail(f,normalfuse);
        for(int f=0;f<2*N;f++)
        {
            if(f==f0)continue;
            corner__ += count_corner(f);
        }
        assert(__corner==corner__);
        #undef ACTIVE_INNER
        #endif // ACTVIE_INNER
        #undef ACTIVE
        #endif // ACTIVE
    }
    inline void check_corner_gamma5(int f1,int f2)
    {
        #define ACTIVE_INNER
        #ifdef ACTIVE
        #ifdef ACTVIE_INNER
        int corner__ = count_corner_fuse_local(f1,f2,gamma5fuse);
        for(int f=0;f<2*N;f++)
        {
            if(f==f1)continue;
            if(f==f2)continue;
            corner__ += count_corner(f);
        }
        assert(__corner==corner__);
        #undef ACTIVE_INNER
        #endif // ACTVIE_INNER
        #undef ACTIVE
        #endif // ACTIVE
    }
    void refresh_n()
    {
        for(int f=0;f<=2*N;f++)
            n[f] = 0;
        for(int i=0;i<X;i++)
        {
            for(int j=0;j<Y;j++)
            {
                int occu=0;
                for(int f=0;f<2*N;f++)
                {
                    if(link[f][i][j])
                    {
                        occu ++;
                    }
                }
                assert(occu == occupying[i][j]);
                n[occu]++;
            }
        }
    }
    #define ACTIVE
    inline void start_check_occupy()
    {
        #ifdef ACTIVE
        refresh_n();
        #endif
    }
    void print_occupy()
    {
        //#define PRINT_OCCUPY
        #ifdef PRINT_OCCUPY
        fprintf(stderr,"n:");
        for(int i=0;i<=2*N;i++)
        {
            fprintf(stderr,"%d ",n[i]);
        }
        fprintf(stderr,"\n");
        #undef PRINT_OCCUPY
        #endif // PRINT_OCCUPY
    }
    void change_occupy(const int occupy,const int doccupy)
    {
        n[occupy]         --;
        n[occupy+doccupy] ++;
    }
    void end_check_occupy()
    {
        #ifdef ACTIVE
        int n__[2*N+1];
        for(int f=0;f<=2*N;f++)
            n__[f] = 0;
        for(int i=0;i<X;i++)
        {
            for(int j=0;j<Y;j++)
            {
                int occu=0;
                for(int f=0;f<2*N;f++)
                {
                    if(link[f][i][j])
                    {
                        occu ++;
                    }
                }
                if(occu != occupying[i][j])
                {
                    for(int f=0;f<2*N;f++)
                    {
                        printLink(link,f,0,0,0,0);
                    }
                    fprintf(stderr,"processed %d\n",count);
                    assert(occu == occupying[i][j]);
                }
                n__[occu]++;
            }
        }
        for(int f=0;f<=2*N;f++)
            if(n[f]!=n__[f])
            {
                for(int f=0;f<2*N;f++)
                {
                    printLink(link,f,0,0,0,0);
                }
                fprintf(stderr,"processed %d\n",count);
                assert(n[f]==n__[f]);
            }
        #undef ACTIVE
        #endif
    }
public:
    unsigned int quasipdf_reallocate(unsigned long long quasipdf[][X][Y][length])
    {
        if(debug_quasi_pdf)
        {
            fprintf(stderr,">> Entered quasi-PDF reallocate\n");
        }
        unsigned char f1 = randomizer()%(2*N);
        unsigned char f2 = (f1+1+randomizer()%(2*N-1))%(2*N);
        loctype tailx = randomizer()%X;
        loctype taily = randomizer()%Y;
        unsigned char _d = randomizer();
        unsigned char d = _d&3;
        unsigned char tail1 = link[f1][tailx][taily];
        unsigned char tail2 = link[f2][tailx][taily];
        unsigned char _tail1 = tail1^(1<<d);
        unsigned char _tail2 = tail2^(1<<d);
        if(needbreak&(1<<_tail1))return 0;
        if(needbreak&(1<<_tail2))return 0;
        if((_tail1<<2)==(_tail2)||(_tail1)==(_tail2<<2))return 0;
        loctype headx = movex(tailx,d);
        loctype heady = movey(taily,d);
        unsigned char head1 = link[f1][headx][heady];
        unsigned char head2 = link[f2][headx][heady];
        char dntail = (tail1==0)+(tail2==0);
        char dnhead  = ( head1==0)+( head2==0);
        unsigned char _head1 = head1^(1<<(d^2));
        unsigned char _head2 = head2^(1<<(d^2));
        unsigned char tailn = occupying[tailx][taily];
        if((!tail1)&&(!tail2))
        {
            if(decision(2,0,0,0,0,tailn,2))
            do
            {
                pcount[2][0][0][0] += 1;
            }
            while(!decision(2,0,0,0,0,tailn+2,-2));
        }
        if((needbreak&(1<<_head1))&&(needbreak&(1<<_head2)))
        {
            branch("break relink both");
            unsigned char _d;
            if(head1==head2)
            {
                _d = breakopt[(_d>>3)&1][head1];
            }
            else _d = link2direct[head1&head2];
            assert(_d>=0&&_d<4);
            loctype nextx = movex(headx,_d);
            loctype nexty = movey(heady,_d);
            unsigned char _head1 = head1^(1<<_d)^(1<<(d^2));
            unsigned char _head2 = head2^(1<<_d)^(1<<(d^2));
            unsigned char next1 = link[f1][nextx][nexty];
            unsigned char next2 = link[f2][nextx][nexty];
            unsigned char _next1= next1^(1<<(_d^2));
            unsigned char _next2= next2^(1<<(_d^2));
            // break relink both
            int dcorner =
                +(gamma5fuse[link2direct[_next1]][link2direct[_next2]])
                +((iscorner&(1<<_head1))!=0)
                +((iscorner&(1<<_head2))!=0)
                -((iscorner&(1<< head1))!=0)
                -((iscorner&(1<< head2))!=0)
                -((iscorner&(1<< next1))!=0)
                -((iscorner&(1<< next2))!=0);
            assert(dcorner%2==0);
            if(decision(1,dcorner,2,0,0,tailn,2))
            {
                accept();
                begin_reallocate("gamma5");
                start_check_corner();
                start_check_occupy();
                change_corner(dcorner);
                print_occupy();
                change_occupy(tailn,2);
                print_occupy();
                assert(dntail==2);
                occupying[tailx][taily]=tailn + 2;
                link[f1][tailx][taily] = _tail1;
                link[f2][tailx][taily] = _tail2;
                link[f1][headx][heady] = _head1;
                link[f2][headx][heady] = _head2;
                link[f1][nextx][nexty] = _next1;
                link[f2][nextx][nexty] = _next2;
                unsigned char state = psibargamma5psi_state<true>(f1,f2,tailx,taily);
                end_reallocate(state,tailx,taily,headx,heady,f1,f2);
                quasipdf_move(state,f1,f2,_next1,_next2,tailx,taily,nextx,nexty,quasipdf);
                check_is_loop();
                windings = (windings&(~(3<<(2*f1))))|(countWinding(link[f1])<<(2*f1));
                windings = (windings&(~(3<<(2*f2))))|(countWinding(link[f2])<<(2*f2));
                return 1;
            }
            else
            {
                if(debug_quasi_pdf)fprintf(stderr,"denied quasi-PDF move.\n");
                denied();
                return 0;
            }
        }
        else
        {
            branch("other");
            int dcorner =
                +(((needbreak&(1<<_head1))||(needbreak&(1<<_head2)))?2:gamma5fuse[link2direct[_head1]][link2direct[_head2]])
                +((iscorner&(1<<(tail1^tail2)))!=0)
                -((iscorner&(1<<tail1))!=0)
                -((iscorner&(1<<tail2))!=0)
                -((iscorner&(1<< head1))!=0)
                -((iscorner&(1<< head2))!=0);
            unsigned char headn = occupying[headx][heady];
            if(decision(1,dcorner,1,tailn,dntail,headn,dnhead))
            {
                accept();
                begin_reallocate("gamma5");
                start_check_corner();
                start_check_occupy();
                change_corner(dcorner);
                print_occupy();
                change_occupy(tailn,dntail);
                change_occupy(headn,dnhead);
                print_occupy();
                occupying[tailx][taily] = tailn + dntail;
                occupying[headx][heady] = headn + dnhead;
                link[f1][tailx][taily] = _tail1;
                link[f2][tailx][taily] = _tail2;
                link[f1][headx][heady] = _head1;
                link[f2][headx][heady] = _head2;
                unsigned char state = psibargamma5psi_state<true>(f1,f2,tailx,taily);
                end_reallocate(state,tailx,taily,headx,heady,f1,f2);
                quasipdf_move(state,f1,f2,_head1,_head2,tailx,taily,headx,heady,quasipdf);
                check_is_loop();
                windings = (windings&(~(3<<(2*f1))))|(countWinding(link[f1])<<(2*f1));
                windings = (windings&(~(3<<(2*f2))))|(countWinding(link[f2])<<(2*f2));
                return 1;
            }
            else
            {
                if(debug_quasi_pdf)fprintf(stderr,"denied quasi-PDF move.\n");
                denied();
                return 0;
            }
        }
        check_corner_gamma5(f1,f2);
        #pragma omp critical
        if(true&&runned[0]>X*Y)
        {
            fprintf(stderr,"runned %lld steps\n",runned[1]);
        }

    }
    template<bool enable_phase = true>
    unsigned int psi_reallocate()
    {
        unsigned char f=randomizer()%(2*N);
        loctype tailx = randomizer()%X;
        loctype taily = randomizer()%Y;
        unsigned char tail = link[f][tailx][taily];
        unsigned char _d= randomizer();
        unsigned char d = _d&3;
        unsigned char _tail = tail ^ (1<<d);
        if(needbreak&(1<<_tail))return 0;
        loctype headx = movex(tailx,d);
        loctype heady = movey(taily,d);
        unsigned char head = link[f][headx][heady];
        unsigned char _head= head^(1<<(d^2));
        unsigned char tailn = occupying[tailx][taily];
        if(!tail)
        {
            if(decision(0,0,0,0,0,tailn,1))
            do
            {
                pcount[0][0][0][0] += 2;
            }
            while(!decision(0,0,0,0,0,tailn+1,-1));
        }
        if(needbreak&(1<<_head))
        {
            branch("break relink");
            unsigned char _d = breakopt[(_d>>3)&1][head];
            loctype nextx = movex(headx,_d);
            loctype nexty = movey(heady,_d);
            unsigned char next = link[f][nextx][nexty];
            _head = head^(1<<( d^2))^(1<<_d);
            unsigned char _next = next^(1<<(_d^2));
            int dcorner =
                + normalfuse[link2direct[_tail]][link2direct[_next]]
                + ((iscorner&(1<<_head))!=0)
                - ((iscorner&(1<< head))!=0)
                - ((iscorner&(1<< next))!=0)
                ;
            if(decision(0,dcorner,2,0,0,tailn,1))
            {
                accept();
                begin_reallocate("normal");
                start_check_corner();
                start_check_occupy();
                change_corner(dcorner);
                print_occupy();
                change_occupy(tailn,1);
                print_occupy();
                link[f][tailx][taily] = _tail;
                link[f][headx][heady] = _head;
                link[f][nextx][nexty] = _next;
                occupying[tailx][taily] = tailn + 1;
                unsigned char state = refresh_state<enable_phase>(f,tailx,taily);
                end_reallocate(state,tailx,taily,headx,heady,f);
                psi_move<enable_phase>(state,f,_next,tailx,taily,nextx,nexty);
                windings = (windings&(~(3<<(2*f))))|(countWinding(link[f])<<(2*f));
                return 1;
            }
            else
            {
                denied();
                return 0;
            }
        }
        else
        {
            char dcorner =
                +normalfuse[link2direct[_tail]][link2direct[_head]]
                -((iscorner&(1<<tail))!=0)
                -((iscorner&(1<<head))!=0)
                ;
            char dnhead =
                (_head!=0)-(head!=0);
            char dntail =
                ( tail==0);
            unsigned char headn = occupying[headx][heady];
            if(decision(0,dcorner,1,tailn,dntail,headn,dnhead))
            {
                accept();
                begin_reallocate("normal");
                start_check_corner();
                start_check_occupy();
                change_corner(dcorner);
                print_occupy();
                change_occupy(tailn,dntail);
                change_occupy(headn,dnhead);
                print_occupy();
                occupying[tailx][taily] = tailn+dntail;
                occupying[headx][heady] = headn+dnhead;
                link[f][tailx][taily] = _tail;
                link[f][headx][heady] = _head;
                unsigned char state = refresh_state<enable_phase>(f,tailx,taily);
                end_reallocate(state,tailx,taily,headx,heady,f);
                psi_move<enable_phase>(state,f,_head,tailx,taily,headx,heady);
                windings = (windings&(~(3<<(2*f))))|(countWinding(link[f])<<(2*f));
                return 1;
            }
            else
            {
                denied();
                return 0;
            }
        }
        #pragma omp critical
        if(true&&runned[0]>X*Y)
        {
            fprintf(stderr,"runned %lld steps\n",runned[0]);
        }
    }
    template<bool enable_phase = true>
    unsigned int psibargamma5psi_reallocate()
    {
        unsigned char f1 = randomizer()%(2*N);
        unsigned char f2 = (f1+1+randomizer()%(2*N-1))%(2*N);
        loctype tailx = randomizer()%X;
        loctype taily = randomizer()%Y;
        unsigned char _d = randomizer();
        unsigned char d = _d&3;
        unsigned char tail1 = link[f1][tailx][taily];
        unsigned char tail2 = link[f2][tailx][taily];
        unsigned char _tail1 = tail1^(1<<d);
        unsigned char _tail2 = tail2^(1<<d);
        if(needbreak&(1<<_tail1))return 0;
        if(needbreak&(1<<_tail2))return 0;
        if((_tail1<<2)==(_tail2)||(_tail1)==(_tail2<<2))return 0;
        loctype headx = movex(tailx,d);
        loctype heady = movey(taily,d);
        unsigned char head1 = link[f1][headx][heady];
        unsigned char head2 = link[f2][headx][heady];
        char dntail = (tail1==0)+(tail2==0);
        char dnhead  = ( head1==0)+( head2==0);
        unsigned char _head1 = head1^(1<<(d^2));
        unsigned char _head2 = head2^(1<<(d^2));
        unsigned char tailn = occupying[tailx][taily];
        if((!tail1)&&(!tail2))
        {
            if(decision(1,0,0,0,0,tailn,2))
            do
            {
                pcount[1][0][0][0] += 2;
            }
            while(!decision(1,0,0,0,0,tailn+2,-2));
        }
        if((needbreak&(1<<_head1))&&(needbreak&(1<<_head2)))
        {
            branch("break relink both");
            unsigned char _d;
            if(head1==head2)
            {
                _d = breakopt[(_d>>3)&1][head1];
            }
            else _d = link2direct[head1&head2];
            assert(_d>=0&&_d<4);
            loctype nextx = movex(headx,_d);
            loctype nexty = movey(heady,_d);
            unsigned char _head1 = head1^(1<<_d)^(1<<(d^2));
            unsigned char _head2 = head2^(1<<_d)^(1<<(d^2));
            unsigned char next1 = link[f1][nextx][nexty];
            unsigned char next2 = link[f2][nextx][nexty];
            unsigned char _next1= next1^(1<<(_d^2));
            unsigned char _next2= next2^(1<<(_d^2));
            // break relink both
            int dcorner =
                +(gamma5fuse[link2direct[_next1]][link2direct[_next2]])
                +((iscorner&(1<<_head1))!=0)
                +((iscorner&(1<<_head2))!=0)
                -((iscorner&(1<< head1))!=0)
                -((iscorner&(1<< head2))!=0)
                -((iscorner&(1<< next1))!=0)
                -((iscorner&(1<< next2))!=0);
            assert(dcorner%2==0);
            if(decision(1,dcorner,2,0,0,tailn,2))
            {
                accept();
                begin_reallocate("gamma5");
                start_check_corner();
                start_check_occupy();
                change_corner(dcorner);
                print_occupy();
                change_occupy(tailn,2);
                print_occupy();
                assert(dntail==2);
                occupying[tailx][taily]=tailn + 2;
                link[f1][tailx][taily] = _tail1;
                link[f2][tailx][taily] = _tail2;
                link[f1][headx][heady] = _head1;
                link[f2][headx][heady] = _head2;
                link[f1][nextx][nexty] = _next1;
                link[f2][nextx][nexty] = _next2;
                unsigned char state = psibargamma5psi_state<enable_phase>(f1,f2,tailx,taily);
                end_reallocate(state,tailx,taily,headx,heady,f1,f2);
                psibargamma5psi_move<enable_phase>(state,f1,f2,_next1,_next2,tailx,taily,nextx,nexty);
                windings = (windings&(~(3<<(2*f1))))|(countWinding(link[f1])<<(2*f1));
                windings = (windings&(~(3<<(2*f2))))|(countWinding(link[f2])<<(2*f2));
                return 1;
            }
            else
            {
                denied();
                return 0;
            }
        }
        else
        {
            branch("other");
            int dcorner =
                +(((needbreak&(1<<_head1))||(needbreak&(1<<_head2)))?2:gamma5fuse[link2direct[_head1]][link2direct[_head2]])
                +((iscorner&(1<<(tail1^tail2)))!=0)
                -((iscorner&(1<<tail1))!=0)
                -((iscorner&(1<<tail2))!=0)
                -((iscorner&(1<< head1))!=0)
                -((iscorner&(1<< head2))!=0);
            unsigned char headn = occupying[headx][heady];
            if(decision(1,dcorner,1,tailn,dntail,headn,dnhead))
            {
                accept();
                begin_reallocate("gamma5");
                start_check_corner();
                start_check_occupy();
                change_corner(dcorner);
                print_occupy();
                change_occupy(tailn,dntail);
                change_occupy(headn,dnhead);
                print_occupy();
                occupying[tailx][taily] = tailn + dntail;
                occupying[headx][heady] = headn + dnhead;
                link[f1][tailx][taily] = _tail1;
                link[f2][tailx][taily] = _tail2;
                link[f1][headx][heady] = _head1;
                link[f2][headx][heady] = _head2;
                unsigned char state = psibargamma5psi_state<enable_phase>(f1,f2,tailx,taily);
                end_reallocate(state,tailx,taily,headx,heady,f1,f2);
                psibargamma5psi_move<enable_phase>(state,f1,f2,_head1,_head2,tailx,taily,headx,heady);
                windings = (windings&(~(3<<(2*f1))))|(countWinding(link[f1])<<(2*f1));
                windings = (windings&(~(3<<(2*f2))))|(countWinding(link[f2])<<(2*f2));
                return 1;
            }
            else
            {
                denied();
                return 0;
            }
        }
        check_corner_gamma5(f1,f2);
        #pragma omp critical
        if(true&&runned[0]>X*Y)
        {
            fprintf(stderr,"runned %lld steps\n",runned[1]);
        }
    }
private:
    inline void psi_contribute(const unsigned char state,const int f,const loctype tailx,const loctype taily,const loctype headx,const loctype heady)
    {
        assert(state>=0&&state<=7);
        //#define ACTIVE
        #ifdef ACTIVE
        unsigned char state__ = refresh_state(f,tailx,taily);
        //#define PRINT_CONTRIBUTE
        #ifdef PRINT_CONTRIBUTE
        fprintf(stderr,"state:%d(should %d)",state,state__);
        #undef PRINT_CONTRIBUTE
        #endif // PRINT_CONTRIBUTE
        assert(state==state__);
        #undef ACTVIE
        #endif // ACTIVE
        pcount[0][state][relativex(tailx,headx)][relativey(taily,heady)] ++;
    }
    inline void psibargamma5psi_contribute(const unsigned char state,const int f1,const int f2,const loctype tailx,const loctype taily,const loctype headx,const loctype heady)
    {
        assert(state>=0&&state<=7);
        //#define ACTIVE
        #ifdef ACTIVE
        unsigned char state__ = psibargamma5psi_state(f1,f2,tailx,taily);
        //#define PRINT_CONTRIBUTE
        #ifdef PRINT_CONTRIBUTE
        fprintf(stderr,"state:%d(should %d)",state,state__);
        #undef PRINT_CONTRIBUTE
        #endif // PRINT_CONTRIBUTE
        assert(state==state__);
        #undef ACTVIE
        #endif // ACTIVE
        pcount[1][state][relativex(tailx,headx)][relativey(taily,heady)] ++;
    }
    inline void quasipdf_contribute(const unsigned char state,const int f1,const int f2,
                                           const loctype tailx,const loctype taily,
                                           const loctype headx,const loctype heady,
                                           unsigned long long quasipdf[][X][Y][length])
    {
        if(debug_quasi_pdf)fprintf(stderr,">> Entered quasi-PDF contribute.\n");
        if(0)system("pause");
        psibargamma5psi_contribute(state,f1,f2,tailx,taily,headx,heady);
        switch(randomizer()&2)
        {
        case 0:
            quasipdf_scan(state,f1,f2,tailx,taily,headx,heady,quasipdf);
            break;
        case 1:
            quasipdf_scan(state,f2,f1,tailx,taily,headx,heady,quasipdf);
            break;
        case 2:
            quasipdf_scan(state,f1,f2,headx,heady,tailx,taily,quasipdf);
            break;
        default://case 3:
            quasipdf_scan(state,f2,f1,headx,heady,tailx,taily,quasipdf);
            break;
        }
    }
    template<bool enable_phase = true>
    void psi_move(unsigned char state,const int f,unsigned char head
                   ,const loctype tailx,const loctype taily,loctype headx,loctype heady)
    {
        runned[0] = 0;
        unsigned char d;
        unsigned char sel;
        unsigned char _d;
        const unsigned char taild = link2direct[link[f][tailx][taily]];
        unsigned char temp;
        unsigned char next;
        unsigned char _head;
        unsigned char _temp;
        unsigned char _next;
        unsigned char _state;
        loctype tempx;
        loctype tempy;
        while(tailx!=headx||taily!=heady)
        {
            psi_contribute(state,f,tailx,taily,headx,heady);
            _d = randomizer();
             d = _d    &3;
            sel=(_d>>3)&1;
            tempx = movex(headx,d);
            tempy = movey(heady,d);
            temp = link[f][tempx][tempy];
            _head=head^(1<< d   );
            _temp=temp^(1<<(d^2));
            {
                if(needbreak&(1<<_temp))
                {
                    branch("break relink");
                    _d = breakopt[sel][temp];
                    loctype nextx = movex(tempx,_d);
                    loctype nexty = movey(tempy,_d);
                    next = link[f][nextx][nexty];
                    _temp = temp^(1<<(d^2))^(1<< _d   );
                    _next = next           ^(1<<(_d^2));
                    if(X<=2||Y<=2)
                    {
                        if(nextx==headx&&nexty==heady)
                        {
                            _head= _next = next^(1<<d)^(1<<(_d^2));
                        }
                    }
                    int dcorner =
                        + ((isopen&(1<<_next))?normalfuse[taild][link2direct[_next]]:0)
                        + ((iscorner&(1<<_head))!=0)
                        + ((iscorner&(1<<_temp))!=0)
                        - normalfuse[taild][link2direct[head]]
                        - ((iscorner&(1<< temp))!=0)
                        - ((iscorner&(1<< next))!=0);
                    char dnnext = -(_next==0);
                    char ddistance = distancex(tailx,nextx)+distancey(taily,nexty)-distancex(tailx,headx)-distancey(taily,heady);
                    unsigned char nextn = occupying[nextx][nexty];
                    if(decision(0,dcorner,ddistance,0,0,nextn,dnnext))
                    {
                        accept();
                        change_corner(dcorner);
                        print_occupy();
                        change_occupy(occupying[nextx][nexty],dnnext);
                        print_occupy();
                        link[f][headx][heady] = _head;
                        link[f][tempx][tempy] = _temp;
                        link[f][nextx][nexty] = _next;
                        state = ( state + turnphase[d][head] - turnphase[d^2][temp^(1<<_d)]
                        + breakphase[_d][temp][next])&7;
                        headx = nextx;
                        heady = nexty;
                        occupying[nextx][nexty] += dnnext;
                        head = _next;
                        if(enable_phase&&(tailx!=headx||taily!=heady))state = refresh_state<enable_phase>(f,tailx,taily);
                        end_step(state,runned,tailx,taily,headx,heady,f);
                        runned[0]++;
                    }
                    else denied();
                }
                else
                {
                    branch("retrieve or occupy or remove");
                    char dntemp= (_temp!=0)-(temp!=0);
                    char dnhead =-(_head==0);
                    char dcorner =
                        +((isopen&(1<<_temp))?normalfuse[taild][link2direct[_temp]]:((iscorner&(1<<_temp))!=0))
                        +((iscorner&(1<<_head))!=0)
                        -normalfuse[taild][link2direct[head]]
                        -((iscorner&(1<<temp))!=0);
                    char ddistance =
                        +distancex(tailx,tempx)
                        +distancey(taily,tempy)
                        -distancex(tailx,headx)
                        -distancey(taily,heady);
                    unsigned char headn = occupying[headx][heady];
                    unsigned char tempn = occupying[tempx][tempy];
                    if(decision(0,dcorner,ddistance,headn,dnhead,tempn,dntemp))
                    {
                        accept();
                        change_corner(dcorner);
                        print_occupy();
                        change_occupy(occupying[tempx][tempy],dntemp);
                        change_occupy(occupying[headx][heady],dnhead);
                        print_occupy();
                        occupying[tempx][tempy] = tempn+dntemp;
                        occupying[headx][heady] = headn+dnhead;
                        link[f][headx][heady] = _head;
                        link[f][tempx][tempy] = _temp;
                        state = (state +(((1<<d)==head)?turnphase[d][temp]:turnphase[d][head]))&7;
                        headx = tempx;
                        heady = tempy;
                        head = _temp;
                        end_step(state,runned,tailx,taily,headx,heady,f);
                        runned[0]++;
                    }
                    else denied();
                }
            }
            check_corner_head_tail(f,headx,heady,tailx,taily);
        }
        check_is_loop();
        end_check_corner();
        end_check_occupy();
    }
    template<bool enable_phase = true>
    void psibargamma5psi_move(unsigned char state,const unsigned char f1,const unsigned char f2
        ,unsigned char cur1,unsigned char cur2
        ,const loctype tailx,const loctype taily,loctype headx,loctype heady)
    {
        runned[1] = 0;
        while(tailx!=headx||taily!=heady)
        {
            unsigned int d = randomizer()&3;
            if((!(needbreak&(1<<cur1)))&&(!(needbreak&(1<<cur2))))
                psibargamma5psi_contribute(state,f1,f2,tailx,taily,headx,heady);
            else
            {
                if(((needbreak&(1<<cur1))&&(!(cur1&(1<<d))))
                || ((needbreak&(1<<cur2))&&(!(cur2&(1<<d)))))
                {
                    continue;
                }
            }
            loctype tempx = movex(headx,d);
            loctype tempy = movey(heady,d);
            unsigned char _cur1 = cur1^(1<<d);
            unsigned char _cur2 = cur2^(1<<d);
            unsigned char temp1 = link[f1][tempx][tempy];
            unsigned char temp2 = link[f2][tempx][tempy];
            unsigned char _temp1 = temp1^(1<<(d^2));
            unsigned char _temp2 = temp2^(1<<(d^2));
            char dnhead =-( _cur1==0)-( _cur2==0);
            char dntemp =+(_temp1!=0)+(_temp2!=0)-(temp1!=0)-(temp2!=0);
            if(isopen&(1<<temp1))
            {
                assert(isopen&(1<<temp2));
                int dcorner =
                    +((iscorner&(1<<_temp1))!=0)
                    +((iscorner&(1<<_temp2))!=0)
                    +((iscorner&(1<< _cur1))!=0)
                    +((iscorner&(1<< _cur2))!=0)
                    -gamma5fuse[link2direct[temp1]][link2direct[temp2]]
                    -(((needbreak&(1<<cur1))||(needbreak&(1<<cur2)))?2:gamma5fuse[link2direct[ cur1]][link2direct[ cur2]])
                    ;
                char ddistance =
                     -distancex(tailx,headx)
                     -distancey(taily,heady);
                unsigned char headn = occupying[headx][heady];
                unsigned char tempn = occupying[tempx][tempy];
                if(decision(1,dcorner,ddistance,headn,dnhead,tempn,dntemp))
                {
                    accept();
                    change_corner(dcorner);
                    print_occupy();
                    change_occupy(occupying[headx][heady],dnhead);
                    change_occupy(occupying[tempx][tempy],dntemp);
                    print_occupy();
                    occupying[headx][heady] = headn+dnhead;
                    occupying[tempx][tempy] = tempn+dntemp;
                    link[f1][headx][heady] = _cur1;
                    link[f2][headx][heady] = _cur2;
                    link[f1][tempx][tempy] = _temp1;
                    link[f2][tempx][tempy] = _temp2;
                    break;
                }
            }
            else if((needbreak&(1<<_temp1))&&(needbreak&(1<<_temp2)))
            {
                unsigned char _d;
                if(temp1==temp2)
                {
                    _d = breakopt[randomizer()&1][temp1];
                }
                else _d = link2direct[temp1&temp2];
                loctype nextx = movex(tempx,_d);
                loctype nexty = movey(tempy,_d);
                _temp1 = temp1 ^(1<<(d^2))^(1<<_d);
                _temp2 = temp2 ^(1<<(d^2))^(1<<_d);
                unsigned char next1 = link[f1][nextx][nexty];
                unsigned char next2 = link[f2][nextx][nexty];
                unsigned char _next1= next1^(1<<(_d^2));
                unsigned char _next2= next2^(1<<(_d^2));
                char dntemp =+(_temp1!=0)+(_temp2!=0)-(temp1!=0)-(temp2!=0);
                char dnnext =+(_next1!=0)+(_next2!=0)-(next1!=0)-(next2!=0);
                assert(dntemp==0);
                int dcorner =
                    +(gamma5fuse[link2direct[_next1]][link2direct[_next2]])
                    +((iscorner&(1<<_temp1))!=0)
                    +((iscorner&(1<<_temp2))!=0)
                    +((iscorner&(1<< _cur1))!=0)
                    +((iscorner&(1<< _cur2))!=0)
                    -(gamma5fuse[link2direct[  cur1]][link2direct[  cur2]])
                    -((iscorner&(1<< temp1))!=0)
                    -((iscorner&(1<< temp2))!=0)
                    -((iscorner&(1<< next1))!=0)
                    -((iscorner&(1<< next2))!=0);
                char ddistance =
                    +distancex(tailx,nextx)
                    +distancey(taily,nexty)
                    -distancex(tailx,headx)
                    -distancey(taily,heady);
                unsigned char headn = occupying[headx][heady];
                unsigned char nextn = occupying[nextx][nexty];
                if(decision(1,dcorner,ddistance,headn,dnhead,nextn,dnnext))
                {
                    accept();
                    change_corner(dcorner);
                    print_occupy();
                    change_occupy(occupying[headx][heady],dnhead);
                    change_occupy(occupying[nextx][nexty],dnnext);
                    print_occupy();
                    link[f1][headx][heady] = _cur1;
                    link[f2][headx][heady] = _cur2;
                    link[f1][tempx][tempy] = _temp1;
                    link[f2][tempx][tempy] = _temp2;
                    link[f1][nextx][nexty] = _next1;
                    link[f2][nextx][nexty] = _next2;
                    occupying[headx][heady] += dnhead;
                    occupying[nextx][nexty] += dnnext;
                    headx = nextx;
                    heady = nexty;
                    cur1 = _next1;
                    cur2 = _next2;
                    if(tailx!=headx||taily!=heady)state = psibargamma5psi_state<enable_phase>(f1,f2,tailx,taily);
                    end_step(state,runned,tailx,taily,headx,heady,f1,f2);
                    runned[1]++;
                }
                else denied();
            }
            else
            {
                branch("other");
                int dcorner =
                    +(((needbreak&(1<<_temp1))||(needbreak&(1<<_temp2)))?2:gamma5fuse[link2direct[_temp1]][link2direct[_temp2]])
                    +((iscorner&(1<<_cur1))!=0)
                    +((iscorner&(1<<_cur2))!=0)
                    -(((needbreak&(1<<  cur1))||(needbreak&(1<<  cur2)))?2:gamma5fuse[link2direct[  cur1]][link2direct[  cur2]])
                    -((iscorner&(1<<temp1))!=0)
                    -((iscorner&(1<<temp2))!=0);
                char ddistance =
                    +distancex(tailx,tempx)
                    +distancey(taily,tempy)
                    -distancex(tailx,headx)
                    -distancey(taily,heady);
                unsigned char headn = occupying[headx][heady];
                unsigned char tempn = occupying[tempx][tempy];
                if(decision(1,dcorner,ddistance,headn,dnhead,tempn,dntemp))
                {
                    accept();
                    change_corner(dcorner);
                    print_occupy();
                    change_occupy(occupying[headx][heady],dnhead);
                    change_occupy(occupying[tempx][tempy],dntemp);
                    print_occupy();
                    link[f1][headx][heady] = _cur1;
                    link[f2][headx][heady] = _cur2;
                    link[f1][tempx][tempy] =_temp1;
                    link[f2][tempx][tempy] =_temp2;
                    occupying[headx][heady] = headn+dnhead;
                    occupying[tempx][tempy] = tempn+dntemp;
                    if(!enable_phase) state = 0;
                    else if((!(needbreak&(1<<_temp1)))&&(!(needbreak&(1<<_temp2))))
                    {
                        if((needbreak&(1<<  cur1))||(needbreak&(1<<  cur2)))
                        {
                            state = psibargamma5psi_state<enable_phase>(f1,f2,tailx,taily);
                        }
                        else
                        {
                            state = (state
                                +((1<<d)==cur1?turnphase[d][temp1]:turnphase[d][cur1])
                                -((1<<d)==cur2?turnphase[d][temp2]:turnphase[d][cur2]))&7;
                        }
                    }
                    cur1 = _temp1;
                    cur2 = _temp2;
                    headx = tempx;
                    heady = tempy;
                    end_step(state,runned,tailx,taily,headx,heady,f1,f2);
                    runned[1]++;
                }
                else denied();
            }
        }
        check_is_loop();
        end_check_corner();
        end_check_occupy();
        return;
    }

    void quasipdf_move(unsigned char state,const unsigned char f1,const unsigned char f2,
        unsigned char cur1,unsigned char cur2,
        const loctype tailx,const loctype taily,
        loctype headx,loctype heady,
        unsigned long long quasipdf[][X][Y][length])
    {
        if(debug_quasi_pdf)
        {
            fprintf(stderr,">> Entered quasi-PDF move.\n");
        }
        runned[1] = 0;
        while(tailx!=headx||taily!=heady)
        {
            unsigned int d = randomizer()&3;
            if((!(needbreak&(1<<cur1)))&&(!(needbreak&(1<<cur2))))
            {
                quasipdf_contribute(state,f1,f2,tailx,taily,headx,heady,quasipdf);
                cur1 = link[f1][headx][heady];
                cur2 = link[f2][headx][heady];
            }
            else
            {
                if(((needbreak&(1<<cur1))&&(!(cur1&(1<<d))))
                || ((needbreak&(1<<cur2))&&(!(cur2&(1<<d)))))
                {
                    continue;
                }
            }
            loctype tempx = movex(headx,d);
            loctype tempy = movey(heady,d);
            unsigned char _cur1 = cur1^(1<<d);
            unsigned char _cur2 = cur2^(1<<d);
            unsigned char temp1 = link[f1][tempx][tempy];
            unsigned char temp2 = link[f2][tempx][tempy];
            unsigned char _temp1 = temp1^(1<<(d^2));
            unsigned char _temp2 = temp2^(1<<(d^2));
            char dnhead =-( _cur1==0)-( _cur2==0);
            char dntemp =+(_temp1!=0)+(_temp2!=0)-(temp1!=0)-(temp2!=0);
            if(isopen&(1<<temp1))
            {
                assert(isopen&(1<<temp2));
                int dcorner =
                    +((iscorner&(1<<_temp1))!=0)
                    +((iscorner&(1<<_temp2))!=0)
                    +((iscorner&(1<< _cur1))!=0)
                    +((iscorner&(1<< _cur2))!=0)
                    -gamma5fuse[link2direct[temp1]][link2direct[temp2]]
                    -(((needbreak&(1<<cur1))||(needbreak&(1<<cur2)))?2:gamma5fuse[link2direct[ cur1]][link2direct[ cur2]])
                    ;
                char ddistance =
                     -distancex(tailx,headx)
                     -distancey(taily,heady);
                unsigned char headn = occupying[headx][heady];
                unsigned char tempn = occupying[tempx][tempy];
                if(decision(1,dcorner,ddistance,headn,dnhead,tempn,dntemp))
                {
                    accept();
                    change_corner(dcorner);
                    print_occupy();
                    change_occupy(occupying[headx][heady],dnhead);
                    change_occupy(occupying[tempx][tempy],dntemp);
                    print_occupy();
                    occupying[headx][heady] = headn+dnhead;
                    occupying[tempx][tempy] = tempn+dntemp;
                    link[f1][headx][heady] = _cur1;
                    link[f2][headx][heady] = _cur2;
                    link[f1][tempx][tempy] = _temp1;
                    link[f2][tempx][tempy] = _temp2;
                    break;
                }
            }
            else if((needbreak&(1<<_temp1))&&(needbreak&(1<<_temp2)))
            {
                unsigned char _d;
                if(temp1==temp2)
                {
                    _d = breakopt[randomizer()&1][temp1];
                }
                else _d = link2direct[temp1&temp2];
                loctype nextx = movex(tempx,_d);
                loctype nexty = movey(tempy,_d);
                _temp1 = temp1 ^(1<<(d^2))^(1<<_d);
                _temp2 = temp2 ^(1<<(d^2))^(1<<_d);
                unsigned char next1 = link[f1][nextx][nexty];
                unsigned char next2 = link[f2][nextx][nexty];
                unsigned char _next1= next1^(1<<(_d^2));
                unsigned char _next2= next2^(1<<(_d^2));
                char dntemp =+(_temp1!=0)+(_temp2!=0)-(temp1!=0)-(temp2!=0);
                char dnnext =+(_next1!=0)+(_next2!=0)-(next1!=0)-(next2!=0);
                assert(dntemp==0);
                int dcorner =
                    +(gamma5fuse[link2direct[_next1]][link2direct[_next2]])
                    +((iscorner&(1<<_temp1))!=0)
                    +((iscorner&(1<<_temp2))!=0)
                    +((iscorner&(1<< _cur1))!=0)
                    +((iscorner&(1<< _cur2))!=0)
                    -(gamma5fuse[link2direct[  cur1]][link2direct[  cur2]])
                    -((iscorner&(1<< temp1))!=0)
                    -((iscorner&(1<< temp2))!=0)
                    -((iscorner&(1<< next1))!=0)
                    -((iscorner&(1<< next2))!=0);
                char ddistance =
                    +distancex(tailx,nextx)
                    +distancey(taily,nexty)
                    -distancex(tailx,headx)
                    -distancey(taily,heady);
                unsigned char headn = occupying[headx][heady];
                unsigned char nextn = occupying[nextx][nexty];
                if(decision(1,dcorner,ddistance,headn,dnhead,nextn,dnnext))
                {
                    accept();
                    change_corner(dcorner);
                    print_occupy();
                    change_occupy(occupying[headx][heady],dnhead);
                    change_occupy(occupying[nextx][nexty],dnnext);
                    print_occupy();
                    link[f1][headx][heady] = _cur1;
                    link[f2][headx][heady] = _cur2;
                    link[f1][tempx][tempy] = _temp1;
                    link[f2][tempx][tempy] = _temp2;
                    link[f1][nextx][nexty] = _next1;
                    link[f2][nextx][nexty] = _next2;
                    occupying[headx][heady] += dnhead;
                    occupying[nextx][nexty] += dnnext;
                    headx = nextx;
                    heady = nexty;
                    cur1 = _next1;
                    cur2 = _next2;
                    if(tailx!=headx||taily!=heady)state = psibargamma5psi_state<true>(f1,f2,tailx,taily);
                    end_step(state,runned,tailx,taily,headx,heady,f1,f2);
                    runned[1]++;
                }
                else denied();
            }
            else
            {
                branch("other");
                int dcorner =
                    +(((needbreak&(1<<_temp1))||(needbreak&(1<<_temp2)))?2:gamma5fuse[link2direct[_temp1]][link2direct[_temp2]])
                    +((iscorner&(1<<_cur1))!=0)
                    +((iscorner&(1<<_cur2))!=0)
                    -(((needbreak&(1<<  cur1))||(needbreak&(1<<  cur2)))?2:gamma5fuse[link2direct[  cur1]][link2direct[  cur2]])
                    -((iscorner&(1<<temp1))!=0)
                    -((iscorner&(1<<temp2))!=0);
                char ddistance =
                    +distancex(tailx,tempx)
                    +distancey(taily,tempy)
                    -distancex(tailx,headx)
                    -distancey(taily,heady);
                unsigned char headn = occupying[headx][heady];
                unsigned char tempn = occupying[tempx][tempy];
                if(decision(1,dcorner,ddistance,headn,dnhead,tempn,dntemp))
                {
                    accept();
                    change_corner(dcorner);
                    print_occupy();
                    change_occupy(occupying[headx][heady],dnhead);
                    change_occupy(occupying[tempx][tempy],dntemp);
                    print_occupy();
                    link[f1][headx][heady] = _cur1;
                    link[f2][headx][heady] = _cur2;
                    link[f1][tempx][tempy] =_temp1;
                    link[f2][tempx][tempy] =_temp2;
                    occupying[headx][heady] = headn+dnhead;
                    occupying[tempx][tempy] = tempn+dntemp;
                    if((!(needbreak&(1<<_temp1)))&&(!(needbreak&(1<<_temp2))))
                    {
                        if((needbreak&(1<<  cur1))||(needbreak&(1<<  cur2)))
                        {
                            state = psibargamma5psi_state<true>(f1,f2,tailx,taily);
                        }
                        else
                        {
                            state = (state
                                +((1<<d)==cur1?turnphase[d][temp1]:turnphase[d][cur1])
                                -((1<<d)==cur2?turnphase[d][temp2]:turnphase[d][cur2]))&7;
                        }
                    }
                    cur1 = _temp1;
                    cur2 = _temp2;
                    headx = tempx;
                    heady = tempy;
                    end_step(state,runned,tailx,taily,headx,heady,f1,f2);
                    runned[1]++;
                }
                else denied();
            }
        }
        check_is_loop();
        end_check_corner();
        end_check_occupy();
        return;
    }
    long long windings_sign()
    {
        long long sign;
        if(windings&1)
        {
            sign = -1;
        }
        else sign = 1;
        for(int l=0;l<2;l++)
        {
            if((windings>>(2*l))&3)
            {
                sign = -sign;
            }
        }
        return sign;
    }
public:
    void skip(int runs)
    {
        for(int i=0;i<runs;i++)
        {
            psi_reallocate<false>();
        }
        reset_statistic();
    }
    void sample_single(int runs,bool enable_phase)
    {
        for(int i=0;i<runs;i++)
        {
            if(enable_phase)psi_reallocate<true>();
            else psi_reallocate<false>();
            record();
        }
    }
    void sample_single(int runs)
    {
        for(int i=0;i<runs;i++)
        {
            psi_reallocate<true>();
            record();
        }
    }
    void sample_gamma5(int runs)
    {
        for(int i=0;i<runs;i++)
        {
            psi_reallocate<true>();
            psibargamma5psi_reallocate<true>();
            record();
        }
    }
    void record()
    {
        Ztopo += windings_sign();
        count += 1;
    }
    void get_reweight_observer(int i,double __reach[distancex.maxd+1])
    {
        unsigned long long reach[X];
        for(int u=0;u<=distancex.maxd;u++)
        {
            __reach[u] = 0;
        }
        for(int u=0;u<X;u++)
        {
            reach[u] = 0;
            for(int v=0;v<Y;v++)
            {
                for(int s=0;s<8;s++)
                    reach[u] += pcount[i][s][u][v];
            }
            __reach[distancex(0,u)] += reach[u];
        }
        __reach[0] *= 2;
        if(X%2==1)__reach[distancex.maxd]*=2;
        double normalize = __reach[0];
        for(int u=0;u<=distancex.maxd;u++)
        {
            __reach[u] /= normalize;
        }
        return ;
    }
    bool set_reweight_by_observer(int i,double reach[distancex.maxd+1])
    {
        int maxnz = 1;
        int bound = distancex.maxd/2+1;
        for(int u=2;u<=bound;u++)
        {
            if((reach[u]>0.01)&&maxnz<u)maxnz = u;
        }
        double r;
        if(maxnz<bound)
        {
            r = exp(log(10)/bound);
            reweight[i] *= r;
            return false;
        }
        else
        {
            r = 1.0;
            return true;
        }
        return maxnz==bound;
    }
    void smart_reweight_by_table(int target)
    {
        if(target==0)
        {
            reweight[target] = 1.0;
        }
        else if(target == 1)
        {
            double predict[2][2] = {{1.0,1.0},{1.0,exp(0.3)}};
            double m_div[1] = {-0.3};
            double g_div[1] = {1.0};
            int m_id = 0;
            int g_id = 0;
            while(m>m_div[m_id]&&m_id<1)m_id++;
            while(g>g_div[g_id]&&g_id<1)g_id++;
            reweight[target] = predict[m_id][g_id];
        }
    }
    void smart_reweight()
    {
        //bool need = true;
        for(int need = 0;need<2;need+=1)
        {
            double reach[2][distancex.maxd+1];
            for(int i=0;i<10000;i++)
            {
                psi_reallocate();
                psibargamma5psi_reallocate();
                record();
            }
            get_reweight_observer(0,reach[0]);
            get_reweight_observer(1,reach[1]);
            bool ready[2];
            ready[0] = set_reweight_by_observer(0,reach[0]);
            ready[1] = set_reweight_by_observer(1,reach[1]);
            #pragma omp critical
            {
                fprintf(stderr,"set re-weight: %d:%f,%d,%f\n",0,reweight[0],1,reweight[1]);
            }
            refresh_reweight(0);
            refresh_reweight(1);
            reset_statistic();
            if(!ready[0]||!ready[1])need--;
        }
    }
    void reset_statistic()
    {
        count = 0;
        acceptCount = 0;
        randomCount = 0;
        Ztopo = 0;
        ndarray_zero(ncount,2*N+1);
        ndarray_zero(runned,maxpropagators);
        ndarray_zero(pcount,maxpropagators,8,X,Y);
    }
    void save_G_to(int id,int p,FILE* fp)
    {
        fprintf(fp,"[");
        for(int i=0;i<X;i++)
        {
            if(i)fprintf(fp,",");
            fprintf(fp,"%e",G[id][i][p]);
        }
        fprintf(fp,"]");
    }
    void summary(int i)
    {
        for(int u=0;u<X;u++)
        {
            for(int v=0;v<Y;v++)
            {
                propagator[i][u][v] = 0;
                for(int s=0;s<8;s++)
                {
                    propagator[i][u][v] += statephase[s]*pcount[i][s][u][v];
                }
                propagator[i][u][v] *= exp(reweight[i],-distancex(0,u)-distancey(0,v))/count;
            }
            if(upperY==Y)fft_real_even(propagator[i][u],G[i][u],log2Y,false);
            else
            {
                constexpr double pi = 4*atan(1);
                for(int p=0;p<Y;p++)
                {
                    G[i][u][p] = 0;
                    for(int v=0;v<Y;v++)
                    {
                        G[i][u][p] += propagator[i][u][v]*cos(p*v*2*pi/Y);
                    }
                }
            }
        }
    }
    #undef RANDOM
    #undef RANDOM64
};
#define WORM_PHASE
#ifdef WORM_PHASE
#define WORM_PHASE_ARRAY {1,1,0,-1,-1,-1,0,1}
#else
#define WORM_PHASE_ARRAY {1,1,0,1,1,1,0,1}
#endif
template<int X,int Y,int N>
unsigned char fermionworm<X,Y,N>::statecorner[8] = {0,1,2,1,0,1,2,1};
template<int X,int Y,int N>
unsigned char fermionworm<X,Y,N>::vproducts[4][4] =
        {{0,1,0,7}
        ,{7,0,1,0}
        ,{0,7,0,1}
        ,{1,0,7,0}};
template<int X,int Y,int N>
unsigned char fermionworm<X,Y,N>::gamma5fuse[4][4] =
        {{0,1,2,1}
        ,{1,0,1,2}
        ,{2,1,0,1}
        ,{1,2,1,0}};
template<int X,int Y,int N>
unsigned char fermionworm<X,Y,N>::normalfuse[4][4] =
        {{2,1,0,1}
        ,{1,2,1,0}
        ,{0,1,2,1}
        ,{1,0,1,2}};
template<int X,int Y,int N>
char fermionworm<X,Y,N>::statephase [8] = WORM_PHASE_ARRAY;
template<int X,int Y,int N>
unsigned char fermionworm<X,Y,N>::arbiter[4][16] = WORM_ARBITER_ARRAY_4;
template<int X,int Y,int N>
unsigned char fermionworm<X,Y,N>::ones[512] = ONES_ARRAY_8;
template<int X,int Y,int N>
unsigned char fermionworm<X,Y,N>::link2direct[16] = LINK2DIRECT_ARRAY;
template<int X,int Y,int N>
unsigned char fermionworm<X,Y,N>::breakphase[4][16][16] = WORM_BREAK_PHASE_ARRAY;
template<int X,int Y,int N>
unsigned char fermionworm<X,Y,N>::turnphase[4][16] = WORM_TURN_PHASE_ARRAY;
template<int X,int Y,int N>
unsigned char fermionworm<X,Y,N>::breakopt[2][16] = REALLOC_BREAK_OPT_ARRAY;
template<int X,int Y,int N>
relative<X,typename fermionworm<X,Y,N>::loctype> fermionworm<X,Y,N>::relativex;
template<int X,int Y,int N>
relative<Y,typename fermionworm<X,Y,N>::loctype> fermionworm<X,Y,N>::relativey;
template<int X,int Y,int N>
distance<X,typename fermionworm<X,Y,N>::loctype> fermionworm<X,Y,N>::distancex;
template<int X,int Y,int N>
distance<Y,typename fermionworm<X,Y,N>::loctype> fermionworm<X,Y,N>::distancey;
template<int X,int Y,int N>
movement<X,4,8,2,typename fermionworm<X,Y,N>::loctype> fermionworm<X,Y,N>::movex;
template<int X,int Y,int N>
movement<Y,4,1,4,typename fermionworm<X,Y,N>::loctype> fermionworm<X,Y,N>::movey;
template<int X,int Y,int N>
int fermionworm<X,Y,N>::isopen
        = 1<<(1)
        | 1<<(2)
        | 1<<(4)
        | 1<<(8);
template<int X,int Y,int N>
int fermionworm<X,Y,N>::iscorner
        = 1<<(1|2)
        | 1<<(2|4)
        | 1<<(4|8)
        | 1<<(8|1);
template<int X,int Y,int N>
int fermionworm<X,Y,N>::needbreak
        = 1<<(2|4|8)
        | 1<<(4|8|1)
        | 1<<(8|1|2)
        | 1<<(1|2|4);
#undef STAT_ACC_RATE
#endif // FERMION_WORM_H
