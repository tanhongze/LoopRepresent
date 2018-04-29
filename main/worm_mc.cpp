#include "../version/worm_version.h"
#include "../common.h"
#include "../loop.h"
#include "../worm.h"
#include "../utils/timer.h"
#define PARALLEL
int firstc= 0;
int firstm= 0;
void pushversionforward()
{
    FILE* fp = fopen("../version/worm_version.h","wt");
    fprintf(fp,"#define WORM_VERSION %d",WORM_VERSION+1);
    fclose(fp);
}
template<int X,int Y,int N>
void getzxx(auto g,int cases,int resultid,FILE *fp)
{
    int topt = 0;
    unsigned int seeds[cases];
    for(int i=0;i<cases;i++)
    {
        seeds[i] = mtrand();
    }
    #ifdef PARALLEL
    #pragma omp parallel for
    #endif // PARALLEL
    for(int i=0;i<cases;i++)
    {
        int task;
        #pragma omp critical
        {
            task = topt++;
        }
        double step = 0.4;
        double mc   = 1.0;
        int seedid = 0;
        fermionworm<X,Y,N> WR(g(task),mc,seeds[task]);
        WR.zeroinitlattice();
        WR.skip(40000);
        double mcr,mcl = -2.0;
        double Z;
        while(true)
        {
            WR.reset_statistic();
            WR.sample_single(40000);
            Z = WR.Ztopo*1.0/WR.count;
            if(Z>0.1)break;
            else
            {
                mc+=step;
                WR.reset_coeff(g(task),mc);
                WR.skip(40000);
            }
        }
        mcr = mc;
        mc -= step;
        double sign = 1;
        #ifndef PARALLEL
        fprintf(stderr,"start from %f\n",mc);
        #endif // PARALLEL
        while(mcr-mcl>0.001)
        {
            #pragma omp critical
            {
                fprintf(stderr,"task %d:%f<%f<%f\n",task,mcl,mc,mcr);
            }
            WR.reset_coeff(g(task),mc);
            WR.skip(40000);
            WR.sample_single(40000);
            Z = sign*WR.Ztopo*1.0/WR.count;
            WR.reset_statistic();
            WR.sample_single(40000);
            Z = max(Z,sign*WR.Ztopo*1.0/WR.count);
            if(Z<0)
            {
                if(sign>0)
                {
                    mcl = mc;
                    sign = -1.1;
                }
                else
                {
                    mcr = mc;
                    sign = +0.9;
                }
                step = step/2;
            }
            else
            {
                if(sign>0)
                {
                    mcr = mc;
                }
                else mcl = mc;
                if((mc<mcl)||(mc>mcr))
                {
                    mc = mcl + mcr;
                    break;
                }
            }
            mc -= sign *step;
        }
        #pragma omp critical
        {
            if(firstm)
            {
                fprintf(fp,"\n,");
            }
            firstm++;
            fprintf(fp,"(%d,%d,%d,%d,%d):(%f,%f,%f)",resultid,task,X,Y,N,g(task),(mcr+mcl)/2,mcr-mcl);
        }
    }
}
int main()
{
    using namespace std;
    pushversionforward();
    char filename [256];
    sprintf(filename,"../result/worm/mc_%d.py",WORM_VERSION);
    FILE *mcfp = fopen(filename,"wt");
    //FILE *mcfp = stderr;
    fprintf(mcfp,"#X Y N  m  g  Ztopo\n""observers = {");
    for(int i=0;i<10;i++)
    {
        getzxx<32,32,1>([](int i)->double{return 0.1*i;},12,i,mcfp);
        getzxx<32,32,2>([](int i)->double{return 0.1*i;},12,i,mcfp);
        getzxx<32,32,4>([](int i)->double{return 0.1*i;},12,i,mcfp);
        getzxx<32,32,8>([](int i)->double{return 0.1*i;},12,i,mcfp);
    }
    fprintf(mcfp,"}\n");
    fclose(mcfp);
    return 0;
}
