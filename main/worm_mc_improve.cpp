#include "../version/worm_version.h"
#include "../common.h"
#include "../loop.h"
#include "../worm.h"
#include "../utils/timer.h"
#include <algorithm>
using std::sort;
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
    unsigned int seed[cases];
    std::mt19937 *rand_thread[omp_get_max_threads()];
    unsigned int max_threads = omp_get_max_threads();
    for(int i=0;i<cases;i++)
    {
        seed[i] = mtrand();
    }
    for(int i=0;i<max_threads;i++)
    {
        rand_thread[i] = nullptr;
    }
    fprintf(stderr,"launch tasks with maximum %d threads.\n",max_threads);
    #ifdef PARALLEL
    #pragma omp parallel for
    #endif
    for(int i=0;i<cases;i++)
    {
        int thread_id = omp_get_thread_num();
        if(rand_thread[thread_id]==nullptr)
        {
            rand_thread[thread_id] = new std::mt19937(seed[thread_id]);
        }
        int task;
        #pragma omp critical
        {
            task = topt++;
        }
        const double mc_min = -1.0;
        double mc   = 0.1;
        double mc_last = 0.1;
        double mcr = 1.0,mcl = mc_min;
        double Zr = 1.0;
        double Zl = -0.25;
        int seedid = 0;
        fermionworm<X,Y,N> WR(g(task),mc,(*rand_thread[thread_id])());
        WR.zeroinitlattice();
        WR.skip(10000);
        double alpha = 0.5;
        bool flag = false;
        long long Ztopo[200];
        long long Ztopo_temp[200];
        long long Ztopo_jacknife[200];
        int records = 0;
        const int maxrecords = 10;
        double m_rec[maxrecords];
        double Z_avr[maxrecords];
        double Z_err[maxrecords];
        double Z_max[maxrecords];
        double Z_min[maxrecords];
        while(records<maxrecords)
        {
            WR.reset_coeff(g(task),mc);
            WR.skip(20000);
            for(int j=0;j<200;j++)
            {
                WR.sample_single(200,false);
                Ztopo[j] = WR.Ztopo;
                WR.count = 0;
                WR.Ztopo = 0;
            }
            Ztopo_temp[0] = 0;
            Ztopo_jacknife[199] = 0;
            for(int j=1;j<200;j++)
            {
                Ztopo_temp[j] = Ztopo_temp[j-1]+Ztopo[j-1];
                Ztopo_jacknife[199-j] = Ztopo_jacknife[200-j]+Ztopo[200-j];
            }
            for(int j=0;j<200;j++)
            {
                Ztopo_jacknife[j] += Ztopo_temp[j];
            }
            sort(Ztopo_jacknife,Ztopo_jacknife+200);
            double Zmin = Ztopo_jacknife[15] *1.0/(199*200);
            double Zmax = Ztopo_jacknife[184]*1.0/(199*200);
            m_rec[records] = mc;
            Z_avr[records] = 0.0;
            Z_err[records] = 0.0;
            for(int j=0;j<200;j++)
            {
                Z_avr[records]+= Ztopo_jacknife[j];
            }
            Z_avr[records]/= 200;
            for(int j=0;j<200;j++)
            {
                Z_err[records] += (Z_avr[records]-Ztopo_jacknife[j])*(Z_avr[records]-Ztopo_jacknife[j]);
            }
            Z_err[records]*= 199*1.0/200;
            Z_err[records] = sqrt(Z_err[records]);
            Z_avr[records]/=199*200;
            Z_err[records]/=199*200;
            Z_max[records] = Zmax;
            Z_min[records] = Zmin;
            #pragma omp critical
            {
                if(false)
                for(int j=0;j<200;j++)
                {
                    fprintf(stderr,"%lld ",Ztopo_jacknife[j]);
                }
                fprintf(stderr,"model<%d,%d,%d> task %d:%f<%f<%f(%f<Z<%f)[err:%e]\n",X,Y,N,task,mcl,mc,mcr,Zmin,Zmax,Z_err[records]);
                if((Zmin<0&&Zmax>0)||(mcr<=mcl))
                {
                    fprintf(stderr,"get records %d\n",records);
                }
                if(false)system("pause");
            }
            if(Zmax<0)
            {
                mcl = mc;
                if(mcr==mc)
                {
                    mcr = mc - Zmin*alpha;
                }
                Zl = Zmin;
                mc -= Zmin*alpha;
                if(flag == false)
                {
                    alpha = (alpha+(mcr-mcl>0?(mcr-mcl)/(Zr-Zl):0.00))*0.5;
                    flag = true;
                }
            }
            else if(Zmin>0)
            {
                mcr = mc;
                if(mcl==mc)
                {
                    mcl = mc - Zmax*alpha;
                }
                Zr = Zmax;
                mc -= Zmax*alpha;
                if(flag == true)
                {
                    alpha = (alpha+(mcr-mcl>0?(mcr-mcl)/(Zr-Zl):0.00))*0.5;
                    flag = false;
                }
            }
            else //if((Zmin<=0&&0<=Zmax)||(mcr<=mcl))
            {
                double k = ((*rand_thread[thread_id])()*0.5)/(1u<<31);
                double step = (1.0+0.2*k)*Zmax+(1.0-0.2*k)*Zmin;
                if(step>0)flag = false;
                else flag = true;
                mc -= alpha*step;
            }
            if((Zmin<=0&&0<=Zmax)||(mcr<=mcl))
            {
                records++;
            }
            if(mc>mcr)
            {
                mc = mcr;
            }
            if(mc<mcl)
            {
                mc = mcl;
            }
        }
        #pragma omp critical
        {
            for(int j=0;j<records;j++)
            {
                if(firstm)
                {
                    fprintf(fp,"\n,");
                }
                firstm++;
                fprintf(fp,"(%d,%d,%d,%d,%d,%d):(%f,%f,%e,%e,%e,%e)",resultid,task,X,Y,N,j,g(task),m_rec[j],Z_avr[j],Z_err[j],Z_max[j],Z_min[j]);
            }
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
    mtrand.seed(mtrand());
    for(int i=0;i<1;i++)
    {
        //getzxx<32,32,1>([](int i)->double{return 0.1*i;},12,i,mcfp);
        getzxx<32,32,4>([](int i)->double{return 0.1*i;},12,i,mcfp);
        //getzxx<32,32,4>([](int i)->double{return 0.1*i;},12,i,mcfp);
        //getzxx<32,32,8>([](int i)->double{return 0.1*i;},12,i,mcfp);
    }
    fprintf(mcfp,"}\n");
    fclose(mcfp);
    return 0;
}
