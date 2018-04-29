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
void getzxx(auto m,auto g,int cases_m,int cases_g,int resultid,FILE *fp)
{
    int topt = 0;
    unsigned int seed[omp_get_max_threads()];
    std::mt19937 *rand_thread[omp_get_max_threads()];
    unsigned int max_threads = omp_get_max_threads();
    for(int i=0;i<max_threads;i++)
    {
        seed[i] = mtrand();
    }
    for(int i=0;i<max_threads;i++)
    {
        rand_thread[i] = nullptr;
    }
    fprintf(stderr,"launch tasks with maximum %d threads.\n",max_threads);
    int cases = cases_g*cases_m;
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
        int task_m = task/cases_g;
        int task_g = task%cases_g;
        fermionworm<X,Y,N> WR(g(task_g),m(task_m),(*rand_thread[thread_id])());
        WR.zeroinitlattice();
        WR.skip(30000);
        const int jackknife_groups = 800;
        const int jackknife_elems  =1250;
        long long Ztopo[jackknife_groups];
        long long Ztopo_temp[jackknife_groups];
        long long Ztopo_jacknife[jackknife_groups];
        double Z_avr;
        double Z_err;
        double Z_max;
        double Z_min;
        WR.skip(20000);
        for(int j=0;j<jackknife_groups;j++)
        {
            WR.sample_single(jackknife_elems,false);
            Ztopo[j] = WR.Ztopo;
            WR.count = 0;
            WR.Ztopo = 0;
        }
        Ztopo_temp[0] = 0;
        Ztopo_jacknife[jackknife_groups-1] = 0;
        for(int j=1;j<jackknife_groups;j++)
        {
            Ztopo_temp[j] = Ztopo_temp[j-1]+Ztopo[j-1];
            Ztopo_jacknife[jackknife_groups-1-j] = Ztopo_jacknife[jackknife_groups-j]+Ztopo[jackknife_groups-j];
        }
        for(int j=0;j<jackknife_groups;j++)
        {
            Ztopo_jacknife[j] += Ztopo_temp[j];
        }
        sort(Ztopo_jacknife,Ztopo_jacknife+jackknife_groups);
        Z_min = Ztopo_jacknife[(jackknife_groups* 4)/25]*1.0/((jackknife_groups-1)*jackknife_elems);
        Z_max = Ztopo_jacknife[(jackknife_groups*21)/25]*1.0/((jackknife_groups-1)*jackknife_elems);
        Z_avr = 0.0;
        Z_err = 0.0;
        for(int j=0;j<jackknife_groups;j++)
        {
            Z_avr+= Ztopo_jacknife[j];
        }
        Z_avr/= jackknife_groups;
        for(int j=0;j<200;j++)
        {
            Z_err += (Z_avr-Ztopo_jacknife[j])*(Z_avr-Ztopo_jacknife[j]);
        }
        Z_err*= (jackknife_groups-1)*1.0/jackknife_groups;
        Z_err = sqrt(Z_err);
        Z_avr/=(jackknife_groups-1)*jackknife_elems;
        Z_err/=(jackknife_groups-1)*jackknife_elems;
        #pragma omp critical
        {
            fprintf(stderr,"model<%d,%d,%d> task %d:(%f,%f):(%f<Z<%f)[avr:%e,err:%e]\n"
                    ,X,Y,N,task,g(task_g),m(task_m),Z_min,Z_max,Z_avr,Z_err);
            if(firstm)
            {
                fprintf(fp,"\n,");
            }
            firstm++;
            fprintf(fp,"(%d,%d,%d,%d,%d,%d):(%f,%f,%e,%e)",resultid,task_g,X,Y,N,task_m,g(task_g),m(task_m),Z_avr,Z_err);
        }
    }
}
int main()
{
    using namespace std;
    pushversionforward();
    char filename [256];
    sprintf(filename,"../result/worm/Ztopo_%d.py",WORM_VERSION);
    FILE *mcfp = fopen(filename,"wt");
    //FILE *mcfp = stderr;
    fprintf(mcfp,"#X Y N  m  g  Ztopo\n""observers = {");
    mtrand.seed(mtrand());

    for(int i=0;i<1;i++)
    {
        getzxx< 8, 8,1>([](int i)->double{return 0.2-0.02*i;},[](int i)->double{return 0.6+0.2*i;},50,1,i,mcfp);
        getzxx<16,16,1>([](int i)->double{return 0.2-0.02*i;},[](int i)->double{return 0.6+0.2*i;},50,1,i,mcfp);
        getzxx<32,32,1>([](int i)->double{return 0.2-0.02*i;},[](int i)->double{return 0.6+0.2*i;},50,1,i,mcfp);
        getzxx< 8, 8,2>([](int i)->double{return 0.2-0.02*i;},[](int i)->double{return 0.6+0.2*i;},50,1,i,mcfp);
        getzxx<16,16,2>([](int i)->double{return 0.2-0.02*i;},[](int i)->double{return 0.6+0.2*i;},50,1,i,mcfp);
        getzxx<32,32,2>([](int i)->double{return 0.2-0.02*i;},[](int i)->double{return 0.6+0.2*i;},50,1,i,mcfp);
    }
    fprintf(mcfp,"}\n");
    fclose(mcfp);
    return 0;
}
