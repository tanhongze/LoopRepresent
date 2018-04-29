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
void getzxxpsi(auto m,auto g,auto p,auto r,int cases_m,int cases_g,int cases_p,int repeats,int runs,
               int resultid,FILE* cfp,FILE* mfp)
{
    const int cases = cases_m*cases_g;
    const unsigned int max_threads = omp_get_max_threads();
    std::mt19937 *rand_thread[max_threads];
    fermionworm<X,Y,N> *WR[max_threads];
    mytimer timer[max_threads]={};
    static const int Xh = (X+31)>>5;
    static const int Yh = (Y+31)>>5;
    static const int Xl =  X&31;
    static const int Yl =  Y&31;
    unsigned int seeds[max_threads];
    for(int i=0;i<max_threads;i++)
    {
        seeds[i] = mtrand();
        rand_thread[i] = nullptr;
        WR[i] = nullptr;
    }
    int top = 0;
    #ifdef PARALLEL
    #pragma omp parallel for
    #endif // PARALLEL
    for(int i=0;i<cases;i++)
    {
        const unsigned int thread_id = omp_get_thread_num();
        int task = 0;
        #pragma omp critical
        {
            task = top;
            top = top +1;
            fprintf(stderr,"launch task %d \n",task);
        }
        int task_g = task%cases_g;
        int task_m = task/cases_g;
        if(rand_thread[thread_id]==nullptr)
        {
            rand_thread[thread_id] = new std::mt19937(seeds[thread_id]);
        }
        std::mt19937 &randomizer = *rand_thread[thread_id];
        if(WR[thread_id]==nullptr)
        {
            WR[thread_id] = new fermionworm<X,Y,N>(g(task_g),m(task_m),randomizer());
        }
        else
        {
            WR[thread_id]->reset_coeff(g(task_g),m(task_m));
        }
        if(r(g(task_g),m(task_m))==0)
        {
            WR[thread_id] -> smart_reweight_by_table(0);
            WR[thread_id] -> smart_reweight_by_table(1);
        }
        else if(r(g(task_g),m(task_m))>0)
        {
            WR[thread_id]->set_reweight(0,sqrt(r(g(task_g),m(task_m))));
            WR[thread_id]->set_reweight(1,r(g(task_g),m(task_m)));
        }
        WR[thread_id] -> zeroinitlattice();
        WR[thread_id] -> skip(20000);
        #pragma omp critical
        {
            fprintf(stderr,"#with configure:\n(g,m)=(%f,%f)\n",g(task_g),m(task_m));
            fprintf(stderr,"re-weight:\n");
            fprintf(stderr,"%d:%f  ,",0,WR[thread_id] -> get_reweight(0));
            fprintf(stderr,"%d:%f \n",1,WR[thread_id] -> get_reweight(1));
        }
        for(int repeat=0;repeat<repeats;repeat++)
        {
            timer[thread_id].begin();
            WR[thread_id] -> skip(20000);
            WR[thread_id] -> reset_statistic();
            WR[thread_id] -> sample_gamma5(runs);
            WR[thread_id] -> summary(0);
            WR[thread_id] -> summary(1);
            timer[thread_id].end();
            #pragma omp critical
            {
                if(repeat%1==0)
                {
                    fprintf(stderr,"Task %d runs %d, accept rate:%f\n",task,repeat,
                        WR[thread_id]->acceptCount*1.0/WR[thread_id]->randomCount);
                }
                for(int task_p = 0;task_p<cases_p;task_p++)
                {
                    if(firstc)fprintf(cfp,"\n,");
                    firstc ++;
                    fprintf(cfp,"(%d,%d,%d,%d,%d,%d,%d,\"gamma5\",%d):",X,Y,N,resultid,task_m,task_g,repeat,task_p);
                    fprintf(cfp,"(%d,%d,%lld,%f,%f,%f,%f,",
                            timer[thread_id].s,timer[thread_id].ms,WR[thread_id]->count,
                            WR[thread_id]->get_reweight(1),WR[thread_id]->acceptCount*1.0/WR[thread_id]->randomCount,
                            WR[thread_id]->m,WR[thread_id]->g);
                    WR[thread_id]->save_G_to(1,p(task_p),cfp);
                    fprintf(cfp,")");
                    if(firstc)fprintf(cfp,"\n,");
                    firstc ++;
                    fprintf(cfp,"(%d,%d,%d,%d,%d,%d,%d,\"single\",%d):",X,Y,N,resultid,task_m,task_g,repeat,task_p);
                    fprintf(cfp,"(%d,%d,%lld,%f,%f,%f,%f,",
                            timer[thread_id].s,timer[thread_id].ms,WR[thread_id]->count,
                            WR[thread_id]->get_reweight(0),WR[thread_id]->acceptCount*1.0/WR[thread_id]->randomCount,
                            WR[thread_id]->m,WR[thread_id]->g);
                    WR[thread_id]->save_G_to(0,p(task_p),cfp);
                    fprintf(cfp,")");
                }
            }
            #pragma omp critical
            {
                if(firstm)fprintf(mfp,"\n,");
                firstm ++;
                fprintf(mfp,"(%d,%d,%d,%d,%d,%d):",X,Y,N,resultid,task_m,task_g,repeat);
                fprintf(mfp,"(%.4f,%.4f,%.4e,%.4e,%.4e,%.4e,%.4e)",
                    WR[thread_id]->m,WR[thread_id]->g,WR[thread_id]->Ztopo*1.0/WR[thread_id]->count);
            }
        }
        #pragma omp critical
        {
            fprintf(stderr,"Task %d end, accept rate:%f\n",task);
        }
    }
    for(int i=0;i<max_threads;i++)
    {
        if(WR[i]!=nullptr)delete WR[i];
        if(rand_thread[i]!=nullptr)delete rand_thread[i];
    }
}
int main()
{
    pushversionforward();
    char filename [256];
    sprintf(filename,"../result/worm/worm_propagator_%d.py",WORM_VERSION);
    FILE *cfp = fopen(filename,"wt");
    sprintf(filename,"../result/worm/mc_meff_%d.py",WORM_VERSION);
    FILE *mfp = fopen(filename,"wt");
    fprintf(cfp,"propagators = {");
    fprintf(mfp,"#X Y N  m  g  Ztopo  meff_b   fit_error_b  meff_f   fit_error_f\n"
                "observers = {");
    //getzxxpsi<5,5,1>([](int i)->double{return -0.0-0.02*i;},[](int i)->double{return 1.0;},1,777,cfp,mfp);
    //getzxxpsi<16,16,1>([](int i)->double{return 0.3-0.02*i;},[](int i)->double{return 1.1;},51,4,cfp,mfp);

    //getzxxpsi<32,32,1>([](int i)->double{return 0.3-0.01*i;},[](int i)->double{return 1.1;},90,0,cfp,mfp);
    //getzxxpsi<32,32,1>([](int i)->double{return 0.3-0.01*i;},[](int i)->double{return 1.1;},90,1,cfp,mfp);
    //getzxxpsi<32,32,1>([](int i)->double{return 0.3-0.01*i;},[](int i)->double{return 1.1;},90,2,cfp,mfp);
    #pragma omp parallel for
    for(int r=0;r<=3;r++)
    getzxxpsi<128,128,1>(
        [](int i)->double{return 0.38-0.04*i;},
        [](int i)->double{return 1.1;},
        [](int i)->int{return 4*i;},
        [r](double g,double m)->double{return exp(0.01+0.35*r);},
        1,1,1,40,2000000,r,cfp,mfp);
    if(false)
    getzxxpsi<128,128,1>(
        [](int i)->double{return 0.4-0.4*i;},
        [](int i)->double{return 1.1;},
        [](int i)->int{return 4*i;},
        [](double g,double m)->double{return 1.0;},
        3,1,4,4,2000000,0,cfp,mfp);
    //getzxxpsi<128,128,1>([](int i)->double{return 0.29-0.02*i;},[](int i)->double{return 1.1;},44,4,cfp,mfp);
    //getzxxpsi<128,128,1>([](int i)->double{return 0.29-0.02*i;},[](int i)->double{return 1.1;},44,5,cfp,mfp);
    /*
    getzxxpsi<32,32,1>([](int i)->double{return 0.3-0.01*i;},[](int i)->double{return 0.0;},90,6,cfp,mfp);
    getzxxpsi<32,32,2>([](int i)->double{return 0.3-0.01*i;},[](int i)->double{return 1.1;},90,7,cfp,mfp);
    getzxxpsi<32,32,3>([](int i)->double{return 0.3-0.01*i;},[](int i)->double{return 1.1;},90,8,cfp,mfp);
    getzxxpsi<128,128,1>([](int i)->double{return 0.30-0.01*i;},[](int i)->double{return 0.0;},90,9,cfp,mfp);
    getzxxpsi<128,128,2>([](int i)->double{return 0.30-0.01*i;},[](int i)->double{return 1.1;},90,10,cfp,mfp);
    getzxxpsi<128,128,3>([](int i)->double{return 0.30-0.01*i;},[](int i)->double{return 1.1;},90,11,cfp,mfp);
    getzxxpsi<128,128,1>([](int i)->double{return 0.0;},[](int i)->double{return 0.0+0.015*i;},81,12,cfp,mfp);
    getzxxpsi<128,128,2>([](int i)->double{return 0.0;},[](int i)->double{return 0.0+0.015*i;},81,13,cfp,mfp);
    getzxxpsi<128,128,3>([](int i)->double{return 0.0;},[](int i)->double{return 0.0+0.015*i;},81,14,cfp,mfp);
    */
    fprintf(cfp,"}\n");
    fprintf(mfp,"}\n");
    fclose(cfp);
    fclose(mfp);
    return 0;
}
