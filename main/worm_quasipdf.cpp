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
void get_quasipdf(auto m,auto g,auto p,auto t,auto r,
                  int cases_m,int cases_g,int cases_p,int cases_t,int repeats,int runs,
                int resultid,FILE* cfp)
{
    const int cases = cases_m*cases_g;
    static const int length = X>Y?X:Y;
    static const int max_threads = 4;
    assert(omp_get_max_threads()<=max_threads);
    std::mt19937 *rand_thread[max_threads];
    fermionworm<X,Y,N> *WR[max_threads];
    unsigned long long (*quasipdf)[max_threads][2][length][length][length];
    quasipdf = (decltype(quasipdf))malloc(max_threads*2*length*length*length*sizeof(unsigned long long));
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
        for(int j=0;j<length;j++)
        for(int k=0;k<length;k++)
        for(int u=0;u<length;u++)
        {
            (*quasipdf)[i][0][j][k][u] = 0;
            (*quasipdf)[i][1][j][k][u] = 0;
        }
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
        WR[thread_id]->debug_quasi_pdf = false;
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
            WR[thread_id] -> skip(2000);
            WR[thread_id] -> reset_statistic();
            for(int ii=0;ii<runs;ii++)
            {

                if(runs%2==2)
                {
                    fprintf(stderr,"##");
                    //system("pause");
                }
                //if(ii>=100)WR[thread_id]->debug_quasi_pdf = true;
                WR[thread_id] -> record();
                WR[thread_id] -> psi_reallocate();
                if(false)system("pause");
                WR[thread_id] -> quasipdf_reallocate((*quasipdf)[thread_id]);
            }
            timer[thread_id].end();
            #pragma omp critical
            {
                if(repeat%1==0)
                {
                    fprintf(stderr,"Task %d runs %d, accept rate:%f\n",task,repeat,
                        WR[thread_id]->acceptCount*1.0/WR[thread_id]->randomCount);
                }
                double qpdf[length][length];
                double qpdf_pz[length][length];
                WR[thread_id]->summary(1);
                for(int task_t = 0;task_t<cases_t;task_t++)
                {
                    for(int z=0;z<X;z++)
                    {
                        for(int pz=0;pz<X;pz++)
                        {
                            qpdf[z][pz] =1.0*(*quasipdf)[thread_id][0][z][pz][task_t];
                            qpdf[z][pz]-=1.0*(*quasipdf)[thread_id][1][z][pz][task_t];
                        }
                        if((1<<log2_upper(Y))==Y)fft_real_even(qpdf[z],qpdf_pz[z],log2_upper(Y),false);
                        else
                        {
                            constexpr double pi = 4*atan(1);
                            for(int p=0;p<Y;p++)
                            {
                                qpdf_pz[z][p] = 0;
                                for(int v=0;v<Y;v++)
                                {
                                    qpdf_pz[z][p] += qpdf[z][v]*cos(p*v*2*pi/Y);
                                }
                            }
                        }
                    }
                    for(int task_p = 0;task_p<cases_p;task_p++)
                    {
                        if(firstc)fprintf(cfp,"\n,");
                        firstc ++;
                        fprintf(cfp,"(%d,%d,%d,%d,%d,%d,%d,\"gamma5_pdfs\",%d,%d):",X,Y,N,
                            resultid,task_m,task_g,repeat,task_p,task_t);
                        fprintf(cfp,"(%d,%d,%lld,%f,%f,%d,%d,%e,[",
                                timer[thread_id].s,timer[thread_id].ms,WR[thread_id]->count,
                                WR[thread_id]->m,WR[thread_id]->g,
                                t(task_t),p(task_p),
                                WR[thread_id]->G[1][task_t][task_p]);
                        for(int z=0;z<length;z++)
                        {
                            if(z)fprintf(cfp,",");
                            fprintf(cfp,"%e",qpdf_pz[z][p(task_p)]);
                        }
                        fprintf(cfp,"])");
                    }
                }
            }
            for(int j=0;j<length;j++)
            for(int k=0;k<length;k++)
            for(int u=0;u<length;u++)
            {
                (*quasipdf)[thread_id][0][j][k][u] = 0;
                (*quasipdf)[thread_id][1][j][k][u] = 0;
            }
        }
        #pragma omp critical
        {
            fprintf(stderr,"Task %d end\n",task);
        }
    }
    for(int i=0;i<max_threads;i++)
    {
        if(WR[i]!=nullptr)delete WR[i];
        if(rand_thread[i]!=nullptr)delete rand_thread[i];
    }
    delete quasipdf;
}
int main()
{
    pushversionforward();
    char filename [256];
    sprintf(filename,"../result/worm/worm_quasi_pdf_%d.py",WORM_VERSION);
    FILE *cfp = fopen(filename,"wt");
    fprintf(cfp,"pdfs = {");
    for(int r=0;r<=0;r++)
    get_quasipdf<32,32,1>(
        [](int i)->double{return 0.2-0.1*i;},
        [](int i)->double{return 1.1;},
        [](int i)->int{return 4*i;},
        [](int i)->int{return 4*i;},
        [r](double g,double m)->double{return 1.0;},
        1,1,5,5,30,1000000,r,cfp);
    fprintf(cfp,"}\n");
    fclose(cfp);
    return 0;
}
