#include "../version/worm_version.h"
#include "../common.h"
#include "../loop.h"
#include "../worm.h"
#define PARALLEL
void pushversionforward()
{
    FILE* fp = fopen("../version/worm_version.h","wt");
    fprintf(fp,"#define WORM_VERSION %d",WORM_VERSION+1);
    fclose(fp);
}
template<int X,int Y,int N>
void getzxxpsi(auto m,auto g,int cases,int resultid)
{
    #define RECORD_TYPE   std::tuple<int,double,double,double,double,double>
    #define RECORD_TITLE  "#task  m  g  Ztopo  meff   fit_error"
    #define RECORD_FORMAT "%d:(%.3f,%.3f,%.3f,%.3f,%.3f)"
    static const int Xh = (X+31)>>5;
    static const int Yh = (Y+31)>>5;
    static const int Xl =  X&31;
    static const int Yl =  Y&31;
    unsigned int seeds[cases];
    fermionworm<X,Y,N> *WR[cases];
    char filename [256];
    RECORD_TYPE *records[cases];
    for(int i=0;i<cases;i++)
    {
        seeds[i] = mtrand();
    }
    int top = 0;
    int first = 0;
    sprintf(filename,"../result/worm/worm_propagator_%d_%d.py",WORM_VERSION,resultid);
    FILE * cfp = fopen(filename,"wt");
    fprintf(cfp,"propagators = {");
    #ifdef PARALLEL
    #pragma omp parallel for
    #endif // PARALLEL
    for(int i=0;i<cases;i++)
    {
        int task = 0;
        #pragma omp critical
        {
            task = top;
            top = top +1;
            fprintf(stderr,"launch task %d \n",task);
            fprintf(stderr,"#with configure:\n(g,m)=(%f,%f)\n",g(task),m(task));
        }
        WR[task] = new fermionworm<X,Y,N>(g(task),m(task),seeds[task]);
        WR[task] -> zeroinitlattice();
        WR[task] -> skip(0);
        WR[task]->reset_statistic();
        WR[task] -> sample(1000000);
        int len = WR[task] -> summary_psi();
        WR[task] -> fitting_meff(len,1);
        #pragma omp critical
        {
            if(first)
            {
                fprintf(cfp,"\n,");
            }
            first ++;
            fprintf(cfp,"(%d,%d,%d,\"psi\",%f,%f):",X,Y,N,WR[task]->m,WR[task]->g);
            WR[task]->save_G_to(len,cfp);
        }
        records[task] = new RECORD_TYPE(task,WR[task]->m,WR[task]->g,WR[task]->Ztopo*1.0/WR[task]->count,WR[task]->meff,WR[task]->fit_error);
        delete WR[task];
    }
    fprintf(cfp,"}\n");
    fclose(cfp);
    sprintf(filename,"../result/worm/Ztopo_meff_%d_%d.py",WORM_VERSION,resultid);
    FILE* fp = fopen(filename,"wt");
    fprintf(fp,"#Worm Lattice\n(X,Y,N) = (%d,%d,%d)\n",X,Y,N);
    fprintf(fp,RECORD_TITLE"\n");
    fprintf(fp,"observers = {");
    for(int i=0;i<cases;i++)
    {
        if(i)fprintf(fp,"\n,");
        fprintf(fp,RECORD_FORMAT
                ,std::get<0>(*records[i]),std::get<1>(*records[i]),std::get<2>(*records[i])
                ,std::get<3>(*records[i]),std::get<4>(*records[i]),std::get<5>(*records[i]));
    }
    fprintf(fp,"}");
    fclose(fp);
    for(int i=0;i<cases;i++)
    {
        delete records[i];
    }
}
int main()
{
    pushversionforward();
    getzxxpsi<3,3,1>([](int i)->double{return 0.0-0.02*i;},[](int i)->double{return 1.1;},1,0);
    //getzxxpsi<16,16,1>([](int i)->double{return 0.1-0.02*i;},[](int i)->double{return 0.9;},41,1);
    //getzxxpsi<32,32,1>([](int i)->double{return 0.1-0.02*i;},[](int i)->double{return 0.9;},41,2);
    return 0;
}
