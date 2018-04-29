#include "../version/worm_version.h"
#include "../common.h"
#include "../loop.h"
#include "../worm.h"
//#define PARALLEL
int firstc= 0;
int firstm= 0;
void pushversionforward()
{
    FILE* fp = fopen("../version/worm_version.h","wt");
    fprintf(fp,"#define WORM_VERSION %d",WORM_VERSION+1);
    fclose(fp);
}
template<int logX,int logY,int N>
void getzxxpsi(auto m,auto g,int cases,int resultid)
{
    static_assert(logX>1,"should:X>2");
    static_assert(logY>1,"should:Y>2");
    static const int X = 1<<logX;
    static const int Y = 1<<logY;
    fermionworm<X,Y,N> *WR[cases];
    static const int repeats = 16;
    float format[4]= {Y,1,X,repeats};
    float indexs[X];
    for(int i=0;i<X;i++)
    {
        indexs[i] = i;
    }
    unsigned int seeds[cases];
    for(int i=0;i<cases;i++)
    {
        seeds[i] = mtrand();
    }
    int top = 0;
    fprintf(stderr,"enter worm algorithm\n");
    #ifdef PARALLEL
    #pragma omp parallel for
    #endif // PARALLEL
    for(int i=0;i<cases;i++)
    {
        int task = 0;
        float results[repeats][X][Y];
        double buffer[Y];
        double result[Y];
        #pragma omp critical
        {
            task = top;
            top = top +1;
            fprintf(stderr,"launch task %d \n",task);
            fprintf(stderr,"#with configure:\n(g,m)=(%f,%f)\n",g(task),m(task));
        }
        WR[task] = new fermionworm<X,Y,N>(g(task),m(task),seeds[task]);
        WR[task] -> zeroinitlattice();
        WR[task] -> skip(20000);
        WR[task] -> smart_reweight();
        #pragma omp critical
        {
            fprintf(stderr,"initialize end\n");
        }
        for(int j=0;j<repeats;j++)
        {
            WR[task] -> reset_statistic();
            WR[task] -> sample_gamma5(400000);
            WR[task] -> summary(1);
            for(int u=0;u<X;u++)
            {
                for(int v=0;v<Y;v++)
                {
                    results[j][u][v] = WR[task]->G[1][u][v];
                }
            }
            #pragma omp critical
            {
                fprintf(stderr,"task %d, complete %d of %d runs\n",task,j+1,repeats);
            }
        }
        char filename[512];
        sprintf(filename,"../result/worm/propagator_%d_%d_%d.mbf",WORM_VERSION,resultid,task);
        FILE* mbf = fopen(filename,"wb");
        fwrite(format,sizeof(float),4,mbf);
        fwrite(indexs,sizeof(float),X,mbf);
        fwrite(results,sizeof(float),X*Y*repeats,mbf);
        fclose(mbf);
        delete WR[task];
    }
}
int main()
{
    pushversionforward();
    //getzxxpsi<5,5,1>([](int i)->double{return -0.0-0.02*i;},[](int i)->double{return 1.0;},1,777,cfp,mfp);
    //getzxxpsi<16,16,1>([](int i)->double{return 0.3-0.02*i;},[](int i)->double{return 1.1;},51,4,cfp,mfp);

    //getzxxpsi<32,32,1>([](int i)->double{return 0.3-0.01*i;},[](int i)->double{return 1.1;},90,0,cfp,mfp);
    //getzxxpsi<32,32,1>([](int i)->double{return 0.3-0.01*i;},[](int i)->double{return 1.1;},90,1,cfp,mfp);
    //getzxxpsi<32,32,1>([](int i)->double{return 0.3-0.01*i;},[](int i)->double{return 1.1;},90,2,cfp,mfp);

    getzxxpsi<7,7,1>([](int i)->double{return 0.36-0.060*i;},[](int i)->double{return 1.1;},16,3);
    //getzxxpsi<128,128,1>([](int i)->double{return 0.29-0.02*i;},[](int i)->double{return 1.1;},44,4,cfp,mfp);
    //getzxxpsi<128,128,1>([](int i)->double{return 0.29-0.02*i;},[](int i)->double{return 1.1;},44,5,cfp,mfp);
    return 0;
}
