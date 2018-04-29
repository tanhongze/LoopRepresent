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
template<int X,int Y,int N>
void getzxxpsi(auto m,auto g,int cases,int resultid)
{
    fermionworm<X,Y,N> *WR[cases];
    static const int repeats = 1;
    float format[4]= {Y,1,X,repeats};
    float indexs[X];
    for(int i=0;i<X;i++)
    {
        indexs[i] = i;
    }
    mtrand.seed(103);
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
        long long phase[8] = WORM_PHASE_ARRAY;
        #pragma omp critical
        {
            task = top;
            top = top +1;
            fprintf(stderr,"launch task %d \n",task);
            fprintf(stderr,"#with configure:\n(g,m)=(%f,%f)\n",g(task),m(task));
        }
        WR[task] = new fermionworm<X,Y,N>(g(task),m(task),seeds[task]);
        WR[task] -> zeroinitlattice();
        if(true)
        {
            #pragma omp critical
            {
                fprintf(stderr,"#model factors:\n f=[");
                for(int i=0;i<=2*N;i++)
                {
                    fprintf(stderr,"%f,",WR[task]->get_factor(i,true));
                }
                fprintf(stderr,"]\n");
            }
        }
        //WR[task] -> set_reweight(0,2.0);
        //WR[task] -> set_reweight(1,2.0);
        WR[task] -> skip(100000);
        WR[task] -> smart_reweight();
        #pragma omp critical
        {
            fprintf(stderr,"initialize end\n");
            fprintf(stderr,"reweight:%f %f\n",WR[task]->get_reweight(0),WR[task]->get_reweight(1));
        }
        WR[task] -> reset_statistic();
        WR[task] -> sample_gamma5(1000000);
        #pragma omp critical
        {
            fprintf(stderr,"100000 runs end\n");
        }
        WR[task] -> sample_gamma5(1000000);
        #pragma omp critical
        {
            fprintf(stderr,"100000 runs end\n");
        }
        WR[task]->summary(0);
        WR[task]->summary(1);
        /*
        for(int u=0;u<X;u++)
        {
            for(int v=0;v<Y;v++)
            {
                long long sum = 0;
                for(int s=0;s<8;s++)
                {
                    sum += WR[task]->pcount[0][s][u][v]*phase[s];
                }
                fprintf(stderr,"%lld ",sum);
            }
            fprintf(stderr,"\n");
        }
        */
        /*
        for(int u=0;u<X;u++)
        {
            for(int v=0;v<Y;v++)
            {
                long long sum = 0;
                for(int s=0;s<8;s++)
                {
                    sum += WR[task]->pcount[1][s][u][v]*phase[s];
                }
                fprintf(stderr,"%lld ",sum);
            }
            fprintf(stderr,"\n");
        }
        */
        for(int u=0;u<X;u++)
        {
            for(int v=0;v<Y;v++)
            {
                fprintf(stderr,"%.3e ",WR[task]->propagator[0][u][v]);
            }
            fprintf(stderr,"\n");
        }
        system("pause");
        char filename[512];
        sprintf(filename,"../result/worm/observers_%d_%d.py",WORM_VERSION,resultid);
        FILE* py = fopen(filename,"wb");
        #pragma omp critical
        {
            fprintf(py,"m%d = %f\n",task,m(task));
            fprintf(py,"g%d = %f\n",task,g(task));
            fprintf(py,"x%d = [",task);
            for(int u=0;u<X;u++)
            {
                if(u)fprintf(py,"\n,");
                fprintf(py,"[");
                for(int v=0;v<Y;v++)
                {
                    if(v)fprintf(py,",");
                    fprintf(py,"%d",u);
                }
                fprintf(py,"]");
            }
            fprintf(py,"]\n");
            fprintf(py,"y%d = [",task);
            for(int u=0;u<X;u++)
            {
                if(u)fprintf(py,"\n,");
                fprintf(py,"[");
                for(int v=0;v<Y;v++)
                {
                    if(v)fprintf(py,",");
                    fprintf(py,"%d",v);
                }
                fprintf(py,"]");
            }
            fprintf(py,"]\n");
            fprintf(py,"p%d = [",task);
            for(int u=0;u<X;u++)
            {
                if(u)fprintf(py,"\n,");
                fprintf(py,"[");
                for(int v=0;v<Y;v++)
                {
                    if(v)fprintf(py,",");
                    fprintf(py,"%e",WR[task]->propagator[0][u][v]);
                }
                fprintf(py,"]");
            }
            fprintf(py,"]\n");
            fprintf(py,"q%d = [",task);
            for(int u=0;u<X;u++)
            {
                if(u)fprintf(py,"\n,");
                fprintf(py,"[");
                for(int v=0;v<Y;v++)
                {
                    if(v)fprintf(py,",");
                    fprintf(py,"%e",WR[task]->propagator[1][u][v]);
                }
                fprintf(py,"]");
            }
            fprintf(py,"]\n");
        }
        fclose(py);
        delete WR[i];
    }
}
int main()
{
    pushversionforward();

    getzxxpsi<4,4,1>([](int i)->double{return 0.5;},[](int i)->double{return 0.3;},1,0);
    return 0;
}
