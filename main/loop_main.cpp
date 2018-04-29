#include "../version/loop_version.h"
#include "../common.h"
#include "../loop.h"
//#include "worm.h"
#define PARALLEL
void pushversionforward()
{
    FILE* fp = fopen("../version/loop_version.h","wt");
    fprintf(fp,"#define LOOP_REPRESENT_VERSION %d",LOOP_REPRESENT_VERSION+1);
    fclose(fp);
}
template<int X,int Y,int N>
void getCxChi(auto m,auto g,int cases,int resultid)
{
    #define RECORDTYPE std::tuple<int,double,double,double,double,int,double,double,double,double,double,double>
    #define RECORDTITLE "task g m Cx Chi count n0avr n1avr n2avr n0var n1var n2var"
    #define RECORDFORMAT "%d  %f  %f  %f  %f  %d  %f  %f  %f  %f  %f  %f  "
    fermionloop<X,Y,N> *LR[cases];
    unsigned int seeds[cases];
    RECORDTYPE *records[cases];
    for(int i=0;i<cases;i++)
    {
        seeds[i] = mtrand();
    }
    int top = 0;
    #ifdef PARALLEL
    #pragma omp parallel for
    #endif
    for(int i=0;i<cases;i++)
    {
        int task = 0;
        #pragma omp critical
        {
            task = top;
            top = top +1;
            fprintf(stderr,"launch task %d \n",task);
            fprintf(stderr,"with configure:(%f,%f)\n",g(task),m(task));
        }

        LR[task] = new fermionloop<X,Y,N>(g(task),m(task),seeds[task]);
        #pragma omp critical
        {
            fprintf(stderr,"task %d: prepare end \n",task);
        }
        LR[task]->zeroinitlattice();
        LR[task]->skip(10000);
        LR[task]->sample(50000,10);
        LR[task]->getMoreStatistic();
        #pragma omp critical
        {
            fprintf(stderr,"task %d: compute end \n",task);
            fprintf(stderr,"Accumulate to: %llu,%llu,%llu \n",LR[task]->nSum[0],LR[task]->nSum[1],LR[task]->nSum[2]);
            fprintf(stderr,"acceptRate:%f \n",(double)LR[task]->acceptCount*1.0/LR[task]->randomCount);
        }
        records[task] = new RECORDTYPE(task,g(task),m(task),LR[task]->Cx,LR[task]->Chi,LR[task]->count
            ,LR[task]->nAvr[0],LR[task]->nAvr[1],LR[task]->nAvr[2]
            ,LR[task]->nVar[0],LR[task]->nVar[1],LR[task]->nVar[2]);
    }
    fprintf(stderr,"start result processing \n");
    char filename[256];
    sprintf(filename,"../result/loop/Cxchi_%d_%d.txt",LOOP_REPRESENT_VERSION,resultid);
    fprintf(stderr,"result file: %s \n",filename);
    FILE* fp = fopen(filename,"wt");
    if(fp==nullptr)
    {
        fprintf(stderr,"Error:failed on open %s",filename);
        return ;
    }
    fprintf(fp,"Lattice (X,Y,N) = (%d,%d,%d)\n",X,Y,N);
    fprintf(fp,"Supporting multi-flavor \n");
    fprintf(fp,RECORDTITLE"\n");
    double mx[cases];
    double cx[cases];
    double chi[cases];
    for(int i=0;i<cases;i++)
    {
        fprintf(fp,RECORDFORMAT"\n"
            ,std::get<0>(*records[i]),std::get<1>(*records[i]),std::get<2>(*records[i])
            ,std::get<3>(*records[i]),std::get<4>(*records[i]),std::get<5>(*records[i])
            ,std::get<6>(*records[i]),std::get<7>(*records[i]),std::get<8>(*records[i])
            ,std::get<9>(*records[i]),std::get<10>(*records[i]),std::get<11>(*records[i])
                );
        mx[i]   = m(i);
        cx[i]  = std::get<3>(*records[i]);
        chi[i] = std::get<4>(*records[i]);
        delete records[i];
        delete LR[i];
    }
    fclose(fp);

    sprintf(filename,"../graph/loop/m-Cchi_%d_%d.asy",LOOP_REPRESENT_VERSION,resultid);
    FILE* asyfp = fopen(filename,"wt");
    if(asyfp==nullptr)
    {
        fprintf(stderr,"Error:failed on open %s",filename);
        return ;
    }
    fprintf(asyfp,"size(200);\n");
    fprintf(asyfp,"label(\"$C_\\chi $\",(0,0.11));\n");
    fprintf(asyfp,"label(\"$m$\",(0.5,0));\n");
    printASYaxis(asyfp,-0.5,0.5,-1.0,0.1);
    printASYdefcross(asyfp,0.01);
    printASYmarks<double,double>(asyfp,cases,mx,cx,"cross",false,"black");
    fclose(asyfp);

    sprintf(filename,"../graph/loop/m-chi_%d_%d.asy",LOOP_REPRESENT_VERSION,resultid);
    asyfp = fopen(filename,"wt");
    if(asyfp==nullptr)
    {
        fprintf(stderr,"Error:failed on open %s",filename);
        return ;
    }
    fprintf(asyfp,"size(200);\n");
    fprintf(asyfp,"label(\"$\\chi $\",(0,0.11));\n");
    fprintf(asyfp,"label(\"$m$\",(0.5,0));\n");
    printASYaxis(asyfp,-0.5,0.5,-1.0,0.1);
    printASYdefcross(asyfp,0.01);
    printASYmarks<double,double>(asyfp,cases,mx,chi,"cross",false,"black");
    fclose(asyfp);
    #undef RECORDTYPE
    #undef RECORDTITLE
    #undef RECORDFORMAT
}
int main()
{
    pushversionforward();
    getCxChi<3,3,1>([](int i)->double{return 0.3-0.02*i;},[](int i)->double{return 1.1;},21,0);
    //getCxChi<32,32,2>([](int i)->double{return 0.3-0.02*i;},[](int i)->double{return 1.1;},21,2);
    //getCxChi<128,128,1>([](int i)->double{return 0.3-0.02*i;},[](int i)->double{return 1.1;},21,1);
    return 0;
}
