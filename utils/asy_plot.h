#ifndef UTILS_ASY_H
#define UTILS_ASY_H
void printASYaxis(FILE* asyfile,double x1,double x2,double y1,double y2)
{
    fprintf(asyfile,"draw(");
    fprintf(asyfile,"(%.3f,0)--(%.3f,0)",x1,x2);
    fprintf(asyfile,",EndArrow);");
    fprintf(asyfile,"draw(");
    fprintf(asyfile,"(0,%.3f)--(0,%.3f)",y1,y2);
    fprintf(asyfile,",EndArrow);\n");
}
template<typename T1,typename T2>
void printASYline(FILE* asyfile,int n,T1 x[],T2 y[],bool isloop,bool putarrow,bool colored,char* color)
{
    fprintf(asyfile,"draw(");
    for(int i=0;i<n;i++)
    {
        if(n%8==0)fprintf(asyfile,"\n");
        if(n)fprintf(asyfile,"--");
        fprintf(asyfile,"(%.3f,%.3f)",(double)x[i],(double)y[i]);
    }
    if(isloop)fprintf(asyfile,"--cycle");
    if(putarrow)fprintf(asyfile,",EndArrow");
    if(colored)fprintf(asyfile,"%c%s",putarrow?'+':',',color);
    fprintf(asyfile,");");
}
template<typename T1,typename T2>
void printASYmarks(FILE* asyfile,int n,T1 x[],T2 y[],char *name,bool colored,char* color)
{
    for(int i=0;i<n;i++)
    {
        fprintf(asyfile,"draw(shift(%.3f,%.3f)*%s%c%s);\n",(double)x[i],(double)y[i],name,colored?',':' ',colored?color:"");
    }
}
void printASYdefcross(FILE* asyfile,double r)
{
    fprintf(asyfile,"path[] cross = (-%.3f,0)--(%.3f,0)^^(0,-%.3f)--(0,%.3f);\n",r,r,r,r);
}
#endif // UTILS_ASY_H
