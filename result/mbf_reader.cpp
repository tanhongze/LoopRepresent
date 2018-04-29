#include<cstdio>
#include<iostream>
using namespace std;
void print_mbf(int K,int V,int M,int N,float* index,float* propagator,FILE* target)
{

    fprintf(target,"%d\n%d\n%d\n%d\n",K,V,M,N);
    for(int i=0;i<M;i++)
    {
        fprintf(target,"%d ",i);
        for(int j=0;j<V;j++)
        {
            fprintf(target,"%f%c",index[i*V+j],(j==V-1)?'\n':' ');
        }
    }
    for(int n=0;n<N;n++)
    {
        for(int m=0;m<M;m++)
        {
            fprintf(target,"%d %d ",n,m);
            for(int k=0;k<K;k++)
            {
                fprintf(target,"%f%c",propagator[n*M*K+m*K+k],(k==K-1)?'\n':' ');
            }
        }
    }
}
int main()
{
    char filename[512];
    int end = 0;
    fprintf(stderr,
            "MBF Reader Lite\n"
            "Enter exit to exit program)\n");
    do
    {
        int i=0;
        fprintf(stderr,"Please input .mbf file:\n");
        while(i<511)
        {
            scanf("%c",&filename[i]);
            if(filename[i]==' ')break;
            if(filename[i]=='\n')break;
            if(filename[i]=='\t')break;
            i++;
        }
        if(!i){fprintf(stderr,"(Enter exit to exit program)\n");}
        filename[i] = 0;
        if(filename[0]=='e'
        && filename[1]=='x'
        && filename[2]=='i'
        && filename[3]=='t'
        && filename[4]== 0 )
        {
            break;
        }
        else if(
           filename[i-4]=='.'
        && filename[i-3]=='m'
        && filename[i-2]=='b'
        && filename[i-1]=='f'
        && filename[i  ]== 0 )
        {
            FILE* mbf = fopen(filename,"rb");
            if(mbf!=nullptr)
            {
                float format[4];
                fread(format,sizeof(float),4,mbf);
                int K,V,M,N;
                K = format[0]+0.5;
                V = format[1]+0.5;
                M = format[2]+0.5;
                N = format[3]+0.5;
                float index[M][V];
                fread(index,sizeof(float),M*V,mbf);
                float propagator[N][M][K];
                fread(propagator,sizeof(float),M*N*K,mbf);
                print_mbf(K,V,M,N,&index[0][0],&propagator[0][0][0],stdout);
                filename[0] = 0;
                fprintf(stderr,"output to .txt file ?\n(Yes/No)[Default No]\n%s\n",filename);
                filename[i-3]=='t';
                filename[i-2]=='x';
                filename[i-1]=='t';
                int j=0;
                char sel[16];
                while(j<15)
                {
                    scanf("%c",&sel[j]);
                    if(sel[j]==' ')break;
                    if(sel[j]=='\n')break;
                    if(sel[j]=='\t')break;
                    j++;
                }
                sel[j] = 0;
                if((sel[0]=='Y'||sel[0]=='y')
                && (sel[1]=='E'||sel[1]=='e')
                && (sel[2]=='S'||sel[2]=='s')
                &&  sel[3]== 0)
                {
                    FILE* sav = fopen(filename,"wt");
                    print_mbf(K,V,M,N,&index[0][0],&propagator[0][0][0],sav);
                }
            }
        }
    }
    while(!end);
    return 0;
}
