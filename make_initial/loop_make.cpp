#include<cstdio>
void printOccupationArray(FILE* fp)
{
    unsigned char occupation[512];
    for(int i=0;i<512;i++)
    {
        char result = 0;
        {
            unsigned int u;
            u = (i   )&9;
            if(u==1||u==8)result|= 1;
            u = (i>>1)&9;
            if(u==1||u==8)result|= 3;
            u = (i>>2)&9;
            if(u==1||u==8)result|= 2;
            u = (i>>3)&9;
            if(u==1||u==8)result|= 4;
            u = (i>>4)&9;
            if(u==1||u==8)result|=12;
            u = (i>>5)&9;
            if(u==1||u==8)result|= 8;
            u = (i   )&3;
            if(u==1||u==2)result|= 1;
            u = (i>>1)&3;
            if(u==1||u==2)result|= 2;
            u = (i>>3)&3;
            if(u==1||u==2)result|= 5;
            u = (i>>4)&3;
            if(u==1||u==2)result|= 10;
            u = (i>>6)&3;
            if(u==1||u==2)result|= 4;
            u = (i>>7)&3;
            if(u==1||u==2)result|= 8;
        }
        occupation[i] = result;
    }
    fprintf(fp,"#define LOOP_OCCUPATION_ARRAY ");
    for(int i=0;i<512;i++)
    {
        if(i%8==0)fprintf(fp,"\\\n");
        if(i)fprintf(fp,",");
        else fprintf(fp,"{");
        fprintf(fp,"%d",occupation[i]);
    }
    fprintf(fp,"}\n");
}
void printCornerArray(FILE* fp)
{
    unsigned char corners[512];
    for(unsigned int i=0;i<512;i++)
    {
        char result = 0;
        {
            unsigned int v[4] = {(i&27),(i&54)>>1,(i&216)>>3,(i&432)>>4};
            for(int k=0;k<4;k++)
            {
                {
                    if(v[k]== 1||(v[k]^ 1)==27
                    || v[k]== 2||(v[k]^ 2)==27
                    || v[k]== 8||(v[k]^ 8)==27
                    || v[k]==16||(v[k]^16)==27)
                    {
                        result++;
                    }
                }
            }
        }
        corners[i] = result;
    }
    fprintf(fp,"#define LOOP_CORNER_ARRAY ");
    for(int i=0;i<512;i++)
    {
        if(i%8==0)fprintf(fp,"\\\n");
        if(i)fprintf(fp,",");
        else fprintf(fp,"{");
        fprintf(fp,"%d",(int)corners[i]);
    }
    fprintf(fp,"}\n");
}
void printForbiddenArray(FILE* fp)
{
    unsigned int forbidden[16];
    for(int i=0;i<16;i++)
    {
        unsigned int result = 0;
        for(int j=0;j<32;j++)
        {
            unsigned int u = (i<<5)+j;
            unsigned int v;
            v =  (u& 27)    ^16;
            if(v==10||v==17)result|=1<<j;
            v = ((u& 54)>>1)^8;
            if(v==10||v==17)result|=1<<j;
            v = ((u&216)>>3)^2;
            if(v==10||v==17)result|=1<<j;
            v = ((u&432)>>4)^1;
            if(v==10||v==17)result|=1<<j;
        }
        forbidden[i] = result;
    }
    fprintf(fp,"#define LOOP_FORBIDDEN_ARRAY ");
    for(int i=0;i<16;i++)
    {
        if(i%8==0)fprintf(fp,"\\\n");
        if(i)fprintf(fp,",");
        else fprintf(fp,"{");
        fprintf(fp,"%u ",(unsigned int)forbidden[i]);
    }
    fprintf(fp,"}\n");
}
template<int bits>
void printOnes(FILE* fp)
{
    unsigned char ones[1<<(bits+1)];
    ones[0] = 0;
    ones[1] = 1;
    for(int i=1;i<=bits;i++)
    {
        int upper = 1<<(i+1);
        int mask = (1<<i)-1;
        for(int j=1<<i;j<upper;j++)
        {
            ones[j] = ones[j&mask]+1;
        }
    }
    fprintf(fp,"#define ONES_ARRAY_%d ",bits);
    for(int i=0;i<(1<<(bits+1));i++)
    {
        if(i%8==0)fprintf(fp,"\\\n");
        if(i)fprintf(fp,",");
        else fprintf(fp,"{");
        fprintf(fp,"%d",(unsigned int)ones[i]);
    }
    fprintf(fp,"}\n");
}
void printnchangeArray(FILE* fp)
{
    unsigned char nchange[16][16];
    unsigned char nchangeid[16][16];
    for(int i=0;i<16;i++)
    {
        for(int j=0;j<16;j++)
        {
            int inc = j&(~i);
            int dec = i&(~j);
            int info = 0;
            int base = 1;
            int id = 0;
            for(int k=0;k<4;k++,base*=3)
            {
                if(inc&(1<<k))info |= 2<<(2*k);
                else if(dec&(1<<k))info |= 0;
                else info |= 1<<(2*k);
                if(inc&(1<<k))id += 2*base;
                else if(dec&(1<<k)) id += 0;
                else id += base;
            }
            nchange[i][j] = info;
            nchangeid[i][j] = id;
        }
    }
    fprintf(fp,"#define LOOP_NCHANGE_ARRAY ");
    for(int i=0;i<16;i++)
    {
        if(i)fprintf(fp,",");
        else fprintf(fp,"\\\n{");
        for(int j=0;j<16;j++)
        {
            if(j)fprintf(fp,",");
            else fprintf(fp,"{");
            fprintf(fp,"%d",nchange[i][j]);
        }
        fprintf(fp,"}\\\n");

    }
    fprintf(fp,"}\n");
    fprintf(fp,"#define LOOP_NCHANGE_ID_ARRAY ");
    for(int i=0;i<16;i++)
    {
        if(i)fprintf(fp,",");
        else fprintf(fp,"\\\n{");
        for(int j=0;j<16;j++)
        {
            if(j)fprintf(fp,",");
            else fprintf(fp,"{");
            fprintf(fp,"%d",nchangeid[i][j]);
        }
        fprintf(fp,"}\\\n");

    }
    fprintf(fp,"}\n");
}
int main()
{
    FILE* fp = fopen("../initial/loop_initial.h","wt");
    fprintf(fp,"#ifndef INITIAL_H\n");
    fprintf(fp,"#define INITIAL_H\n");
    printOccupationArray(fp);
    printCornerArray(fp);
    printForbiddenArray(fp);
    printOnes<8>(fp);
    printnchangeArray(fp);
    //fprintf(fp,"#define SPIN2LINK_ARRAY {0,6,3,5,12,10,15,9,9,15,10,12,5,3,6,0}");
    fprintf(fp,"#endif");
    fclose(fp);
    return 0;
}
