#include<cstdio>
void worm_realloc_break(FILE* fp)
{
    unsigned char dopt[2][16] = {};
    for(int di=0;di<4;di++)
    {
        for(int dj=di+1;dj<4;dj++)
        {
            dopt[0][(1<<di)|(1<<dj)] = di;
            dopt[1][(1<<di)|(1<<dj)] = dj;
        }
    }
    fprintf(fp,"#define REALLOC_BREAK_OPT_ARRAY \\\n");
    fprintf(fp,"{{%d",dopt[0][0]);
    for(int i=1;i<16;i++)
    {
        fprintf(fp,",%d",dopt[0][i]);
    }
    fprintf(fp,"}\\\n,{%d",dopt[1][0]);
    for(int i=1;i<16;i++)
    {
        fprintf(fp,",%d",dopt[1][i]);
    }
    fprintf(fp,"}}\n");
}
void worm_break_phase(FILE* fp)
{
    unsigned char phase[4][16][16] = {};
    phase[0][3][ 5] = 5;
    phase[0][3][ 6] = 6;
    phase[0][3][12] = 4;
    phase[0][5][ 5] = 4;
    phase[0][5][ 6] = 5;
    phase[0][5][12] = 3;
    phase[0][9][ 5] = 3;
    phase[0][9][ 6] = 4;
    phase[0][9][12] = 2;
    for(int d = 1;d<4;d++)
    {
        for(int i=0;i<16;i++)
        {
            for(int j=0;j<16;j++)
            {
                phase[d][((i<<d)|(i>>(4-d)))&0xf][((j<<d)|(j>>(4-d)))&0xf] = phase[0][i][j];
            }
        }
    }
    fprintf(fp,"#define WORM_BREAK_PHASE_ARRAY \\\n");
    fprintf(fp,"{");
    for(int d = 0;d<4;d++)
    {
        if(d)fprintf(fp,",");
        fprintf(fp,"{");
        for(int i=0;i<16;i++)
        {
            if(i)fprintf(fp,",");
            fprintf(fp,"{%d",phase[d][i][0]);
            for(int j=1;j<16;j++)
            {
                fprintf(fp,",%d",phase[d][i][j]);
            }
            fprintf(fp,"}\\\n");
        }
        fprintf(fp,"}");
    }
    fprintf(fp,"}\n");
}
void worm_turn_phase(FILE* fp)
{
    unsigned char phase[4][16] = {};
    phase[0][1] = 0;
    phase[0][2] = 1;
    phase[0][4] = 0;
    phase[0][8] = 7;
    for(int d = 1;d<4;d++)
    {
        for(int i=0;i<16;i++)
        {
            phase[d][((i<<d)|(i>>(4-d)))&0xf] = phase[0][i];
        }
    }
    phase[0][ 6] = (-phase[2][2])&7;
    phase[0][12] = (-phase[2][8])&7;
    for(int d = 1;d<4;d++)
    {
        for(int i=0;i<16;i++)
        {
            phase[d][((i<<d)|(i>>(4-d)))&0xf] = phase[0][i];
        }
    }
    fprintf(fp,"#define WORM_TURN_PHASE_ARRAY \\\n");
    for(int d = 0;d<4;d++)
    {
        if(d)fprintf(fp,",");
        else fprintf(fp,"{");
        for(int i=0;i<16;i++)
        {
            if(i)fprintf(fp,",");
            else fprintf(fp,"{");
            fprintf(fp,"%d",phase[d][i]);
        }
        fprintf(fp,"}\\\n");
    }
    fprintf(fp,"}\n");
}
void print_worm_arbiter(FILE* fp)
{
    int size = 4;
    unsigned char arbiter[size][1<<size] = {};
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<(1<<size);j++)
        {
            int num = i+1;
            for(int k=0;k<size;k++)
            {
                if(j&(1<<k))num--;
                if(!num)
                {
                    arbiter[i][j] = k;
                    break;
                }
            }
        }
    }
    fprintf(fp,"#define WORM_ARBITER_ARRAY_%d \\\n",size);
    for(int i=0;i<size;i++)
    {
        if(i)fprintf(fp,",");
        else fprintf(fp,"{");
        for(int j=0;j<(1<<size);j++)
        {
            if(j%8==0)fprintf(fp,"\\\n");
            if(j)fprintf(fp,",");
            else fprintf(fp,"{");
            fprintf(fp,"%d",arbiter[i][j]);
        }
        fprintf(fp,"}\\\n");
    }
    fprintf(fp,"}\n");
}
int main()
{
    FILE* fp = fopen("../initial/worm_initial.h","wt");
    fprintf(fp,"#ifndef WORM_INITIAL_H\n");
    fprintf(fp,"#define WORM_INITIAL_H\n");
    fprintf(fp,"#define LOOP2LINK_ARRAY {0,6,3,5,12,10,15,9,9,15,10,12,5,3,6,0}\n");
    fprintf(fp,"#define LINK2DIRECT_ARRAY {4,0,1,4,2,4,4,4,3,4,4,4,4,4,4,4}\n");
    worm_realloc_break(fp);
    worm_break_phase(fp);
    worm_turn_phase(fp);
    print_worm_arbiter(fp);
    fprintf(fp,"#endif //WORM_INITIAL_H");
    fclose(fp);
    return 0;
}
