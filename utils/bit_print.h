#ifndef UTILS_BIT_PRINT_H
#define UTILS_BIT_PRINT_H
template<int x,int y,bool end,typename T>
void printBits(T b)
{
    for(int i=0;i<x*y;i+=y)
    {
        for(int j=0;j<y;j++)
        {
            std::cout<<((b>>(i+j))&1);
        }
        if(end||!(i+y==x*y))std::cout<<std::endl;
    }
}
template<typename ... T>
void printContextLine(const unsigned char i,const unsigned int ctx1,T ... other)
{
    printContextLine(i,ctx1);
    printContextLine(i,other...);
}
void printContextLine(const unsigned char i,const unsigned int ctx)
{
    for(int j=0;j<9;j+=3)
        std::cout<<((ctx>>(i+j))&1);
    std::cout<<"  ";
}
template<bool tags,typename ... T>
void printContexts(const unsigned int ctx1,T ... other)
{
    if(tags)std::cout<<"----Context Begin----"<<std::endl;
    for(unsigned char i=0;i<3;i+=1)
    {
        printContextLine(i,ctx1,other...);
        std::cout<<std::endl;
    }
    if(tags)std::cout<<"----Context  End ----"<<std::endl;
    return;
}
template<int x,int y>
void printLinkwithLoop(unsigned char link[][x][y],int loop[][(x+31)>>5][y],int f,int hx,int hy,int tx,int ty)
{
    int u,v;
    char buffer[x][y][3][4];
    for(int i=0;i<x;i++)
        for(int j=0;j<y;j++)
        {
            buffer[i][j][0][3]=0;
            buffer[i][j][1][3]=0;
            buffer[i][j][2][3]=0;
            if(((int)link[f][i][j]&1))buffer[i][j][1][2]='-';
            else buffer[i][j][1][2]=' ';
            if(((int)link[f][i][j]&2))buffer[i][j][0][1]='|';
            else buffer[i][j][0][1]=' ';
            if(((int)link[f][i][j]&4))buffer[i][j][1][0]='-';
            else buffer[i][j][1][0]=' ';
            if(((int)link[f][i][j]&8))buffer[i][j][2][1]='|';
            else buffer[i][j][2][1]=' ';
            if(hx==tx&&i==hx&&hy==ty&&j==hy)buffer[i][j][1][1] = 'X';
            else if(i==hx&&j==hy)buffer[i][j][1][1] = 'H';
            else if(i==tx&&j==ty)buffer[i][j][1][1] = 'T';
            else buffer[i][j][1][1] = 'E';
            u=i;v=j;
            if((loop[f][u>>5][v]>>(u&31))&1)buffer[i][j][0][0] = '1';
            else buffer[i][j][0][0] = '0';
            u=i+1;v=j;
            if(u>=x)u-=x;
            if((loop[f][u>>5][v]>>(u&31))&1)buffer[i][j][2][0] = '1';
            else buffer[i][j][2][0] = '0';
            u=i+1;v=j+1;
            if(u>=x)u-=x;
            if(v>=y)v-=y;
            if((loop[f][u>>5][v]>>(u&31))&1)buffer[i][j][2][2] = '1';
            else buffer[i][j][2][2] = '0';
            u=i;v=j+1;
            if(v>=y)v-=y;
            if((loop[f][u>>5][v]>>(u&31))&1)buffer[i][j][0][2] = '1';
            else buffer[i][j][0][2] = '0';
        }
    for(int i=0;i<x;i++)
    {
        for(int j=0;j<y;j++)
        {
            fprintf(stderr,"%s",buffer[i][j][0]);
        }
        fprintf(stderr,"\n");
        for(int j=0;j<y;j++)
        {
            fprintf(stderr,"%s",buffer[i][j][1]);
        }
        fprintf(stderr,"\n");
        for(int j=0;j<y;j++)
        {
            fprintf(stderr,"%s",buffer[i][j][2]);
        }
        fprintf(stderr,"\n");
    }
    return;
}
template<int x,int y>
void printLink(unsigned char link[][x][y],int f,int hx,int hy,int tx,int ty,int bx,int by,int nx,int ny)
{
    int u,v;
    char buffer[x][y][3][4];
    for(int i=0;i<x;i++)
        for(int j=0;j<y;j++)
        {
            buffer[i][j][0][3]=0;
            buffer[i][j][1][3]=0;
            buffer[i][j][2][3]=0;
            if(((int)link[f][i][j]&1))buffer[i][j][1][2]='-';
            else buffer[i][j][1][2]=' ';
            if(((int)link[f][i][j]&2))buffer[i][j][0][1]='|';
            else buffer[i][j][0][1]=' ';
            if(((int)link[f][i][j]&4))buffer[i][j][1][0]='-';
            else buffer[i][j][1][0]=' ';
            if(((int)link[f][i][j]&8))buffer[i][j][2][1]='|';
            else buffer[i][j][2][1]=' ';
            if(hx==tx&&i==hx&&hy==ty&&j==hy)buffer[i][j][1][1] = 'X';
            else if(i==hx&&j==hy)buffer[i][j][1][1] = 'H';
            else if(i==tx&&j==ty)buffer[i][j][1][1] = 'T';
            else if(i==bx&&j==by)buffer[i][j][1][1] = 'B';
            else if(i==nx&&j==ny)buffer[i][j][1][1] = 'N';
            else buffer[i][j][1][1] = 'E';
            buffer[i][j][0][0] = '0';
            buffer[i][j][2][0] = '0';
            buffer[i][j][2][2] = '0';
            buffer[i][j][0][2] = '0';
        }
    for(int i=0;i<x;i++)
    {
        for(int j=0;j<y;j++)
        {
            fprintf(stderr,"%s",buffer[i][j][0]);
        }
        fprintf(stderr,"\n");
        for(int j=0;j<y;j++)
        {
            fprintf(stderr,"%s",buffer[i][j][1]);
        }
        fprintf(stderr,"\n");
        for(int j=0;j<y;j++)
        {
            fprintf(stderr,"%s",buffer[i][j][2]);
        }
        fprintf(stderr,"\n");
    }
    return;
}
template<int x,int y>
void printLink(unsigned char link[][x][y],int f,int hx,int hy,int tx,int ty)
{
    int u,v;
    char buffer[x][y][3][4];
    for(int i=0;i<x;i++)
        for(int j=0;j<y;j++)
        {
            buffer[i][j][0][3]=0;
            buffer[i][j][1][3]=0;
            buffer[i][j][2][3]=0;
            if(((int)link[f][i][j]&1))buffer[i][j][1][2]='-';
            else buffer[i][j][1][2]=' ';
            if(((int)link[f][i][j]&2))buffer[i][j][0][1]='|';
            else buffer[i][j][0][1]=' ';
            if(((int)link[f][i][j]&4))buffer[i][j][1][0]='-';
            else buffer[i][j][1][0]=' ';
            if(((int)link[f][i][j]&8))buffer[i][j][2][1]='|';
            else buffer[i][j][2][1]=' ';
            if(hx==tx&&i==hx&&hy==ty&&j==hy)buffer[i][j][1][1] = 'X';
            else if(i==hx&&j==hy)buffer[i][j][1][1] = 'H';
            else if(i==tx&&j==ty)buffer[i][j][1][1] = 'T';
            else buffer[i][j][1][1] = 'E';
            buffer[i][j][0][0] = '0';
            buffer[i][j][2][0] = '0';
            buffer[i][j][2][2] = '0';
            buffer[i][j][0][2] = '0';
        }
    for(int i=0;i<x;i++)
    {
        for(int j=0;j<y;j++)
        {
            fprintf(stderr,"%s",buffer[i][j][0]);
        }
        fprintf(stderr,"\n");
        for(int j=0;j<y;j++)
        {
            fprintf(stderr,"%s",buffer[i][j][1]);
        }
        fprintf(stderr,"\n");
        for(int j=0;j<y;j++)
        {
            fprintf(stderr,"%s",buffer[i][j][2]);
        }
        fprintf(stderr,"\n");
    }
    return;
}
#endif // UTILS_BIT_PRINT_H

