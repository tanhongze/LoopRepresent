#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cmath>
#include<omp.h>
#include "../common.h"
using namespace std;
namespace exhaustionLoopRepresent
{
    static const int X = 4;
    static const int Y = 4;
    // TODO: split even odd for X axis
    unsigned long long num[2][2][2][1<<(Y-1)][1<<(Y-1)][2][2][2*X*Y+1][X*Y+1][X*Y+1];
    template<int sizeX,int sizeY>
    class linkStateNumber
    {
      public:
        static const int X = (sizeX<sizeY)?sizeY:sizeX;
        static const int Y = (sizeX<sizeY)?sizeX:sizeY;
        static const int maxcorner = 2*X*Y;
        static const int maxpoints =   X*Y;
        static const int mask = (1<<Y)-1;
        static const int numindex = (1<<(Y-1));
        static const int numlink = (1<<Y);
        int search_ylink[X+1];
        int search_kinds[X+1];
        //static const bool disableFolding = false;
        int ones[1<<8];
        int evenodd[2][1<<7];
        int evenoddindex[1<<8];
        unsigned int shift[1<<(Y-1)];
        unsigned long long evenoddnumber[2][2][2][2][maxcorner+1][maxpoints+1][maxpoints+1];
        unsigned int singlelayer[2][2][maxpoints+1][maxpoints+1];
        unsigned int graphnumber[2][2];
        unsigned int occupy[2][2][1<<(Y-1)][1<<(Y-1)];
        unsigned int kinds[2][numindex][numindex];
        unsigned int oddx[2][2][1<<(Y-1)][1<<(Y-1)];
        #define CORNER(i,j) (ones[i^j])
        void initSingleLayer()
        {
            for(int ex=0;ex<2;ex++)
                for(int ey=0;ey<2;ey++)
                    graphnumber[ex][ey] = 0;
            for(int ex=0;ex<2;ex++)
                for(int ey=0;ey<2;ey++)
                    for(int corner = 0;corner<=maxpoints;corner++)
                        for(int n1 = 0;n1<=maxpoints;n1++)
                                singlelayer[ex][ey][corner][n1] = 0;
            for(int i=0;i<(1<<(Y-1));i++)
            {
                getsinglegraphnumber(0,i,i,X);
                getsinglegraphnumber(1,i,i,X);
            }
            for(int ex=0;ex<2;ex++)
                for(int ey=0;ey<2;ey++)
                    for(int corner = 0;corner<=maxpoints;corner++)
                        for(int n1 = 0;n1<=maxpoints;n1++)
                                graphnumber[ex][ey] += singlelayer[ex][ey][corner][n1];
        }
        void initBitCount()
        {
            ones[0] = 0;
            ones[1] = 1;
            for(int i=2;i<(1<<8);i++)
            {
                ones[i] = ones[i>>1]+(i%2);
            }
            int evenoddnum[2] = {};
            for(int i=0;i<(1<<8);i++)
            {
                evenodd[ones[i]%2][evenoddnum[ones[i]%2]++] = i;
            }
            for(int i=0;i<(1<<8);i++)
            {
                evenoddindex[evenodd[i%2][i>>1]] = i>>1;
            }
        }
        void initBitShift()
        {
            for(int i=0;i<(1<<(Y-1));i++)
            {
                unsigned differ = evenodd[0][i];
                shift[i] = 0;
                int haslink = 0;
                for(int j=0;j<Y;j++)
                {
                    if(differ&(1<<j))haslink ^= 1;
                    if(haslink)shift[i] |= 1<<j;
                }
            }
        }
        void initBitOccupy()
        {
            for(int i=0;i<numindex;i++)
            {
                for(int j=0;j<numindex;j++)
                {
                    for(int e=0;e<2;e++)
                    {
                        unsigned int up   = evenodd[e][i];
                        unsigned int down = evenodd[e][j];
                        unsigned int differ = up^down;
                        unsigned int index = evenoddindex[differ];
                        unsigned int shifted[2] = {shift[index],shift[index]^mask};
                        unsigned int transvers[2] = {};
                        unsigned int occupying[2] = {};
                        transvers[0] = (shifted[0]|(shifted[0]<<1)|(shifted[0]>>(Y-1)))&mask;
                        transvers[1] = (shifted[1]|(shifted[1]<<1)|(shifted[1]>>(Y-1)))&mask;
                        occupying[0] = (up&down)|transvers[0];
                        occupying[1] = (up&down)|transvers[1];
                        if(up&down)
                        {
                            if(up&down&transvers[0])
                            {
                                if(up&down&transvers[1])
                                {
                                    kinds[e][i][j] = 0;
                                    occupy[0][e][i][j] = 0;
                                    occupy[1][e][i][j] = 0;
                                    oddx[0][e][i][j] = 0;
                                    oddx[1][e][i][j] = 0;
                                }
                                else
                                {
                                    kinds[e][i][j] = 1;
                                    occupy[0][e][i][j] = occupying[1];
                                    occupy[1][e][i][j] = occupying[1];
                                    oddx[0][e][i][j] = shifted[1]&1;
                                    oddx[1][e][i][j] = shifted[1]&1;
                                }
                            }
                            else
                            {
                                kinds[e][i][j] = 1;
                                occupy[0][e][i][j] = occupying[0];
                                occupy[1][e][i][j] = occupying[0];
                                oddx[0][e][i][j] = shifted[0]&1;
                                oddx[1][e][i][j] = shifted[0]&1;
                            }
                        }
                        else
                        {
                            kinds[e][i][j] = 2;
                            occupy[0][e][i][j] = occupying[0];
                            occupy[1][e][i][j] = occupying[1];
                            oddx[0][e][i][j] = shifted[0]&1;
                            oddx[1][e][i][j] = shifted[1]&1;
                        }
                    }
                }
            }
        }
        linkStateNumber<sizeX,sizeY>()
        {
            fprintf(stderr,"initialize start\n");
            initBitCount();
            fprintf(stderr,"bit count  end\n");
            initBitShift();
            fprintf(stderr,"bit shift  end\n");
            initBitOccupy();
            fprintf(stderr,"bit occupy end\n");
            initSingleLayer();
            fprintf(stderr,"initialize end\n");
        }
        int searchoprint(int ex,int ey,int i,int j,int depth)
        {
            int num = 0;
            if(depth<=1)
            {
                if(kinds[ey][i][j]==0)return 0;
                search_ylink[depth] = j;
                int eo;
                for(int s=0;s<kinds[ey][i][j];s++)
                {
                    eo = 0;
                    for(int k=X;k>1;k--)
                    {
                        eo ^= oddx[search_kinds[k]][ey][search_ylink[k]][search_ylink[k-1]];
                    }
                    eo ^= oddx[s][ey][j][i];
                    if(eo!=ex)continue;
                    num+=1;
                    cerr<<"link:"<<endl;
                    for(int k=X;k>1;k--)
                    {
                        cerr<<search_kinds[k]<<" "<<search_ylink[k]<<endl;
                    }
                    cerr<<s<<" "<<j<<" "<<endl;
                    for(int k=X;k>1;k--)
                    {
                        printLink(ey,search_ylink[k],search_ylink[k-1],search_kinds[k]);
                    }
                    printLink(ey,j,i,s);
                    cerr<<endl;
                }
                return num;
            }
            else
            {
                search_ylink[depth] = j;
                for(int k=0;k<(1<<(Y-1));k++)
                {
                    if(kinds[ey][j][k]==0)continue;
                    for(int s = 0;s<kinds[ey][j][k];s++)
                    {
                        search_kinds[depth] = s;
                        num+=searchoprint(ex,ey,i,k,depth-1);
                    }
                }
                return num;
            }
        }
        void searchoprint(int e,int i,int j,int depth)
        {
            cerr<<"depth:"<<depth<<endl;
            if(depth<=1)
            {
                if(kinds[e][i][j]==0)return;
                search_ylink[depth] = j;
                for(int s=0;s<kinds[e][i][j];s++)
                {
                    cerr<<"link:"<<endl;
                    for(int k=X;k>1;k--)
                    {
                        printLink(e,search_ylink[k],search_ylink[k-1],search_kinds[k]);
                    }
                    printLink(e,j,i,s);
                    cerr<<endl;
                }
                system("pause");
                return;
            }
            else
            {
                search_ylink[depth] = j;
                for(int k=0;k<(1<<(Y-1));k++)
                {
                    if(kinds[e][j][k]==0)continue;
                    for(int s = 0;s<kinds[e][j][k];s++)
                    {
                        search_kinds[depth] = s;
                        printLink(e,j,k,s);
                        searchoprint(e,i,k,depth-1);
                    }
                }
                return;
            }
        }
        void searchoprint(int e)
        {
            for(int i=0;i<(1<<(Y-1));i++)
            {
                searchoprint(e,i,i,X);
            }
        }
        void searchoprint(int ex,int ey)
        {
            for(int i=0;i<(1<<(Y-1));i++)
            {
                cerr<<"get "<<searchoprint(ex,ey,i,i,X)<<" in search for "<<i<<endl;
            }
            if(true)system("pause");
        }
        void getsinglegraphnumber(int ey,int i,int j,int depth)
        {
            if(depth<=1)
            {
                if(kinds[ey][i][j]==0)return;
                int eo;
                int corner;
                int n1;
                search_ylink[depth] = j;
                int num = 0;
                for(int s=0;s<kinds[ey][i][j];s++)
                {
                    eo = 0;
                    corner = 0;
                    n1 = 0;
                    for(int k=X;k>1;k--)
                    {
                        eo ^= oddx[search_kinds[k]][ey][search_ylink[k]][search_ylink[k-1]];
                        corner += CORNER(evenodd[ey][search_ylink[k]],evenodd[ey][search_ylink[k-1]]);
                        n1 += ones[occupy[search_kinds[k]][ey][search_ylink[k]][search_ylink[k-1]]];
                    }
                    eo ^= oddx[s][ey][j][i];
                    corner += CORNER(evenodd[ey][j],evenodd[ey][i]);
                    n1 += ones[occupy[s][ey][j][i]];
                    singlelayer[eo][ey][corner][n1] ++;
                }
                return;
            }
            else
            {
                search_ylink[depth] = j;
                for(int k=0;k<(1<<(Y-1));k++)
                {
                    if(kinds[ey][j][k]==0)continue;
                    for(int s = 0;s<kinds[ey][j][k];s++)
                    {
                        search_kinds[depth] = s;
                        getsinglegraphnumber(ey,i,k,depth-1);
                    }
                }
                return;
            }
        }

        void count()
        {
            //init record board
            {
                for(int e1y=0;e1y<2;e1y++)
                  for(int e2y=0;e2y<2;e2y++)
                    for(int e1x=0;e1x<2;e1x++)
                      for(int e2x=0;e2x<2;e2x++)
                        for(int corner = 0;corner<=maxcorner;corner++)
                          for(int n1 = 0;n1<=maxpoints;n1++)
                            for(int n2 = 0;n2<=maxpoints;n2++)
                              evenoddnumber[e1y][e2y][e1x][e2x][corner][n1][n2] = 0;

            }
            // compute for each start
            {
                #pragma omp parallel for
                for(int Ey=0;Ey<4;Ey++)
                {
                    int e1y = Ey&1;
                    int e2y = (Ey>>1)&1;
                    for(int i=0;i<(1<<(Y-1));i++)
                      for(int j=0;j<(1<<(Y-1));j++)
                      {
                          // init the first line
                          {
                              for(int k=0;k<(1<<(Y-1));k++)
                                for(int l=0;l<(1<<(Y-1));l++)
                                  for(int Ex=0;Ex<4;Ex++)
                                  {
                                    int e1x = Ex&1;
                                    int e2x = (Ex>>1)&1;
                                    for(int corner=0;corner<=maxcorner;corner++)
                                      for(int n1=0;n1<=maxpoints;n1++)
                                        for(int n2=0;n2<=maxpoints;n2++)
                                        {
                                            num[0][e1y][e2y][k][l][e1x][e2x][corner][n1][n2] = 0;
                                        }
                                  }
                              num[0][e1y][e2y][i][j][0][0][0][0][0] = 1;
                          }
                          // go forward until last line
                          for(int x = 0;x<X;x++)
                          {
                              int cur  =  x   &1;
                              int next = (x+1)&1;
                              for(int k=0;k<(1<<(Y-1));k++)
                                for(int l=0;l<(1<<(Y-1));l++)
                                  for(int e1x = 0;e1x<2;e1x++)
                                    for(int e2x = 0;e2x<2;e2x++)
                                      for(int corner = 0;corner<=maxcorner;corner++)
                                        for(int n1 = 0;n1<=maxpoints;n1++)
                                          for(int n2 = 0;n2<=maxpoints;n2++)
                                            num[next][e1y][e2y][k][l][e1x][e2x][corner][n1][n2] = 0;
                              for(int k=0;k<(1<<(Y-1));k++)
                                for(int l=0;l<(1<<(Y-1));l++)
                                  for(int p=0;p<(1<<(Y-1));p++)
                                    for(int q=0;q<(1<<(Y-1));q++)
                                    {
                                        unsigned int rconnect = evenodd[e1y][p];
                                        unsigned int bconnect = evenodd[e2y][q];
                                        unsigned int rtop = evenodd[e1y][k];
                                        unsigned int btop = evenodd[e2y][l];
                                        int rS = kinds[e1y][p][k];
                                        int bS = kinds[e2y][q][l];
                                        for(int s1=0;s1<rS;s1++)
                                          for(int s2=0;s2<bS;s2++)
                                          {
                                              unsigned int roccupied = occupy[s1][e1y][p][k];
                                              unsigned int boccupied = occupy[s2][e2y][q][l];
                                              unsigned int dcorner = CORNER(rconnect,rtop)+ CORNER(bconnect,btop);
                                              unsigned int dn1 = ones[roccupied^boccupied];
                                              unsigned int dn2 = ones[roccupied&boccupied];
                                              unsigned int deo1 = oddx[s1][e1y][p][k];
                                              unsigned int deo2 = oddx[s2][e2y][q][l];
                                              for(int Ex = 0;Ex<4;Ex++)
                                              {
                                                  int e1x = Ex&1;
                                                  int e2x = (Ex>>1)&1;
                                                  for(int corner = dcorner;corner<=maxcorner;corner++)
                                                    for(int n1 = dn1;n1<=maxpoints;n1++)
                                                      for(int n2 = dn2;n2<=maxpoints;n2++)
                                                      {
                                                          num[next][e1y][e2y][k][l][e1x][e2x][corner][n1][n2] +=
                                                          num[cur ][e1y][e2y][p][q][e1x^deo1][e2x^deo2][corner-dcorner][n1-dn1][n2-dn2];
                                                      }
                                              }
                                          }
                                    }
                          }
                          //the last line return to the first
                          {
                              int last = X&1;
                              for(int corner = 0;corner<=maxcorner;corner++)
                                for(int n1 = 0;n1<=maxpoints;n1++)
                                  for(int n2 = 0;n2<=maxpoints;n2++)
                                    for(int e1x = 0;e1x<2;e1x++)
                                      for(int e2x = 0;e2x<2;e2x++)
                                        evenoddnumber[e1y][e2y][e1x][e2x][corner][n1][n2] += num[last][e1y][e2y][i][j][e1x][e2x][corner][n1][n2];
                          }
                      }
                }
            }
        }
        void printInitShift()
        {
            for(int i=0;i<(1<<(Y-1));i++)
            {
                cerr<<"state ";
                printBits<1,Y,false>((int)evenodd[0][i]);
                cerr<<":";
                printBits<1,Y,true>((int)shift[i]);
            }
        }
        void printLink()
        {
            int connectcount = 0;
            for(int i=0;i<(1<<(Y-1));i++)
            {
                for(int j=0;j<(1<<(Y-1));j++)
                {
                    for(int e=0;e<2;e++)
                    {
                        if(e)cerr<<"odd :";
                        else cerr<<"even:";
                        connectcount += kinds[e][i][j];
                        cerr<<"kinds:"<<kinds[e][i][j]<<endl;
                        for(int s=0;s<kinds[e][i][j];s++)
                        {
                            printLink(e,i,j,s);
                        }
                    }
                }
            }
        }
        void connectcount()
        {
            int connectcount[2][1<<(Y-1)][1<<(Y-1)] = {};
            int totalcount[X][2][1<<(Y-1)][1<<(Y-1)] = {};
            for(int i=0;i<(1<<(Y-1));i++)
            {
                for(int j=0;j<(1<<(Y-1));j++)
                {
                    for(int s = 0;s<kinds[0][i][j];s++)
                    {
                        connectcount[oddx[s][0][i][j]][i][j]++;
                    }
                }
            }
            for(int e=0;e<2;e++)
            {
                for(int i=0;i<(1<<(Y-1));i++)
                {
                    for(int j=0;j<(1<<(Y-1));j++)
                    {
                        cerr<<connectcount[e][i][j]<<" ";
                    }
                    cerr<<endl;
                }
                cerr<<endl;
            }
            for(int e1=0;e1<2;e1++)
              for(int u=0;u<(1<<(Y-1));u++)
                for(int v=0;v<(1<<(Y-1));v++)
                {
                    totalcount[0][e1][u][v] = connectcount[e1][u][v];
                }
            for(int i=1;i<X;i++)
              for(int e1=0;e1<2;e1++)
                for(int u=0;u<(1<<(Y-1));u++)
                  for(int v=0;v<(1<<(Y-1));v++)
                  {
                    totalcount[i][e1][u][v] = 0;
                    for(int e2=0;e2<2;e2++)
                      for(int w=0;w<(1<<(Y-1));w++)
                        totalcount[i][e1][u][v] += totalcount[i-1][e2][u][w] * connectcount[e1^e2][w][v];
                  }
            for(int u=0;u<(1<<(Y-1));u++)
                cerr<<totalcount[X-1][0][u][u]<<endl;
        }
        void printLink(int e,int i,int j,int s)
        {
            char buffer[3][3*10+1];
            buffer[0][3*Y]=0;
            buffer[1][3*Y]=0;
            buffer[2][3*Y]=0;
            unsigned int up   = evenodd[e][i];
            unsigned int down = evenodd[e][j];
            unsigned int mid = occupy[s][e][i][j];
            unsigned int left = shift[evenoddindex[up^down]];
            for(int u=0;u<Y;u++)
            {
                buffer[0][3*u] = ' ';
                if(up&(1<<u)) buffer[0][3*u+1] = '|';
                else buffer[0][3*u+1] = ' ';
                buffer[0][3*u+2] = ' ';
            }
            for(int u=0;u<Y;u++)
            {
                buffer[2][3*u] = ' ';
                if(down&(1<<u)) buffer[2][3*u+1] = '|';
                else buffer[2][3*u+1] = ' ';
                buffer[2][3*u+2] = ' ';
            }
            if(up&down&left)left^=mask;
            else if((up&down)==0&&s==1)left^=mask;
            //unsigned int right;
            for(int u=0;u<Y;u++)
            {
                if(left&(1<<u))buffer[1][3*u+2] = '-';
                else buffer[1][3*u+2] = ' ';
            }
            for(int u=0;u<Y;u++)
            {
                if(left&(1<<((u+Y-1)%Y)))buffer[1][3*u] = '-';
                else buffer[1][3*u] = ' ';
            }
            for(int u=0;u<Y;u++)
            {
                if(mid&(1<<u)) buffer[1][3*u+1] = '+';
                else buffer[1][3*u+1] = ' ';
            }
            cerr<<buffer[0]<<endl;
            cerr<<buffer[1]<<endl;
            cerr<<buffer[2]<<endl;
        }
        void outputSingle(ofstream& fileout)
        {
            fileout<<"statenum = "<<graphnumber[0][0]+graphnumber[0][1]+graphnumber[1][0]+graphnumber[1][1]<<endl;
            fileout<<"state_evenodd_num = ["
                  " [ "<<graphnumber[0][0]<<","<<graphnumber[0][1]<<"]"
                  ",[ "<<graphnumber[1][0]<<","<<graphnumber[1][1]<<"] ]"<<endl;
            fileout<<"states = {";
            bool first = true;
            for(int ex=0;ex<2;ex++)
            {
                for(int ey=0;ey<2;ey++)
                {
                    for(int corner = 0;corner<=maxpoints;corner++)
                    {
                        for(int n1 = 0;n1<=maxpoints;n1++)
                        {
                            if(singlelayer[ex][ey][corner][n1]==0)continue;
                            if(!first)
                            {
                                fileout<<",";
                            }
                            else first = false;
                            fileout<<"("<<ex<<","<<ey<<","<<corner<<","<<n1<<"):"<<singlelayer[ex][ey][corner][n1]<<endl;
                        }
                    }
                }
            }
            fileout<<"}"<<endl;
        }
        void saveSingleTo(char filename[])
        {
            ofstream  fileout(filename);
            if(fileout.bad())
            {
                cerr<<"Error on opening "<<filename<<endl;
                return ;
            }
            outputSingle(fileout);
        }
        void output(ofstream& fileout)
        {
            int total = 0;
            fileout<<"lenx = "<<X<<endl;
            fileout<<"leny = "<<Y<<endl;
            fileout<<"states = {";
            bool first = true;
            for(int E=0;E<16;E++)
            {
                int e1y =  E    &1;
                int e2y = (E>>1)&1;
                int e1x = (E>>2)&1;
                int e2x = (E>>3)&1;
                for(int cn=0;cn<=maxcorner;cn++)
                  for(int n1=0;n1<=maxpoints;n1++)
                    for(int n2=0;n2<=maxpoints;n2++)
                    {
                        if(evenoddnumber[e1y][e2y][e1x][e2x][cn][n1][n2]==0)continue;
                        if(!first)
                        {
                            fileout<<",";
                        }
                        else first = false;
                        fileout<<"("<<e1y<<","<<e2y<<","<<e1x<<","<<e2x<<","<<cn<<","<<n1<<","<<n2<<"):";
                        fileout<<evenoddnumber[e1y][e2y][e1x][e2x][cn][n1][n2]<<endl;
                        total += evenoddnumber[e1y][e2y][e1x][e2x][cn][n1][n2];
                    }
            }
            fileout<<"}"<<endl;
            fileout<<"#total state number:"<<total<<endl;
            double each = sqrt(total);
            int eachint = each + 0.5;
            fileout<<"#sqrt:"<<each<<endl;
            fileout<<"#recover:"<<eachint*eachint<<endl;
        }
        void saveTo(char filename[])
        {
            ofstream  fileout(filename);
            if(fileout.bad())
            {
                cerr<<"Error on opening "<<filename<<endl;
                return ;
            }
            output(fileout);
        }
        void print()
        {
            #define fileout cout
            int total = 0;
            fileout<<"lenx = "<<X<<endl;
            fileout<<"leny = "<<Y<<endl;
            fileout<<"states = {";
            bool first = true;
            for(int E=0;E<16;E++)
            {
                int e1y =  E    &1;
                int e2y = (E>>1)&1;
                int e1x = (E>>2)&1;
                int e2x = (E>>3)&1;
                for(int cn=0;cn<=maxcorner;cn++)
                  for(int n1=0;n1<=maxpoints;n1++)
                    for(int n2=0;n2<=maxpoints;n2++)
                    {
                        if(evenoddnumber[e1y][e2y][e1x][e2x][cn][n1][n2]==0)continue;
                        if(!first)
                        {
                            fileout<<",";
                        }
                        else first = false;
                        fileout<<"("<<e1y<<","<<e2y<<","<<e2x<<","<<e2x<<","<<cn<<","<<n1<<","<<n2<<"):";
                        fileout<<evenoddnumber[e1y][e2y][e1x][e2x][cn][n1][n2]<<endl;
                        total += evenoddnumber[e1y][e2y][e1x][e2x][cn][n1][n2];
                    }
            }
            fileout<<"}"<<endl;
            fileout<<"#total state number:"<<total<<endl;
            double each = sqrt(total);
            int eachint = each + 0.5;
            fileout<<"#sqrt:"<<each<<endl;
            fileout<<"#recover:"<<eachint*eachint<<endl;
            #undef fileout
        }
        #undef CORNER
    };
    linkStateNumber<5,5> counter;
}
using namespace exhaustionLoopRepresent;
int main()
{
    fprintf(stderr,"space need: %d\n",sizeof(linkStateNumber<5,5>));
    char filename[256];
    sprintf(filename,"../result/exhaustion/model_%d_%d_2/data.py",counter.X,counter.Y);
    if(false)
    {
        counter.printLink();
        counter.connectcount();
        //counter.printInitShift();
    }
    if(false)
    {
        counter.count();
        counter.saveTo(filename);
        counter.print();
    }
    if(true)
    {
        sprintf(filename,"../result/exhaustion/model_%d_%d_1/check_data.py",counter.X,counter.Y);
        counter.saveSingleTo(filename);
        //
    }
    if(false)
    {
        counter.searchoprint(0,0);
    }
    if(false)
    {
        cerr<<counter.kinds[0][2][2];
    }
    return 0;
}
