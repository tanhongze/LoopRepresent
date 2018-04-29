#ifndef UTILS_MOVE_H
#define UTILS_MOVE_H
#include "statistic.h"
template<int X,int D,int inc,int dec,typename T>
class movement
{
private:
    T move[X][D];
public:
    movement()
    {
        for(int d=0;d<D;d++)
            for(int i=0;i<X;i++)
            {
                if((1<<d)&inc)     move[i][d] = ((i==X-1)?0  :i+1);
                else if((1<<d)&dec)move[i][d] = ((i==0  )?X-1:i-1);
                else move[i][d] = i;
            }
    }
    inline T operator ()(T x,T d)
    {
        return move[x][d];
    }
};
template<int X,typename T>
class distance
{
private:
    T d[X][X];
    #define POS(i,j) ((i>j)?(i-j):(j-i))
public:
    static const T maxd = X/2;
    distance()
    {
        for(int i=0;i<X;i++)
        {
            for(int j=0;j<X;j++)
            {
                d[i][j] = min(POS(i,j),POS(i,j+X),POS(i+X,j));
                assert(maxd>=d[i][j]);
            }
        }
        #undef POS
    }
    inline T operator() (T x1,T x2)
    {
        return d[x1][x2];
    }
};
template<int X,typename T>
class relative
{
private:
    T d[X][X];
    #define POS(i,j) ((i>=j)?(i-j):(i-j+X))
public:
    static const T maxd = X;
    relative()
    {
        for(int i=0;i<X;i++)
        {
            for(int j=0;j<X;j++)
            {
                d[i][j] = POS(i,j);
                assert(maxd>=d[i][j]);
            }
        }
        #undef POS
    }
    inline T operator() (T x1,T x2)
    {
        return d[x1][x2];
    }
};
#endif // UTILS_MOVE_H

