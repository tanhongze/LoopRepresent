#ifndef TIMER_H
#define TIMER_H
#include <sys/time.h>
#include <cstdio>
class mytimer
{
    struct timeval t[2];
public:
    int  s;
    int ms;
    int us;
    void refresh()
    {
        us = t[1].tv_usec-t[0].tv_usec;
        s  = t[1].tv_sec -t[0].tv_sec ;
        if(us<0)
        {
            s -=1;
            us+=1000000;
        }
        ms = us/1000;
    }
    void begin()
    {
        gettimeofday(&t[0], NULL);
    }
    void end()
    {
        gettimeofday(&t[1], NULL);
        refresh();
    }
    void printsm()
    {
        printf("%d s %d ms ",s,ms);
    }
    void printsec()
    {
        printf("%d s",s);
    }
    void print()
    {
        printf("%d s %d ms %d us",s,ms,us);
    }
};
#endif // TIMER_H
