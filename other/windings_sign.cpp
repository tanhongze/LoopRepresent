#include<cstdio>
int windings_sign(int windings)
    {
        int sign;
        if(windings&1)
        {
            sign = -1;
        }
        else sign = 1;
        for(int l=0;l<2;l++)
        {
            if((windings>>(2*l))&3)
            {
                sign = -sign;
            }
        }
        return sign;
    }

    void statisticZu(long long Z[1<<(2*2)],long long Zu[1<<(2*2)])
    {
        {
            for(int i=0;i<4;i++)
            {
                for(int j=0;j<4;j++)
                {
                    Zu[i*4+j] = 4ll*Z[0];
                    for(int u=0;u<4;u++)
                    {
                        if((i&u)==0||(i&u)==3)Zu[i*4+j]-=2ll*Z[u*4];
                        else Zu[i*4+j]+=2*Z[u*4];
                        if((j&u)==0||(j&u)==3)Zu[i*4+j]-=2ll*Z[u];
                        else Zu[i*4+j]+=2*Z[u  ];
                        for(int v=0;v<4;v++)
                        {
                            int si = (u&i)^(v&j);
                            if(si==0||si==3)Zu[i*4+j]+=Z[u*4+v];
                            else Zu[i*4+j]-=Z[u*4+v];
                        }
                    }
                }
            }
        }
    }
int main()
{
    for(int i=0;i<16;i++)
    {
        long long Z[16] = {};
        long long Zu[16] = {};
        Z[i] = 1;
        statisticZu(Z,Zu);
        printf("%d %d ",i,windings_sign(i));
        printf("%lld \n",Zu[1]);
    }
    return 0;
}
