#include <iostream>
#include <random>
using namespace std;



int main()
{
    default_random_engine e; 
    uniform_real_distribution<double> u(0, 1);
    for(int i=0;i<10;i++)
    cout<<u(e)<<endl;
    
    double z[10],sum=0;
    int n=10;
    for(int i=0;i<n;i++)
    {
        z[i]=0;
        sum+=z[i];
    }
    z[1]=1;
    sum=1;
        for(int i=0;i<n;i++)
    {
        z[i]/=sum;
    }
    double x[10],y[10];
    double lowb=-1+2*z[0],upb=1;
    double lb,ub;
    for(int i=0;i<n-1;i++)
    {
        lb=lowb>0?lowb:0;
        ub=upb>2*z[i]?2*z[i]:upb;
        x[i]=u(e)*(ub-lb)+lb;
        y[i]=2*z[i]-x[i];
        if(x[i]<0||y[i]<0||x[i]>1||y[i]>1)
        cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        lowb+=(2*z[i+1]-z[i]);
        upb-=x[i];
        cout<<z[i]<<","<<x[i]<<","<<y[i]<<","<<2*z[i]-x[i]-y[i]<<endl;
    }
    
    double sx=0,sy=0;
    for(int i=0;i<n-1;i++)
    {
        sx+=x[i];
        sy+=y[i];
    }
    x[n-1]=1-sx;
    y[n-1]=1-sy;
    cout<<z[n-1]<<","<<x[n-1]<<","<<y[n-1]<<","<<2*z[n-1]-x[n-1]-y[n-1]<<endl;

}