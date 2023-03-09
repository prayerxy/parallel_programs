#include<iostream>
#include<windows.h>
#include<stdlib.h>
#pragma GCC optimize(3)
using namespace std;

const int N=2000;

double b[N][N],a[N],sum[N];

void init(int n) //对N*N数组赋初值
{
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++)
            b[i][j]=i+j;
    for(int i=0;i<N;i++)
        a[i]=i;
}
void pattern1()//平凡算法
{
    for(int i=0;i<N;i++){
        sum[i]=0.0;
        for(int j=0;j<N;j++)
            sum[i]+=b[j][i]*a[j];
    }
}
void pattern2()//Cache优化算法
{
    for(int i=0;i<N;i++){
        sum[i]=0.0;
    }
    for(int j=0;j<N;j++){
        for(int i=0;i<N;i++)
            sum[i]+=b[j][i]*a[j];
    }
}
int main()
{
    long long head,tail,freq; //timers
    long double tmp=0.0;
    init(N);
    //记时间
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    for(int k=0;k<100;k++)
    {
        //start time
      QueryPerformanceCounter((LARGE_INTEGER*)&head);
      pattern1();
      //end time
      QueryPerformanceCounter((LARGE_INTEGER*)&tail);
      tmp+=(tail-head)*1000.0/freq;
      for(int m=0;m<N;m++)//重置
        sum[m]=0.0;
    }

    cout<<"算法用时:"<<tmp/100<<"ms"<<endl;
    system("pause");
}
