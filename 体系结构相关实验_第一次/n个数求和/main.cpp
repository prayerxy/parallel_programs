#include<iostream>
#include<windows.h>
#include<stdlib.h>
//#pragma GCC optimize(3)
using namespace std;

const long long int N=200000000;

double a[N],sum,sum1,sum2,sum3,sum4,sum5;
void init()
{
    for(int i=0;i<N;i++)
        a[i]=i;
    sum=0;
}
void add1(double a[N])
{
    for(int i=0;i<N;i++)
        sum+=a[i];
}
void add2(double a[N])
{
    int i;sum1=0;sum2=0;sum3=0;sum4=0;sum5=0;
    for(i=0;i<N-4;i+=5){
      sum1+=a[i];
      sum2+=a[i+1];
      sum3+=a[i+2];
      sum4+=a[i+3];
      sum5+=a[i+4];
    }
    for(;i<N;i++)
        sum1+=a[i];
    sum=sum1+sum2+sum3+sum4+sum5;
}

int main()
{
    long long head,tail,freq; //timers
    long double tmp=0.0;
    init();
    //��ʱ��
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    for(int k=0;k<100;k++)
    {
    //start time
       QueryPerformanceCounter((LARGE_INTEGER*)&head);
       add1(a);
    //end time
       QueryPerformanceCounter((LARGE_INTEGER*)&tail);
       sum=0;
       tmp+=(tail-head)*1000.0/freq;

    }

    cout<<"�㷨��ʱ:"<<tmp/100<<"ms"<<endl;
    system("pause");
}
