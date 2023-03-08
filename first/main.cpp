#include<iostream>
#include<windows.h>
#include<stdlib.h>

using namespace std;

const int N=10240;

double b[N][N],a[N],sum[N];

void init(int n) //��N*N���鸳��ֵ
{
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++)
            b[i][j]=i+j;
    for(int i=0;i<N;i++)
        a[i]=i;
}

int main()
{
    long long head,tail,freq; //timers
    init(N);
    //��ʱ��
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    //start time
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    for(int i=0;i<N;i++){
        sum[i]=0.0;
    }
    for(int j=0;j<N;j++){
        for(int i=0;i<N;i++)
            sum[i]+=b[j][i]*a[j];
    }
    //end time
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout<<"ƽ���㷨��ʱ:"<<(tail-head)*1000.0/freq<<"ms"<<endl;
    system("pause");
}
