#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include<omp.h>
using namespace std;
#define thread_nums 8
const int N = 2000;
float m[N][N];
void m_reset()
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < i; j++)
            m[i][j] = 0;
        m[i][i] = 1.0;
        for (int j = i + 1; j < N; j++)
        {
            m[i][j] = rand();
        }
    }
    for (int k = 0; k < N; k++)
        for (int i = k + 1; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                if (k && k % 8 == 0)
                    m[i][j] -= m[k][j];
                else
                    m[i][j] += m[k][j];
            }
}

int main()
{

    m_reset();
    struct timeval start;
    struct timeval end;//clock
    float timecount;

    gettimeofday(&start, NULL);
    int i, j, k;
    float tmp;
#pragma omp parallel num_threads(thread_nums),private(i,j,k,tmp)
    for (k = 0; k < N; k++)
    {
#pragma omp single
        {
            tmp = m[k][k];
            for (int j = k + 1; j < N; j++)
                m[k][j] = m[k][j] / m[k][k];
            m[k][k] = 1.0;
        }
#pragma omp for
        for (i = k + 1; i < N; i++)
        {
            tmp = m[i][k];
            for (j = k + 1; j < N; j++)
                m[i][j] = m[i][j] - tmp * m[k][j];
            m[i][k] = 0;
        }
    }
    gettimeofday(&end, NULL);
    timecount = (end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec;
    cout << "������ʱ:" << timecount << "΢��" << endl;
    return 0;
}
