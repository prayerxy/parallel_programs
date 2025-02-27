#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <arm_neon.h>
#include <sys/time.h>
using namespace std;

const int N = 100;
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
                m[i][j] += m[k][j];
}
int main() {
    struct timeval start;
    struct timeval end;//clock
    float timecount;
    m_reset();
    gettimeofday(&start, NULL);
    for (int times = 0; times < 50; times++) { //50次重复计时
        for (int k = 0; k < N; k++) {
            for (int j = k + 1; j < N; j++)
                m[k][j] = m[k][j] / m[k][k];
            m[k][k] = 1.0;
            for (int i = k + 1; i < N; i++) {
                float32x4_t vaik = vmovq_n_f32(m[i][k]);
                int j;
                for (j = k + 1; j + 4 <= N; j += 4) {//内层展开，对j2
                    float32x4_t vakj = vld1q_f32(m[k] + j);
                    float32x4_t vaij = vld1q_f32(m[i] + j);
                    vakj = vmulq_f32(vaik, vakj);
                    vaij = vsubq_f32(vaij, vakj);
                    vst1q_f32(m[i] + j, vaij);
                }
                for (j; j < N; j++) {
                    m[i][j] = m[i][j] - m[i][k] * m[k][j];
                }
                m[i][k] = 0;
            }
        }
    }
    gettimeofday(&end, NULL);
    timecount = (end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec;
    cout << "程序用时:" << timecount / 50 << "微秒" << endl;
    return 0;
}
