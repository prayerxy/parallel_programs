#include <iostream>
#include<xmmintrin.h>
#include<windows.h>
using namespace std;

const int N = 500;
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
    m_reset();
    long long head, tail, freq; //timers
    //记时间
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    //start time
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    m_reset();
    for (int times = 0; times < 50; times++) { //50次重复计时
        for (int k = 0; k < N; k++) {
            __m128 vt = _mm_set1_ps(m[k][k]);//用m[k][k]初始化向量vt
            int j = k + 1;
            if (j % 4 != 0) {
                while (j % 4 != 0) {
                    m[k][j] = m[k][j] / m[k][k];
                    j++;
                }
            }
            for (j; j + 4 <= N; j += 4) {
                __m128 va = _mm_load_ps(m[k] + j);
                va = _mm_div_ps(va, vt);
                _mm_store_ps(m[k] + j, va);//从向量寄存器储存到内存
            }
            for (j; j < N; j++) {
                m[k][j] = m[k][j] / m[k][k];
            }
            m[k][k] = 1.0;
            for (int i = k + 1; i < N; i++) {
                __m128 vaik = _mm_set1_ps(m[i][k]);
                j = k + 1;
                if (j % 4 != 0) {
                    while (j % 4 != 0) {
                        m[i][j] = m[i][j] - m[i][k] * m[k][j];
                        j++;
                    }
                }
                for (j; j + 4 <= N; j += 4) {//内层展开，对j2
                    __m128 vakj = _mm_load_ps(m[k] + j);
                    __m128 vaij = _mm_load_ps(m[i] + j);
                    vakj = _mm_mul_ps(vaik, vakj);
                    vaij = _mm_sub_ps(vaij, vakj);
                    _mm_store_ps(m[i] + j, vaij);
                }
                for (j; j < N; j++) {
                    m[i][j] = m[i][j] - m[i][k] * m[k][j];
                }
                m[i][k] = 0;
            }
        }
    }
    //end time
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout << "程序用时:" << ((tail - head) * 1000000.0 / freq) / 50 << "微秒" << endl;
    return 0;
}
