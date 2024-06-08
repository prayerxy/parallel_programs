#include <iostream>
#include<immintrin.h>
#include<windows.h>
using namespace std;
const int N = 1000;
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
   //��ʱ��
   QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
   //start time
   QueryPerformanceCounter((LARGE_INTEGER*)&head);
   for (int times = 0; times < 50; times++) { //50���ظ���ʱ
       for (int k = 0; k < N; k++) {
           __m256 vt = _mm256_set1_ps(m[k][k]);//��m[k][k]��ʼ������vt
           int j = k + 1;
           for (j; j + 8 <= N; j += 8) {
               __m256 va = _mm256_loadu_ps(m[k] + j);
               va = _mm256_div_ps(va, vt);
               _mm256_storeu_ps(m[k] + j, va);//�������Ĵ������浽�ڴ�
           }
           for (j; j < N; j++) {
               m[k][j] = m[k][j] / m[k][k];
           }
           m[k][k] = 1.0;
           for (int i = k + 1; i < N; i++) {
               __m256 vaik = _mm256_set1_ps(m[i][k]);
               j = k + 1;
               for (j; j + 8 <= N; j += 8) {//�ڲ�չ������j2
                   __m256 vakj = _mm256_loadu_ps(m[k] + j);
                   __m256 vaij = _mm256_loadu_ps(m[i] + j);
                   vakj = _mm256_mul_ps(vaik, vakj);
                   vaij = _mm256_sub_ps(vaij, vakj);
                   _mm256_storeu_ps(m[i] + j, vaij);
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
   cout << "������ʱ:" << ((tail - head) * 1000000.0 / freq) / 50 << "΢��" << endl;
   return 0;
}
