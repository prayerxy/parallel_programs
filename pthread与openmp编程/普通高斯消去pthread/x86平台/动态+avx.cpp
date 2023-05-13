#include<iostream>
#include<pthread.h>
#include<immintrin.h>
#pragma comment(lib, "pthreadVC2.lib")  //静态链接
#include<Windows.h>

#define count 9   //创建的线程数，由于线程数不是越多越好
using namespace std;

const int N = 1000;
float m[N][N];
//使用随机函数产生测试数据集
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

typedef struct {
    int k; //消去的轮次
    int t_id; // 线程 id


}threadParam_t;

//每个线程的执行函数
void* threadFunc(void* param) {
    threadParam_t* p = (threadParam_t*)param;
    int k = p->k;//消去的轮次
    int t_id = p->t_id;//线程的编号
    int cur = k + t_id + 1;//获取自己的任务
    for (int i = k + t_id + 1; i < N; i += count) {
        __m256 vaik = _mm256_set1_ps(m[i][k]);
        int j = k + 1;
        for (j; j + 8 <= N; j += 8) {//内层展开，对j2
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
    // cout << "finished:" << p->t_id << endl;
    pthread_exit(NULL);
    return NULL;
}

int main()
{
    m_reset();

    int worker_count = count;
    pthread_t* handles = new pthread_t[worker_count];
    threadParam_t* param = new threadParam_t[worker_count];
    long long head, tail, freq; //timers
  //记时间
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    //start time
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    //创建工作线程，进行消去操作
    for (int k = 0; k < N; k++) {
        __m256 vt = _mm256_set1_ps(m[k][k]);//用m[k][k]初始化向量vt
        int j = k + 1;
        for (j; j + 8 <= N; j += 8) {
            __m256 va = _mm256_loadu_ps(m[k] + j);
            va = _mm256_div_ps(va, vt);
            _mm256_storeu_ps(m[k] + j, va);//从向量寄存器储存到内存
        }
        for (j; j < N; j++) {
            m[k][j] = m[k][j] / m[k][k];
        }
        m[k][k] = 1.0;
        //分配任务
        for (int t_id = 0; t_id < worker_count; t_id++) {
            param[t_id].k = k;
            param[t_id].t_id = t_id;
        }
        //创建线程
        for (int t_id = 0; t_id < worker_count; t_id++) {
            pthread_create(&handles[t_id], NULL, threadFunc, (&param[t_id]));
        }
        for (int t_id = 0; t_id < worker_count; t_id++) {
            pthread_join(handles[t_id], NULL);
        }
        //  cout << "------------" << endl;
    }
    delete[] handles;
    delete[]param;
    cout << endl;
    /*for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            cout << m[i][j] << " ";
        cout << endl;
    }*/
    //end time
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout << "程序用时:" << ((tail - head) * 1000000.0 / freq) << "微秒" << endl;
}
