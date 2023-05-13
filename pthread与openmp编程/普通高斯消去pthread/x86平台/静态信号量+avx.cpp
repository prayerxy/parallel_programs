#include<iostream>
#include<pthread.h>
#pragma comment(lib, "pthreadVC2.lib")  //静态链接
#include<Windows.h>
#include<semaphore.h>
#include<immintrin.h>
#define NUM_THREADS 9   //创建的线程数
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
    int t_id;   //线程id
}threadParam_t;

//信号量定义
sem_t sem_main;
sem_t sem_workerstart[NUM_THREADS]; //每个线程有自己专属的信号量
sem_t sem_workerend[NUM_THREADS];

//线程函数定义
void* threadFunc(void* param) {
    threadParam_t* p = (threadParam_t*)param;

    int t_id = p->t_id;//线程的编号

    for (int k = 0; k < N; k++) {
        sem_wait(&sem_workerstart[t_id]);//阻塞，等待主线程完成除法操作（操作专属自己的信号量）

        for (int i = k + t_id + 1; i < N; i += NUM_THREADS) {
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
        sem_post(&sem_main);//唤醒主线程
        sem_wait(&sem_workerend[t_id]);//阻塞，等待主线程唤醒
    }

    pthread_exit(NULL);
    return NULL;
}

int main()
{
    m_reset();
    //初始化信号量
    sem_init(&sem_main, 0, 0);
    for (int i = 0; i < NUM_THREADS; i++) {
        sem_init(&sem_workerstart[i], 0, 0);
        sem_init(&sem_workerend[i], 0, 0);
    }
    //创建线程
    pthread_t* handles = new pthread_t[NUM_THREADS];
    threadParam_t* param = new threadParam_t[NUM_THREADS];

    long long head, tail, freq; //timers
  //记时间
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    //start time
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    //创建工作线程
    for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
        param[t_id].t_id = t_id;
        pthread_create(&handles[t_id], NULL, threadFunc, (&param[t_id]));
    }
    for (int k = 0; k < N; k++) {
        //主线程做除法操作

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
        //开始唤醒工作线程
        for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
            sem_post(&sem_workerstart[t_id]);
        }
        //主线程睡眠
        for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
            sem_wait(&sem_main);
        }

        //主线程再次唤醒工作线程进入下一轮次的消去任务,发出end信号
        for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
            sem_post(&sem_workerend[t_id]);
        }
    }
    //销毁进程
    for (int t_id = 0; t_id < NUM_THREADS; t_id++)
        pthread_join(handles[t_id], NULL);
    sem_destroy(&sem_main);
    for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
        sem_destroy(&sem_workerstart[t_id]);
        sem_destroy(&sem_workerend[t_id]);
    }

    delete[] handles;
    delete[]param;

    //end time
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout << "程序用时:" << ((tail - head) * 1000000.0 / freq) << "微秒" << endl;
}
