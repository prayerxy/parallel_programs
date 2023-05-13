#include<iostream>
#include<pthread.h>
#pragma comment(lib, "pthreadVC2.lib")  //静态链接
#include<Windows.h>
#include<semaphore.h>
#define NUM_THREADS 10   //创建的线程数
using namespace std;

const int N = 100;
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
sem_t sem_leader;
sem_t sem_Divsion[NUM_THREADS - 1]; //对于线程0不需要此信号
sem_t sem_Elimination[NUM_THREADS - 1];

//线程函数定义
void* threadFunc(void* param) {
    threadParam_t* p = (threadParam_t*)param;

    int t_id = p->t_id;//线程的编号

    for (int k = 0; k < N; k++) {
        //t_id为0的线程做除法操作，其他线程先等待
        //这里只采用一个工作线程负责除法操作，而所有线程都做消去
        if (t_id == 0) {
            for (int j = k + 1; j < N; j++)
                m[k][j] = m[k][j] / m[k][k];
            m[k][k] = 1.0;
        }
        else {
            sem_wait(&sem_Divsion[t_id - 1]);
        }

        if (t_id == 0) {
            for (int i = 0; i < NUM_THREADS - 1; i++) {
                sem_post(&sem_Divsion[i]);
            }
        }
        //t_id为0的也进行消去操作,从k+1行开始消元
        for (int i = k + 1 + t_id; i < N; i += NUM_THREADS) {
            //消去
            for (int j = k + 1; j < N; j++) {
                m[i][j] = m[i][j] - m[i][k] * m[k][j];
            }
            m[i][k] = 0.0;
        }
        if (t_id == 0) {
            for (int i = 0; i < NUM_THREADS - 1; i++)
                sem_wait(&sem_leader); //等待其他worker完成消去

            for (int i = 0; i < NUM_THREADS - 1; i++)
                sem_post(&sem_Elimination[i]); //通知进入下一轮

        }
        else {
            sem_post(&sem_leader);
            sem_wait(&sem_Elimination[t_id - 1]);//等待下一轮通知
        }
    }
    pthread_exit(NULL);
    return NULL;
}

int main()
{
    m_reset();
    //初始化信号量
    sem_init(&sem_leader, 0, 0);
    for (int i = 0; i < NUM_THREADS - 1; i++) { //有一个子线程是leader，主线程只进行子线程创建销毁
        sem_init(&sem_Divsion[i], 0, 0);
        sem_init(&sem_Elimination[i], 0, 0);
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

    //销毁进程
    for (int t_id = 0; t_id < NUM_THREADS; t_id++)
        pthread_join(handles[t_id], NULL);
    sem_destroy(&sem_leader);
    for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
        sem_destroy(&sem_Divsion[t_id]);
        sem_destroy(&sem_Elimination[t_id]);
    }

    delete[] handles;
    delete[]param;
    /*
    cout<<endl;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            cout << m[i][j] << " ";
        cout << endl;
    }*/
    //end time
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout << "程序用时:" << ((tail - head) * 1000000.0 / freq) << "微秒" << endl;
}
