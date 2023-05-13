#include<iostream>
#include<pthread.h>
#include <arm_neon.h>
#include <sys/time.h>
#include<semaphore.h>
#define NUM_THREADS 2   //创建的线程数
using namespace std;

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



typedef struct {
    int t_id;   //线程id
}threadParam_t;


sem_t sem_leader;
sem_t sem_Divsion[NUM_THREADS - 1]; //对于线程0不需要此信号
sem_t sem_Elimination[NUM_THREADS - 1];


void* threadFunc(void* param) {
    threadParam_t* p = (threadParam_t*)param;

    int t_id = p->t_id;//线程的编号

    for (int k = 0; k < N; k++) {

        if (t_id == 0) {
            float32x4_t vt = vmovq_n_f32(m[k][k]);//用m[k][k]初始化向量vt
            int j;
            for (j = k + 1; j + 4 <= N; j += 4) {
                float32x4_t va = vld1q_f32(m[k] + j);
                va = vdivq_f32(va, vt);
                vst1q_f32(m[k] + j, va);//从向量寄存器储存到内存
            }
            for (j; j < N; j++) {
                m[k][j] = m[k][j] / m[k][k];
            }
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
        for (int i = k + 1 + t_id; i < N; i += NUM_THREADS) {
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

    sem_init(&sem_leader, 0, 0);
    for (int i = 0; i < NUM_THREADS - 1; i++) { //有一个子线程是leader，主线程只进行子线程创建销毁
        sem_init(&sem_Divsion[i], 0, 0);
        sem_init(&sem_Elimination[i], 0, 0);
    }

    pthread_t* handles = new pthread_t[NUM_THREADS];
    threadParam_t* param = new threadParam_t[NUM_THREADS];

    struct timeval start;
    struct timeval end;//clock
    float timecount;

    gettimeofday(&start, NULL);

    for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
        param[t_id].t_id = t_id;
        pthread_create(&handles[t_id], NULL, threadFunc, (&param[t_id]));
    }


    for (int t_id = 0; t_id < NUM_THREADS; t_id++)
        pthread_join(handles[t_id], NULL);
    sem_destroy(&sem_leader);
    for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
        sem_destroy(&sem_Divsion[t_id]);
        sem_destroy(&sem_Elimination[t_id]);
    }

    delete[] handles;
    delete[]param;

    gettimeofday(&end, NULL);
    timecount = (end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec;
    cout << "程序用时:" << timecount << "微秒" << endl;
    return 0;
}
