#include<iostream>
#include<pthread.h>
#include <sys/time.h>
#include <arm_neon.h>
#define count 6   //创建的线程数
using namespace std;

const int N = 400;
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
    int k; //消去的轮次
    int t_id; // 线程 id


}threadParam_t;

void* threadFunc(void* param) {
    threadParam_t* p = (threadParam_t*)param;
    int k = p->k;//消去的轮次
    int t_id = p->t_id;//线程的编号
    int cur = k + t_id + 1;
    for (int i = cur; i < N; i += count) {
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
    pthread_exit(NULL);
    return NULL;
}

int main()
{
    m_reset();

    int worker_count = count;
    pthread_t* handles = new pthread_t[worker_count];
    threadParam_t* param = new threadParam_t[worker_count];
    struct timeval start;
    struct timeval end;//clock
    float timecount;

    gettimeofday(&start, NULL);

    for (int k = 0; k < N; k++) {

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

        for (int t_id = 0; t_id < worker_count; t_id++) {
            param[t_id].k = k;
            param[t_id].t_id = t_id;
        }

        for (int t_id = 0; t_id < worker_count; t_id++) {
            pthread_create(&handles[t_id], NULL, threadFunc, (&param[t_id]));
        }
        for (int t_id = 0; t_id < worker_count; t_id++) {
            pthread_join(handles[t_id], NULL);
        }

    }
    delete[] handles;
    delete[]param;
    cout << endl;

    gettimeofday(&end, NULL);
    timecount = (end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec;
    cout << "程序用时:" << timecount << "微秒" << endl;
    return 0;
}
