#include<iostream>
#include<pthread.h>

#include <arm_neon.h>
#include <sys/time.h>
#include<semaphore.h>
#define NUM_THREADS 8   //�������߳���
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
            {
                if (k && k % 8 == 0)
                    m[i][j] -= m[k][j];
                else
                    m[i][j] += m[k][j];
            }
}



typedef struct {
    int t_id;   //�߳�id
}threadParam_t;


sem_t sem_main;
sem_t sem_workerstart[NUM_THREADS]; //ÿ���߳����Լ�ר�����ź���
sem_t sem_workerend[NUM_THREADS];


void* threadFunc(void* param) {
    threadParam_t* p = (threadParam_t*)param;

    int t_id = p->t_id;//�̵߳ı��

    for (int k = 0; k < N; k++) {
        sem_wait(&sem_workerstart[t_id]);//�������ȴ����߳���ɳ�������������ר���Լ����ź�����


        for (int i = k + 1; i < N; i++) {
            for (int j = k + 1 + t_id; j < N; j += NUM_THREADS)
                m[i][j] = m[i][j] - m[i][k] * m[k][j];
            m[i][k] = 0.0;
        }
        sem_post(&sem_main);//�������߳�
        sem_wait(&sem_workerend[t_id]);//�������ȴ����̻߳���
    }
    pthread_exit(NULL);
    return NULL;
}

int main()
{
    m_reset();

    sem_init(&sem_main, 0, 0);
    for (int i = 0; i < NUM_THREADS; i++) {
        sem_init(&sem_workerstart[i], 0, 0);
        sem_init(&sem_workerend[i], 0, 0);
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
    for (int k = 0; k < N; k++) {

        for (int j = k + 1; j < N; j++)
            m[k][j] = m[k][j] * 1.0 / m[k][k];
        m[k][k] = 1.0;

        for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
            sem_post(&sem_workerstart[t_id]);
        }

        for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
            sem_wait(&sem_main);
        }


        for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
            sem_post(&sem_workerend[t_id]);
        }
    }

    for (int t_id = 0; t_id < NUM_THREADS; t_id++)
        pthread_join(handles[t_id], NULL);
    sem_destroy(&sem_main);
    for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
        sem_destroy(&sem_workerstart[t_id]);
        sem_destroy(&sem_workerend[t_id]);
    }

    delete[] handles;
    delete[]param;

    gettimeofday(&end, NULL);
    timecount = (end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec;
    cout << "������ʱ:" << timecount << "΢��" << endl;
    return 0;
}
