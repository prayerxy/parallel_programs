#include<iostream>
#include<pthread.h>
#pragma comment(lib, "pthreadVC2.lib")  //��̬����
#include<Windows.h>
#include<semaphore.h>
#include<immintrin.h>
#define NUM_THREADS 10   //�������߳���
using namespace std;

const int N = 100;
float m[N][N];
//ʹ��������������������ݼ�
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

//�ź�������
sem_t sem_leader;
sem_t sem_Divsion[NUM_THREADS - 1]; //�����߳�0����Ҫ���ź�
sem_t sem_Elimination[NUM_THREADS - 1];

//�̺߳�������
void* threadFunc(void* param) {
    threadParam_t* p = (threadParam_t*)param;

    int t_id = p->t_id;//�̵߳ı��

    for (int k = 0; k < N; k++) {
        //t_idΪ0���߳������������������߳��ȵȴ�
        //����ֻ����һ�������̸߳�������������������̶߳�����ȥ
        if (t_id == 0) {
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
        }
        else {
            sem_wait(&sem_Divsion[t_id - 1]);
        }

        if (t_id == 0) {
            for (int i = 0; i < NUM_THREADS - 1; i++) {
                sem_post(&sem_Divsion[i]);
            }
        }
        //t_idΪ0��Ҳ������ȥ����,��k+1�п�ʼ��Ԫ

        for (int i = k + t_id + 1; i < N; i += NUM_THREADS) {
            __m256 vaik = _mm256_set1_ps(m[i][k]);
            int j = k + 1;
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
        if (t_id == 0) {
            for (int i = 0; i < NUM_THREADS - 1; i++)
                sem_wait(&sem_leader); //�ȴ�����worker�����ȥ

            for (int i = 0; i < NUM_THREADS - 1; i++)
                sem_post(&sem_Elimination[i]); //֪ͨ������һ��

        }
        else {
            sem_post(&sem_leader);
            sem_wait(&sem_Elimination[t_id - 1]);//�ȴ���һ��֪ͨ
        }
    }
    pthread_exit(NULL);
    return NULL;
}

int main()
{
    m_reset();
    //��ʼ���ź���
    sem_init(&sem_leader, 0, 0);
    for (int i = 0; i < NUM_THREADS - 1; i++) { //��һ�����߳���leader�����߳�ֻ�������̴߳�������
        sem_init(&sem_Divsion[i], 0, 0);
        sem_init(&sem_Elimination[i], 0, 0);
    }
    //�����߳�
    pthread_t* handles = new pthread_t[NUM_THREADS];
    threadParam_t* param = new threadParam_t[NUM_THREADS];

    long long head, tail, freq; //timers
  //��ʱ��
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    //start time
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    //���������߳�
    for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
        param[t_id].t_id = t_id;
        pthread_create(&handles[t_id], NULL, threadFunc, (&param[t_id]));
    }

    //���ٽ���
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
    cout << "������ʱ:" << ((tail - head) * 1000000.0 / freq) << "΢��" << endl;
}
