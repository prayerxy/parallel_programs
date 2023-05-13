#include<iostream>
#include<pthread.h>
#pragma comment(lib, "pthreadVC2.lib")  //��̬����
#include<Windows.h>
#include<semaphore.h>
#define NUM_THREADS 4  //�������߳���
using namespace std;

const int N = 1000;
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

//barrier����
pthread_barrier_t barrier_Divsion;
pthread_barrier_t barrier_Elimination;

//�̺߳�������
void* threadFunc(void* param) {
    threadParam_t* p = (threadParam_t*)param;

    int t_id = p->t_id;//�̵߳ı��

    for (int k = 0; k < N; k++) {
        //t_idΪ0���߳����������������������߳��ȵȴ�
        if (t_id == 0) {
            for (int j = k + 1; j < N; j++)
                m[k][j] = m[k][j] / m[k][k];
            m[k][k] = 1.0;
        }
        //��һ��ͬ����
        pthread_barrier_wait(&barrier_Divsion);

        //ѭ����������,Ҳ���Գ��Զ������񻮷ֵķ�ʽ
        for (int i = k + 1 + t_id; i < N; i += NUM_THREADS) {
            //��ȥ
            for (int j = k + 1; j < N; j++)
                m[i][j] = m[i][j] - m[i][k] * m[k][j];
            m[i][k] = 0.0;
        }

        //�ڶ���ͬ����
        pthread_barrier_wait(&barrier_Elimination);
    }

    pthread_exit(NULL);
    return NULL;
}

int main()
{
    m_reset();
    //��ʼ��barrier
    pthread_barrier_init(&barrier_Divsion, NULL, NUM_THREADS);
    pthread_barrier_init(&barrier_Elimination, NULL, NUM_THREADS);

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
    pthread_barrier_destroy(&barrier_Divsion);
    pthread_barrier_destroy(&barrier_Elimination);
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
