#include<iostream>
#include<pthread.h>
#pragma comment(lib, "pthreadVC2.lib")  //��̬����
#include<Windows.h>

#define count 8   //�������߳����������߳�������Խ��Խ��
using namespace std;

const int N = 800;
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
    int k; //��ȥ���ִ�
    int t_id; // �߳� id


}threadParam_t;

//ÿ���̵߳�ִ�к���
void* threadFunc(void* param) {
    threadParam_t* p = (threadParam_t*)param;
    int k = p->k;//��ȥ���ִ�
    int t_id = p->t_id;//�̵߳ı��
    int cur = k + t_id + 1;//��ȡ�Լ�������

    for (int i = k + 1; i < N; i++) {
        for (int j = cur; j < N; j += count)  //��ֱ����
            m[i][j] = m[i][j] - m[i][k] * m[k][j];
        m[i][k] = 0.0;
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
  //��ʱ��
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    //start time
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    //���������̣߳�������ȥ����
    for (int k = 0; k < N; k++) {
        //���߳�����������
        for (int j = k + 1; j < N; j++)
            m[k][j] = m[k][j] * 1.0 / m[k][k];
        m[k][k] = 1.0;
        //��������
        for (int t_id = 0; t_id < worker_count; t_id++) {
            param[t_id].k = k;
            param[t_id].t_id = t_id;
        }
        //�����߳�
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

    //end time
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout << "������ʱ:" << ((tail - head) * 1000000.0 / freq) << "΢��" << endl;
}
