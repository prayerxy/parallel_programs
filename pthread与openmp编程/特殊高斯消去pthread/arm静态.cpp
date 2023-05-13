#include<iostream>
#include<fstream>
#include<sstream>
#include<string.h>
#include <arm_neon.h>
#include <sys/time.h>
#include<pthread.h>

#include<semaphore.h>

using namespace std;
#define E_LineN 14291 //����Ԫ���� line number
#define RN 37960//�������� row number
#define NUM_THREADS 6

char** eliminer = new char* [RN];

struct line {

    bool ifUprade;

    int num;

    char* bit = new char[RN / 8 + 1];
};
line eline[E_LineN];//����Ԫ��

void init() {
    for (int i = 0; i < RN; i++)
        eliminer[i] = new char[RN / 8 + 1];
    for (int i = 0; i < E_LineN; i++) {
        for (int j = RN / 8; j >= 0; j--)
            eline[i].bit[j] = 0;
    }
    for (int i = 0; i < RN; i++) {
        for (int j = RN / 8; j >= 0; j--)
            eliminer[i][j] = 0;
    }
}
bool eline_ifnull(int i) { //��i��
    bool flag = true;
    for (int j = RN / 8; j >= 0; j--) {
        if (!eline[i].bit[j] == 0)  //ֻҪ��һ����Ϊ�գ���return false;
        {
            flag = false;
            return flag;
        }
    }
    return flag;
}
bool eliminer_ifnull(int i) {
    if (eliminer[i][i / 8] == 0)
        return true;
    return false;
}
typedef struct {
    int t_id;   //�߳�id
}threadParam_t;


pthread_barrier_t barrier_eliminerNULL;
pthread_barrier_t barrier_AfterEliminerNull;

int nextstart = 0;
int flagUpgrade = 0;  //ȫ�ֱ�����������Ϊ�̺߳����ֲ�����


void* threadFunc(void* param) {
    threadParam_t* p = (threadParam_t*)param;

    int t_id = p->t_id;//�̵߳ı��
    for (int i = RN - 1; i >= 0; i--) {
        if (!eliminer_ifnull(i)) {
            for (int j = t_id; j < E_LineN; j += NUM_THREADS) { //�Ե�j������Ԫ�н�����Ԫ
                if (eline[j].num == i) {
                    int ss = 0;
                    for (ss; ss <= i / 8; ss++)
                        eline[j].bit[ss] ^= eliminer[i][ss];
                    if (eline_ifnull(j)) {
                        eline[j].ifUprade = true;
                        eline[j].num = -2;
                    }
                    else {
                        for (int n = eline[j].num; n >= 0; n--) {//�ҵ��׸�1
                            if (eline[j].bit[n / 8] & (1 << (n % 8))) {  //n%8Ϊ0-7��λ�ã�n/8Ϊchar��λ��
                                eline[j].num = n;
                                break;
                            }
                        }
                    }
                }

            }
        }
        else {
            pthread_barrier_wait(&barrier_eliminerNULL);

            if (t_id == 0) {
                for (int j = 0; j < E_LineN; j++) {
                    if (eline[j].num == i) {
                        for (int ss = eline[j].num / 8; ss >= 0; ss--)
                            eliminer[eline[j].num][ss] = eline[j].bit[ss];
                        eline[j].ifUprade = true;
                        flagUpgrade = 1;
                        nextstart = j;
                        eline[j].num = -2;
                        break;
                    }
                }
            }
            pthread_barrier_wait(&barrier_AfterEliminerNull);
            if (flagUpgrade == 0)
                continue;
            int newid = t_id;
            while (newid <= nextstart) {
                newid += NUM_THREADS;
            }
            for (int j = newid; j < E_LineN; j += NUM_THREADS) {
                if (eline[j].num == i) {
                    int ss = 0;
                    for (ss; ss <= i / 8; ss++)
                        eline[j].bit[ss] ^= eliminer[i][ss];
                    if (eline_ifnull(j)) {
                        eline[j].ifUprade = true;
                        eline[j].num = -2;
                    }
                    else {//����num
                        for (int n = eline[j].num; n >= 0; n--) {//�ҵ��׸�1
                            if (eline[j].bit[n / 8] & (1 << (n % 8))) {  //n%8Ϊ0-7��λ�ã�n/8Ϊchar��λ��
                                eline[j].num = n;
                                break;
                            }
                        }
                    }
                }
            }
        }

    }


    pthread_exit(NULL);
    return NULL;
}

int main() {
    init();
    ifstream file_eliminer;
    file_eliminer.open("/home/data/Groebner/1_130_22_8/1.txt", ios_base::in);
    string s1;
    if (!file_eliminer.is_open())
        cout << " ��ʧ��";
    while (getline(file_eliminer, s1)) { //���ж���
        stringstream st;//���ַ���ת������
        st << s1;
        int tmp;
        int row;//��Ԫ�ӵ����������һ��1��λ�����
        bool first = true;
        while (st >> tmp) {
            if (first) {
                row = tmp;
                first = false;
            }
            int i = tmp / 8; //���еĵڼ���char��
            int j = tmp % 8;//�����char�ĵڼ���λ����
            eliminer[row][i] |= (1 << j);
        }
    }
    file_eliminer.close();

    ifstream file_eline;
    file_eline.open("/home/data/Groebner/1_130_22_8/2.txt", ios_base::in);
    string s2;
    if (!file_eline.is_open())
        cout << " ��ʧ��";
    int x = 0;//����Ԫ�е��±�
    while (getline(file_eline, s2)) { //���ж���
        stringstream st;//���ַ���ת������
        st << s2;
        int tmp;
        bool first = true; //Ҫ��¼����Ԫ�е���1��λ��
        while (st >> tmp) {
            if (first) {
                eline[x].num = tmp;
                first = false;
            }
            int i = tmp / 8;
            int j = tmp % 8;
            eline[x].bit[i] |= (1 << j);
        }
        x++;
    }
    file_eline.close();

    pthread_barrier_init(&barrier_eliminerNULL, NULL, NUM_THREADS);
    pthread_barrier_init(&barrier_AfterEliminerNull, NULL, NUM_THREADS);


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
    pthread_barrier_destroy(&barrier_eliminerNULL);
    pthread_barrier_destroy(&barrier_AfterEliminerNull);
    delete[] handles;
    delete[]param;
    gettimeofday(&end, NULL);
    timecount = (end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec;
    cout << "������ʱ:" << timecount << "΢��" << endl;

}
