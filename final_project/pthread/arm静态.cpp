#include<iostream>
#include<fstream>
#include<sstream>
#include<string.h>
#include <arm_neon.h>
#include <sys/time.h>
#include<pthread.h>

#include<semaphore.h>

using namespace std;
#define E_LineN 14291 //被消元行数 line number
#define RN 37960//矩阵列数 row number
#define NUM_THREADS 6

char** eliminer = new char* [RN];

struct line {

    bool ifUprade;

    int num;

    char* bit = new char[RN / 8 + 1];
};
line eline[E_LineN];//被消元行

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
bool eline_ifnull(int i) { //第i行
    bool flag = true;
    for (int j = RN / 8; j >= 0; j--) {
        if (!eline[i].bit[j] == 0)  //只要有一个不为空，就return false;
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
    int t_id;   //线程id
}threadParam_t;


pthread_barrier_t barrier_eliminerNULL;
pthread_barrier_t barrier_AfterEliminerNull;

int nextstart = 0;
int flagUpgrade = 0;  //全局变量，不能作为线程函数局部变量


void* threadFunc(void* param) {
    threadParam_t* p = (threadParam_t*)param;

    int t_id = p->t_id;//线程的编号
    for (int i = RN - 1; i >= 0; i--) {
        if (!eliminer_ifnull(i)) {
            for (int j = t_id; j < E_LineN; j += NUM_THREADS) { //对第j个被消元行进行消元
                if (eline[j].num == i) {
                    int ss = 0;
                    for (ss; ss <= i / 8; ss++)
                        eline[j].bit[ss] ^= eliminer[i][ss];
                    if (eline_ifnull(j)) {
                        eline[j].ifUprade = true;
                        eline[j].num = -2;
                    }
                    else {
                        for (int n = eline[j].num; n >= 0; n--) {//找到首个1
                            if (eline[j].bit[n / 8] & (1 << (n % 8))) {  //n%8为0-7的位置，n/8为char的位置
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
                    else {//重置num
                        for (int n = eline[j].num; n >= 0; n--) {//找到首个1
                            if (eline[j].bit[n / 8] & (1 << (n % 8))) {  //n%8为0-7的位置，n/8为char的位置
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
        cout << " 打开失败";
    while (getline(file_eliminer, s1)) { //逐行读入
        stringstream st;//从字符串转到数字
        st << s1;
        int tmp;
        int row;//消元子的行数与其第一个1的位置相等
        bool first = true;
        while (st >> tmp) {
            if (first) {
                row = tmp;
                first = false;
            }
            int i = tmp / 8; //这行的第几个char中
            int j = tmp % 8;//在这个char的第几个位置上
            eliminer[row][i] |= (1 << j);
        }
    }
    file_eliminer.close();

    ifstream file_eline;
    file_eline.open("/home/data/Groebner/1_130_22_8/2.txt", ios_base::in);
    string s2;
    if (!file_eline.is_open())
        cout << " 打开失败";
    int x = 0;//被消元行的下标
    while (getline(file_eline, s2)) { //逐行读入
        stringstream st;//从字符串转到数字
        st << s2;
        int tmp;
        bool first = true; //要记录被消元行的首1的位置
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
    cout << "程序用时:" << timecount << "微秒" << endl;

}
