#include<iostream>
#include<fstream>
#include<sstream>
#include<string.h>
#include <sys/time.h>
#include <arm_neon.h>
#include<omp.h>
using namespace std;
#define E_LineN 4535 //����Ԫ���� line number
#define RN 8399//�������� row number
#define thread_nums 6
char** eliminer = new char* [RN]; //��Ϊ����Ԫ����RN�У�������Ԫ������ʱ���RN��

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
    if (eliminer[i][i / 8] == 0)  //��Ԫ�ӵ�i�еĵ�iλ��һ����1�����û��˵��Ϊ��
        return true;
    return false;
}
int nextstart = 0;
int flagUpgrade = 0;  //ȫ�ֱ�����������Ϊ�̺߳����ֲ�����
void elimination() {
    int i, j, ss, n, newid;
    uint8x16_t ve;
    uint8x16_t vr;
    uint8x16_t tmp;
#pragma omp parallel num_threads(thread_nums),private(i,j,ss,n,newid,ve,vr,tmp)
    {
        for (i = RN - 1; i >= 0; i--) { //��Ԫ��
            if (!eliminer_ifnull(i)) {
#pragma omp for 
                for (j = 0; j < E_LineN; j++) { //�Ե�j������Ԫ�н�����Ԫ
                    if (eline[j].num == i) {
                        ss = 0;
                        for (ss; ss + 16 <= i / 8; ss += 16) { //32���ֽڣ�8��int
                            ve = vld1q_u8((const uint8_t*)&eline[j].bit[ss]);
                            vr = vld1q_u8((const uint8_t*)&eliminer[i][ss]);


                            tmp = veorq_u8(ve, vr);
                            vst1q_u8((uint8_t*)&eline[j].bit[ss], tmp);
                        }
                        for (ss; ss <= i / 8; ss++)  //��i/8�����ܶԲ���Ҫ������λ���������
                            eline[j].bit[ss] ^= eliminer[i][ss];
                        if (eline_ifnull(j)) {
                            eline[j].ifUprade = true;
                            eline[j].num = -2;
                        }
                        else {
                            for (n = eline[j].num; n >= 0; n--) {//�ҵ��׸�1
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

#pragma omp barrier

#pragma omp single 
                {
                    flagUpgrade = 0;
                    for (j = 0; j < E_LineN; j++) {
                        if (eline[j].num == i) {
                            for (ss = eline[j].num / 8; ss >= 0; ss--)  //��ֵ����
                                eliminer[eline[j].num][ss] = eline[j].bit[ss];
                            eline[j].ifUprade = true;
                            flagUpgrade = 1;
                            nextstart = j;//֮ǰ�ı���Ԫ�е���λ������i,������Ԫ
                            eline[j].num = -2;
                            break;
                        }
                    }
                }
#pragma omp barrier
                if (flagUpgrade == 0)
                    continue;
                newid = nextstart + 1;

#pragma omp for
                for (j = newid; j < E_LineN; j++) {
                    if (eline[j].num == i) {
                        ss = 0;
                        for (ss; ss + 16 <= i / 8; ss += 16) { //32���ֽڣ�8��int
                            ve = vld1q_u8((const uint8_t*)&eline[j].bit[ss]);
                            vr = vld1q_u8((const uint8_t*)&eliminer[i][ss]);


                            tmp = veorq_u8(ve, vr);
                            vst1q_u8((uint8_t*)&eline[j].bit[ss], tmp);
                        }
                        for (ss; ss <= i / 8; ss++)  //��i/8�����ܶԲ���Ҫ������λ���������
                            eline[j].bit[ss] ^= eliminer[i][ss];
                        if (eline_ifnull(j)) {
                            eline[j].ifUprade = true;
                            eline[j].num = -2;
                        }
                        else {
                            for (n = eline[j].num; n >= 0; n--) {//�ҵ��׸�1
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
    }
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
    struct timeval start;
    struct timeval end;//clock
    float timecount;

    gettimeofday(&start, NULL);
    elimination();
    gettimeofday(&end, NULL);
    timecount = (end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec;
    cout << "������ʱ:" << timecount << "΢��" << endl;

}