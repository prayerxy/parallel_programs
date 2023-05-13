#include<iostream>
#include<fstream>
#include<sstream>
#include<string.h>
#include<Windows.h>  

#include<vector>
#include<immintrin.h>

using namespace std;
#define E_LineN 14291 //����Ԫ���� line number
#define RN 37960//�������� row number


char** eliminer = new char* [RN]; //��Ϊ����Ԫ����RN�У�������Ԫ������ʱ���RN��

struct line {  //������Ԫ��,���������׸�1�ı������Ա������Ԫ����
    //�Ƿ�����
    bool ifUprade;
    //�׸�1��λ��
    int num;
    //λͼ��ʾһ��
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
    for (int i = RN - 1; i >= 0; i--) { //��Ԫ��
        if (!eliminer_ifnull(i)) {
            for (int j = 0; j < E_LineN; j++) { //�Ե�j������Ԫ�н�����Ԫ
                if (eline[j].num == i) {
                    int ss = 0;
                    for (ss; ss + 32 <= i / 8; ss += 32) { //32���ֽڣ�8��int
                        __m256i ve = _mm256_loadu_epi8(eline[j].bit + ss);
                        __m256i vr = _mm256_loadu_epi8(eliminer[i] + ss);
                        __m256i tmp = _mm256_xor_si256(ve, vr);
                        _mm256_storeu_epi8(eline[j].bit + ss, tmp);
                    }
                    for (ss; ss <= i / 8; ss++)  //��i/8�����ܶԲ���Ҫ������λ���������
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
        else {


            for (int j = 0; j < E_LineN; j++) {
                if (eline[j].num == i) {
                    for (int ss = eline[j].num / 8; ss >= 0; ss--)  //��ֵ����
                        eliminer[eline[j].num][ss] = eline[j].bit[ss];
                    eline[j].ifUprade = true;
                    flagUpgrade = 1;
                    nextstart = j;//֮ǰ�ı���Ԫ�е���λ������i,������Ԫ
                    eline[j].num = -2;
                    break;
                }
            }

            if (flagUpgrade == 0)
                continue;
            int newid = nextstart + 1;
            for (int j = newid; j < E_LineN; j++) {
                if (eline[j].num == i) {
                    int ss = 0;
                    for (ss; ss + 32 <= i / 8; ss += 32) { //32���ֽڣ�8��int
                        __m256i ve = _mm256_loadu_epi8(eline[j].bit + ss);
                        __m256i vr = _mm256_loadu_epi8(eliminer[i] + ss);
                        __m256i tmp = _mm256_xor_si256(ve, vr);
                        _mm256_storeu_epi8(eline[j].bit + ss, tmp);
                    }
                    for (ss; ss <= i / 8; ss++)  //��i/8�����ܶԲ���Ҫ������λ���������
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
}
int main() {
    init();
    //��Ԫ�Ӷ���
    ifstream file_eliminer;
    file_eliminer.open("��Ԫ��.txt", ios_base::in);
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
    //����Ԫ�еĶ���
    ifstream file_eline;
    file_eline.open("����Ԫ��.txt", ios_base::in);
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
    long long head, tail, freq; //timers
    //��ʱ��
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    //start time
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    elimination();
    //end time
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout << "������ʱ:" << ((tail - head) * 1000000.0 / freq) << "΢��" << endl;
    //��������Ԫ����Ϊ���д����Ԫ�����
    ofstream out;
    out.open("��Ԫ���2.txt", ios::trunc);
    for (int i = 0; i < E_LineN; i++) {
        if (eline_ifnull(i))continue;
        string s;
        string stemp;
        stringstream ss;
        for (int j = RN / 8; j >= 0; j--) {
            for (int k = 7; k >= 0; k--) {
                if (eline[i].bit[j] & (1 << k)) {
                    ss.clear();
                    ss << (j * 8 + k);
                    ss >> stemp;
                    s.append(stemp); s.append(" ");
                }
            }

        }
        out << s << "\n";
        s.clear();
    }
}
