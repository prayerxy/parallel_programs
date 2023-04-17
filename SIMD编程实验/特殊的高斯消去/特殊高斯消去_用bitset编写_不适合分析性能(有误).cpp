#include<iostream>
#include<fstream>
#include<sstream>
#include<string.h>
#include<Windows.h>  //�߾��ȼ�ʱ
#include<bitset>   //������stl����λͼ��ʾ�����ÿһ�У������������
using namespace std;
#define E_LineN 453 //����Ԫ���� line number
#define RN 2362  //�������� row number


bitset<RN> eliminer[RN]; //��Ϊ����Ԫ����NUM�У�������Ԫ������ʱ���NUM��

struct line {  //������Ԫ��,���������׸�1�ı������Ա������Ԫ����
    //�Ƿ�����
    bool ifUprade;
    //�׸�1��λ��
    int num;
    //λͼ��ʾһ��
    bitset<RN>bit;
};
line eline[E_LineN];//����Ԫ��

void elimination() {
    int o = 5;  //��Ԫ��һ��ĸ���
    for (int i = RN - 1; i >= 0; i -= o) { //��Ԫ��5��һ��������Ԫ
        if (i - o < 0)o = i + 1;//���һ��ֱ��ȫ��һ����Ԫ��Ȼ�����ѭ��
        for (int j = 0; j < E_LineN; j++) { //�Ե�j������Ԫ�н�����Ԫ
            if (eline[j].ifUprade == 0) {//û�б����񣬽����������
                for (int k = 0; k < o; k++) {//����Ԫ��1��ǰ�Ŀ�ʼ��Ԫ������Ϊ����
                    if (eline[j].bit[i - k] == 1) {//�������Ԫ��i-kλ��Ϊ1���������eliminer[i-k]��Ԫ
                        eline[j].bit ^= eliminer[i - k];//���
                    }
                    //�Ƚ������������ټ�鱻��Ԫ�д�ʱΪ�ջ��߱�����Ϊ��Ԫ��
                    if (eline[j].bit.none()) {//�˱���Ԫ��Ϊ��ֵ����Ԫ��ϣ��˳�����ѭ��,������һ������Ԫ��
                        eline[j].ifUprade = true;
                        break;
                    }
                    else {//����Ƿ�����
                        for (int n = eline[j].num; n >= 0; n--) {//�ҵ��׸�1
                            if (eline[j].bit[n] == 1) {
                                eline[j].num = n;
                                break;
                            }
                        }
                        if (eliminer[eline[j].num].none()) {
                            if (eline[j].num >= i - 4 && eline[j].num <= i)
                            {
                                eliminer[eline[j].num] = eline[j].bit;
                                eline[j].ifUprade = true;
                                break;//��һ������Ԫ��
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
}
int main() {
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
            eliminer[row].set(tmp);
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
            eline[x].bit.set(tmp);
        }
        x++;
    }
    file_eline.close();
    elimination();

    //��������Ԫ����Ϊ���д����Ԫ�����
    ofstream out;
    out.open("��Ԫ���2.txt", ios::trunc);
    for (int i = 0; i < E_LineN; i++) {
        if (eline[i].bit.none()) { continue; }
        cout << i << endl;
        string s;
        string stemp;
        stringstream ss;
        for (int j = RN - 1; j >= 0; j--) {
            if (eline[i].bit[j] == 1)
            {
                ss.clear();
                ss << j;
                ss >> stemp;
                s.append(stemp); s.append(" ");
            }
        }
        out << s << "\n";
        s.clear();
    }
}
