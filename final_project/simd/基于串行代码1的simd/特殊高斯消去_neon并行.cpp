#include<iostream>
#include<fstream>
#include<sstream>
#include<string.h>
#include <arm_neon.h>
#include <sys/time.h>
using namespace std;
#define E_LineN 8 //被消元行数 line number
#define RN 130  //矩阵列数 row number


char** eliminer = new char* [RN];

struct line {

    bool ifUprade;

    int num;

    char* bit = new char[RN / 8 + 1];
};
line eline[E_LineN];

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
bool eline_ifnull(int i) {
    bool flag = true;
    for (int j = RN / 8; j >= 0; j--) {
        if (!eline[i].bit[j] == 0)
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
void elimination() {
    int o = 5;
    for (int i = RN - 1; i >= 0; i -= o) {
        if (i - o < 0)o = i + 1;
        for (int j = 0; j < E_LineN; j++) {
            if (eline[j].ifUprade == 0) {
                for (int k = 0; k < o; k++) {
                    if (eline[j].bit[(i - k) / 8] & (1 << (i - k) % 8)) {

                        int ss = 0;
                        for (ss; ss + 16 <= (i - k) / 8; ss += 16) { //32个字节，8个int
                            uint8x16_t ve = vld1q_u8((const uint8_t*)&eline[j].bit[ss]);
                            uint8x16_t vr = vld1q_u8((const uint8_t*)&eliminer[i - k][ss]);


                            uint8x16_t tmp = veorq_u8(ve, vr);
                            vst1q_u8((uint8_t*)&eline[j].bit[ss], tmp);
                        }
                        for (ss; ss <= (i - k) / 8; ss++)
                            eline[j].bit[ss] ^= eliminer[i - k][ss];


                    }
                    if (eline_ifnull(j)) {
                        eline[j].ifUprade = true;
                        break;
                    }
                    else {
                        for (int n = eline[j].num; n >= 0; n--) {
                            if (eline[j].bit[n / 8] & (1 << (n % 8))) {
                                eline[j].num = n;
                                break;
                            }
                        }
                        if (eliminer_ifnull(eline[j].num)) {
                            if (eline[j].num >= i - 4 && eline[j].num <= i)
                            {
                                for (int ss = eline[j].num / 8; ss >= 0; ss--)  //赋值过来
                                    eliminer[eline[j].num][ss] = eline[j].bit[ss];
                                //eliminer[eline[j].num] ^= eline[j].bit;
                                eline[j].ifUprade = true;
                                break;//下一个被消元行
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
    init();

    ifstream file_eliminer;
    file_eliminer.open("/home/data/Groebner/1_130_22_8/1.txt", ios_base::in);
    string s1;
    if (!file_eliminer.is_open())
    {
        return 1;
        cout << " 打开失败";
    }
    while (getline(file_eliminer, s1)) { //逐行读入
        stringstream st;//从字符串转到数字
        st << s1;
        int tmp;
        int row;
        bool first = true;
        while (st >> tmp) {
            if (first) {
                row = tmp;
                first = false;
            }
            int i = tmp / 8;
            int j = tmp % 8;
            eliminer[row][i] |= (1 << j);
        }
    }
    file_eliminer.close();
    ifstream file_eline;
    file_eline.open("/home/data/Groebner/1_130_22_8/2.txt", ios_base::in);
    string s2;
    if (!file_eline.is_open())
    {
        return 1;
        cout << " 打开失败";
    }
    int x = 0;
    while (getline(file_eline, s2)) { //逐行读入
        stringstream st;
        st << s2;
        int tmp;
        bool first = true;
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
    cout << "程序用时:" << timecount << "微秒" << endl;

}
