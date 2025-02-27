#include<iostream>
#include<fstream>
#include<sstream>
#include<string.h>
#include<Windows.h>  
#include<immintrin.h>


using namespace std;
#define E_LineN 756 //被消元行数 line number
#define RN 85401//矩阵列数 row number


char** eliminer = new char* [RN]; //因为被消元行有RN列，所以消元子最终时最多RN个

struct line {  //定义消元行,增加升格，首个1的变量，以便等下消元处理
    //是否升格
    bool ifUprade;
    //首个1的位置
    int num;
    //位图表示一行
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
    if (eliminer[i][i / 8] == 0)  //消元子第i行的第i位置一定有1，如果没有说明为空
        return true;
    return false;
}
int nextstart = 0;
int flagUpgrade = 0;  //全局变量，不能作为线程函数局部变量
void elimination() {
    for (int i = RN - 1; i >= 0; i--) { //消元子
        if (!eliminer_ifnull(i)) {
            for (int j = 0; j < E_LineN; j++) { //对第j个被消元行进行消元
                if (eline[j].num == i) {
                    int ss = 0;

                    for (ss; ss + 64 <= i / 8; ss += 64) { //32个字节，8个int
                        __m512i ve = _mm512_loadu_epi8(eline[j].bit + ss);
                        __m512i vr = _mm512_loadu_epi8(eliminer[i] + ss);
                        __m512i tmp = _mm512_xor_si512(ve, vr);
                        _mm512_storeu_epi8(eline[j].bit + ss, tmp);
                    }

                    for (ss; ss <= i / 8; ss++)  //第i/8个可能对不需要操作的位进行了异或
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
        else {


            for (int j = 0; j < E_LineN; j++) {
                if (eline[j].num == i) {
                    for (int ss = eline[j].num / 8; ss >= 0; ss--)  //赋值过来
                        eliminer[eline[j].num][ss] = eline[j].bit[ss];
                    eline[j].ifUprade = true;
                    flagUpgrade = 1;
                    nextstart = j;//之前的被消元行的首位不等于i,不用消元
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

                    for (ss; ss + 64 <= i / 8; ss += 64) { //32个字节，8个int
                        __m512i ve = _mm512_loadu_epi8(eline[j].bit + ss);
                        __m512i vr = _mm512_loadu_epi8(eliminer[i] + ss);
                        __m512i tmp = _mm512_xor_si512(ve, vr);
                        _mm512_storeu_epi8(eline[j].bit + ss, tmp);
                    }

                    for (ss; ss + 16 <= i / 8; ss += 16) { //32个字节，8个int
                        __m128i ve = _mm_loadu_epi8(eline[j].bit + ss);
                        __m128i vr = _mm_loadu_epi8(eliminer[i] + ss);
                        __m128i tmp = _mm_xor_si128(ve, vr);
                        _mm_storeu_epi8(eline[j].bit + ss, tmp);

                    }

                    for (ss; ss <= i / 8; ss++)  //第i/8个可能对不需要操作的位进行了异或
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
}
int main() {
    init();
    //消元子读入
    ifstream file_eliminer;
    //file_eliminer.open("特殊高斯消去的data//测试样例1 矩阵列数130，非零消元子22，被消元行8//消元子.txt", ios_base::in);
    //file_eliminer.open("特殊高斯消去的data//测试样例2 矩阵列数254，非零消元子106，被消元行53//消元子.txt", ios_base::in);
    //file_eliminer.open("特殊高斯消去的data//测试样例3 矩阵列数562，非零消元子170，被消元行53//消元子.txt", ios_base::in);
     //file_eliminer.open("特殊高斯消去的data//测试样例4 矩阵列数1011，非零消元子539，被消元行263//消元子.txt", ios_base::in);
     //file_eliminer.open("特殊高斯消去的data//测试样例5 矩阵列数2362，非零消元子1226，被消元行453//消元子.txt", ios_base::in);
    //file_eliminer.open("特殊高斯消去的data//测试样例6 矩阵列数3799，非零消元子2759，被消元行1953//消元子.txt", ios_base::in);
    //file_eliminer.open("特殊高斯消去的data//测试样例7 矩阵列数8399，非零消元子6375，被消元行4535//消元子.txt", ios_base::in);
    //file_eliminer.open("特殊高斯消去的data//测试样例8 矩阵列数23075，非零消元子18748，被消元行14325//消元子.txt", ios_base::in);
    //file_eliminer.open("特殊高斯消去的data//测试样例9 矩阵列数37960，非零消元子29304，被消元行14291//消元子.txt", ios_base::in);
    file_eliminer.open("特殊高斯消去的data//测试样例11 矩阵列数85401，非零消元子5724，被消元行756//消元子.txt", ios_base::in);
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
    //被消元行的读入
    ifstream file_eline;
     //file_eline.open("特殊高斯消去的data//测试样例1 矩阵列数130，非零消元子22，被消元行8//被消元行.txt", ios_base::in);
     //file_eline.open("特殊高斯消去的data//测试样例2 矩阵列数254，非零消元子106，被消元行53//被消元行.txt", ios_base::in);
    //file_eline.open("特殊高斯消去的data//测试样例3 矩阵列数562，非零消元子170，被消元行53//被消元行.txt", ios_base::in);
    //file_eline.open("特殊高斯消去的data//测试样例4 矩阵列数1011，非零消元子539，被消元行263//被消元行.txt", ios_base::in);
     //file_eline.open("特殊高斯消去的data//测试样例5 矩阵列数2362，非零消元子1226，被消元行453//被消元行.txt", ios_base::in);
    //file_eline.open("特殊高斯消去的data//测试样例6 矩阵列数3799，非零消元子2759，被消元行1953//被消元行.txt", ios_base::in);
    //file_eline.open("特殊高斯消去的data//测试样例7 矩阵列数8399，非零消元子6375，被消元行4535//被消元行.txt", ios_base::in);
    //file_eline.open("特殊高斯消去的data//测试样例8 矩阵列数23075，非零消元子18748，被消元行14325//被消元行.txt", ios_base::in);
    //file_eline.open("特殊高斯消去的data//测试样例9 矩阵列数37960，非零消元子29304，被消元行14291//被消元行.txt", ios_base::in);
    file_eline.open("特殊高斯消去的data//测试样例11 矩阵列数85401，非零消元子5724，被消元行756//被消元行.txt", ios_base::in);

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
    long long head, tail, freq; //timers
    //记时间
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    //start time
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    elimination();
    //end time
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout << "程序用时:" << ((tail - head) * 1000000.0 / freq) << "微秒" << endl;
    //把最后的消元子作为结果写入消元结果中

}
