#include<iostream>
#include<fstream>
#include<sstream>
#include<string.h>
#include<Windows.h>  

#include<vector>
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
void elimination() {
    int o = 5;  //消元子一组的个数
    for (int i = RN - 1; i >= 0; i -= o) { //消元子5个一组依次消元
        if (i - o < 0)o = i + 1;//最后一次直接全部一起消元，然后结束循环
        for (int j = 0; j < E_LineN; j++) { //对第j个被消元行进行消元
            if (eline[j].ifUprade == 0) {//没有被升格或者不为空，进行下面操作
                for (int k = 0; k < o; k++) {//从消元子1靠前的开始消元，所以为倒序
                    if (eline[j].bit[(i - k) / 8] & (1 << (i - k) % 8)) {//如果被消元行i-k位置为1，则可以用eliminer[i-k]消元

                        int ss = 0;
                        for (ss; ss + 32 <= (i - k) / 8; ss += 32) { //32个字节，8个int
                            __m256i ve = _mm256_loadu_epi8(eline[j].bit + ss);
                            __m256i vr = _mm256_loadu_epi8(eliminer[i - k] + ss);
                            __m256i tmp = _mm256_xor_si256(ve, vr);
                            _mm256_storeu_epi8(eline[j].bit + ss, tmp);
                        }
                        for (ss; ss <= (i - k) / 8; ss++)  //第(i-k)/8个可能对不需要操作的位进行了异或
                            eline[j].bit[ss] ^= eliminer[i - k][ss];


                    }
                    //先进行异或操作，再检查被消元行此时为空或者被升格为消元子
                    if (eline_ifnull(j)) {//此被消元行j为空值，消元完毕，退出此轮循环,进行下一个被消元行
                        eline[j].ifUprade = true;
                        break;
                    }
                    else {//检查是否升格
                        for (int n = eline[j].num; n >= 0; n--) {//找到首个1
                            if (eline[j].bit[n / 8] & (1 << (n % 8))) {  //n%8为0-7的位置，n/8为char的位置
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
    //消元子读入
    ifstream file_eliminer;
    file_eliminer.open("消元子.txt", ios_base::in);
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
    file_eline.open("被消元行.txt", ios_base::in);
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
    ofstream out;
    out.open("消元结果1.txt", ios::trunc);
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
