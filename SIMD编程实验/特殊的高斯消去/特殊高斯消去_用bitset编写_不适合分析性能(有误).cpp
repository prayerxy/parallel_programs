#include<iostream>
#include<fstream>
#include<sstream>
#include<string.h>
#include<Windows.h>  //高精度计时
#include<bitset>   //我们用stl内置位图表示矩阵的每一行，方便异或运算
using namespace std;
#define E_LineN 453 //被消元行数 line number
#define RN 2362  //矩阵列数 row number


bitset<RN> eliminer[RN]; //因为被消元行有NUM列，所以消元子最终时最多NUM个

struct line {  //定义消元行,增加升格，首个1的变量，以便等下消元处理
    //是否升格
    bool ifUprade;
    //首个1的位置
    int num;
    //位图表示一行
    bitset<RN>bit;
};
line eline[E_LineN];//被消元行

void elimination() {
    int o = 5;  //消元子一组的个数
    for (int i = RN - 1; i >= 0; i -= o) { //消元子5个一组依次消元
        if (i - o < 0)o = i + 1;//最后一次直接全部一起消元，然后结束循环
        for (int j = 0; j < E_LineN; j++) { //对第j个被消元行进行消元
            if (eline[j].ifUprade == 0) {//没有被升格，进行下面操作
                for (int k = 0; k < o; k++) {//从消元子1靠前的开始消元，所以为倒序
                    if (eline[j].bit[i - k] == 1) {//如果被消元行i-k位置为1，则可以用eliminer[i-k]消元
                        eline[j].bit ^= eliminer[i - k];//异或
                    }
                    //先进行异或操作，再检查被消元行此时为空或者被升格为消元子
                    if (eline[j].bit.none()) {//此被消元行为空值，消元完毕，退出此轮循环,进行下一个被消元行
                        eline[j].ifUprade = true;
                        break;
                    }
                    else {//检查是否升格
                        for (int n = eline[j].num; n >= 0; n--) {//找到首个1
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
            eliminer[row].set(tmp);
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
            eline[x].bit.set(tmp);
        }
        x++;
    }
    file_eline.close();
    elimination();

    //把最后的消元子作为结果写入消元结果中
    ofstream out;
    out.open("消元结果2.txt", ios::trunc);
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
