#include<iostream>
#include<fstream>
#include<sstream>
#include<string.h>
#include<Windows.h>
#include"cuda_runtime.h"
#include"device_launch_parameters.h"
#include"device_functions.h"
#define CUDA_CHECK(call) \
do { \
    cudaError_t result = call; \
    if (result != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d code=%d(%s) \"%s\" \n", \
                __FILE__, __LINE__, result, cudaGetErrorString(result), #call); \
        exit(EXIT_FAILURE); \
    } \
} while(0)

using namespace std;
#define E_LineN 756 //被消元行数 line number
#define RN 85401 //矩阵列数 row number

char* eliminer;
struct line {
    bool ifUprade;
    int num;
    char bit[RN / 8 + 1];
};
line* eline;
int flagUpgrade;
int nextstart;
void init() {

    for (int i = 0; i < E_LineN; i++) {
        for (int j = RN / 8; j >= 0; j--)
            eline[i].bit[j] = 0;
    }
    for (int i = 0; i < RN; i++) {
        for (int j = RN / 8; j >= 0; j--)
            eliminer[i * (RN / 8 + 1) + j] = 0;
    }
}
//在当前消元子不为空，传入消元子下标i，GPU开始消元工作，每个线程得到整体的索引t_id，stride是所有线程数量，作为步长
__global__ void eliminate1_kernel(char* eliminer, line* eline, int i) {
    int t_id = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;
    for (int j = t_id; j < E_LineN; j += stride) {
        if (eline[j].num == i) {
            int ss = 0;
            for (ss; ss <= i / 8; ss++)
                eline[j].bit[ss] ^= eliminer[i * (RN / 8 + 1) + ss];
            bool flag2 = true;
            for (int j2 = RN / 8; j2 >= 0; j2--) {
                if (!(eline[j].bit[j2] == 0)) {
                    flag2 = false;
                    break;
                }
            }
            if (flag2) {
                eline[j].ifUprade = true;
                eline[j].num = -2;
            }
            else {
                for (int n = eline[j].num; n >= 0; n--) {
                    if (eline[j].bit[n / 8] & (1 << (n % 8))) {
                        eline[j].num = n;
                        break;
                    }
                }
            }
        }
    }

}
//在当前消元子为空时，cpu首先完成升格操作，传入当前是否升格成功标志flagUpgrade
//在成功升格后，需要gpu完成剩余的被消元行的消去工作，其中nextstart是开始的被消元行的下标
__global__ void eliminate2_kernel(char* eliminer, line* eline, int i, int flagUpgrade, int nextstart) {
    int t_id = threadIdx.x + blockIdx.x * blockDim.x;//全局索引
    int stride = blockDim.x * gridDim.x;  //步长
    int newid = nextstart + 1;
    for (int j = newid + t_id; j < E_LineN; j += stride) {
        if (eline[j].num == i) {
            int ss = 0;
            for (ss; ss <= i / 8; ss++)  //第i/8个可能对不需要操作的位进行了异或
                eline[j].bit[ss] ^= eliminer[i * (RN / 8 + 1) + ss];
            bool flag2 = true;
            for (int j2 = RN / 8; j2 >= 0; j2--) {
                if (!(eline[j].bit[j2] == 0)) {
                    flag2 = false;
                    break;
                }
            }
            if (flag2) {
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

int main() {

    cudaMallocManaged(&eliminer, sizeof(char) * RN * (RN / 8 + 1));
    cudaMallocManaged((void**)&eline, sizeof(line) * E_LineN);
   
    init();
    // 消元子读入
    ifstream file_eliminer;
    //file_eliminer.open("特殊高斯消去的data//测试样例1 矩阵列数130，非零消元子22，被消元行8//消元子.txt", ios_base::in);
    //file_eliminer.open("特殊高斯消去的data//测试样例2 矩阵列数254，非零消元子106，被消元行53//消元子.txt", ios_base::in);
    //file_eliminer.open("特殊高斯消去的data//测试样例3 矩阵列数562，非零消元子170，被消元行53//消元子.txt", ios_base::in);
    //file_eliminer.open("特殊高斯消去的data//测试样例4 矩阵列数1011，非零消元子539，被消元行263//消元子.txt", ios_base::in);
   // file_eliminer.open("特殊高斯消去的data//测试样例5 矩阵列数2362，非零消元子1226，被消元行453//消元子.txt", ios_base::in);
    //file_eliminer.open("特殊高斯消去的data//测试样例6 矩阵列数3799，非零消元子2759，被消元行1953//消元子.txt", ios_base::in);
    //file_eliminer.open("特殊高斯消去的data//测试样例7 矩阵列数8399，非零消元子6375，被消元行4535//消元子.txt", ios_base::in);
   // file_eliminer.open("特殊高斯消去的data//测试样例8 矩阵列数23075，非零消元子18748，被消元行14325//消元子.txt", ios_base::in);
    //file_eliminer.open("特殊高斯消去的data//测试样例9 矩阵列数37960，非零消元子29304，被消元行14291//消元子.txt", ios_base::in);
    file_eliminer.open("特殊高斯消去的data//测试样例11 矩阵列数85401，非零消元子5724，被消元行756//消元子.txt", ios_base::in);
    string s1;
    if (!file_eliminer.is_open())
        cout << " 打开失败";
    while (getline(file_eliminer, s1)) {
        stringstream st;
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
            eliminer[row * (RN / 8 + 1) + i] |= (1 << j);
        }
    }
    file_eliminer.close();
    // 被消元行的读入
    ifstream file_eline;
    
    //file_eline.open("特殊高斯消去的data//测试样例1 矩阵列数130，非零消元子22，被消元行8//被消元行.txt", ios_base::in);
    //file_eline.open("特殊高斯消去的data//测试样例3 矩阵列数562，非零消元子170，被消元行53//被消元行.txt", ios_base::in);
    //file_eline.open("特殊高斯消去的data//测试样例4 矩阵列数1011，非零消元子539，被消元行263//被消元行.txt", ios_base::in);
   // file_eline.open("特殊高斯消去的data//测试样例5 矩阵列数2362，非零消元子1226，被消元行453//被消元行.txt", ios_base::in);
   // file_eline.open("特殊高斯消去的data//测试样例6 矩阵列数3799，非零消元子2759，被消元行1953//被消元行.txt", ios_base::in);
    //file_eline.open("特殊高斯消去的data//测试样例7 矩阵列数8399，非零消元子6375，被消元行4535//被消元行.txt", ios_base::in);
    //file_eline.open("特殊高斯消去的data//测试样例8 矩阵列数23075，非零消元子18748，被消元行14325//被消元行.txt", ios_base::in);
   // file_eline.open("特殊高斯消去的data//测试样例9 矩阵列数37960，非零消元子29304，被消元行14291//被消元行.txt", ios_base::in);
    file_eline.open("特殊高斯消去的data//测试样例11 矩阵列数85401，非零消元子5724，被消元行756//被消元行.txt", ios_base::in);
    string s2;
    if (!file_eline.is_open())
        cout << " 打开失败";
    int x = 0;
    while (getline(file_eline, s2)) {
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

    // 创建CUDA事件
    cudaEvent_t start, stop;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));
    CUDA_CHECK(cudaEventRecord(start, 0));

    for (int i = RN - 1; i >= 0; i--) {
        cudaDeviceSynchronize();
        bool flag1 = (eliminer[i * (RN / 8 + 1) + i / 8] == 0) ? true : false;//判断当前消元子是否为空

        if (!flag1) { //不为空时
            eliminate1_kernel << <128, 1024 >> > (eliminer, eline, i);
        }
        else {//为空时，需要cpu升格
            flagUpgrade = 0;
            nextstart = 0;
            for (int j = 0; j < E_LineN; j++) {
                if (eline[j].num == i) {
                    for (int ss = eline[j].num / 8; ss >= 0; ss--)  //赋值过来
                        eliminer[eline[j].num * (RN / 8 + 1) + ss] = eline[j].bit[ss];
                    eline[j].ifUprade = true;
                    flagUpgrade = 1;
                    nextstart = j;//之前的被消元行的首位不等于i,不用消元
                    eline[j].num = -2;
                    break;
                }
            }
            if (flagUpgrade == 0)
                continue;
            eliminate2_kernel << <32, 512 >> > (eliminer, eline, i, flagUpgrade, nextstart);

        }
    }
    // 销毁CUDA事件
    CUDA_CHECK(cudaEventRecord(stop, 0));
    CUDA_CHECK(cudaEventSynchronize(stop));
    float elapsedTime;
    CUDA_CHECK(cudaEventElapsedTime(&elapsedTime, start, stop));
    //记录时间
    cout << "程序用时：" << elapsedTime << "毫秒" << endl;

    // 把最后的消元子作为结果写入消元结果中
    ofstream out;
    out.open("消元结果1.txt", ios::trunc);
    for (int i = 0; i < E_LineN; i++) {
        string s;
        string stemp;
        stringstream ss;
        for (int j = RN / 8; j >= 0; j--) {
            for (int k = 7; k >= 0; k--) {
                if (eline[i].bit[j] & (1 << k)) {
                    ss.clear();
                    ss << (j * 8 + k);
                    ss >> stemp;
                    s.append(stemp);
                    s.append(" ");
                }
            }
        }
        out << s << "\n";
        s.clear();
    }

    // 释放内存
    cudaFree(eliminer);
    cudaFree(eline);
    
    return 0;
}
