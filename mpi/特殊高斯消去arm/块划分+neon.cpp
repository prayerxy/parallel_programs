#include<iostream>
#include <mpi.h>
#include<fstream>
#include<sstream>
#include<string.h>
#include <arm_neon.h>
#include <sys/time.h>


using namespace std;
#define E_LineN 8//����Ԫ���� line number
#define RN 130//�������� row number


char eliminer[RN][RN / 8 + 1]; //��Ϊ����Ԫ����RN�У�������Ԫ������ʱ���RN��

struct line {

    bool ifUprade;

    int num;

    char bit[RN / 8 + 1];
};
line eline[E_LineN];//����Ԫ��

void init() {

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
int nextstart = -1;
int flagUpgrade = 0;

int main(int argc, char* argv[]) {
    int myid, numprocs;

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
            int j = tmp % 8;
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
    while (getline(file_eline, s2)) {
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

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    struct timeval start;
    struct timeval end;//clock
    float timecount;

    gettimeofday(&start, NULL);


    int r1, r2;
    int range = E_LineN / numprocs;
    r1 = myid * range;
    if (myid != numprocs - 1)
        r2 = r1 + range - 1;
    else
        r2 = E_LineN - 1;
    for (int i = RN - 1; i >= 0; i--) {
        if (!eliminer_ifnull(i)) {
            for (int j = r1; j <= r2; j++) {

                if (eline[j].num == i) {
                    int ss = 0;
                    for (ss; ss + 16 <= i / 8; ss += 16) { //32���ֽڣ�8��int
                        uint8x16_t ve = vld1q_u8((const uint8_t*)&eline[j].bit[ss]);
                        uint8x16_t vr = vld1q_u8((const uint8_t*)&eliminer[i][ss]);


                        uint8x16_t tmp = veorq_u8(ve, vr);
                        vst1q_u8((uint8_t*)&eline[j].bit[ss], tmp);
                    }
                    for (ss; ss <= i / 8; ss++)
                        eline[j].bit[ss] ^= eliminer[i][ss];
                    if (eline_ifnull(j)) {
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
        else {
            if (myid != 0)
                MPI_Send(&eline[r1], (r2 - r1 + 1) * sizeof(line), MPI_BYTE, 0, 1, MPI_COMM_WORLD);
            else {
                for (int tmps = 1; tmps < numprocs; tmps++) {
                    if (tmps != numprocs - 1)
                        MPI_Recv(&eline[tmps * range], range * sizeof(line), MPI_BYTE, tmps, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    else
                        MPI_Recv(&eline[tmps * range], sizeof(line) * (E_LineN - tmps * range), MPI_BYTE, tmps, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            if (myid == 0) {
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

                for (int tmps = 1; tmps < numprocs; tmps++)
                {
                    MPI_Send(&nextstart, 1, MPI_INT, tmps, 1, MPI_COMM_WORLD);
                    MPI_Send(&eliminer[i], RN / 8 + 1, MPI_CHAR, tmps, 1, MPI_COMM_WORLD);
                }
            }
            else {
                MPI_Recv(&nextstart, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&eliminer[i], RN / 8 + 1, MPI_CHAR, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            if (nextstart == -1)
                continue;
            for (int j = nextstart + 1; j < E_LineN; j++) {
                if (r1 <= j && j <= r2 && eline[j].num == i) {
                    int ss = 0;
                    for (ss; ss + 16 <= i / 8; ss += 16) { //32���ֽڣ�8��int
                        uint8x16_t ve = vld1q_u8((const uint8_t*)&eline[j].bit[ss]);
                        uint8x16_t vr = vld1q_u8((const uint8_t*)&eliminer[i][ss]);


                        uint8x16_t tmp = veorq_u8(ve, vr);
                        vst1q_u8((uint8_t*)&eline[j].bit[ss], tmp);
                    }
                    for (ss; ss <= i / 8; ss++)  //��i/8�����ܶԲ���Ҫ������λ���������
                        eline[j].bit[ss] ^= eliminer[i][ss];
                    if (eline_ifnull(j)) {
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
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == 0) {
        gettimeofday(&end, NULL);
        timecount = (end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec;
        cout << "������ʱ:" << timecount << "΢��" << endl;

    }

    MPI_Finalize();

    return 0;
}
