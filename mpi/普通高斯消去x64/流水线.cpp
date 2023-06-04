#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
using namespace std;
const int N = 3200;
float m[N][N];
void m_reset()
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < i; j++)
            m[i][j] = 0;
        m[i][i] = 1.0;
        for (int j = i + 1; j < N; j++)
        {
            m[i][j] = rand();
        }
    }
    for (int k = 0; k < N; k++)
        for (int i = k + 1; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                if (k && k % 8 == 0)
                    m[i][j] -= m[k][j];
                else
                    m[i][j] += m[k][j];
            }
}
int main(int argc, char* argv[]) {
    int myid, numprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    if (myid == 0) {
        m_reset();
        for (int i = 1; i < numprocs; i++)
            MPI_Send(m, N * N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
    }
    else
        MPI_Recv(m, N * N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Barrier(MPI_COMM_WORLD);
    long long head, tail, freq; //timers
    //记时间
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    //start time
    QueryPerformanceCounter((LARGE_INTEGER*)&head);

    //一维行划分
    for (int k = 0; k < N; k++) {
        if (k % numprocs == myid) {
            for (int j = k + 1; j < N; j++)
                m[k][j] = m[k][j] / m[k][k];
            m[k][k] = 1.0;
            MPI_Send(&m[k][0], N, MPI_FLOAT, (myid == numprocs - 1) ? 0 : myid + 1, 2, MPI_COMM_WORLD);
        }
        else {
            MPI_Recv(&m[k][0], N, MPI_FLOAT, (myid == 0) ? numprocs - 1 : myid - 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (k % numprocs != myid + 1 && k % numprocs) {
                MPI_Send(&m[k][0], N, MPI_FLOAT, (myid == numprocs - 1) ? 0 : myid + 1, 2, MPI_COMM_WORLD);
            }
            if (k % numprocs == 0 && myid != numprocs - 1) {
                MPI_Send(&m[k][0], N, MPI_FLOAT, (myid == numprocs - 1) ? 0 : myid + 1, 2, MPI_COMM_WORLD);
            }
        }

        for (int i = k + 1; i < N; i++) {
            if (myid == i % numprocs) {
                for (int j = k + 1; j < N; j++)
                    m[i][j] = m[i][j] - m[i][k] * m[k][j];
                m[i][k] = 0.0;
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (myid == 0) {
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        cout << "程序用时:" << ((tail - head) * 1000000.0 / freq) << "微秒" << endl;

    }

    MPI_Finalize();

    return 0;
}
