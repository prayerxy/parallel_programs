#include <mpi.h>
#include<omp.h>
#include <stdio.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#define thread_nums 6
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
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided < MPI_THREAD_MULTIPLE)
        MPI_Abort(MPI_COMM_WORLD, 1);
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
    int i, j, k, tmp;

#pragma omp parallel num_threads(thread_nums),private(i,j,k,tmp)
    for (k = 0; k < N; k++) {
#pragma omp single
        {
            if (k % numprocs == myid) {
                for (j = k + 1; j < N; j++)
                    m[k][j] = m[k][j] / m[k][k];
                m[k][k] = 1.0;

                for (tmp = 0; tmp < numprocs; tmp++)
                    if (tmp != myid)
                        MPI_Send(&m[k][0], N, MPI_FLOAT, tmp, 2, MPI_COMM_WORLD);
            }
            else {
                MPI_Recv(&m[k][0], N, MPI_FLOAT, k % numprocs, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
#pragma omp for      
        for (i = k + 1; i < N; i++) {
            if (myid == i % numprocs) {
                for (j = k + 1; j < N; j++)
                    m[i][j] = m[i][j] - m[i][k] * m[k][j];
                m[i][k] = 0.0;
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == 0) {
        /*   for (int i = 0; i < N; i++) {
               for (int j = 0; j < N; j++)
                   cout << m[i][j] << " ";
               cout << endl;
           }*/

        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        cout << "程序用时:" << ((tail - head) * 1000000.0 / freq) << "微秒" << endl;
    }

    MPI_Finalize();

    return 0;
}
