#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include<immintrin.h>
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
    //��ʱ��
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    //start time
    QueryPerformanceCounter((LARGE_INTEGER*)&head);

    //һά�л���
    for (int k = 0; k < N; k++) {
        if (k % numprocs == myid) {
            for (int j = k + 1; j < N; j++)
                m[k][j] = m[k][j] / m[k][k];
            m[k][k] = 1.0;

            for (int tmp = 0; tmp < numprocs; tmp++)
                if (tmp != myid)
                    MPI_Send(&m[k][0], N, MPI_FLOAT, tmp, 2, MPI_COMM_WORLD);
        }
        else {
            MPI_Recv(&m[k][0], N, MPI_FLOAT, k % numprocs, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        for (int i = k + 1; i < N; i++) {
            if (myid == i % numprocs) {
                __m256 vaik = _mm256_set1_ps(m[i][k]);
                int j = k + 1;
                for (j; j + 8 <= N; j += 8) {//�ڲ�չ������j2
                    __m256 vakj = _mm256_loadu_ps(m[k] + j);
                    __m256 vaij = _mm256_loadu_ps(m[i] + j);
                    vakj = _mm256_mul_ps(vaik, vakj);
                    vaij = _mm256_sub_ps(vaij, vakj);
                    _mm256_storeu_ps(m[i] + j, vaij);
                }
                for (j; j < N; j++) {
                    m[i][j] = m[i][j] - m[i][k] * m[k][j];
                }
                m[i][k] = 0.0;
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    /*  if (myid != 0) {
          for (int i = myid; i < N; i += numprocs)
              MPI_Send(&m[i][0], N, MPI_FLOAT, 0, i, MPI_COMM_WORLD);

      }
      else {
          for (int i = 1; i < N; i++)
              if (i % numprocs)
                  MPI_Recv(&m[i][0], N, MPI_FLOAT, i % numprocs, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      }*/
    if (myid == 0) {

        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        cout << "������ʱ:" << ((tail - head) * 1000000.0 / freq) << "΢��" << endl;
    }

    MPI_Finalize();

    return 0;
}
