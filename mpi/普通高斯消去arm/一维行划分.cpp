#include <mpi.h>
#include <sys/time.h>
#include <stdio.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
using namespace std;
const int N = 400;
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
            m[i][j] = rand() % 1000;
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
    struct timeval start;
    struct timeval end;//clock
    float timecount;
    gettimeofday(&start, NULL);
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
                for (int j = k + 1; j < N; j++)
                    m[i][j] = m[i][j] - m[i][k] * m[k][j];
                m[i][k] = 0.0;
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == 0) {
        gettimeofday(&end, NULL);
        timecount = (end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec;
        cout << "³ÌÐòÓÃÊ±:" << timecount << "Î¢Ãë" << endl;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++)
                cout << m[i][j] << " ";
            cout << endl;
        }
    }

    MPI_Finalize();

    return 0;
}
