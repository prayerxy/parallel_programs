#include"cuda_runtime.h"
#include"device_launch_parameters.h"
#include<stdio.h>
#include<iostream>
#include"device_functions.h"
const int N = 1024;
float* m;
void reset() {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < i; j++)
			m[i * N + j] = 0;
		m[i * N + i] = 1.0;
		for (int j = i + 1; j < N; j++)
			m[i * N + j] = rand();
	}
		for (int k = 0; k < N; k++)
			for (int i = k + 1; i < N; i++)
				for (int j = 0; j < N; j++)
					m[i * N + j] += m[k * N + j];
	
}
__global__ void division_kernel(float* m, int k) {
	int tid = blockDim.x * blockIdx.x + threadIdx.x;//计算线程索引
	int element = m[k * N + k];
	int temp = m[k * N + tid];
	m[k * N + tid] = (float)temp / element;

}
//除法部分
__global__ void eliminate_kernel(float* m, int k) {
	int tk = blockDim.x * blockIdx.x + threadIdx.x;
	if (tk == 0)
		m[k * N + k] = 1.0;//对角线元素设为1
	int row = k + 1 + blockIdx.x;//每个块负责一行
	while (row < N) {
		int tid = threadIdx.x;
		while (k + 1 + tid < N) {
			int col = k + 1 + tid;
			float temp_1 = m[(row * N) + col];
			float temp_2 = m[(row * N) + k];
			float temp_3 = m[k * N + col];
			m[(row * N) + col] = temp_1 - temp_2 * temp_3;
			tid = tid + blockDim.x;
		}
		__syncthreads();//同步
		if (threadIdx.x == 0)
			m[row * N + k] = 0;
		row += gridDim.x;
	}
	return;
		
}
//消去部分
int main() {
	float timecount;
	size_t size = N * N * sizeof(float);
	cudaMallocManaged(&m, size);
	reset();
	cudaEvent_t start, stop;
	float elapsedTime = 0.0;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);//开始计时
	cudaError_t ret;
	for (int k = 0; k < N; k++) {
		division_kernel<<<1,1024>>>(m, k);
		cudaDeviceSynchronize();
		ret = cudaGetLastError();
		if (ret != cudaSuccess)
			printf("division_kernel failed,%s\n", cudaGetErrorString(ret));
		eliminate_kernel <<<128,1024>>> (m, k);
		cudaDeviceSynchronize();
		ret = cudaGetLastError();
		if(ret!=cudaSuccess)
			printf("division_kernel failed,%s\n", cudaGetErrorString(ret));

	}
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);
	cudaFree(m);
	printf("GPU_LU:%f ms\n", elapsedTime);
	return 0;
}