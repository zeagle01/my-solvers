

#include <cuda_runtime.h>
#include "device_launch_parameters.h"
#include "gpu_les.h"
#include <iostream>
/*
typedef float T;

void GPU_GS_Solver<T>::solve(T * A, T * IA, T * JA, T * X) {

}
*/

__global__ void add_kernel(float * a, float * b, float* c, int N){
	int tid = threadIdx.x + blockIdx.x*blockDim.x;
	if (tid < N){
		c[tid] = a[tid] + b[tid];
	}
	
}


__global__ void test(){
	int tid = blockDim.x*blockIdx.x + threadIdx.x;
	printf("I'm in kernel");
	printf("tid=%d",tid);
}
void GPU_GS_Solver<float>::solve(float * A, float * IA, float * JA, float * X) {
	test << <1, 64 >> >();
	
	std::cout << "hi" << std::endl;
}
