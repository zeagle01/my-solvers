#include <cuda_runtime.h>
#include "device_launch_parameters.h"
#include "gpu_solver.h"

//typedef float T;

void GPU_DiffusionSolver<float>::callGPU(){

}
void GPU_DiffusionSolver<float>::copyH2D(){
	cudaMalloc(&d_X, sizeof(float)*mesh.cellNum);
	cudaMemcpy(d_X, X, sizeof(float)*mesh.cellNum, cudaMemcpyHostToDevice);

	cudaMalloc(&d_A, sizeof(float)*mesh.JA.size());
	cudaMemcpy(d_A, A, sizeof(float)*mesh.JA.size(), cudaMemcpyHostToDevice);
	cudaMalloc(&d_JA, sizeof(float)*mesh.JA.size());
	cudaMemcpy(d_JA, JA, sizeof(float)*mesh.JA.size(), cudaMemcpyHostToDevice);
	cudaMalloc(&d_IA, sizeof(float)*mesh.IA.size());
	cudaMemcpy(d_IA, IA, sizeof(float)*mesh.IA.size(), cudaMemcpyHostToDevice);
}
void GPU_DiffusionSolver<float>::copyD2H(){
	cudaMemcpy(X, d_X, sizeof(float)*mesh.cellNum, cudaMemcpyDeviceToHost);
}