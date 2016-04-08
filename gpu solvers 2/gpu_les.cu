

#include "case.h"
#include "Solver.h"
#include "configReader.h"
#include "linearEquationSolver.h"
#include <cuda_runtime.h>
#include "device_launch_parameters.h"


/*
__global__ void kernel(CSR* eq, CellField* phi, Mesh* mesh){
	int tid = blockDim.x*blockIdx.x + threadIdx.x;
	if (tid < mesh->cellNum){
		printf("numof%d", mesh->X[tid * 2+0]);
		printf("hello??");
	}
}
*/
/*
for (int c = 0; c < mesh.cellNum; c++) {
	double temp = 0;
	for (int i = mesh.IA[c] + 1; i < mesh.IA[c + 1]; i++) {
		temp += eq.A[i] * phi.inner[mesh.JA[i]];
	}
	phi.inner[c] = (eq.b[c] - temp) / eq.A[mesh.IA[c]];
}
*/
__global__ void kernel(double*A, double* b,int* IA, int *JA, double* X, double* preX,int N){
	int tid = blockDim.x*blockIdx.x + threadIdx.x;
	if (tid < N){
			double temp = 0;
			for (int i = IA[tid] + 1; i < IA[tid + 1]; i++){
				temp += A[i] * preX[JA[i]];
			}
			X[tid] = (b[tid] - temp) / A[IA[tid]];
	}
}

void GPU_Jacobi_Solver::solve(CSR eq, CellField& phi, Mesh mesh){

	/*
	CSR* d_eq;
	CellField* d_phi;
	Mesh* d_mesh;
	cudaMalloc(&d_eq, sizeof(CSR));
	cudaMalloc(&d_phi, sizeof(CellField));
	cudaMalloc(&d_mesh, sizeof(Mesh));
	cudaMemcpy(d_mesh, &mesh, sizeof(Mesh), cudaMemcpyHostToDevice);
	*/
	/*
	double *A =new double[eq.A.size()];
	for (int i = 0; i < eq.A.size(); i++){
		A[i] = eq.A[i];
	}
	*/
	double *A = eq.A.data();
	double *b = eq.b.data();
	int *IA = mesh.IA.data();
	int  *JA = mesh.JA.data();
	double* X = phi.inner.data();
	int N = mesh.cellNum;
	double*d_A,*d_b, *d_X, *d_preX;
	int*d_IA, *d_JA;


	
	cudaMalloc(&d_X, sizeof(double)*phi.inner.size());
	cudaMalloc(&d_preX, sizeof(double)*phi.inner.size());
	cudaMalloc(&d_A, sizeof(double)*eq.A.size());
	cudaMalloc(&d_b, sizeof(double)*eq.b.size());
	cudaMalloc(&d_IA, sizeof(int)*mesh.IA.size());
	cudaMalloc(&d_JA, sizeof(int)*mesh.JA.size());

	cudaMemcpy(d_A, A, sizeof(double)*eq.A.size(), cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, b, sizeof(double)*eq.b.size(), cudaMemcpyHostToDevice);
	cudaMemcpy(d_X, X, sizeof(double)*phi.inner.size(), cudaMemcpyHostToDevice);
	cudaMemcpy(d_preX, X, sizeof(double)*phi.inner.size(), cudaMemcpyHostToDevice);
	cudaMemcpy(d_IA, IA, sizeof(int)*mesh.IA.size(), cudaMemcpyHostToDevice);
	cudaMemcpy(d_JA, JA, sizeof(int)*mesh.JA.size(), cudaMemcpyHostToDevice);

	int tpb = 256;
	int bpg =(  N+ 256 - 1) / tpb;

	double error = 0; int step = 0;
	do{
		CellField pre = phi;

		for (int it = 0; it < check_step; it++){
			cudaMemcpy(d_preX, d_X, sizeof(double)*phi.inner.size(), cudaMemcpyDeviceToDevice);
			kernel << <bpg, tpb >> >(d_A, d_b, d_IA, d_JA, d_X, d_preX, N);
		}
		cudaMemcpy(X, d_X, sizeof(double)*phi.inner.size(), cudaMemcpyDeviceToHost);

		if (step%check_step == 0){
			error = VectorMath<double>::rootOfSquareSum(phi.inner, pre.inner);
			cout << step << endl;
		}
		step++;
		/*
		transform(phi->inner.cbegin(), phi->inner.cend(), pre.inner.cbegin(),pre.inner.begin(), minus<T>());
		transform(pre.inner.begin(), pre.inner.end(), pre.inner.begin(), pre.inner.begin(), multiplies<T>());
		error = accumulate(pre.inner.begin(), pre.inner.end(), (T)0);//have to force convert to T type unless the result is always 0
		//error = inner_product(pre.inner.cbegin(), pre.inner.cend(), pre.inner.cbegin(),
		//	1, plus<T>(), plus<T>());
		error = sqrt(error);
		*/
		
	} while (error>converge_threhold&&step<max_step);

	//phi.inner . assign(X, X + N);
	
}



//vecotr data structure

__global__ void vector_kernel(CSR eq,CellField phi, CellField pre_phi,Mesh mesh){
	int tid = blockDim.x*blockIdx.x + threadIdx.x;
	double *A = eq.A.data();
	double *b = eq.b.data();
	int* IA = mesh.IA.data();
	int* JA = mesh.JA.data();
	double* X = phi.inner.data();
	double*preX = pre_phi.inner.data();

	if (tid < mesh.cellNum){
		double temp = 0;
		X[tid];
		/*
		for (int i = IA[tid] + 1; i < IA[tid + 1]; i++){
			temp += A[i] * preX[JA[i]];
		}
		X[tid] = (b[tid] - temp) / A[IA[tid]];
		*/
	}
}


void GPU_Jacobi_Solver_vector::solve(CSR eq, CellField& phi, Mesh mesh){

	
	CSR d_eq;
	CellField d_phi,d_pre_phi;
	Mesh d_mesh;

	
	/*
	double *A =new double[eq.A.size()];
	for (int i = 0; i < eq.A.size(); i++){
	A[i] = eq.A[i];
	}
	*/
	double *A = eq.A.data();
	double *b = eq.b.data();
	int *IA = mesh.IA.data();
	int  *JA = mesh.JA.data();
	double* X = phi.inner.data();
	int N = mesh.cellNum;
	double*d_A, *d_b, *d_X, *d_preX;
	int*d_IA, *d_JA;



	cudaMalloc(&d_X, sizeof(double)*phi.inner.size());
	cudaMalloc(&d_preX, sizeof(double)*phi.inner.size());
	cudaMalloc(&d_A, sizeof(double)*eq.A.size());
	cudaMalloc(&d_b, sizeof(double)*eq.b.size());
	cudaMalloc(&d_IA, sizeof(int)*mesh.IA.size());
	cudaMalloc(&d_JA, sizeof(int)*mesh.JA.size());

	cudaMemcpy(d_A, A, sizeof(double)*eq.A.size(), cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, b, sizeof(double)*eq.b.size(), cudaMemcpyHostToDevice);
	cudaMemcpy(d_X, X, sizeof(double)*phi.inner.size(), cudaMemcpyHostToDevice);
	cudaMemcpy(d_preX, X, sizeof(double)*phi.inner.size(), cudaMemcpyHostToDevice);
	cudaMemcpy(d_IA, IA, sizeof(int)*mesh.IA.size(), cudaMemcpyHostToDevice);
	cudaMemcpy(d_JA, JA, sizeof(int)*mesh.JA.size(), cudaMemcpyHostToDevice);

	d_eq.A.assign(d_A, d_A + eq.A.size() - 1);
	d_eq.b.assign(d_b, d_b + eq.b.size() - 1);
	d_mesh.IA.assign(d_IA, d_IA + mesh.IA.size() - 1);
	d_mesh.JA.assign(d_JA, d_JA + mesh.JA.size() - 1);
	d_phi.inner.assign(d_X, d_X + phi.inner.size() - 1);
	d_pre_phi.inner.assign(d_preX, d_preX + phi.inner.size() - 1);

	int tpb = 256;
	int bpg = (N + 256 - 1) / tpb;

	double error = 0; int step = 0;
	do{
		CellField pre = phi;

		for (int it = 0; it < check_step; it++){
			cudaMemcpy(d_preX, d_X, sizeof(double)*phi.inner.size(), cudaMemcpyDeviceToDevice);
			vector_kernel << <bpg, tpb >> >(d_eq,d_phi,d_pre_phi,d_mesh);
		}
		cudaMemcpy(X, d_X, sizeof(double)*phi.inner.size(), cudaMemcpyDeviceToHost);

		if (step%check_step == 0){
			error = VectorMath<double>::rootOfSquareSum(phi.inner, pre.inner);
			cout << step << endl;
		}
		step++;
		/*
		transform(phi->inner.cbegin(), phi->inner.cend(), pre.inner.cbegin(),pre.inner.begin(), minus<T>());
		transform(pre.inner.begin(), pre.inner.end(), pre.inner.begin(), pre.inner.begin(), multiplies<T>());
		error = accumulate(pre.inner.begin(), pre.inner.end(), (T)0);//have to force convert to T type unless the result is always 0
		//error = inner_product(pre.inner.cbegin(), pre.inner.cend(), pre.inner.cbegin(),
		//	1, plus<T>(), plus<T>());
		error = sqrt(error);
		*/

	} while (error>converge_threhold&&step<max_step);

	//phi.inner . assign(X, X + N);

}