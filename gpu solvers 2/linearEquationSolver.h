#ifndef LINEAREQUATIONSOLVER_H_H
#define LINEAREQUATIONSOLVER_H_H


#include "dataStructure.h"
#include "mesh.h"

class I_LinearEquationSolver :public interface_t{

public:
	double converge_threhold;
	int max_step;
	int check_step;
	virtual void solve(CSR eq, CellField& phi, Mesh mesh) {
		double error = 0; int step = 0;
		do{
			CellField pre = phi;
			iterate(eq, phi, mesh, 1);
			if (step%check_step==0){
				error = VectorMath<double>::rootOfSquareSum(phi.inner, pre.inner);
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
	}
	virtual void iterate(CSR eq, CellField& phi, Mesh mesh, int iterateTimes) = 0;
};


class  GS_Solver : public I_LinearEquationSolver ,public member_t<GS_Solver>{
public:
	virtual void iterate(CSR eq, CellField& phi, Mesh mesh, int iterateTimes) {
		for (int it = 0; it < iterateTimes; ++it){
			//#pragma omp parallel for
			for (int c = 0; c < mesh.cellNum; c++) {
				double temp = 0;
				for (int i = mesh.IA[c] + 1; i < mesh.IA[c + 1]; i++) {
					temp += eq.A[i] * phi.inner[mesh.JA[i]];
				}
				phi.inner[c] = (eq.b[c] - temp) / eq.A[mesh.IA[c]];
			}
		}
	}

};


class  OpenMP_Jacobi_Solver : public I_LinearEquationSolver, public member_t< OpenMP_Jacobi_Solver>{
public:
	virtual void iterate(CSR eq, CellField& phi, Mesh mesh, int iterateTimes) {
		for (int it = 0; it < iterateTimes; ++it){
			vector<double> temp(mesh.cellNum,0);
			#pragma omp parallel for
			for (int c = 0; c < mesh.cellNum; c++) {
				for (int i = mesh.IA[c] + 1; i < mesh.IA[c + 1]; i++) {
					temp[c] += eq.A[i] * phi.inner[mesh.JA[i]];
				}
				phi.inner[c] = (eq.b[c] - temp[c]) / eq.A[mesh.IA[c]];
			}
		}
	}

};


class  GPU_Jacobi_Solver : public I_LinearEquationSolver, public member_t<GPU_Jacobi_Solver>{
public:
	virtual void solve(CSR eq, CellField& phi, Mesh mesh); 
	void iterate(CSR eq, CellField& phi, Mesh mesh, int iterateTimes){};
};


class  GPU_Jacobi_Solver_vector : public I_LinearEquationSolver, public member_t<GPU_Jacobi_Solver_vector>{


public:
	virtual void solve(CSR eq, CellField& phi, Mesh mesh);
	void iterate(CSR eq, CellField& phi, Mesh mesh, int iterateTimes);
};


#endif