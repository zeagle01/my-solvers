#ifndef LAPLACE_SOLVER_H_H
#define LAPLACE_SOLVER_H_H



#include "solver.h"
#include "Laplace.h"
#include "dataStructure.h"
 class DiffusionSolver :public Solver{
public:

	I_Laplace* diffusion;

	double gama;
	

	DiffusionSolver(Case* cas):Solver(cas) {
		
	}

	virtual void solve(){
		cas->config();
		gama = cas->configReader->readDouble("gama");
		diffusion = cas->configReader->readLaplace("Laplace_operator");

		CSR eq = diffusion->apply(gama, cas->phi[0], cas->mesh);
		cas->les[0]->solve(eq, cas->phi[0], cas->mesh);
		cas->phi[0].assignBoundary(cas->mesh);
		cas->printers[0]->print("out.dat", cas->mesh, cas->phi);

		//std::ofstream fout("a.dat");
		//copy(phi[0]->inner.cbegin(), phi[0]->inner.cend(), ostream_iterator<T>(fout, "\n") );
		//for_each(phi[0]->inner.begin(), phi[0]->inner.end(), print<T>());
	}
};



#endif