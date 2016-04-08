#ifndef LAPLACE_H_H
#define LAPLACE_H_H

#include "dataStructure.h"
#include <vector>
#include "serviceLocator.h"
class I_Laplace:public interface_t{
public:
	virtual CSR apply(vector<double> gama_f, CellField phi, Mesh mesh) = 0;
	virtual CSR apply(double  gama_f, CellField phi, Mesh mesh) = 0;
};



class Laplace :public I_Laplace,public member_t<Laplace> {
public:


	virtual CSR apply(double gama_f, CellField phi, Mesh mesh){
		vector<double> gama_field(mesh.faceNum, gama_f);
		return apply(gama_field, phi, mesh);
	}

	virtual CSR apply(vector<double> gama_f, CellField phi, Mesh mesh){
		CSR eq(mesh);
	

		for (int f = 0; f < mesh.iF.size(); f++){
			int inerF = mesh.iF[f];
			int owner = mesh.FC[inerF * 2 + 0];
			int neighbor = mesh.FC[inerF * 2 + 1];
			double S_f = mesh.F_area[inerF];
			double d = mesh.F_d[inerF];
			double diffusion_on = gama_f[inerF] * S_f / d;
			double diffusion_oo = -diffusion_on;
			double diffusion_no = diffusion_on;
			double diffusion_nn = -diffusion_on;

			eq.A[mesh.IA[owner]] += diffusion_oo;
			eq.A[mesh.ON[f * 2 + 0]] += diffusion_on;
			eq.A[mesh.IA[neighbor]] += diffusion_nn;
			eq.A[mesh.ON[f * 2 + 1]] += diffusion_no;
		}

		for (int pat = 0; pat < mesh.bF.size(); pat++) {
			for (int f = 0; f < mesh.bF[pat]->size(); f++) {
				int boundaryF = mesh.bF[pat]->at(f);
				int owner = mesh.FC[boundaryF * 2 + 0];
				double S_f = mesh.F_area[boundaryF];
				double d = mesh.F_d[boundaryF];

				double diffusion_on = gama_f[boundaryF] * S_f / d;
				double diffusion_oo = -diffusion_on;
				phi.boundaryCondition[pat][1];
				double a = phi.boundaryCondition[pat][0];
				double b1 = phi.boundaryCondition[pat][1];
				double c = phi.boundaryCondition[pat][2];
				eq.A[mesh.IA[owner]] += diffusion_oo + a / d / (a / d + b1)*diffusion_on;
				eq.b[owner] -= c / (a / d + b1)*diffusion_on;
			}
		}
		return eq;
	}


};

#endif